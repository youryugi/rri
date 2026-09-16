"""Time integration for the reference port: both a fixed-step RK4 (for an
apples-to-apples comparison against `rri_torch.integrate.integrate_fixed`
at the same dt/substep count) and the reference's real adaptive
Runge-Kutta-Fehlberg 4(5) scheme -- the Cash-Karp coefficients from
RRI_Mod2.f90, used exactly as RRI.f90's main loop does.

Two flavors of state are supported: a 2-D numpy array (slope's `hs`) and a
dict keyed by (i, j) river-cell coordinates (river's `hr`). Each gets its
own small set of functions since sharing one generic implementation would
obscure the (very simple) arithmetic against the Fortran it mirrors.
"""

from __future__ import annotations

import numpy as np

# Cash-Karp RKF45 coefficients, RRI_Mod2.f90::runge_mod.
A2, A3, A4, A5, A6 = 0.2, 0.3, 0.6, 1.0, 0.875
B21 = 0.2
B31, B32 = 3.0 / 40.0, 9.0 / 40.0
B41, B42, B43 = 0.3, -0.9, 1.2
B51, B52, B53, B54 = -11.0 / 54.0, 2.5, -70.0 / 27.0, 35.0 / 27.0
B61, B62, B63, B64, B65 = 1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0
C1, C3, C4, C6 = 37.0 / 378.0, 250.0 / 621.0, 125.0 / 594.0, 512.0 / 1771.0
DC1 = C1 - 2825.0 / 27648.0
DC3 = C3 - 18575.0 / 48384.0
DC4 = C4 - 13525.0 / 55296.0
DC5 = -277.0 / 14336.0
DC6 = C6 - 0.25

SAFETY, PSHRNK = 0.9, -0.25


# ---- grid (numpy array) state: used for slope's hs -------------------

def rk4_step_grid(rhs, y0: np.ndarray, dt: float) -> np.ndarray:
    def c(y):
        return np.clip(y, 0.0, None)
    k1 = rhs(c(y0))
    k2 = rhs(c(y0 + 0.5 * dt * k1))
    k3 = rhs(c(y0 + 0.5 * dt * k2))
    k4 = rhs(c(y0 + dt * k3))
    return c(y0 + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4))


def integrate_fixed_grid(rhs, y0: np.ndarray, dt: float, n_substeps: int) -> np.ndarray:
    sub_dt = dt / n_substeps
    y = y0
    for _ in range(n_substeps):
        y = rk4_step_grid(rhs, y, sub_dt)
    return y


def integrate_adaptive_grid(rhs, y0: np.ndarray, dt_outer: float, domain: np.ndarray, eps: float, ddt_min: float):
    """RRI.f90's slope adaptive-RK loop (`funcs` stages 1..6), restricted
    to the fixed outer window `dt_outer` (i.e. one `t` in the main loop).
    """
    time = 0.0
    y = y0
    ddt = dt_outer
    while True:
        if time + ddt > dt_outer:
            ddt = dt_outer - time
        while True:
            k1 = rhs(y)
            y2 = np.clip(y + B21 * ddt * k1, 0.0, None)
            k2 = rhs(y2)
            y3 = np.clip(y + ddt * (B31 * k1 + B32 * k2), 0.0, None)
            k3 = rhs(y3)
            y4 = np.clip(y + ddt * (B41 * k1 + B42 * k2 + B43 * k3), 0.0, None)
            k4 = rhs(y4)
            y5 = np.clip(y + ddt * (B51 * k1 + B52 * k2 + B53 * k3 + B54 * k4), 0.0, None)
            k5 = rhs(y5)
            y6 = np.clip(y + ddt * (B61 * k1 + B62 * k2 + B63 * k3 + B64 * k4 + B65 * k5), 0.0, None)
            k6 = rhs(y6)

            y_new = np.clip(y + ddt * (C1 * k1 + C3 * k3 + C4 * k4 + C6 * k6), 0.0, None)
            err = ddt * (DC1 * k1 + DC3 * k3 + DC4 * k4 + DC5 * k5 + DC6 * k6)
            err = np.where(domain, err, 0.0)
            errmax = np.max(np.abs(err)) / eps

            if errmax > 1.0 and ddt >= ddt_min:
                ddt = max(SAFETY * ddt * errmax ** PSHRNK, 0.5 * ddt)
                continue
            else:
                if time + ddt > dt_outer:
                    ddt = dt_outer - time
                time += ddt
                y = y_new
                break
        if time >= dt_outer:
            break
    return y


# ---- dict (river, keyed by (i, j)) state ------------------------------

def _dict_op(fn, *ds):
    keys = ds[0].keys()
    return {k: fn(*(d[k] for d in ds)) for k in keys}


def rk4_step_dict(rhs, y0: dict, dt: float) -> dict:
    def c(y):
        return {k: max(v, 0.0) for k, v in y.items()}
    k1 = rhs(c(y0))
    k2 = rhs(c(_dict_op(lambda a, b: a + 0.5 * dt * b, y0, k1)))
    k3 = rhs(c(_dict_op(lambda a, b: a + 0.5 * dt * b, y0, k2)))
    k4 = rhs(c(_dict_op(lambda a, b: a + dt * b, y0, k3)))
    y1 = {k: y0[k] + (dt / 6.0) * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) for k in y0}
    return c(y1)


def integrate_fixed_dict(rhs, y0: dict, dt: float, n_substeps: int) -> dict:
    sub_dt = dt / n_substeps
    y = y0
    for _ in range(n_substeps):
        y = rk4_step_dict(rhs, y, sub_dt)
    return y


def integrate_adaptive_dict(rhs, y0: dict, dt_outer: float, eps: float, ddt_min: float):
    time = 0.0
    y = y0
    ddt = dt_outer
    keys = list(y0.keys())

    def clip0(d):
        return {k: max(v, 0.0) for k, v in d.items()}

    while True:
        if time + ddt > dt_outer:
            ddt = dt_outer - time
        while True:
            k1 = rhs(y)
            y2 = clip0({k: y[k] + B21 * ddt * k1[k] for k in keys})
            k2 = rhs(y2)
            y3 = clip0({k: y[k] + ddt * (B31 * k1[k] + B32 * k2[k]) for k in keys})
            k3 = rhs(y3)
            y4 = clip0({k: y[k] + ddt * (B41 * k1[k] + B42 * k2[k] + B43 * k3[k]) for k in keys})
            k4 = rhs(y4)
            y5 = clip0({k: y[k] + ddt * (B51 * k1[k] + B52 * k2[k] + B53 * k3[k] + B54 * k4[k]) for k in keys})
            k5 = rhs(y5)
            y6 = clip0({k: y[k] + ddt * (B61 * k1[k] + B62 * k2[k] + B63 * k3[k] + B64 * k4[k] + B65 * k5[k]) for k in keys})
            k6 = rhs(y6)

            y_new = clip0({k: y[k] + ddt * (C1 * k1[k] + C3 * k3[k] + C4 * k4[k] + C6 * k6[k]) for k in keys})
            err = {k: ddt * (DC1 * k1[k] + DC3 * k3[k] + DC4 * k4[k] + DC5 * k5[k] + DC6 * k6[k]) for k in keys}
            errmax = max(abs(v) for v in err.values()) / eps if keys else 0.0

            if errmax > 1.0 and ddt >= ddt_min:
                ddt = max(SAFETY * ddt * errmax ** PSHRNK, 0.5 * ddt)
                continue
            else:
                if time + ddt > dt_outer:
                    ddt = dt_outer - time
                time += ddt
                y = y_new
                break
        if time >= dt_outer:
            break
    return y

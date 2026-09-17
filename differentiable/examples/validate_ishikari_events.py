"""Calibrate on one Ishikari event and evaluate an independent event.

This experiment addresses two hydrologic omissions in the first Ishikari
smoke run without pretending the currently implemented physics can produce
groundwater baseflow:

* MERV-Jp PET is applied through RRI's existing ET process.
* Green-Ampt infiltration is enabled with the same finite soil-storage cap
  construction used by RRI (``soildepth * gammaa``).
* A constant observed pre-event baseflow is added to simulated *stormflow*
  only for hydrograph comparison.  It is deliberately not inserted into
  river storage, where it would drain away without a groundwater source.

The defaults were selected on the 1991 event, then frozen before running
the 1996 event.  The independent result is intentionally retained even
though it is poor: it demonstrates that a fixed effective infiltration
rate does not generalize when a distributed intensity-sensitive model is
driven by daily basin-mean precipitation.
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import sys
import time
from dataclasses import dataclass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from examples.run_ishikari import CATCHMENT_AREA_KM2, PROJECT_DIR, build_ishikari_scenario
from rri_torch.model import RRIModel


@dataclass(frozen=True)
class Event:
    name: str
    start: str
    days: int
    baseflow_days: int
    role: str


EVENTS = {
    "calibration": Event("calibration", "1991-08-25", 15, 10, "parameter-selection event"),
    "validation": Event("validation", "1996-09-10", 27, 10, "independent validation event"),
}


def to_device(obj, device):
    if torch.is_tensor(obj):
        return obj.to(device)
    if dataclasses.is_dataclass(obj):
        return dataclasses.replace(obj, **{
            field.name: to_device(getattr(obj, field.name), device)
            for field in dataclasses.fields(obj)
        })
    return obj


def load_event(event: Event):
    path = os.path.join(PROJECT_DIR, "obs", "varssim005_ver2_1.csv")
    frame = pd.read_csv(path)
    frame["date"] = pd.to_datetime(frame[["Year", "Month", "Day"]])
    frame = frame.set_index("date").loc[event.start:].iloc[:event.days]
    if len(frame) != event.days:
        raise ValueError(f"{event.name}: requested {event.days} days, found {len(frame)}")
    observed_m3s = frame["Obs flow"].to_numpy() * CATCHMENT_AREA_KM2 * 1000.0 / 86400.0
    return (
        frame.index,
        frame["Precip"].to_numpy(),
        frame["PET"].to_numpy(),
        observed_m3s,
    )


def configure_params(params, ns_slope: float, ns_river: float, ksv_mm_day: float):
    params.slope.ns = torch.full_like(params.slope.ns, ns_slope)
    params.river.ns_river = torch.full_like(params.river.ns_river, ns_river)
    params.ksv = torch.full_like(params.ksv, ksv_mm_day * 1e-3 / 86400.0)
    if ksv_mm_day > 0.0:
        params.infilt_limit = params.slope.soildepth * params.slope.gammaa
    else:
        params.infilt_limit = torch.full_like(params.infilt_limit, -1.0)
    params.evp_switch = 1


def metrics(simulated: np.ndarray, observed: np.ndarray) -> tuple[float, float]:
    residual = simulated - observed
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    denominator = float(np.sum((observed - observed.mean()) ** 2))
    nse = 1.0 - float(np.sum(residual ** 2)) / denominator
    return nse, rmse


def simulate_event(event: Event, args) -> dict:
    device = torch.device(args.device)
    dtype = torch.float64
    grid, params, outlet = build_ishikari_scenario(dtype=dtype)
    dates, precipitation, pet, observed = load_event(event)

    grid = grid.to(device)
    params.slope = to_device(params.slope, device)
    params.river = to_device(params.river, device)
    params.ksv = params.ksv.to(device)
    params.faif = params.faif.to(device)
    params.infilt_limit = params.infilt_limit.to(device)
    configure_params(params, args.ns_slope, args.ns_river, args.ksv_mm_day)

    model = RRIModel(grid, adaptive=True, track_qr_avg=True)
    hs, gampt_ff, hr = model.initial_state(dtype=dtype)
    rain = torch.zeros((grid.ny, grid.nx), dtype=dtype, device=device)
    pet_grid = torch.zeros_like(rain)
    dt = 600.0
    steps_per_day = 144

    stormflow = []
    hs_max = []
    hr_max = []
    infiltration_mm = []
    actual_et_mm = []
    started = time.time()
    with torch.no_grad():
        for day in range(event.days):
            rain.zero_()
            rain[grid.domain] = precipitation[day] * 1e-3 / 86400.0
            pet_grid.zero_()
            pet_grid[grid.domain] = pet[day] * 1e-3 / 86400.0
            q_day = []
            hs_day = []
            hr_day = []
            infiltration_day = 0.0
            et_day = 0.0
            for _ in range(steps_per_day):
                hs, gampt_ff, hr, diag = model.step(
                    hs, gampt_ff, hr, params, rain, pet_grid, dt,
                )
                q_day.append(float(diag.qr_avg[outlet]))
                hs_day.append(float(hs[grid.domain].max()))
                hr_day.append(float(hr.max()))
                infiltration_day += float((diag.infilt_rate[grid.domain] * dt).mean())
                et_day += float((diag.aevp[grid.domain] * dt).mean())
            if device.type == "cuda":
                torch.cuda.synchronize()
            stormflow.append(float(np.mean(q_day)))
            hs_max.append(max(hs_day))
            hr_max.append(max(hr_day))
            infiltration_mm.append(infiltration_day * 1000.0)
            actual_et_mm.append(et_day * 1000.0)
            print(
                f"[{event.name}] {dates[day].date()} P={precipitation[day]:5.1f}mm "
                f"Qstorm={stormflow[-1]:7.1f}m3/s hs={hs_max[-1]:.3f}m "
                f"hr={hr_max[-1]:.3f}m",
                flush=True,
            )

    stormflow = np.asarray(stormflow)
    baseflow = float(np.median(observed[:event.baseflow_days]))
    simulated = stormflow + baseflow
    nse, rmse = metrics(simulated, observed)
    elapsed = time.time() - started
    result = {
        "dates": np.asarray(dates.astype(str), dtype="U10"),
        "precip_mm_day": precipitation,
        "pet_mm_day": pet,
        "observed_m3s": observed,
        "stormflow_m3s": stormflow,
        "simulated_m3s": simulated,
        "baseflow_m3s": np.asarray(baseflow),
        "hs_max_m": np.asarray(hs_max),
        "hr_max_m": np.asarray(hr_max),
        "infiltration_mm": np.asarray(infiltration_mm),
        "actual_et_mm": np.asarray(actual_et_mm),
        "nse": np.asarray(nse),
        "rmse_m3s": np.asarray(rmse),
        "elapsed_seconds": np.asarray(elapsed),
        "ns_slope": np.asarray(args.ns_slope),
        "ns_river": np.asarray(args.ns_river),
        "ksv_mm_day": np.asarray(args.ksv_mm_day),
    }
    print(
        f"[{event.name}] baseflow={baseflow:.1f}m3/s NSE={nse:.3f} "
        f"RMSE={rmse:.1f}m3/s peak={simulated.max():.1f}m3/s "
        f"on {dates[int(simulated.argmax())].date()} ({elapsed:.1f}s)"
    )
    return result


def plot_results(results: dict[str, dict], output_path: str) -> None:
    fig, axes = plt.subplots(len(results), 1, figsize=(10, 4.2 * len(results)), squeeze=False)
    for ax, (name, result) in zip(axes[:, 0], results.items()):
        days = np.arange(len(result["dates"])) + 0.5
        ax.plot(days, result["observed_m3s"], "o-", color="#d62728", label="Observed")
        ax.plot(days, result["simulated_m3s"], "s-", color="#1f77b4",
                label="RRI stormflow + constant baseflow")
        ax.axhline(float(result["baseflow_m3s"]), color="#1f77b4", ls=":", lw=1,
                   label="Estimated constant baseflow")
        rain_ax = ax.twinx()
        rain_ax.bar(days, result["precip_mm_day"], width=0.8, color="0.7", alpha=0.35)
        rain_ax.invert_yaxis()
        rain_ax.set_ylabel("Precipitation [mm/day]")
        ax.set_ylabel(r"Discharge [m$^3$/s]")
        ax.set_xlabel(f"Days since {result['dates'][0]}")
        ax.set_title(
            f"{name.capitalize()}: NSE={float(result['nse']):.3f}, "
            f"RMSE={float(result['rmse_m3s']):.1f} m3/s"
        )
        lines, labels = ax.get_legend_handles_labels()
        ax.legend(lines, labels, loc="upper left", fontsize=9)
    fig.suptitle("Ishikari event calibration and independent validation", fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--event", choices=("calibration", "validation", "both"), default="both")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--ns-slope", type=float, default=0.2)
    parser.add_argument("--ns-river", type=float, default=0.03)
    parser.add_argument("--ksv-mm-day", type=float, default=1.0)
    parser.add_argument("--output-dir", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results",
    ))
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    selected = EVENTS.values() if args.event == "both" else (EVENTS[args.event],)
    results = {}
    for event in selected:
        result = simulate_event(event, args)
        results[event.name] = result
        np.savez_compressed(
            os.path.join(args.output_dir, f"ishikari_{event.name}_event.npz"), **result,
        )
    plot_results(results, os.path.join(args.output_dir, "fig_ishikari_event_validation.png"))


if __name__ == "__main__":
    main()

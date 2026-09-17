"""Derive RRI rectangular-channel geometry from Ishikari cross sections.

The Hokkaido Development Bureau survey bundle contains two different kinds
of tables that are easy to confuse:

* ``WZAA2001.CSV`` is a *longitudinal* table.  Its second column is the
  along-channel distance from the preceding survey station, not channel
  width.
* ``WZAA3001.XLS`` and the corresponding ``WZAA4xxx.CSV`` files contain
  cross-section metadata and (offset, elevation) point clouds.

For each ordinary (non-auxiliary) cross section this script uses the lower
of the surveyed left/right floodplain elevations as bankfull level.  It
selects the connected below-bankfull part containing the thalweg, integrates
its area, and constructs an equivalent RRI rectangle with:

    depth = bankfull elevation - thalweg elevation
    width = cross-sectional area / depth

Thus the simplified rectangle preserves surveyed bankfull area.  Auxiliary
sections (bridges and other special stations) are excluded because they are
not representative reach geometry.  Later survey years overwrite older
ones where their KP ranges overlap.

Example::

    python build_channel_geometry.py \
      "/path/to/01.石狩川/③定期横断測量測量成果（数値データ）" \
      --output kp_channel_geometry.csv
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


SURVEY_YEARS = ("H24①", "H24②", "H26", "R03")


def _crossing(x1: float, z1: float, x2: float, z2: float, level: float) -> float:
    if z2 == z1:
        return 0.5 * (x1 + x2)
    return x1 + (level - z1) * (x2 - x1) / (z2 - z1)


def equivalent_rectangle(points: np.ndarray, bank_elevation: float) -> dict[str, float] | None:
    """Return geometry of the thalweg-connected wetted cross section."""
    order = np.argsort(points[:, 0], kind="stable")
    x = points[order, 0]
    z = points[order, 1]
    thalweg_index = int(np.argmin(z))
    thalweg = float(z[thalweg_index])
    depth = bank_elevation - thalweg
    if depth <= 0.0:
        return None

    left = thalweg_index
    while left > 0 and z[left - 1] < bank_elevation:
        left -= 1
    right = thalweg_index
    while right + 1 < len(z) and z[right + 1] < bank_elevation:
        right += 1

    x_left = (
        _crossing(x[left - 1], z[left - 1], x[left], z[left], bank_elevation)
        if left > 0 else float(x[left])
    )
    x_right = (
        _crossing(x[right], z[right], x[right + 1], z[right + 1], bank_elevation)
        if right + 1 < len(z) else float(x[right])
    )
    section_x = np.r_[x_left, x[left:right + 1], x_right]
    section_z = np.r_[bank_elevation, z[left:right + 1], bank_elevation]
    area = float(np.trapezoid(np.maximum(bank_elevation - section_z, 0.0), section_x))
    return {
        "depth_m": depth,
        "top_width_m": x_right - x_left,
        "area_m2": area,
        "equivalent_width_m": area / depth,
        "thalweg_m": thalweg,
        "bank_elevation_m": bank_elevation,
    }


def _metadata_by_kp(workbook: str) -> dict[float, pd.Series]:
    sheet = pd.read_excel(workbook, header=None)
    out: dict[float, pd.Series] = {}
    for _, row in sheet.iloc[11:].iterrows():
        numeric = pd.to_numeric(row, errors="coerce")
        if pd.isna(numeric.iloc[0]):
            continue
        out[round(float(numeric.iloc[0]), 3)] = numeric
    return out


def _read_profile(path: str) -> tuple[float, int, np.ndarray] | None:
    with open(path, encoding="cp932", errors="replace") as stream:
        lines = stream.readlines()
    header = lines[0].strip().split(",")
    try:
        kp = round(float(header[0]), 3)
        auxiliary = int(float(header[11]))
    except (ValueError, IndexError):
        return None

    points = []
    for line in lines[1:]:
        fields = line.strip().split(",")
        if len(fields) < 3:
            continue
        try:
            points.append((float(fields[1]), float(fields[2])))
        except ValueError:
            continue
    if not points:
        return None
    return kp, auxiliary, np.asarray(points, dtype=float)


def build_geometry_table(survey_root: str | os.PathLike[str]) -> pd.DataFrame:
    by_kp: dict[float, dict] = {}
    for year in SURVEY_YEARS:
        year_dir = os.path.join(survey_root, year)
        workbook = os.path.join(year_dir, "WZAA3001.XLS")
        if not os.path.exists(workbook):
            continue
        metadata = _metadata_by_kp(workbook)
        for profile_path in sorted(glob.glob(os.path.join(year_dir, "WZAA4*.CSV"))):
            profile = _read_profile(profile_path)
            if profile is None:
                continue
            kp, auxiliary, points = profile
            if auxiliary != 0 or kp not in metadata:
                continue

            row = metadata[kp]
            floodplain = [float(v) for v in (row.iloc[10], row.iloc[11]) if pd.notna(v)]
            if not floodplain:
                continue
            geometry = equivalent_rectangle(points, min(floodplain))
            if geometry is None:
                continue

            # The point-cloud minimum and workbook thalweg column must
            # describe the same section.  Treat a disagreement as bad input.
            workbook_thalweg = float(row.iloc[8])
            if not np.isclose(geometry["thalweg_m"], workbook_thalweg, atol=1e-6):
                raise ValueError(
                    f"Thalweg mismatch at KP {kp}: profile={geometry['thalweg_m']}, "
                    f"workbook={workbook_thalweg}"
                )
            by_kp[kp] = {
                "kp": kp,
                **geometry,
                "survey_year": year,
                "source_file": os.path.basename(profile_path),
            }

    columns = [
        "kp", "depth_m", "equivalent_width_m", "top_width_m", "area_m2",
        "thalweg_m", "bank_elevation_m", "survey_year", "source_file",
    ]
    return pd.DataFrame(by_kp.values(), columns=columns).sort_values("kp").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("survey_root", help="Directory containing H24①/H24②/H26/R03")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    table = build_geometry_table(args.survey_root)
    if args.output:
        table.to_csv(args.output, index=False)
        print(f"wrote {len(table)} sections to {args.output}")
    else:
        print(table.to_csv(index=False), end="")


if __name__ == "__main__":
    main()

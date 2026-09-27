#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import shutil
import glob
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# Try to import readline for tab-completion (standard on Unix, requires pyreadline on Windows)
try:
    import readline
    READLINE_AVAILABLE = True
except ImportError:
    READLINE_AVAILABLE = False


# Default to Final.csv in the Current Working Directory
DEFAULT_FINAL_CSV = Path("Final.csv")

DEFAULT_DISTANCE_KPC = 5.3
DEFAULT_INCLINATION_DEG = 10.0

# Diskbb correction constants:
# Shimura & Takahara-style spectral hardening factor.
DEFAULT_DISK_FCOL = 1.7
DEFAULT_DISK_BOUNDARY_CORRECTION = 0.412

# Representative color-correction factor for helium-rich NS atmospheres.
DEFAULT_BBODY_FCOL = 1.40


def path_completer(text, state):
    """Auto-completes file paths for the input prompt in the CWD."""
    matches = glob.glob(text + "*")
    return matches[state] if state < len(matches) else None


def ask_float(label: str, default: float) -> float:
    """Prompt for a float while allowing Enter to keep the default value."""
    while True:
        text = input(f"{label} [{default:g}]: ").strip()
        if not text:
            return default
        try:
            return float(text)
        except ValueError:
            print("Please enter a number, or press Enter to keep the default.")


def ask_path(label: str, default: Path) -> Path:
    """Prompt for a path while allowing Enter to keep the default value."""
    if READLINE_AVAILABLE:
        # Enable tab completion for paths
        readline.set_completer_delims('\t\n')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(path_completer)

    text = input(f"{label} [{default}]: ").strip()

    if READLINE_AVAILABLE:
        # Disable completion for subsequent prompts (like floats)
        readline.set_completer(None)

    return Path(text).expanduser() if text else default


def parse_triplet(cell: object) -> tuple[float, float, float]:
    """Parse a Final.csv cell into value, negative error, positive error."""
    text = str(cell).strip()
    if not text:
        return np.nan, np.nan, np.nan

    parts = [part.strip() for part in text.split(",")]
    values: list[float] = []

    for part in parts[:3]:
        try:
            values.append(float(part))
        except ValueError:
            values.append(np.nan)

    while len(values) < 3:
        values.append(np.nan)

    value, neg_error, pos_error = values
    return value, abs(neg_error), abs(pos_error)


def format_triplet(value: float, neg_error: float, pos_error: float) -> str:
    """Write output cells in the same compact style as Final.csv."""
    if not np.isfinite(value):
        return ""
    return f"{value:.8g}, {neg_error:.8g}, {pos_error:.8g}"


def find_row(table: pd.DataFrame, candidates: list[str]) -> str:
    for row in candidates:
        if row in table.index:
            return row
    raise KeyError(f"None of these rows were found: {', '.join(candidates)}")


def radius_from_norm(
    norm: float,
    norm_err_minus: float,
    norm_err_plus: float,
    scale: float,
) -> tuple[float, float, float]:
    """Propagate asymmetric normalization errors exactly through R=scale*sqrt(K)."""
    if not np.isfinite(norm) or norm < 0:
        return np.nan, np.nan, np.nan

    err_minus = 0.0 if not np.isfinite(norm_err_minus) else max(norm_err_minus, 0.0)
    err_plus = 0.0 if not np.isfinite(norm_err_plus) else max(norm_err_plus, 0.0)

    radius = scale * math.sqrt(norm)
    lower_norm = max(norm - err_minus, 0.0)
    upper_norm = norm + err_plus

    radius_lower = scale * math.sqrt(lower_norm)
    radius_upper = scale * math.sqrt(upper_norm)

    return radius, radius - radius_lower, radius_upper - radius


def scaled_triplet(
    values: tuple[float, float, float],
    factor: float,
) -> tuple[float, float, float]:
    value, neg_error, pos_error = values
    return value * factor, neg_error * factor, pos_error * factor


def build_radius_rows(
    table: pd.DataFrame,
    distance_kpc: float,
    inclination_deg: float,
    disk_fcol: float,
    disk_boundary_correction: float,
    bbody_fcol: float,
) -> dict[str, list[str]]:
    """Return new Final.csv rows as {row_name: list_of_cells}."""
    bbody_row = find_row(table, ["bbodyrad_norm", "norm"])
    disk_row = find_row(table, ["diskbb_norm", "norm_1"])

    d10 = distance_kpc / 10.0
    cos_i = math.cos(math.radians(inclination_deg))

    if d10 <= 0:
        raise ValueError("distance_kpc must be positive")
    if cos_i <= 0:
        raise ValueError("inclination_deg must be less than 90 degrees")
    if disk_fcol <= 0 or disk_boundary_correction <= 0 or bbody_fcol <= 0:
        raise ValueError("Correction factors must be positive")

    # XSPEC bbodyrad: K = (R_bb,app / D_10)^2.
    bbody_scale = d10
    bbody_color_factor = bbody_fcol

    # XSPEC diskbb: K = (R_in,app / D_10)^2 cos(i).
    disk_app_scale = d10 / math.sqrt(cos_i)
    disk_correction = disk_fcol**2 * disk_boundary_correction

    rows = {
        "bbodyrad_R_app_km": [],
        "bbodyrad_R_colorcorr_km": [],
        "diskbb_Rin_app_km": [],
        "diskbb_Rin_corr_km": [],
    }

    for obsid in table.columns:
        bbody_norm = parse_triplet(table.loc[bbody_row, obsid])
        disk_norm = parse_triplet(table.loc[disk_row, obsid])

        bbody_app = radius_from_norm(*bbody_norm, scale=bbody_scale)
        bbody_corr = scaled_triplet(bbody_app, bbody_color_factor)

        disk_app = radius_from_norm(*disk_norm, scale=disk_app_scale)
        disk_corr = scaled_triplet(disk_app, disk_correction)

        rows["bbodyrad_R_app_km"].append(format_triplet(*bbody_app))
        rows["bbodyrad_R_colorcorr_km"].append(format_triplet(*bbody_corr))
        rows["diskbb_Rin_app_km"].append(format_triplet(*disk_app))
        rows["diskbb_Rin_corr_km"].append(format_triplet(*disk_corr))

    return rows


def update_final_csv(
    final_csv: Path,
    distance_kpc: float,
    inclination_deg: float,
    disk_fcol: float,
    disk_boundary_correction: float,
    bbody_fcol: float,
    no_backup: bool,
) -> Path | None:
    table = pd.read_csv(final_csv, sep="\t", index_col=0, dtype=str, keep_default_na=False)
    table.index = table.index.astype(str)
    table.columns = table.columns.astype(str)

    new_rows = build_radius_rows(
        table=table,
        distance_kpc=distance_kpc,
        inclination_deg=inclination_deg,
        disk_fcol=disk_fcol,
        disk_boundary_correction=disk_boundary_correction,
        bbody_fcol=bbody_fcol,
    )

    for row_name, cells in new_rows.items():
        table.loc[row_name, :] = cells

    backup_path = None
    if not no_backup:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = final_csv.with_name(f"{final_csv.stem}.backup_before_radii_{timestamp}{final_csv.suffix}")
        shutil.copy2(final_csv, backup_path)

    table.to_csv(final_csv, sep="\t", index=True, lineterminator="\n")
    return backup_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-csv", type=Path, default=DEFAULT_FINAL_CSV)
    parser.add_argument("--distance-kpc", type=float, default=DEFAULT_DISTANCE_KPC)
    parser.add_argument("--inclination-deg", type=float, default=DEFAULT_INCLINATION_DEG)
    parser.add_argument("--disk-fcol", type=float, default=DEFAULT_DISK_FCOL)
    parser.add_argument(
        "--disk-boundary-correction",
        type=float,
        default=DEFAULT_DISK_BOUNDARY_CORRECTION,
    )
    parser.add_argument("--bbody-fcol", type=float, default=DEFAULT_BBODY_FCOL)
    parser.add_argument(
        "--ask",
        action="store_true",
        help="Ask for advanced colour factors interactively.",
    )
    parser.add_argument("--no-backup", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Always ask for the CSV file name upfront
    args.final_csv = ask_path("CSV file name", args.final_csv)

    # Always ask for Distance and Inclination upfront
    args.distance_kpc = ask_float("Distance in kpc", args.distance_kpc)
    args.inclination_deg = ask_float("Inclination in degrees", args.inclination_deg)

    # Only ask for advanced parameters if --ask flag is provided
    if args.ask:
        args.disk_fcol = ask_float("Disk colour factor", args.disk_fcol)
        args.disk_boundary_correction = ask_float(
            "Disk inner-boundary correction",
            args.disk_boundary_correction,
        )
        args.bbody_fcol = ask_float("Blackbody colour factor", args.bbody_fcol)

    backup_path = update_final_csv(
        final_csv=args.final_csv,
        distance_kpc=args.distance_kpc,
        inclination_deg=args.inclination_deg,
        disk_fcol=args.disk_fcol,
        disk_boundary_correction=args.disk_boundary_correction,
        bbody_fcol=args.bbody_fcol,
        no_backup=args.no_backup,
    )

    print(f"\nUpdated: {args.final_csv.absolute()}")
    if backup_path is not None:
        print(f"Backup:  {backup_path.absolute()}")
    print("Added/updated rows:")
    print("  bbodyrad_R_app_km")
    print("  bbodyrad_R_colorcorr_km")
    print("  diskbb_Rin_app_km")
    print("  diskbb_Rin_corr_km")
    print(
        f"\nAssumptions: "
        f"D={args.distance_kpc:g} kpc, i={args.inclination_deg:g} deg, "
        f"disk_fcol={args.disk_fcol:g}, "
        f"disk_boundary={args.disk_boundary_correction:g}, "
        f"bbody_fcol={args.bbody_fcol:g}"
    )


if __name__ == "__main__":
    main()

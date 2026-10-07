"""
Verify that the lateral boundaries are interpolated onto the pressure of each boundary
time, and that the interface pressure uses the hybrid vertical coordinate.

Before this, every boundary time was remapped onto the interface pressure of the
wrfinput passed on the command line. With one wrfbdy per wrfinput (one per cycle) that
is only wrong for the end of the interval; with a wrfbdy spanning many days (WRF
restart cycling) every record after the first used the wrong pressure. The interface
pressure also used the terrain-following formula ZNW·mu + p_top, which is off at upper
levels and over high terrain with hybrid_opt = 2.

Needs real files: a wrfbdy spanning several boundary times (LONG_DIR, from real.exe run
over the whole period) and per-cycle wrfinput/wrfbdy for the same period, domain and
meteorology (CYCLE_DIR). Checks:

1. MU_B* of the long wrfbdy at record k equals MU on the boundary of the cycle k
   wrfinput, and the reconstructed MU at the end of the last record equals the next
   cycle's wrfinput.
2. The interface pressure rebuilt from those equals the cycle k wrfinput's.
3. How much the hybrid formula differs from the old ZNW formula.
4. End to end: interpolating the chemistry onto the long wrfbdy (with the cycle 0
   wrfinput) gives the same boundaries as interpolating each per-cycle wrfbdy with its
   own wrfinput, values and tendencies.

Run with the project venv:
    uv run python exploration/verify_boundary_pressure.py LONG_DIR CYCLE_DIR CAMS_DIR SPECIES_MAP
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import netCDF4 as nc
import numpy as np

from interpolator_for_wrfchem import utils
from interpolator_for_wrfchem.wrf import WRFBoundary, WRFInput, hybrid_pressure

BOUNDARIES = ["BXS", "BXE", "BYS", "BYE"]


def rel(a, b) -> float:
    a = np.asarray(a, "f8")
    b = np.asarray(b, "f8")
    return float(np.abs(a - b).max() / max(np.abs(b).max(), 1e-30))


def boundary_row(field: np.ndarray, bdy: str) -> np.ndarray:
    """Outermost row of a (..., south_north, west_east) field, like wrfbdy row 0"""
    return {
        "BXS": field[..., :, 0],
        "BXE": field[..., :, -1],
        "BYS": field[..., 0, :],
        "BYE": field[..., -1, :],
    }[bdy]


def check_pressure(long_dir: Path, cycle_dir: Path) -> None:
    bdy = WRFBoundary(long_dir / "wrfbdy_d01")
    wrf0 = WRFInput(long_dir / "wrfinput_d01_cycle_0", read_only=True)
    ds0 = wrf0.get_dataset()
    print(f"long wrfbdy: {len(bdy.times) - 1} records + end time, {bdy.times[0]} .. {bdy.times[-1]}")

    worst_mu, worst_p = 0.0, 0.0
    for t_idx in range(len(bdy.times)):
        path = cycle_dir / f"wrfinput_d01_cycle_{t_idx}"
        if not path.exists():
            print(f"  no {path.name}, stopping at time {t_idx}")
            break
        wrf_k = WRFInput(path, read_only=True)
        assert wrf_k.time == bdy.times[t_idx], (wrf_k.time, bdy.times[t_idx])
        ds_k = wrf_k.get_dataset()
        mu_k = wrf_k.nc_file["MU"][0]
        for b in BOUNDARIES:
            mu = bdy.boundary_mu(t_idx, b)
            worst_mu = max(worst_mu, float(np.abs(mu - boundary_row(mu_k, b)).max()))
            profile = utils.get_boundary_profile(ds0, b)
            pres_hf = utils.boundary_pressure_hf(profile, mu)
            expected = utils.get_boundary_profile(ds_k, b)["pres_hf"]
            worst_p = max(worst_p, float(np.abs(pres_hf - expected).max()))
        kind = "end of last record" if t_idx == len(bdy.times) - 1 else f"record {t_idx}"
        print(f"  time {t_idx} ({kind}): max |MU - cycle wrfinput| so far {worst_mu:.3g} Pa, "
              f"max |pres_hf - cycle wrfinput| so far {worst_p:.3g} hPa")
        wrf_k.close()

    # How much the hybrid formula moves the interfaces compared to ZNW·mu + p_top
    mu_total = wrf0.nc_file["MU"][0] + wrf0.nc_file["MUB"][0]
    znw = wrf0.nc_file["ZNW"][0]
    old = hybrid_pressure(znw, np.zeros_like(znw), mu_total, ds0.attrs["P_TOP"])
    diff = np.abs(ds0["pres_hf"].to_numpy() - old)
    k = np.unravel_index(diff.argmax(), diff.shape)
    hgt = wrf0.nc_file["HGT"][0][k[1:]]
    print(f"hybrid vs ZNW interface pressure: max {diff.max():.1f} hPa at interface {k[0]} "
          f"({ds0['pres_hf'].to_numpy()[k]:.0f} hPa, terrain {hgt:.0f} m), "
          f"median {np.median(diff):.2f} hPa")
    wrf0.close()
    bdy.close()


def interpolate(cams: Path, species_map: Path, wrfinput: Path, wrfbdy: Path, no_ic: bool):
    args = ["interpolator-for-wrfchem", "cams_global_forecasts_pl", str(cams), str(species_map),
            str(wrfinput), f"--wrfbdy={wrfbdy}"]
    if no_ic:
        args.append("--no-ic")
    res = subprocess.run(args, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(res.stdout + res.stderr)


def check_end_to_end(long_dir: Path, cycle_dir: Path, cams: Path, species_map: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="verify_bdy_pres") as tmp:
        tmp = Path(tmp)
        shutil.copy(long_dir / "wrfbdy_d01", tmp / "wrfbdy_long")
        shutil.copy(long_dir / "wrfinput_d01_cycle_0", tmp / "wrfinput_0")
        interpolate(cams, species_map, tmp / "wrfinput_0", tmp / "wrfbdy_long", no_ic=True)

        with nc.Dataset(tmp / "wrfbdy_long") as long_bdy:
            long_bdy.set_auto_mask(False)
            n_records = long_bdy.dimensions["Time"].size
            chem = species(species_map)
            worst_value, worst_tend = {}, {}
            for k in range(n_records):
                shutil.copy(cycle_dir / f"wrfbdy_d01_cycle_{k}", tmp / f"wrfbdy_{k}")
                shutil.copy(cycle_dir / f"wrfinput_d01_cycle_{k}", tmp / f"wrfinput_{k}")
                interpolate(cams, species_map, tmp / f"wrfinput_{k}", tmp / f"wrfbdy_{k}", no_ic=True)
                with nc.Dataset(tmp / f"wrfbdy_{k}") as cyc:
                    cyc.set_auto_mask(False)
                    for v in (v for v in cyc.variables if v.split("_B")[0] in chem):
                        target = worst_tend if "_BT" in v else worst_value
                        d = rel(long_bdy[v][k], cyc[v][0])
                        if d >= target.get(v, (0.0, k))[0]:
                            target[v] = (d, k)
        for label, worst in (("values", worst_value), ("tendencies", worst_tend)):
            name, (d, k) = max(worst.items(), key=lambda x: x[1][0])
            print(f"long vs per-cycle chemistry boundary {label}: {len(worst)} variables, "
                  f"max relative difference {d:.2e} ({name}, record {k})")


def species(species_map: Path) -> set[str]:
    import tomllib

    return set(tomllib.loads(species_map.read_text())["species_map"])


def main():
    long_dir, cycle_dir, cams, species_map = (Path(a) for a in sys.argv[1:5])
    check_pressure(long_dir, cycle_dir)
    check_end_to_end(long_dir, cycle_dir, cams, species_map)


if __name__ == "__main__":
    main()

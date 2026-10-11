"""Regenerate the pictures used in the README.  Run:  python docs/make_readme_figures.py

This needs no internet and no Skyfield. It runs the real scripts (Ground_Track.py, ISS.py,
TLE.py) but swaps in tests/fake_skyfield.py, a simplified circular-orbit model, and orbital
elements taken from the original project's TLEs. So these pictures show what the output
LOOKS like; they are not real predictions, and the map has no coastlines. Run the scripts
themselves (with Skyfield and Cartopy installed) for real, current results.
"""

import contextlib
import io
import math
import pathlib
import sys
from datetime import datetime, timezone
from unittest import mock

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fake_skyfield

fake_skyfield.install()

import Ground_Track
import ISS
import TLE
from fake_skyfield import build_tle
from Orbit_Tools import EARTH_MU_KM3_S2, EARTH_RADIUS_KM

OUT = ROOT / "docs" / "images"
OUT.mkdir(parents=True, exist_ok=True)
START = "2026-01-01T00:00"
EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)

SATS = {
    "iss": build_tle("ISS (ZARYA)", 25544, 51.6404, 259.2792, 0.0003751, 161.8425, 37.3500, 15.50019787, EPOCH),
    "noaa15": build_tle("NOAA 15", 25338, 98.7406, 210.1389, 0.0011530, 75.7009, 284.4304, 14.25982377, EPOCH),
    "gps": build_tle("GPS (PRN 10)", 40105, 55.1505, 35.0374, 0.0107909, 324.8964, 34.9501, 2.00567658, EPOCH),
}


def run(module, argv):
    with mock.patch.object(module, "get_tle", side_effect=lambda key, *a, **k: SATS[key.lower()]), \
            contextlib.redirect_stdout(io.StringIO()) as captured:
        code = module.main(argv)
    assert code == 0
    return captured.getvalue()


def orbit_scale_picture():
    """The orbits of the satellites in this project, drawn to the same scale as the Earth."""
    def radius_from_period(seconds):
        return (EARTH_MU_KM3_S2 * (seconds / (2 * math.pi)) ** 2) ** (1 / 3)

    orbits = [("ISS", EARTH_RADIUS_KM + 420, "#d1495b", "93 min"),
              ("NOAA polar orbiters", EARTH_RADIUS_KM + 820, "#2e86ab", "101 min"),
              ("GPS", radius_from_period(43082), "#3c8d5a", "11 h 58 min"),
              ("Geostationary", radius_from_period(86164), "#888888", "23 h 56 min")]
    fig, ax = plt.subplots(figsize=(7.5, 8.2))
    ax.add_patch(plt.Circle((0, 0), EARTH_RADIUS_KM, color="#9cc3e6", zorder=3))
    ax.text(0, 0, "Earth", ha="center", va="center", fontsize=9, zorder=4)
    handles = []
    for name, radius, color, period in orbits:
        ax.add_patch(plt.Circle((0, 0), radius, fill=False, color=color, linewidth=2, zorder=2))
        handles.append(plt.Line2D([], [], color=color, linewidth=3,
                                  label=f"{name}: {radius - EARTH_RADIUS_KM:,.0f} km up, {period} per orbit"))
    ax.set_xlim(-44000, 44000)
    ax.set_ylim(-44000, 44000)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("The orbits in this project, drawn to scale", fontsize=13)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.02), frameon=False, fontsize=10)
    fig.savefig(OUT / "orbit_scale.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    run(Ground_Track, ["--start", START, "--minutes", "270", "--no-cartopy", "--no-show",
                       "--save", str(OUT / "ground_tracks.png")])
    run(TLE, ["--start", START, "--minutes", "90", "--no-cartopy", "--no-show",
              "--save", str(OUT / "iss_90_minutes.png")])
    table = run(ISS, ["--lat", "30.0444", "--lon", "31.2357", "--start", START, "--hours", "48",
                      "--tz", "Africa/Cairo", "--no-show", "--save", str(OUT / "iss_passes.png")])
    (OUT / "iss_passes_example.txt").write_text(table)
    orbit_scale_picture()
    print(f"Saved pictures to {OUT}")

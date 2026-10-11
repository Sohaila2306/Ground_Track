"""TLE.py - look at one satellite's orbit data and where it goes in the next little while.

    python TLE.py                         # the ISS over the next 90 minutes
    python TLE.py --sat noaa15 --minutes 200
    python TLE.py --sat 25544 --no-plot   # just print the orbit summary

The "TLE" is the two-line element set that describes a satellite's orbit. It is
downloaded fresh from CelesTrak the first time and then reused for a few hours.
"""

from __future__ import annotations

import argparse
import sys

from Orbit_Tools import (COLORS, DEFAULT_GPS_PRN, MapCanvas, OrbitToolsError, SatelliteTracker,
                         add_common_arguments, cartopy_note, compute_ground_track, describe_period,
                         finish_figure, get_tle, parse_start, plot_ground_track, prepare_matplotlib,
                         stale_warning, summarize_orbit)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Show a satellite's orbit data and plot its ground track.")
    parser.add_argument("--sat", default="iss",
                        help="iss (default), noaa15, noaa18, noaa19, gps, or a NORAD catalogue number")
    parser.add_argument("--minutes", type=float, default=90, help="how far ahead to plot (default: 90)")
    parser.add_argument("--gps-prn", type=int, default=DEFAULT_GPS_PRN, help="which GPS satellite for --sat gps")
    parser.add_argument("--no-plot", action="store_true", help="print the orbit summary and stop")
    parser.add_argument("--no-cartopy", action="store_true", help="plain grid instead of a coastline map")
    add_common_arguments(parser)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.minutes <= 0:
        build_parser().error("--minutes must be positive")
    prepare_matplotlib(show=not (args.no_show or args.no_plot))

    try:
        start = parse_start(args.start)                  # check the input before touching the network
        tle = get_tle(args.sat, args.offline, args.gps_prn)
        summary = summarize_orbit(tle)

        print(f"{tle.name}  (NORAD {tle.catalog_number})")
        print(f"  orbit data : {tle.source}")
        print(f"  measured   : {tle.epoch:%Y-%m-%d %H:%M} UTC  ({abs((start - tle.epoch).days)} days from the start time)")
        print(f"  period     : {summary.period_minutes:.1f} min  ({summary.orbits_per_day:.2f} orbits per day)")
        print(f"  altitude   : {summary.perigee_altitude_km:,.0f} to {summary.apogee_altitude_km:,.0f} km")
        print(f"  inclination: {summary.inclination_deg:.2f} degrees")
        warning = stale_warning(tle, start)
        if warning:
            print("\n" + warning)
        if args.no_plot:
            return 0

        track = compute_ground_track(SatelliteTracker(tle), start, args.minutes)
        canvas = MapCanvas(f"{tle.name} ground track: {args.minutes:g} minutes from {start:%Y-%m-%d %H:%M} UTC "
                           "(dot = start, square = end)", use_cartopy=False if args.no_cartopy else None)
        cartopy_note(canvas)
        key = args.sat.lower()
        plot_ground_track(canvas, track, COLORS.get(key, "#d1495b"),
                          f"{tle.name} ({describe_period(summary.period_minutes)} orbit)",
                          mark_every_minutes=15 if args.minutes <= 240 else None)
        canvas.legend()
        finish_figure(canvas.fig, args.save, show=not args.no_show)
    except OrbitToolsError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

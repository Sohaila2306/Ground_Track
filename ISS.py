"""ISS.py - when can I see the ISS (or another satellite) from where I am?

    python ISS.py                                   # asks for your latitude and longitude
    python ISS.py --lat 30.0444 --lon 31.2357
    python ISS.py --lat 51.5 --lon -0.13 --tz Europe/London --hours 48
    python ISS.py --lat 40.7 --lon -74 --sat noaa15 --min-elevation 20

It prints every pass above the minimum elevation and draws two pictures: how high the
satellite is over time, and the path of the best pass across your sky.

"Above the horizon" is not the same as "visible to the eye". To actually see the ISS the
sky has to be dark where you are while the station is still in sunlight, which happens
in the hour or two after sunset or before sunrise.
"""

from __future__ import annotations

import argparse
import sys
from datetime import timedelta
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from Orbit_Tools import (COLORS, DEFAULT_GPS_PRN, UTC, OrbitToolsError, SatelliteTracker, add_common_arguments,
                         find_passes, finish_figure, format_pass_table, get_tle, latitude_arg, longitude_arg,
                         parse_start, plot_elevation, plot_sky_track, prepare_matplotlib, stale_warning,
                         time_grid)

STEP_SECONDS = 20


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="List and plot the passes of a satellite over your location.")
    parser.add_argument("--lat", type=latitude_arg, help="your latitude in degrees, north is positive")
    parser.add_argument("--lon", type=longitude_arg, help="your longitude in degrees, east is positive")
    parser.add_argument("--elevation", type=float, default=0.0, metavar="METRES",
                        help="your height above sea level (default: 0)")
    parser.add_argument("--sat", default="iss",
                        help="iss (default), noaa15, noaa18, noaa19, gps, or a NORAD catalogue number")
    parser.add_argument("--gps-prn", type=int, default=DEFAULT_GPS_PRN, help="which GPS satellite for --sat gps")
    parser.add_argument("--hours", type=float, default=24, help="how far ahead to look (default: 24)")
    parser.add_argument("--min-elevation", type=float, default=10, metavar="DEGREES",
                        help="count a pass only above this elevation (default: 10)")
    parser.add_argument("--tz", help="show times in this time zone, e.g. Africa/Cairo (default: UTC)")
    add_common_arguments(parser)
    return parser


def ask_number(prompt: str, low: float, high: float) -> float:
    while True:
        try:
            value = float(input(prompt).strip())
        except ValueError:
            print("  Please type a number, for example 30.0444")
            continue
        if low <= value <= high:
            return value
        print(f"  That needs to be between {low:g} and {high:g}.")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.hours <= 0:
        parser.error("--hours must be positive")
    prepare_matplotlib(show=not args.no_show)

    try:
        if args.lat is None or args.lon is None:
            print("Where are you? (north and east are positive, e.g. Cairo is 30.0444, 31.2357)")
        lat = args.lat if args.lat is not None else ask_number("Latitude  (degrees): ", -90, 90)
        lon = args.lon if args.lon is not None else ask_number("Longitude (degrees): ", -180, 180)
    except (KeyboardInterrupt, EOFError):
        print("\nCancelled.")
        return 1

    try:
        tz = ZoneInfo(args.tz) if args.tz else UTC
    except (ZoneInfoNotFoundError, ValueError):
        parser.error(f"unknown time zone {args.tz!r}; use a name like Africa/Cairo")

    try:
        start = parse_start(args.start)                  # check the input before touching the network
        tle = get_tle(args.sat, args.offline, args.gps_prn)
        tracker = SatelliteTracker(tle)
        moments = time_grid(start, args.hours * 60, STEP_SECONDS)
        elevation, azimuth, _ = tracker.look_angles(lat, lon, args.elevation, moments)
        passes = find_passes(moments, elevation, azimuth, args.min_elevation)

        zone_name = args.tz or "UTC"
        print(f"{tle.name} over {lat:.4f}, {lon:.4f}: passes above {args.min_elevation:g} degrees "
              f"in the next {args.hours:g} hours (times in {zone_name})")
        print(f"Orbit data: {tle.source}; measured {tle.epoch:%Y-%m-%d %H:%M} UTC\n")
        warning = stale_warning(tle, start)
        if warning:
            print(warning + "\n")
        if passes:
            print(format_pass_table(passes, tz))
            print("\nAbove the horizon is not the same as visible: you also need a dark sky and a sunlit satellite.")
        else:
            print("No passes in that window. Try a longer one with --hours, or a lower --min-elevation.")

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(13, 5))
        grid = fig.add_gridspec(1, 2, width_ratios=[2.1, 1])
        ax_time = fig.add_subplot(grid[0])
        color = COLORS.get(args.sat.lower(), COLORS["iss"])
        plot_elevation(ax_time, moments, elevation, args.min_elevation, passes, tz, color)
        ax_time.set_xlabel(f"Time ({zone_name})")
        ax_time.set_title(f"{tle.name}: elevation over {lat:.2f}, {lon:.2f}")

        ax_sky = fig.add_subplot(grid[1], projection="polar")
        if passes:
            best = max(passes, key=lambda p: p.peak_elevation)
            t0 = best.rise_time or moments[0]
            t1 = best.set_time or moments[-1]
            sky_times = time_grid(t0, (t1 - t0).total_seconds() / 60, 5)
            sky_elev, sky_az, _ = tracker.look_angles(lat, lon, args.elevation, sky_times)
            plot_sky_track(ax_sky, sky_az, sky_elev, color)
            ax_sky.set_title(f"Best pass: {best.peak_time.astimezone(tz):%a %H:%M}, peak {best.peak_elevation:.0f}\N{DEGREE SIGN}",
                             pad=18)
        else:
            ax_sky.set_title("No pass in this window", pad=18)
        fig.tight_layout()
        finish_figure(fig, args.save, show=not args.no_show)
    except OrbitToolsError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

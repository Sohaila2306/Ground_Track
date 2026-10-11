"""Ground_Track.py - where on Earth are the ISS, a NOAA weather satellite and a GPS satellite?

    python Ground_Track.py                           # the three satellites, next 4.5 hours
    python Ground_Track.py --minutes 720 --save tracks.png
    python Ground_Track.py --sats iss noaa19 noaa18  # pick your own
    python Ground_Track.py --gps-prn 24

The line on the map is the ground track: the point on Earth directly beneath the
satellite. A dot marks where it starts and a square where it ends.
"""

from __future__ import annotations

import argparse
import sys

from Orbit_Tools import (COLORS, DEFAULT_GPS_PRN, MapCanvas, MissingDependency, OrbitToolsError,
                         SatelliteTracker,
                         add_common_arguments, cartopy_note, compute_ground_track, describe_period,
                         finish_figure, get_tle, parse_start, plot_ground_track, prepare_matplotlib,
                         stale_warning, summarize_orbit)

DEFAULT_SATELLITES = ["iss", "noaa15", "gps"]
SPARE_COLORS = ["#8d6a9f", "#edae49", "#444444", "#e07b39"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot the ground tracks of several satellites on one world map.")
    parser.add_argument("--sats", nargs="+", default=DEFAULT_SATELLITES, metavar="SAT",
                        help="satellites to plot: iss, noaa15, noaa18, noaa19, gps, or NORAD catalogue "
                             f"numbers (default: {' '.join(DEFAULT_SATELLITES)})")
    parser.add_argument("--minutes", type=float, default=270, help="how far ahead to plot (default: 270)")
    parser.add_argument("--step", type=float, default=30, metavar="SECONDS", help="time between points (default: 30)")
    parser.add_argument("--gps-prn", type=int, default=DEFAULT_GPS_PRN, help="which GPS satellite for 'gps'")
    parser.add_argument("--no-cartopy", action="store_true", help="plain grid instead of a coastline map")
    add_common_arguments(parser)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.minutes <= 0 or args.step <= 0:
        parser.error("--minutes and --step must be positive")
    prepare_matplotlib(show=not args.no_show)

    try:
        start = parse_start(args.start)
    except OrbitToolsError as error:
        parser.error(str(error))

    tracks, skipped = [], []
    for key in args.sats:
        try:
            tle = get_tle(key, args.offline, args.gps_prn)
            track = compute_ground_track(SatelliteTracker(tle), start, args.minutes, args.step)
        except MissingDependency as error:               # no point trying the other satellites
            print(f"error: {error}", file=sys.stderr)
            return 1
        except OrbitToolsError as error:
            print(f"Skipping {key}: {error}", file=sys.stderr)
            skipped.append(key)
            continue
        warning = stale_warning(tle, start)
        if warning:
            print(warning, file=sys.stderr)
        tracks.append((key, tle, track))

    if not tracks:
        print("error: none of the satellites could be plotted.", file=sys.stderr)
        return 1

    canvas = MapCanvas(f"Ground tracks for {args.minutes / 60:g} hours from {start:%Y-%m-%d %H:%M} UTC "
                       "(dot = start, square = end)", use_cartopy=False if args.no_cartopy else None)
    cartopy_note(canvas)
    spare = iter(SPARE_COLORS)
    for key, tle, track in tracks:
        color = COLORS.get(key.lower()) or next(spare, "#444444")
        period = summarize_orbit(tle).period_minutes
        plot_ground_track(canvas, track, color, f"{tle.name} ({describe_period(period)} orbit)")
    canvas.legend()
    finish_figure(canvas.fig, args.save, show=not args.no_show)
    return 0


if __name__ == "__main__":
    sys.exit(main())

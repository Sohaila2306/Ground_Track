"""Orbit_Tools.py - the code shared by TLE.py, ISS.py and Ground_Track.py.

It has four jobs:

  1. Get fresh orbital data (TLEs) from CelesTrak, cache it, and check that it is valid.
  2. Turn a TLE into positions with Skyfield (ground track, or look angles from your location).
  3. Find the passes of a satellite over an observer.
  4. Draw the maps and charts.

Skyfield, Matplotlib and Cartopy are only imported when they are first needed, so
the TLE parsing and pass-finding logic can be used (and tested) without them.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import re
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np

UTC = timezone.utc

EARTH_RADIUS_KM = 6378.137
EARTH_MU_KM3_S2 = 398600.4418

CELESTRAK_URL = "https://celestrak.org/NORAD/elements/gp.php"
CACHE_DIR = Path(__file__).resolve().parent / ".tle_cache"
CACHE_MAX_AGE_HOURS = 6        # CelesTrak asks people not to download the same data more than every couple of hours
STALE_AFTER_DAYS = 14          # beyond this, a low-orbit prediction can be off by tens of kilometres or more

# NORAD catalogue numbers of the satellites this project knows by name.
SATELLITES = {
    "iss": 25544,
    "noaa15": 25338,
    "noaa18": 28654,
    "noaa19": 33591,
}
GPS_GROUP = "gps-ops"
DEFAULT_GPS_PRN = 10

COLORS = {"iss": "#d1495b", "noaa15": "#2e86ab", "noaa18": "#2e86ab", "noaa19": "#2e86ab", "gps": "#3c8d5a"}


class OrbitToolsError(Exception):
    """Anything the scripts should report politely instead of showing a traceback."""


class TLEError(OrbitToolsError):
    """The orbital data is malformed."""


class TLEUnavailable(OrbitToolsError):
    """No orbital data could be downloaded and none was cached."""


class MissingDependency(OrbitToolsError):
    """Skyfield (or another needed package) isn't installed."""


class PropagationError(OrbitToolsError):
    """Skyfield couldn't work out where the satellite is."""


# ---------------------------------------------------------------------------
# TLEs
# ---------------------------------------------------------------------------

def tle_checksum(line: str) -> int:
    """The check digit for a TLE line: digits count as themselves, '-' counts as 1."""
    return sum(int(ch) if ch.isdigit() else (1 if ch == "-" else 0) for ch in line[:68]) % 10


@dataclass(frozen=True)
class TLE:
    """One two-line element set: a satellite's orbit, as published by the US Space Force."""

    name: str
    line1: str
    line2: str
    source: str = ""            # where it came from, for messages

    @classmethod
    def from_lines(cls, name: str, line1: str, line2: str, source: str = "") -> "TLE":
        line1, line2 = line1.rstrip(), line2.rstrip()
        for number, line in ((1, line1), (2, line2)):
            if len(line) != 69:
                raise TLEError(f"TLE line {number} should be 69 characters long, not {len(line)}.")
            if not line.startswith(f"{number} "):
                raise TLEError(f"TLE line {number} should start with '{number} '.")
            if not line[68].isdigit() or tle_checksum(line) != int(line[68]):
                raise TLEError(f"TLE line {number} fails its checksum, so it has been corrupted or edited.")
        if line1[2:7] != line2[2:7]:
            raise TLEError("The two TLE lines belong to different satellites.")
        return cls(name.strip(), line1, line2, source)

    @property
    def catalog_number(self) -> int:
        return int(self.line1[2:7])

    @property
    def epoch(self) -> datetime:
        """The moment the orbit was measured. Predictions are best close to this time."""
        two_digit_year = int(self.line1[18:20])
        year = 2000 + two_digit_year if two_digit_year < 57 else 1900 + two_digit_year
        day_of_year = float(self.line1[20:32])
        return datetime(year, 1, 1, tzinfo=UTC) + timedelta(days=day_of_year - 1)

    @property
    def inclination_deg(self) -> float:
        return float(self.line2[8:16])

    @property
    def raan_deg(self) -> float:
        return float(self.line2[17:25])

    @property
    def eccentricity(self) -> float:
        return float("0." + self.line2[26:33].strip())

    @property
    def arg_perigee_deg(self) -> float:
        return float(self.line2[34:42])

    @property
    def mean_anomaly_deg(self) -> float:
        return float(self.line2[43:51])

    @property
    def mean_motion_rev_per_day(self) -> float:
        return float(self.line2[52:63])


def parse_tle_text(text: str, source: str = "") -> list[TLE]:
    """Read every TLE in a block of text (with or without a name line above each pair)."""
    lines = [line.rstrip() for line in text.splitlines() if line.strip()]
    found, pending_name, i = [], None, 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            name = pending_name or f"SAT {line[2:7].strip()}"
            found.append(TLE.from_lines(name, line, lines[i + 1], source))
            pending_name, i = None, i + 2
        else:
            pending_name = line[2:].strip() if line.startswith("0 ") else line.strip()
            i += 1
    if not found:
        raise TLEError("No orbital data found in the response. "
                       f"The server said: {text.strip()[:100]!r}")
    return found


@dataclass(frozen=True)
class OrbitSummary:
    period_minutes: float
    orbits_per_day: float
    semi_major_axis_km: float
    perigee_altitude_km: float
    apogee_altitude_km: float
    inclination_deg: float


def summarize_orbit(tle: TLE) -> OrbitSummary:
    """Period, size and shape of the orbit, straight from the TLE numbers (Kepler's third law)."""
    n = tle.mean_motion_rev_per_day
    period_s = 86400.0 / n
    a = (EARTH_MU_KM3_S2 * (period_s / (2 * math.pi)) ** 2) ** (1 / 3)
    e = tle.eccentricity
    return OrbitSummary(period_s / 60, n, a, a * (1 - e) - EARTH_RADIUS_KM,
                        a * (1 + e) - EARTH_RADIUS_KM, tle.inclination_deg)


def stale_warning(tle: TLE, when: datetime) -> Optional[str]:
    days = abs((when - tle.epoch).total_seconds()) / 86400
    if days <= STALE_AFTER_DAYS:
        return None
    return (f"Warning: the orbit data for {tle.name} was measured {days:.0f} days from the time you asked "
            f"about (epoch {tle.epoch:%Y-%m-%d}). Predictions drift as that gap grows, "
            "so treat the result as a rough picture.")


# --- getting TLEs ----------------------------------------------------------------

def _download(url: str, timeout: float = 15.0) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "satellite-orbit-visualizer/1.0"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def fetch_text(query: dict, cache_key: str, offline: bool = False,
               max_age_hours: float = CACHE_MAX_AGE_HOURS,
               download: Callable[[str], str] = _download,
               cache_dir: Optional[Path] = None, now: Optional[datetime] = None) -> tuple[str, str]:
    """TLE text for a CelesTrak query. Returns (text, description of where it came from).

    Uses a cached copy if it is recent. Otherwise downloads, and if the download fails
    falls back to an older cached copy (with a note) rather than giving up.
    """
    cache_dir = CACHE_DIR if cache_dir is None else cache_dir
    path = cache_dir / f"{cache_key}.tle"
    now = now or datetime.now(UTC)
    cached_age = None
    if path.exists():
        cached_age = (now - datetime.fromtimestamp(path.stat().st_mtime, UTC)).total_seconds() / 3600

    if cached_age is not None and (offline or cached_age < max_age_hours):
        return path.read_text(), f"cached copy, {cached_age:.1f} h old"
    if offline:
        raise TLEUnavailable("Offline mode, but there is no cached orbit data yet. "
                             "Run once without --offline while connected to the internet.")

    url = f"{CELESTRAK_URL}?{urllib.parse.urlencode({**query, 'FORMAT': 'tle'})}"
    try:
        text = download(url)
    except (urllib.error.URLError, OSError, ValueError) as error:
        reason = getattr(error, "reason", error)            # URLError wraps the useful part
        if cached_age is not None:
            return path.read_text(), f"cached copy, {cached_age:.1f} h old (couldn't reach CelesTrak: {reason})"
        raise TLEUnavailable(
            f"Couldn't download orbit data from CelesTrak ({reason}). "
            "Check your internet connection and try again.") from error

    parse_tle_text(text)                       # make sure it is real before caching it
    cache_dir.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return text, "downloaded from CelesTrak just now"


def get_tle(key: str = "iss", offline: bool = False, gps_prn: int = DEFAULT_GPS_PRN,
            fetch: Callable = fetch_text) -> TLE:
    """Current orbit data for a satellite.

    `key` is "iss", "noaa15", "noaa18", "noaa19", "gps", or any NORAD catalogue number.
    """
    key = key.strip().lower()
    if key == "gps":
        text, source = fetch({"GROUP": GPS_GROUP}, f"group_{GPS_GROUP}", offline)
        candidates = parse_tle_text(text, source)
        pattern = re.compile(rf"\(PRN\s*0*{gps_prn}\)", re.IGNORECASE)
        for tle in candidates:
            if pattern.search(tle.name):
                return tle
        # PRN numbers get reassigned over time; show something rather than nothing
        return replace(candidates[0], source=f"{source}; PRN {gps_prn} wasn't in the list, showing {candidates[0].name}")

    if key in SATELLITES:
        catalog = SATELLITES[key]
    elif key.isdigit():
        catalog = int(key)
    else:
        raise OrbitToolsError(f"Unknown satellite '{key}'. Use one of {', '.join([*SATELLITES, 'gps'])} "
                              "or a NORAD catalogue number.")
    text, source = fetch({"CATNR": str(catalog)}, f"catnr_{catalog}", offline)
    return parse_tle_text(text, source)[0]


# ---------------------------------------------------------------------------
# Times and geometry helpers
# ---------------------------------------------------------------------------

def parse_start(text: Optional[str]) -> datetime:
    """The start time: 'now' if blank, otherwise an ISO date/time (taken as UTC if it has no zone)."""
    if not text:
        return datetime.now(UTC).replace(microsecond=0)
    try:
        moment = datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    except ValueError:
        raise OrbitToolsError(f"Couldn't read the start time {text!r}. Use a form like 2026-10-10T18:30.") from None
    return moment.replace(tzinfo=UTC) if moment.tzinfo is None else moment.astimezone(UTC)


def time_grid(start: datetime, minutes: float, step_seconds: float) -> list[datetime]:
    """Evenly spaced times from `start` to `start + minutes`, both ends included."""
    count = int(round(minutes * 60 / step_seconds))
    return [start + timedelta(seconds=step_seconds * i) for i in range(count + 1)]


def wrap_longitude(lon) -> np.ndarray:
    """Longitudes folded into the range -180..180."""
    return (np.asarray(lon, dtype=float) + 180.0) % 360.0 - 180.0


def split_at_dateline(lons, lats) -> list[tuple[np.ndarray, np.ndarray]]:
    """Cut a ground track into pieces wherever it crosses the 180 degree meridian.

    Plotting a track that jumps from +179 to -179 draws a line straight across the map.
    This ends each piece exactly on the edge of the map (using the interpolated latitude of
    the crossing) and starts the next one on the opposite edge, so nothing is left dangling.
    Missing (NaN) points also end a piece.
    """
    lons = np.asarray(lons, dtype=float)
    lats = np.asarray(lats, dtype=float)
    pieces: list[tuple[list, list]] = []
    cur_lon: list = []
    cur_lat: list = []

    def close():
        nonlocal cur_lon, cur_lat
        if cur_lon:
            pieces.append((cur_lon, cur_lat))
        cur_lon, cur_lat = [], []

    for i in range(len(lons)):
        if not (np.isfinite(lons[i]) and np.isfinite(lats[i])):
            close()
            continue
        if cur_lon and abs(lons[i] - cur_lon[-1]) > 180:
            previous_lon, previous_lat = cur_lon[-1], cur_lat[-1]
            if previous_lon > 0:                      # heading east across +180
                edge_out, edge_in, shifted = 180.0, -180.0, lons[i] + 360.0
            else:                                     # heading west across -180
                edge_out, edge_in, shifted = -180.0, 180.0, lons[i] - 360.0
            fraction = (edge_out - previous_lon) / (shifted - previous_lon)
            cross_lat = previous_lat + fraction * (lats[i] - previous_lat)
            cur_lon.append(edge_out)
            cur_lat.append(cross_lat)
            close()
            cur_lon, cur_lat = [edge_in], [cross_lat]
        cur_lon.append(float(lons[i]))
        cur_lat.append(float(lats[i]))
    close()
    return [(np.array(a), np.array(b)) for a, b in pieces]


# ---------------------------------------------------------------------------
# Skyfield: positions and look angles
# ---------------------------------------------------------------------------

def _skyfield_api():
    try:
        from skyfield import api
    except ImportError:
        raise MissingDependency("Skyfield isn't installed. Run:  pip install -r requirements.txt") from None
    return api


@dataclass
class GroundTrack:
    name: str
    times: list
    latitude: np.ndarray
    longitude: np.ndarray
    altitude_km: np.ndarray


class SatelliteTracker:
    """A satellite you can ask 'where is it?' and 'where is it in my sky?'."""

    def __init__(self, tle: TLE):
        self.tle = tle
        self._api = _skyfield_api()
        self._ts = self._api.load.timescale()
        self._satellite = self._api.EarthSatellite(tle.line1, tle.line2, tle.name, self._ts)

    def _times(self, moments: Sequence[datetime]):
        return self._ts.from_datetimes(list(moments))

    def subpoints(self, moments: Sequence[datetime]):
        """Latitude, longitude (degrees) and altitude (km) of the point directly below the satellite."""
        wgs84 = self._api.wgs84
        try:
            geocentric = self._satellite.at(self._times(moments))
            if hasattr(wgs84, "subpoint_of"):
                sub = wgs84.subpoint_of(geocentric)
            elif hasattr(wgs84, "subpoint"):
                sub = wgs84.subpoint(geocentric)
            else:
                sub = geocentric.subpoint()
            lat = np.asarray(sub.latitude.degrees, dtype=float)
            lon = wrap_longitude(sub.longitude.degrees)
            alt = np.asarray(sub.elevation.km, dtype=float)
        except Exception as error:               # Skyfield/sgp4 raise several kinds of error
            raise PropagationError(f"Couldn't work out the orbit of {self.tle.name}: {error}") from error
        if not np.isfinite(lat).any():
            raise PropagationError(
                f"{self.tle.name} came out as 'not a number' at every time. "
                "The orbit data may be too old for these dates, or the satellite may have decayed.")
        return lat, lon, alt

    def look_angles(self, latitude: float, longitude: float, elevation_m: float,
                    moments: Sequence[datetime]):
        """Elevation (degrees above the horizon), azimuth (degrees clockwise from north) and range (km)."""
        try:
            observer = self._api.wgs84.latlon(latitude, longitude, elevation_m=elevation_m)
            altitude, azimuth, distance = (self._satellite - observer).at(self._times(moments)).altaz()
            return (np.asarray(altitude.degrees, dtype=float),
                    np.asarray(azimuth.degrees, dtype=float),
                    np.asarray(distance.km, dtype=float))
        except Exception as error:
            raise PropagationError(f"Couldn't work out the orbit of {self.tle.name}: {error}") from error


def compute_ground_track(tracker: SatelliteTracker, start: datetime, minutes: float,
                         step_seconds: float = 30.0) -> GroundTrack:
    moments = time_grid(start, minutes, step_seconds)
    lat, lon, alt = tracker.subpoints(moments)
    return GroundTrack(tracker.tle.name, moments, lat, lon, alt)


# ---------------------------------------------------------------------------
# Passes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Pass:
    """One trip of a satellite across the sky above the minimum elevation."""

    rise_time: Optional[datetime]         # None if it was already up when the window began
    peak_time: datetime
    peak_elevation: float
    set_time: Optional[datetime]          # None if it was still up when the window ended
    rise_azimuth: float
    peak_azimuth: float
    set_azimuth: float

    @property
    def duration(self) -> Optional[timedelta]:
        if self.rise_time is None or self.set_time is None:
            return None
        return self.set_time - self.rise_time


def _crossing(times, elevation, i, j, threshold):
    """When the elevation crosses `threshold` between samples i and j (straight-line estimate)."""
    e0, e1 = elevation[i], elevation[j]
    fraction = 0.0 if e1 == e0 else min(max((threshold - e0) / (e1 - e0), 0.0), 1.0)
    return times[i] + (times[j] - times[i]) * fraction


def find_passes(times: Sequence[datetime], elevation, azimuth, threshold: float = 10.0) -> list[Pass]:
    """Find every stretch where the elevation is at least `threshold` degrees.

    Rise and set times are interpolated between samples, so with samples every 20 seconds
    they are good to a few seconds. The peak is the highest sample.
    """
    elevation = np.asarray(elevation, dtype=float)
    azimuth = np.asarray(azimuth, dtype=float)
    above = elevation >= threshold                # NaN compares as False
    n, i, passes = len(times), 0, []
    while i < n:
        if not above[i]:
            i += 1
            continue
        first = i
        while i + 1 < n and above[i + 1]:
            i += 1
        last = i
        i += 1
        rise = None if first == 0 else _crossing(times, elevation, first - 1, first, threshold)
        fall = None if last == n - 1 else _crossing(times, elevation, last, last + 1, threshold)
        peak = first + int(np.argmax(elevation[first:last + 1]))
        passes.append(Pass(rise, times[peak], float(elevation[peak]), fall,
                           float(azimuth[first]), float(azimuth[peak]), float(azimuth[last])))
    return passes


_COMPASS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]


def compass_point(azimuth: float) -> str:
    return _COMPASS[int(((azimuth % 360) + 22.5) // 45) % 8]


def format_duration(delta: Optional[timedelta]) -> str:
    if delta is None:
        return "-"
    seconds = int(round(delta.total_seconds()))
    return f"{seconds // 60}m {seconds % 60:02d}s"


def format_pass_table(passes: Sequence[Pass], tz=UTC) -> str:
    """A plain-text table of passes, in the time zone `tz`."""
    def stamp(moment, with_date):
        return moment.astimezone(tz).strftime("%a %d %b %H:%M:%S" if with_date else "%H:%M:%S")

    rows = []
    for number, p in enumerate(passes, start=1):
        rise = "already up" if p.rise_time is None else f"{stamp(p.rise_time, True)}  {compass_point(p.rise_azimuth)}"
        peak = f"{stamp(p.peak_time, False)}  {p.peak_elevation:.0f}\N{DEGREE SIGN}"
        sets = "still up" if p.set_time is None else f"{stamp(p.set_time, False)}  {compass_point(p.set_azimuth)}"
        rows.append([str(number), rise, peak, sets, format_duration(p.duration)])
    headers = ["#", "Rises", "Highest", "Sets", "Length"]
    widths = [max(len(h), *(len(r[c]) for r in rows)) for c, h in enumerate(headers)]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    lines = [fmt.format(*headers), "  ".join("-" * w for w in widths)] + [fmt.format(*r) for r in rows]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def has_cartopy() -> bool:
    return importlib.util.find_spec("cartopy") is not None


def prepare_matplotlib(show: bool) -> None:
    """Call before any plotting. If no window will be shown, use the file-only backend,
    which also works on machines with no display."""
    if not show:
        import matplotlib
        matplotlib.use("Agg")


class MapCanvas:
    """A world map on which tracks can be drawn.

    With Cartopy installed you get land, sea and coastlines (Cartopy downloads the
    outlines the first time). Without it you still get a clean latitude/longitude grid.
    """

    def __init__(self, title: str, use_cartopy: Optional[bool] = None, figsize=(14, 7)):
        import matplotlib.pyplot as plt

        self.uses_cartopy = has_cartopy() if use_cartopy is None else use_cartopy
        self.fig = plt.figure(figsize=figsize)
        if self.uses_cartopy:
            import cartopy.crs as ccrs
            import cartopy.feature as cfeature

            self.transform = ccrs.PlateCarree()
            self.ax = self.fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
            self.ax.set_global()
            self.ax.add_feature(cfeature.OCEAN, facecolor="#dbe9f4")
            self.ax.add_feature(cfeature.LAND, facecolor="#f1ede4")
            self.ax.coastlines(linewidth=0.6, color="#555555")
            grid = self.ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
            for side in ("top_labels", "right_labels"):
                if hasattr(grid, side):
                    setattr(grid, side, False)
        else:
            self.transform = None
            self.ax = self.fig.add_subplot(1, 1, 1)
            self.ax.set_facecolor("#eef3f8")
            self.ax.set_xlim(-180, 180)
            self.ax.set_ylim(-90, 90)
            self.ax.set_aspect("equal")
            self.ax.set_xticks(range(-180, 181, 30))
            self.ax.set_yticks(range(-90, 91, 30))
            self.ax.set_xticklabels([f"{abs(v)}\N{DEGREE SIGN}{'E' if v > 0 else 'W' if v < 0 else ''}"
                                     for v in range(-180, 181, 30)], fontsize=8)
            self.ax.set_yticklabels([f"{abs(v)}\N{DEGREE SIGN}{'N' if v > 0 else 'S' if v < 0 else ''}"
                                     for v in range(-90, 91, 30)], fontsize=8)
            self.ax.grid(True, color="#b8c4d0", linewidth=0.6, linestyle="--")
            self.ax.axhline(0, color="#7a8794", linewidth=1.0)
        self.ax.set_title(title, fontsize=13)

    def _kw(self, kwargs):
        if self.transform is not None:
            kwargs["transform"] = self.transform
        return kwargs

    def line(self, lons, lats, **kwargs):
        return self.ax.plot(lons, lats, **self._kw(kwargs))

    def point(self, lon, lat, **kwargs):
        return self.ax.plot([lon], [lat], linestyle="none", **self._kw(kwargs))

    def label(self, lon, lat, text, **kwargs):
        return self.ax.text(lon, lat, text, **self._kw(kwargs))

    def legend(self):
        """Put the legend under the map so it never hides a track."""
        handles, _ = self.ax.get_legend_handles_labels()
        self.ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=max(len(handles), 1),
                       frameon=False, fontsize=10)


def plot_ground_track(canvas: MapCanvas, track: GroundTrack, color: str, label: str,
                      mark_every_minutes: Optional[float] = None) -> None:
    """Draw one track: the line (cut at the date line), a dot at the start and a square at the end."""
    pieces = split_at_dateline(track.longitude, track.latitude)
    for index, (lons, lats) in enumerate(pieces):
        canvas.line(lons, lats, color=color, linewidth=1.8, label=label if index == 0 else None)
    finite = np.isfinite(track.latitude) & np.isfinite(track.longitude)
    if finite.any():
        first, last = np.flatnonzero(finite)[[0, -1]]
        canvas.point(track.longitude[first], track.latitude[first], marker="o", markersize=8,
                     color=color, markeredgecolor="white", zorder=5)
        canvas.point(track.longitude[last], track.latitude[last], marker="s", markersize=7,
                     color=color, markeredgecolor="white", zorder=5)
    if mark_every_minutes:
        origin = track.times[0]
        for k, moment in enumerate(track.times):
            minutes = (moment - origin).total_seconds() / 60
            if k and finite[k] and abs(minutes / mark_every_minutes - round(minutes / mark_every_minutes)) < 1e-6:
                canvas.label(track.longitude[k] + 2, track.latitude[k] + 2, f"+{minutes:g} min",
                             fontsize=7, color=color)


def plot_elevation(ax, times, elevation, threshold: float, passes: Sequence[Pass], tz=UTC,
                   color: str = COLORS["iss"]) -> None:
    """Elevation of the satellite over time, with each pass above the threshold shaded."""
    import matplotlib.dates as mdates

    elevation = np.asarray(elevation, dtype=float)
    ax.plot(times, elevation, color=color, linewidth=1.5)
    ax.fill_between(times, threshold, elevation, where=elevation >= threshold, color=color, alpha=0.3)
    ax.axhline(threshold, color="#7a8794", linestyle="--", linewidth=1,
               label=f"{threshold:g}\N{DEGREE SIGN} minimum elevation")
    for p in passes:
        ax.annotate(f"{p.peak_elevation:.0f}\N{DEGREE SIGN}", (p.peak_time, p.peak_elevation),
                    textcoords="offset points", xytext=(0, 5), ha="center", fontsize=8)
    ax.set_ylim(0, 90)
    ax.set_ylabel("Elevation above horizon")
    ax.yaxis.set_major_formatter(lambda value, _pos: f"{value:g}\N{DEGREE SIGN}")
    long_window = (times[-1] - times[0]) > timedelta(hours=30)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M" if long_window else "%H:%M", tz=tz))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)


def plot_sky_track(ax, azimuth, elevation, color: str = COLORS["iss"]) -> None:
    """A pass as seen looking up: north at the top, east on the right, the zenith in the middle."""
    azimuth, elevation = np.asarray(azimuth, dtype=float), np.asarray(elevation, dtype=float)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.plot(np.radians(azimuth), elevation, color=color, linewidth=2)
    ax.plot(np.radians(azimuth[0]), elevation[0], "o", color=color, markeredgecolor="white", markersize=8)
    ax.plot(np.radians(azimuth[-1]), elevation[-1], "s", color=color, markeredgecolor="white", markersize=7)
    ax.set_ylim(90, 0)
    ax.set_yticks([0, 30, 60])
    ax.set_yticklabels(["horizon", "30\N{DEGREE SIGN}", "60\N{DEGREE SIGN}"], fontsize=7)
    ax.set_xticks(np.radians(np.arange(0, 360, 45)))
    ax.set_xticklabels(_COMPASS, fontsize=9)
    ax.grid(True, alpha=0.4)


def finish_figure(fig, save_path: Optional[str], show: bool) -> None:
    import matplotlib.pyplot as plt

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved {save_path}")
    if show:
        plt.show()
    else:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Command-line helpers shared by the three scripts
# ---------------------------------------------------------------------------

def latitude_arg(text: str) -> float:
    value = float(text)
    if not -90 <= value <= 90:
        raise argparse.ArgumentTypeError("latitude must be between -90 and 90")
    return value


def longitude_arg(text: str) -> float:
    value = float(text)
    if not -180 <= value <= 180:
        raise argparse.ArgumentTypeError("longitude must be between -180 and 180")
    return value


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """The options every script accepts."""
    parser.add_argument("--offline", action="store_true",
                        help="don't download anything; use orbit data cached by an earlier run")
    parser.add_argument("--start", metavar="ISO_TIME",
                        help="start time in UTC, e.g. 2026-10-10T18:30 (default: now)")
    parser.add_argument("--save", metavar="FILE", help="also save the picture to this file (png, pdf, svg...)")
    parser.add_argument("--no-show", action="store_true", help="don't open a window (useful with --save)")


def cartopy_note(canvas: MapCanvas) -> None:
    if not canvas.uses_cartopy:
        print("Note: drawing a plain latitude/longitude grid. Install Cartopy for coastlines:  pip install cartopy")


def describe_period(minutes: float) -> str:
    return f"{minutes:.0f} min" if minutes < 180 else f"{minutes / 60:.1f} h"

"""Tests for the satellite tools. Run with:  python -m pytest   (or:  python tests/test_orbit_tools.py)

Skyfield isn't needed: tests/fake_skyfield.py stands in for it, so these check all of OUR
logic (TLE handling, caching, date-line splitting, pass finding, the three scripts) but not
Skyfield itself. The expected values come from published data and basic orbital mechanics.
"""

import io
import math
import pathlib
import sys
import tempfile
import urllib.error
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timedelta, timezone
from unittest import mock
from zoneinfo import ZoneInfo

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]

import fake_skyfield  # noqa: E402

fake_skyfield.install()

import Ground_Track  # noqa: E402
import ISS  # noqa: E402
import Orbit_Tools as ot  # noqa: E402
import TLE  # noqa: E402
from fake_skyfield import build_tle  # noqa: E402

UTC = timezone.utc
EPOCH = datetime(2026, 1, 1, tzinfo=UTC)

# A real ISS TLE, widely published (epoch 2008-09-20), with valid check digits.
REAL_LINE1 = "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927"
REAL_LINE2 = "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"

# The three TLEs that shipped with the original scripts.
ORIGINAL_TLES = [
    ("1 25544U 98067A   24160.51780556  .00024885  00000+0  41749-3 0  9996",
     "2 25544  51.6404 259.2792 0003751 161.8425  37.3500 15.50019787452273"),
    ("1 28654U 05015A   24160.21513962 -.00000076  00000+0  00000+0 0  9997",
     "2 28654  55.1505  35.0374 0107909 324.8964  34.9501  2.00567658198362"),
    ("1 25338U 98030A   24160.05834304  .00000090  00000+0  77898-4 0  9994",
     "2 25338  98.7406 210.1389 0011530  75.7009 284.4304 14.25982377974135"),
]


def raises(exception, func, *args, **kwargs):
    try:
        func(*args, **kwargs)
    except exception:
        return True
    return False


def iss():
    return build_tle("ISS (ZARYA)", 25544, 51.6404, 259.2792, 0.0003751, 161.8425, 37.35, 15.50019787, EPOCH)


# --- TLEs -------------------------------------------------------------------------

def test_real_tle_parses_and_gives_the_published_epoch():
    tle = ot.TLE.from_lines("ISS (ZARYA)", REAL_LINE1, REAL_LINE2)
    assert tle.catalog_number == 25544
    assert abs((tle.epoch - datetime(2008, 9, 20, 12, 25, 40, tzinfo=UTC)).total_seconds()) < 1
    assert tle.inclination_deg == 51.6416
    assert abs(tle.eccentricity - 0.0006703) < 1e-12
    assert abs(tle.mean_motion_rev_per_day - 15.72125391) < 1e-9


def test_the_tles_in_the_original_scripts_were_corrupted():
    """Every one of them fails its checksum, and one has another satellite's catalogue number."""
    for line1, line2 in ORIGINAL_TLES:
        assert raises(ot.TLEError, ot.TLE.from_lines, "x", line1, line2)


def test_tle_validation_catches_damage():
    assert raises(ot.TLEError, ot.TLE.from_lines, "x", REAL_LINE1[:-1] + "0", REAL_LINE2)      # wrong check digit
    assert raises(ot.TLEError, ot.TLE.from_lines, "x", REAL_LINE1 + " ", REAL_LINE2[:-1])      # wrong length
    assert raises(ot.TLEError, ot.TLE.from_lines, "x", REAL_LINE2, REAL_LINE1)                 # lines swapped
    other = build_tle("other", 99999, 10, 10, 0.001, 10, 10, 14.0, EPOCH)
    assert raises(ot.TLEError, ot.TLE.from_lines, "x", REAL_LINE1, other.line2)                # different satellites


def test_epoch_two_digit_years_are_read_with_the_standard_pivot():
    """TLEs store the year as two digits: 57-99 mean 1957-1999 (Sputnik era), 00-56 mean 2000-2056."""
    for year in (1957, 1998, 1999, 2000, 2026, 2056):
        tle = build_tle("a", 1, 0, 0, 0, 0, 0, 1.0, datetime(year, 6, 15, 12, tzinfo=UTC))
        assert tle.epoch == datetime(year, 6, 15, 12, tzinfo=UTC), year


def test_parse_tle_text_variants():
    three = f"ISS (ZARYA)             \n{REAL_LINE1}\n{REAL_LINE2}\n"
    assert ot.parse_tle_text(three)[0].name == "ISS (ZARYA)"
    assert ot.parse_tle_text(f"{REAL_LINE1}\n{REAL_LINE2}")[0].name == "SAT 25544"
    assert ot.parse_tle_text(f"0 ISS (ZARYA)\n{REAL_LINE1}\n{REAL_LINE2}")[0].name == "ISS (ZARYA)"
    other = build_tle("NOAA 15", 25338, 98.7, 10, 0.001, 75, 284, 14.26, EPOCH)
    two_sets = three + f"{other.name}\n{other.line1}\n{other.line2}\n"
    assert [t.name for t in ot.parse_tle_text(two_sets)] == ["ISS (ZARYA)", "NOAA 15"]
    assert raises(ot.TLEError, ot.parse_tle_text, "No GP data found")
    assert raises(ot.TLEError, ot.parse_tle_text, "")


def test_orbit_summary_matches_keplers_third_law():
    s = ot.summarize_orbit(ot.TLE.from_lines("ISS", REAL_LINE1, REAL_LINE2))
    assert abs(s.period_minutes - 91.60) < 0.01
    assert abs(s.semi_major_axis_km - 6731) < 2
    assert s.perigee_altitude_km < s.apogee_altitude_km
    # a GPS-like orbit: about 12 hours, about 20,200 km up
    gps = ot.summarize_orbit(build_tle("g", 2, 55, 0, 0.0, 0, 0, 2.00567658, EPOCH))
    assert abs(gps.period_minutes / 60 - 11.97) < 0.01
    assert abs(gps.perigee_altitude_km - 20200) < 100


def test_stale_warning():
    tle = iss()
    assert ot.stale_warning(tle, EPOCH + timedelta(days=3)) is None
    assert ot.stale_warning(tle, EPOCH + timedelta(days=60)) is not None
    assert "60 days" in ot.stale_warning(tle, EPOCH + timedelta(days=60))


# --- downloading and caching -----------------------------------------------------------

class FakeNetwork:
    def __init__(self, text=None, fail=False):
        self.calls, self.text, self.fail = [], text, fail

    def __call__(self, url):
        self.calls.append(url)
        if self.fail:
            raise urllib.error.URLError("no route to host")
        return self.text


def test_fetch_downloads_then_reuses_the_cache():
    text = f"ISS (ZARYA)\n{REAL_LINE1}\n{REAL_LINE2}\n"
    with tempfile.TemporaryDirectory() as folder:
        net = FakeNetwork(text)
        first, how = ot.fetch_text({"CATNR": "25544"}, "x", download=net, cache_dir=pathlib.Path(folder))
        assert first == text and "downloaded" in how and len(net.calls) == 1
        assert "CATNR=25544" in net.calls[0] and "FORMAT=tle" in net.calls[0]
        again, how = ot.fetch_text({"CATNR": "25544"}, "x", download=net, cache_dir=pathlib.Path(folder))
        assert again == text and "cached" in how and len(net.calls) == 1        # no second download


def test_fetch_falls_back_to_an_old_cache_when_the_network_is_down():
    text = f"ISS (ZARYA)\n{REAL_LINE1}\n{REAL_LINE2}\n"
    with tempfile.TemporaryDirectory() as folder:
        folder = pathlib.Path(folder)
        ot.fetch_text({"CATNR": "1"}, "x", download=FakeNetwork(text), cache_dir=folder)
        later = datetime.now(UTC) + timedelta(hours=30)
        got, how = ot.fetch_text({"CATNR": "1"}, "x", download=FakeNetwork(fail=True), cache_dir=folder, now=later)
        assert got == text and "couldn't reach CelesTrak" in how


def test_fetch_errors_are_clear():
    with tempfile.TemporaryDirectory() as folder:
        folder = pathlib.Path(folder)
        assert raises(ot.TLEUnavailable, ot.fetch_text, {"CATNR": "1"}, "x",
                      download=FakeNetwork(fail=True), cache_dir=folder)
        assert raises(ot.TLEUnavailable, ot.fetch_text, {"CATNR": "1"}, "x", offline=True,
                      download=FakeNetwork("never used"), cache_dir=folder)
        # a junk response must not be cached
        assert raises(ot.TLEError, ot.fetch_text, {"CATNR": "1"}, "x",
                      download=FakeNetwork("No GP data found"), cache_dir=folder)
        assert not (folder / "x.tle").exists()


def test_offline_uses_the_cache_whatever_its_age():
    text = f"ISS (ZARYA)\n{REAL_LINE1}\n{REAL_LINE2}\n"
    with tempfile.TemporaryDirectory() as folder:
        folder = pathlib.Path(folder)
        ot.fetch_text({"CATNR": "1"}, "x", download=FakeNetwork(text), cache_dir=folder)
        net = FakeNetwork(fail=True)
        got, _ = ot.fetch_text({"CATNR": "1"}, "x", offline=True, download=net,
                               cache_dir=folder, now=datetime.now(UTC) + timedelta(days=30))
        assert got == text and net.calls == []


def test_get_tle_picks_the_right_satellite():
    iss_text = f"ISS (ZARYA)\n{REAL_LINE1}\n{REAL_LINE2}\n"
    seen = []

    def fetch(query, key, offline):
        seen.append((query, key))
        return iss_text, "test"

    assert ot.get_tle("ISS", fetch=fetch).catalog_number == 25544
    assert seen[-1][0] == {"CATNR": "25544"}
    ot.get_tle("noaa15", fetch=fetch)
    assert seen[-1][0] == {"CATNR": "25338"}
    ot.get_tle("33591", fetch=fetch)
    assert seen[-1][0] == {"CATNR": "33591"}
    assert raises(ot.OrbitToolsError, ot.get_tle, "krypton", fetch=fetch)


def test_get_tle_gps_selects_by_prn_and_falls_back():
    sats = [build_tle(f"GPS BIIF-{i} (PRN {prn:02d})", 40000 + i, 55, 0, 0.001, 0, 0, 2.0056, EPOCH)
            for i, prn in enumerate([8, 10, 24])]
    text = "".join(f"{s.name}\n{s.line1}\n{s.line2}\n" for s in sats)
    fetch = lambda query, key, offline: (text, "test")
    assert "PRN 10" in ot.get_tle("gps", gps_prn=10, fetch=fetch).name
    assert "PRN 24" in ot.get_tle("gps", gps_prn=24, fetch=fetch).name
    fallback = ot.get_tle("gps", gps_prn=99, fetch=fetch)
    assert "PRN 08" in fallback.name and "PRN 99 wasn't in the list" in fallback.source


# --- times and the date line ---------------------------------------------------------------

def test_parse_start():
    assert ot.parse_start("2026-10-10T18:30") == datetime(2026, 10, 10, 18, 30, tzinfo=UTC)
    assert ot.parse_start("2026-10-10T18:30Z") == datetime(2026, 10, 10, 18, 30, tzinfo=UTC)
    assert ot.parse_start("2026-10-10T20:30+02:00") == datetime(2026, 10, 10, 18, 30, tzinfo=UTC)
    assert abs((ot.parse_start(None) - datetime.now(UTC)).total_seconds()) < 5
    assert raises(ot.OrbitToolsError, ot.parse_start, "yesterday-ish")


def test_time_grid_includes_both_ends():
    grid = ot.time_grid(EPOCH, 10, 30)
    assert len(grid) == 21 and grid[0] == EPOCH and grid[-1] == EPOCH + timedelta(minutes=10)


def test_wrap_longitude():
    assert list(ot.wrap_longitude([0, 180, 190, -190, 360, 540.0])) == [0, -180, -170, 170, 0, -180]


def test_split_at_dateline_eastward_crossing_is_interpolated():
    # latitude climbs from 10 to 20 while longitude goes 170 -> 190 (= -170)
    pieces = ot.split_at_dateline([170, 175, -175, -170], [10, 12, 18, 20])
    assert len(pieces) == 2
    (lon_a, lat_a), (lon_b, lat_b) = pieces
    assert lon_a[-1] == 180 and lon_b[0] == -180
    assert math.isclose(lat_a[-1], 15) and math.isclose(lat_b[0], 15)          # halfway between 12 and 18
    assert list(lon_a[:-1]) == [170, 175] and list(lon_b[1:]) == [-175, -170]


def test_split_at_dateline_westward_crossing():
    pieces = ot.split_at_dateline([-170, -178, 178, 170], [0, 1, 3, 4])
    assert len(pieces) == 2 and pieces[0][0][-1] == -180 and pieces[1][0][0] == 180


def test_split_at_dateline_threshold_is_half_the_map():
    """A 200 degree jump is really a 160 degree move the other way round the globe, so it is a crossing;
    a 170 degree jump is a genuine 170 degree move and is not."""
    assert len(ot.split_at_dateline([-100, 100], [0, 10])) == 2
    assert len(ot.split_at_dateline([-85, 85], [0, 10])) == 1


def test_split_at_dateline_leaves_ordinary_tracks_alone_and_breaks_on_nan():
    assert len(ot.split_at_dateline([0, 10, 20], [0, 1, 2])) == 1
    assert len(ot.split_at_dateline([0, 10, np.nan, 30, 40], [0, 1, 2, 3, 4])) == 2
    assert ot.split_at_dateline([], []) == []


# --- passes --------------------------------------------------------------------------------

def synthetic(elevation_function, seconds=3000, step=10):
    moments = [EPOCH + timedelta(seconds=s) for s in range(0, seconds + 1, step)]
    t = np.arange(0, seconds + 1, step, dtype=float)
    return moments, elevation_function(t), np.linspace(0, 180, len(moments))


def test_find_passes_interpolates_the_crossing_exactly():
    moments, elevation, azimuth = synthetic(lambda t: t / 10 - 5, seconds=600)      # a straight ramp
    (p,) = ot.find_passes(moments, elevation, azimuth, threshold=10)
    assert p.rise_time == EPOCH + timedelta(seconds=150)       # elevation is 10 exactly at t = 150 s
    assert p.set_time is None                                   # still rising when the window ends


def test_find_passes_two_humps():
    # two identical humps (peaks at 375 s and 1125 s), each rising to 50 degrees from a floor of -10
    moments, elevation, azimuth = synthetic(lambda t: 60 * np.sin(np.pi * t / 750) ** 2 - 10, seconds=1500)
    passes = ot.find_passes(moments, elevation, azimuth, threshold=10)
    assert len(passes) == 2
    for p in passes:
        assert p.rise_time < p.peak_time < p.set_time
        assert abs(p.peak_elevation - 50) < 0.5
        assert p.duration.total_seconds() > 0


def test_find_passes_in_progress_at_both_ends_and_none():
    moments, elevation, azimuth = synthetic(lambda t: np.full_like(t, 30.0))
    (p,) = ot.find_passes(moments, elevation, azimuth, 10)
    assert p.rise_time is None and p.set_time is None and p.duration is None
    assert ot.find_passes(moments, elevation - 100, azimuth, 10) == []
    assert ot.find_passes(moments, np.full_like(elevation, np.nan), azimuth, 10) == []


def test_small_formatting_helpers():
    assert [ot.compass_point(a) for a in (0, 44, 90, 181, 270, 359.9)] == ["N", "NE", "E", "S", "W", "N"]
    assert ot.format_duration(timedelta(seconds=392)) == "6m 32s" and ot.format_duration(None) == "-"
    moments, elevation, azimuth = synthetic(lambda t: 50 * np.sin(np.pi * t / 600) ** 2 - 5, seconds=600)
    table = ot.format_pass_table(ot.find_passes(moments, elevation, azimuth, 10), ZoneInfo("Africa/Cairo"))
    assert "Rises" in table and "02:" in table                  # 00:xx UTC is 02:xx in Cairo in January


# --- tracking with the Skyfield stand-in ----------------------------------------------------------

def test_ground_track_obeys_basic_orbital_mechanics():
    tle = iss()
    track = ot.compute_ground_track(ot.SatelliteTracker(tle), EPOCH, 270)
    assert abs(track.latitude.max() - 51.6404) < 0.1 and abs(track.latitude.min() + 51.6404) < 0.1
    summary = ot.summarize_orbit(tle)
    assert abs(track.altitude_km.mean() - summary.perigee_altitude_km) < 10
    assert (np.abs(track.longitude) <= 180).all()


def test_the_track_shifts_west_each_orbit_because_the_earth_turns():
    tle = iss()
    period = ot.summarize_orbit(tle).period_minutes
    track = ot.compute_ground_track(ot.SatelliteTracker(tle), EPOCH, period, step_seconds=period * 60)
    lat0, lat1 = track.latitude
    shift = (track.longitude[1] - track.longitude[0] + 540) % 360 - 180
    assert abs(lat1 - lat0) < 0.5
    assert abs(shift + 360 * period * 60 / 86164.1) < 0.5        # about 23 degrees west


def test_dateline_pieces_from_a_real_track_are_continuous():
    track = ot.compute_ground_track(ot.SatelliteTracker(iss()), EPOCH, 270)
    pieces = ot.split_at_dateline(track.longitude, track.latitude)
    assert len(pieces) >= 2
    for lons, lats in pieces:
        assert np.abs(np.diff(lons)).max() < 5 and np.abs(np.diff(lats)).max() < 5


def test_overhead_satellite_is_at_90_degrees_and_far_side_is_below_the_horizon():
    tracker = ot.SatelliteTracker(iss())
    when = EPOCH + timedelta(minutes=37)
    lat, lon, _ = tracker.subpoints([when])
    elevation, _, distance = tracker.look_angles(lat[0], lon[0], 0, [when])
    assert elevation[0] > 89.5 and 400 < distance[0] < 440
    elevation, _, _ = tracker.look_angles(-lat[0], lon[0] + 180, 0, [when])
    assert elevation[0] < -20


def test_a_day_of_passes_looks_sensible():
    tracker = ot.SatelliteTracker(iss())
    moments = ot.time_grid(EPOCH, 3 * 24 * 60, 20)
    elevation, azimuth, _ = tracker.look_angles(30.0444, 31.2357, 0, moments)
    passes = ot.find_passes(moments, elevation, azimuth, 10)
    assert len(passes) >= 4
    for p in passes:
        assert p.peak_elevation >= 10
        if p.duration is not None:
            assert timedelta(minutes=1) < p.duration < timedelta(minutes=12)
    starts = [p.peak_time for p in passes]
    assert starts == sorted(starts)
    assert all(b - a > timedelta(minutes=60) for a, b in zip(starts, starts[1:]))


def test_missing_skyfield_is_explained():
    saved = {k: sys.modules.pop(k) for k in ("skyfield", "skyfield.api")}
    try:
        with mock.patch.dict(sys.modules, {"skyfield": None}):
            assert raises(ot.MissingDependency, ot.SatelliteTracker, iss())
    finally:
        sys.modules.update(saved)


# --- the three scripts ------------------------------------------------------------------------------

def fixtures_by_key():
    return {
        "iss": iss(),
        "noaa15": build_tle("NOAA 15", 25338, 98.7406, 210.1389, 0.001153, 75.7, 284.43, 14.25982377, EPOCH),
        "gps": build_tle("GPS BIIF-3 (PRN 10)", 40105, 55.1505, 35.0374, 0.0107, 324.9, 34.95, 2.00567658, EPOCH),
    }


def run(module, argv, typed=()):
    out, err = io.StringIO(), io.StringIO()
    fixtures = fixtures_by_key()
    with redirect_stdout(out), redirect_stderr(err), \
            mock.patch.object(module, "get_tle", side_effect=lambda key, *a, **k: fixtures[key.lower()]), \
            mock.patch("builtins.input", side_effect=list(typed)):
        try:
            code = module.main(argv)
        except SystemExit as exit_:
            code = exit_.code
    return code, out.getvalue(), err.getvalue()


def is_png(path):
    return pathlib.Path(path).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_ground_track_script_draws_a_png():
    with tempfile.TemporaryDirectory() as folder:
        picture = f"{folder}/tracks.png"
        code, out, err = run(Ground_Track, ["--start", "2026-01-01T00:00", "--no-show", "--save", picture])
        assert code == 0 and is_png(picture)
        if not ot.has_cartopy():
            assert "pip install cartopy" in out                        # the plain-grid fallback is announced


def test_ground_track_script_skips_a_missing_satellite_and_carries_on():
    fixtures = fixtures_by_key()

    def lookup(key, *args, **kwargs):
        if key == "noaa99":
            raise ot.TLEUnavailable("no such satellite")
        return fixtures[key]

    out, err = io.StringIO(), io.StringIO()
    with tempfile.TemporaryDirectory() as folder, redirect_stdout(out), redirect_stderr(err), \
            mock.patch.object(Ground_Track, "get_tle", side_effect=lookup):
        picture = f"{folder}/tracks.png"
        code = Ground_Track.main(["--sats", "iss", "noaa99", "--no-show", "--no-cartopy", "--save", picture,
                                  "--start", "2026-01-01T00:00"])
        assert code == 0 and is_png(picture)
    assert "Skipping noaa99: no such satellite" in err.getvalue()


def test_ground_track_script_fails_once_when_skyfield_is_missing():
    saved = {k: sys.modules.pop(k) for k in ("skyfield", "skyfield.api")}
    err = io.StringIO()
    try:
        with redirect_stdout(io.StringIO()), redirect_stderr(err), mock.patch.dict(sys.modules, {"skyfield": None}), \
                mock.patch.object(Ground_Track, "get_tle", side_effect=lambda key, *a, **k: fixtures_by_key()[key]):
            code = Ground_Track.main(["--no-show"])
    finally:
        sys.modules.update(saved)
    assert code == 1 and err.getvalue().count("Skyfield isn't installed") == 1
    assert "Skipping" not in err.getvalue()


def test_network_errors_show_the_reason_not_python_internals():
    with tempfile.TemporaryDirectory() as folder:
        try:
            ot.fetch_text({"CATNR": "1"}, "x", download=FakeNetwork(fail=True), cache_dir=pathlib.Path(folder))
        except ot.TLEUnavailable as error:
            assert "no route to host" in str(error) and "urlopen" not in str(error)
        else:
            raise AssertionError("should have failed")


def test_ground_track_script_reports_when_nothing_works():
    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err), \
            mock.patch.object(Ground_Track, "get_tle", side_effect=ot.TLEUnavailable("no network")):
        code = Ground_Track.main(["--no-show"])
    assert code == 1 and "none of the satellites could be plotted" in err.getvalue() and "no network" in err.getvalue()


def test_iss_script_lists_passes_and_saves_a_picture():
    with tempfile.TemporaryDirectory() as folder:
        picture = f"{folder}/passes.png"
        code, out, _ = run(ISS, ["--lat", "30.0444", "--lon", "31.2357", "--hours", "48", "--tz", "Africa/Cairo",
                                 "--start", "2026-01-01T00:00", "--no-show", "--save", picture])
        assert code == 0 and is_png(picture)
        assert "Rises" in out and "Highest" in out and "Africa/Cairo" in out and "not the same as visible" in out


def test_iss_script_with_no_passes_says_so():
    code, out, _ = run(ISS, ["--lat", "30", "--lon", "31", "--hours", "0.2", "--min-elevation", "89",
                             "--start", "2026-01-01T00:00", "--no-show"])
    assert code == 0 and "No passes in that window" in out


def test_iss_script_asks_for_missing_location_and_retries_on_typos():
    code, out, _ = run(ISS, ["--hours", "6", "--start", "2026-01-01T00:00", "--no-show"],
                       typed=["north", "95", "30.04", "999", "31.23"])
    assert code == 0 and "Please type a number" in out and "needs to be between" in out


def test_iss_script_rejects_bad_arguments():
    assert run(ISS, ["--lat", "95", "--lon", "0"])[0] == 2
    assert run(ISS, ["--lat", "1", "--lon", "1", "--tz", "Mars/Olympus"])[0] == 2
    assert run(ISS, ["--lat", "1", "--lon", "1", "--hours", "-3"])[0] == 2


def test_iss_script_survives_ctrl_c_at_the_prompt():
    code, out, _ = run(ISS, ["--no-show"], typed=[KeyboardInterrupt()])
    assert code == 1 and "Cancelled" in out


def test_tle_script_prints_a_summary_and_plots():
    with tempfile.TemporaryDirectory() as folder:
        picture = f"{folder}/iss.png"
        code, out, _ = run(TLE, ["--start", "2026-01-01T00:00", "--no-show", "--save", picture])
        assert code == 0 and is_png(picture)
        assert "ISS (ZARYA)" in out and "92.9 min" in out and "inclination: 51.64" in out
    code, out, _ = run(TLE, ["--no-plot", "--sat", "noaa15", "--start", "2026-01-01T00:00"])
    assert code == 0 and "NOAA 15" in out and "98.74" in out


def test_tle_script_warns_about_old_orbit_data():
    code, out, _ = run(TLE, ["--no-plot", "--start", "2026-06-01T00:00"])
    assert code == 0 and "Warning" in out and "151 days" in out


def test_a_bad_start_time_is_reported_before_anything_is_downloaded():
    def must_not_be_called(*args, **kwargs):
        raise AssertionError("tried to download before checking the start time")

    for module, argv in ((TLE, ["--start", "someday", "--no-plot"]),
                         (ISS, ["--lat", "1", "--lon", "1", "--start", "someday", "--no-show"])):
        err = io.StringIO()
        with redirect_stdout(io.StringIO()), redirect_stderr(err), \
                mock.patch.object(module, "get_tle", side_effect=must_not_be_called):
            assert module.main(argv) == 1
        assert "Couldn't read the start time" in err.getvalue()


def test_scripts_turn_library_errors_into_messages_not_tracebacks():
    err = io.StringIO()
    with redirect_stdout(io.StringIO()), redirect_stderr(err), \
            mock.patch.object(TLE, "get_tle", side_effect=ot.TLEUnavailable("offline and nothing cached")):
        assert TLE.main(["--no-show"]) == 1
    assert "error: offline and nothing cached" in err.getvalue()


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_")]
    failed = 0
    for name, func in tests:
        try:
            func()
            print(f"ok    {name}")
        except Exception as error:                                   # noqa: BLE001
            failed += 1
            print(f"FAIL  {name}: {type(error).__name__}: {error}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)

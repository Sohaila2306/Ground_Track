"""A stand-in for the small part of Skyfield that this project uses.

It exists for two reasons: the tests run without Skyfield installed or any internet, and
docs/make_readme_figures.py can draw the README pictures offline.

It is NOT accurate. Orbits are perfect circles (with Earth's rotation and the nodal drift
caused by Earth's bulge), the Earth is a sphere, and nothing here knows about drag or the
real SGP4 model. Real runs use Skyfield. Here we only need something that behaves like a
satellite so the surrounding code can be checked.

Use:
    import fake_skyfield
    fake_skyfield.install()      # now `from skyfield import api` finds this module
"""

import math
import sys
import types
from datetime import datetime, timedelta, timezone

import numpy as np

MU = 398600.4418            # km^3 / s^2
EARTH_RADIUS = 6378.137     # km
J2 = 1.08263e-3


# ---------------------------------------------------------------------------
# Building valid TLEs from orbital elements
# ---------------------------------------------------------------------------

def _checksum(line):
    return sum(int(c) if c.isdigit() else (1 if c == "-" else 0) for c in line) % 10


def build_tle(name, catalog, inclination, raan, eccentricity, arg_perigee, mean_anomaly,
              mean_motion, epoch):
    """A TLE with correct formatting and check digits, made from orbital elements.

    Returns an Orbit_Tools.TLE. The numbers are whatever you give it, so it is only as
    real as your inputs.
    """
    from Orbit_Tools import TLE

    day_of_year = (epoch - datetime(epoch.year, 1, 1, tzinfo=timezone.utc)).total_seconds() / 86400 + 1
    line1 = (f"1 {catalog:05d}U 00000A   {epoch.year % 100:02d}{day_of_year:012.8f} "
             f" .00000000  00000-0  00000-0 0  999")
    line2 = (f"2 {catalog:05d} {inclination:8.4f} {raan:8.4f} {int(round(eccentricity * 1e7)):07d} "
             f"{arg_perigee:8.4f} {mean_anomaly:8.4f} {mean_motion:11.8f}00001")
    line1 += str(_checksum(line1))
    line2 += str(_checksum(line2))
    assert len(line1) == 69 and len(line2) == 69, (len(line1), len(line2))
    return TLE.from_lines(name, line1, line2, source="built from orbital elements")


# ---------------------------------------------------------------------------
# The stand-in Skyfield API
# ---------------------------------------------------------------------------

class _Angle:
    def __init__(self, degrees):
        self.degrees = degrees


class _Distance:
    def __init__(self, km):
        self.km = km


class _Time:
    def __init__(self, moments):
        self.moments = list(moments)


class _Timescale:
    def from_datetimes(self, moments):
        moments = list(moments)
        if any(m.tzinfo is None for m in moments):
            raise ValueError("naive datetime: Skyfield wants time-zone-aware datetimes")
        return _Time(moments)


class _Load:
    @staticmethod
    def timescale():
        return _Timescale()


def _gmst_radians(moments):
    jd = np.array([2440587.5 + m.timestamp() / 86400 for m in moments])
    return np.radians((280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360)


class _Subpoint:
    def __init__(self, lat, lon, altitude):
        self.latitude = _Angle(lat)
        self.longitude = _Angle(lon)
        self.elevation = _Distance(altitude)


class _Geocentric:
    def __init__(self, ecef):
        self.ecef = ecef


class _Wgs84:
    @staticmethod
    def subpoint_of(geocentric):
        x, y, z = geocentric.ecef.T
        r = np.sqrt(x * x + y * y + z * z)
        return _Subpoint(np.degrees(np.arcsin(z / r)), np.degrees(np.arctan2(y, x)), r - EARTH_RADIUS)

    @staticmethod
    def latlon(latitude_degrees, longitude_degrees, elevation_m=0.0):
        return _Observer(latitude_degrees, longitude_degrees, elevation_m)


class _Observer:
    def __init__(self, lat, lon, elevation_m):
        self.lat, self.lon = math.radians(lat), math.radians(lon)
        r = EARTH_RADIUS + elevation_m / 1000
        self.position = r * np.array([math.cos(self.lat) * math.cos(self.lon),
                                      math.cos(self.lat) * math.sin(self.lon),
                                      math.sin(self.lat)])


class _Topocentric:
    def __init__(self, ecef, observer):
        self.ecef, self.observer = ecef, observer

    def altaz(self):
        o = self.observer
        d = self.ecef - o.position
        east = np.array([-math.sin(o.lon), math.cos(o.lon), 0.0])
        north = np.array([-math.sin(o.lat) * math.cos(o.lon), -math.sin(o.lat) * math.sin(o.lon), math.cos(o.lat)])
        up = np.array([math.cos(o.lat) * math.cos(o.lon), math.cos(o.lat) * math.sin(o.lon), math.sin(o.lat)])
        distance = np.linalg.norm(d, axis=1)
        altitude = np.degrees(np.arcsin(np.clip(d @ up / distance, -1, 1)))
        azimuth = np.degrees(np.arctan2(d @ east, d @ north)) % 360
        return _Angle(altitude), _Angle(azimuth), _Distance(distance)


class _Difference:
    def __init__(self, satellite, observer):
        self.satellite, self.observer = satellite, observer

    def at(self, t):
        return _Topocentric(self.satellite._ecef(t), self.observer)


class EarthSatellite:
    """A satellite on a circular orbit built from the TLE's own numbers."""

    def __init__(self, line1, line2, name="", ts=None):
        yy, day = int(line1[18:20]), float(line1[20:32])
        self.epoch = datetime(2000 + yy if yy < 57 else 1900 + yy, 1, 1, tzinfo=timezone.utc) \
            + timedelta(days=day - 1)
        self.name = name
        self.inclination = math.radians(float(line2[8:16]))
        self.raan = math.radians(float(line2[17:25]))
        self.arg_perigee = math.radians(float(line2[34:42]))
        self.mean_anomaly = math.radians(float(line2[43:51]))
        self.n = float(line2[52:63]) * 2 * math.pi / 86400          # rad/s
        self.radius = (MU / self.n**2) ** (1 / 3)
        # nodal drift from Earth's bulge: this is what makes sun-synchronous orbits work
        self.raan_rate = -1.5 * J2 * (EARTH_RADIUS / self.radius) ** 2 * self.n * math.cos(self.inclination)

    def __sub__(self, observer):
        return _Difference(self, observer)

    def _ecef(self, t):
        seconds = np.array([(m - self.epoch).total_seconds() for m in t.moments])
        u = self.arg_perigee + self.mean_anomaly + self.n * seconds
        raan = self.raan + self.raan_rate * seconds
        ci, si = math.cos(self.inclination), math.sin(self.inclination)
        x = self.radius * (np.cos(raan) * np.cos(u) - np.sin(raan) * np.sin(u) * ci)
        y = self.radius * (np.sin(raan) * np.cos(u) + np.cos(raan) * np.sin(u) * ci)
        z = self.radius * np.sin(u) * si
        theta = _gmst_radians(t.moments)
        c, s = np.cos(theta), np.sin(theta)
        return np.column_stack([c * x + s * y, -s * x + c * y, z])

    def at(self, t):
        return _Geocentric(self._ecef(t))


def install():
    """Make `from skyfield import api` find this stand-in."""
    api = types.ModuleType("skyfield.api")
    api.EarthSatellite = EarthSatellite
    api.load = _Load()
    api.wgs84 = _Wgs84()
    package = types.ModuleType("skyfield")
    package.api = api
    sys.modules["skyfield"] = package
    sys.modules["skyfield.api"] = api


def uninstall():
    sys.modules.pop("skyfield", None)
    sys.modules.pop("skyfield.api", None)

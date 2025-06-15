import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from skyfield.api import load, EarthSatellite, utc
from datetime import datetime, timedelta

def split_longitude_wraps(lons, lats):
    segments = []
    current_segment_lons = [lons[0]]
    current_segment_lats = [lats[0]]

    for i in range(1, len(lons)):
        if np.isnan(lons[i]) or np.isnan(lats[i]):
            if current_segment_lons and current_segment_lats:
                segments.append((current_segment_lons, current_segment_lats))
            current_segment_lons = []
            current_segment_lats = []
            continue

        if current_segment_lons and abs(lons[i] - lons[i - 1]) > 180:
            segments.append((current_segment_lons, current_segment_lats))
            current_segment_lons = [lons[i]]
            current_segment_lats = [lats[i]]
        else:
            current_segment_lons.append(lons[i])
            current_segment_lats.append(lats[i])

    if current_segment_lons and current_segment_lats:
        segments.append((current_segment_lons, current_segment_lats))

    return segments

def get_satellite_ground_track(name, line1, line2, duration_minutes=270, step_minutes=1):
    satellite = EarthSatellite(line1, line2, name)
    ts = load.timescale()
    start_time = datetime.now(utc)
    times = [ts.utc(start_time + timedelta(minutes=i)) for i in range(0, duration_minutes + 1, step_minutes)]

    latitudes = []
    longitudes = []

    for t in times:
        subpoint = satellite.at(t).subpoint()
        lat = subpoint.latitude.degrees
        lon = subpoint.longitude.degrees
        if lon > 180:
            lon -= 360
        if np.isfinite(lat) and np.isfinite(lon):
            latitudes.append(lat)
            longitudes.append(lon)
        else:
            latitudes.append(np.nan)
            longitudes.append(np.nan)
    return latitudes, longitudes

def main():
    satellites = [
        (
            "ISS (ZARYA)",
            "1 25544U 98067A   24160.51780556  .00024885  00000+0  41749-3 0  9996",
            "2 25544  51.6404 259.2792 0003751 161.8425  37.3500 15.50019787452273"
        ),
        (
            "GPS BIIRM-5 (PRN 10)",
            "1 28654U 05015A   24160.21513962 -.00000076  00000+0  00000+0 0  9997",
            "2 28654  55.1505  35.0374 0107909 324.8964  34.9501  2.00567658198362"
        ),
        (
            "NOAA 15",
            "1 25338U 98030A   24160.05834304  .00000090  00000+0  77898-4 0  9994",
            "2 25338  98.7406 210.1389 0011530  75.7009 284.4304 14.25982377974135"
        ),
    ]

    plt.figure(figsize=(16, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.gridlines(draw_labels=True)

    colors = ['green', 'blue', 'red']

    for i, (name, l1, l2) in enumerate(satellites):
        lat, lon = get_satellite_ground_track(name, l1, l2, duration_minutes=270, step_minutes=1)
        segments = split_longitude_wraps(lon, lat)
        for seg_lons, seg_lats in segments:
            ax.plot(seg_lons, seg_lats, color=colors[i], label=name if seg_lons == segments[0][0] else "", transform=ccrs.PlateCarree())

    plt.title("Ground Tracks for ISS, GPS, and NOAA Satellites (3 Orbits / 270 minutes)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

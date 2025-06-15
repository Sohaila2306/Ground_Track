from skyfield.api import load, EarthSatellite, utc
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Load TLE data (ISS example)
name = "ISS (ZARYA)"
line1 = "1 25544U 98067A   24160.51780556  .00024885  00000+0  41749-3 0  9996"
line2 = "2 25544  51.6404 259.2792 0003751 161.8425  37.3500 15.50019787452273"

satellite = EarthSatellite(line1, line2, name)
ts = load.timescale()

# Make the datetime object UTC-aware
start_time = datetime.now(utc)
times = [ts.utc(start_time + timedelta(minutes=i)) for i in range(0, 90, 2)]

# Get satellite ground track
latitudes = []
longitudes = []

for t in times:
    subpoint = satellite.at(t).subpoint()
    latitudes.append(subpoint.latitude.degrees)
    longitudes.append(subpoint.longitude.degrees)

# Plot
plt.figure(figsize=(10, 5))
plt.plot(longitudes, latitudes, marker='o', linestyle='-', color='blue')
plt.title(f"{name} Ground Track (Next 90 minutes)")
plt.xlabel("Longitude [°]")
plt.ylabel("Latitude [°]")
plt.xlim(-180, 180)
plt.ylim(-90, 90)
plt.grid(True)
plt.show()

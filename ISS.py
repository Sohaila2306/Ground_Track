from skyfield.api import load, EarthSatellite, Topos, utc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# Try to get satellite data from an external file called TLE.py
# If not found, use default data for the ISS (International Space Station)
try:
    from TLE import name, line1, line2
except ImportError:
    name = "ISS (ZARYA)"
    line1 = "1 25544U 98067A   24160.51780556  .00024885  00000+0  41749-3 0  9996"
    line2 = "2 25544  51.6404 259.2792 0003751 161.8425  37.3500 15.50019787452273"

# Ask the user for their current location in degrees
user_lat = float(input("Please enter your latitude (degrees): "))
user_lon = float(input("Please enter your longitude (degrees): "))

# Create a satellite object using the TLE data
satellite = EarthSatellite(line1, line2, name)

# Load a timescale object to work with times in UTC
ts = load.timescale()

# Define where the observer (user) is located on Earth
observer = Topos(latitude_degrees=user_lat, longitude_degrees=user_lon)

# Set the starting time for the observation as the current UTC time
start_time = datetime.now(utc)

# Generate a list of times for the next 90 minutes, one per minute
times = [ts.utc(start_time + timedelta(minutes=i)) for i in range(90)]

# Prepare lists to store the satellite's elevation and the times when it's visible
elevations = []
visible_times = []
plot_times = []

# Calculate the satellite's position relative to the observer at each time
for t in times:
    difference = satellite - observer       # Relative position
    topocentric = difference.at(t)           # Position as seen from observer
    alt, az, distance = topocentric.altaz() # Get altitude, azimuth, distance
    elevations.append(alt.degrees)           # Store elevation angle in degrees
    plot_times.append(t.utc_datetime())      # Store corresponding datetime
    if alt.degrees > 10:                      # Check if satellite is above 10° elevation
        visible_times.append(t.utc_datetime().strftime("%H:%M"))

# Plot the satellite elevation angle over time
plt.figure(figsize=(10, 5))
plt.plot(plot_times, elevations, color="green", marker="o")
plt.axhline(y=10, color='red', linestyle='--', label='Visibility Threshold (10°)')

# Format the x-axis to show time nicely
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
plt.gca().xaxis.set_major_locator(mdates.MinuteLocator(interval=10))
plt.xticks(rotation=45, fontsize=9)

# Label the plot
plt.xlabel("UTC Time")
plt.ylabel("Elevation Angle (°)")
plt.title(f"{name} Elevation as seen from ({user_lat}, {user_lon})")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Print when the satellite will be visible (above 10° elevation)
if visible_times:
    print(f"\nThe {name} will be visible above 10° at these times (UTC):")
    print(", ".join(visible_times))
else:
    print(f"\nThe {name} will NOT be visible above 10° elevation in the next 90 minutes.")

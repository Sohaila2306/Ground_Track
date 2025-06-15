# Ground Track

This repository includes Python scripts that use the [Skyfield](https://rhodesmill.org/skyfield/) library to visualize satellite orbits and ground tracks.

## Contents

| File             | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| `TLE.py`         | Plots the ground track of a satellite (default: ISS) for the next 90 minutes. |
| `ISS.py`         | Asks for your location and shows the ISS elevation above your horizon over 90 minutes. |
| `Ground_Track.py`| Plots the global ground tracks of 3 satellites (ISS, GPS, NOAA) over 3 orbits. |

---

## Installation

Install the required Python packages:

```bash
pip install skyfield matplotlib numpy cartopy

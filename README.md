# Satellite Orbit Visualizer

Three small Python tools for seeing where satellites are and where they are going. Plot the **ground track** of the International Space Station, a NOAA weather satellite and a GPS satellite on one world map, and find out **when the ISS will pass over your house**. It uses live orbit data from [CelesTrak](https://celestrak.org), [Skyfield](https://rhodesmill.org/skyfield/) for the orbit maths, and Matplotlib (and optionally Cartopy) for the pictures.

<p align="center">
  <img src="docs/images/ground_tracks.png" width="900" alt="Ground tracks of the ISS, NOAA 15 and a GPS satellite over 4.5 hours">
</p>

> The pictures in this README were drawn offline by a simplified stand-in for Skyfield, with no coastlines, so you can see what the output looks like. See [Gallery](#about-the-pictures-in-this-readme). Your own runs use real, current orbit data.

## What's in the box

| Script | What it does |
|---|---|
| `Ground_Track.py` | Plots the ground tracks of the ISS, a NOAA satellite and a GPS satellite on one map |
| `ISS.py` | Asks where you are, then lists every pass of the ISS (or any other satellite) overhead and draws it in your sky |
| `TLE.py` | Shows one satellite's orbit data (period, altitude, inclination) and its ground track |
| `Orbit_Tools.py` | The shared code: downloading and checking orbit data, orbit calculations, pass finding, plotting |

## Install

You need Python 3.9 or newer.

```bash
pip install -r requirements.txt
```

For a proper map with land and coastlines, add Cartopy. It's optional, and the tools draw a clean latitude/longitude grid without it:

```bash
pip install cartopy
```

The scripts download the latest orbit data the first time you run them (so you need an internet connection), and reuse it for a few hours after that. Cartopy also downloads its map outlines the first time it draws a map.

## Use it

### See where the satellites are

```bash
python Ground_Track.py
python Ground_Track.py --minutes 720 --save tracks.png      # twelve hours, saved to a file
python Ground_Track.py --sats iss noaa19 noaa18             # pick your own
```

Each line is a ground track: the point on Earth directly beneath the satellite. A dot marks where it starts, a square where it ends.

### Find out when the ISS passes over you

```bash
python ISS.py                                    # it will ask for your latitude and longitude
python ISS.py --lat 30.0444 --lon 31.2357 --tz Africa/Cairo
python ISS.py --lat 51.5 --lon -0.13 --hours 72 --min-elevation 20
```

Latitude is positive north and longitude is positive east, so Cairo is `30.0444, 31.2357` and Los Angeles is `34.05, -118.24`. The output looks like this (illustrative numbers):

```
ISS (ZARYA) over 30.0444, 31.2357: passes above 10 degrees in the next 48 hours (times in Africa/Cairo)

#  Rises                    Highest        Sets          Length
-  -----------------------  -------------  ------------  ------
1  Thu 01 Jan 12:05:50  SW  12:09:20  81°  12:12:31  NE  6m 41s
2  Thu 01 Jan 20:16:10  NW  20:19:40  79°  20:22:51  SE  6m 41s
3  Fri 02 Jan 11:18:36  SW  11:21:40  43°  11:24:58  NE  6m 22s
```

<p align="center">
  <img src="docs/images/iss_passes.png" width="900" alt="Elevation of the ISS over Cairo for two days, and the path of the best pass across the sky">
</p>

On the left is how high the ISS climbs above your horizon over time. On the right is the best pass as you'd see it looking up: north at the top, east on the right, the zenith in the middle. A pass that goes almost through the centre is nearly overhead.

**A pass above the horizon is not the same as a pass you can see.** To see the ISS with your eyes the sky has to be dark where you are while the station is still lit by the Sun, which only happens in the hour or two after sunset or before sunrise. When it does line up, the ISS is one of the brightest things in the sky and crosses it in about six minutes. These tools tell you when it is overhead, not whether it is lit.

### Look at one satellite's orbit

```bash
python TLE.py                              # the ISS, next 90 minutes
python TLE.py --sat noaa15 --minutes 200
python TLE.py --sat gps --gps-prn 24
python TLE.py --sat 25544 --no-plot        # just print the numbers
```

```
ISS (ZARYA)  (NORAD 25544)
  orbit data : downloaded from CelesTrak just now
  measured   : 2026-01-01 00:00 UTC  (0 days from the start time)
  period     : 92.9 min  (15.50 orbits per day)
  altitude   : 414 to 419 km
  inclination: 51.64 degrees
```

<p align="center">
  <img src="docs/images/iss_90_minutes.png" width="800" alt="One ISS orbit over 90 minutes, with 15-minute time marks">
</p>

Notice that the track doesn't join up with its starting point. During the 93 minutes the ISS takes to go once around, the Earth turns about 23 degrees underneath it, so each orbit is shifted west of the last.

### Options common to all three

| Option | What it does |
|---|---|
| `--start 2026-10-10T18:30` | Start at this UTC time instead of now |
| `--save FILE` | Also save the picture (png, pdf, svg...) |
| `--no-show` | Don't open a window (handy with `--save`, and on machines with no screen) |
| `--offline` | Don't download anything; use data saved by an earlier run |
| `--no-cartopy` | (maps) Draw the plain grid even if Cartopy is installed |

Satellites you can name: `iss`, `noaa15`, `noaa18`, `noaa19`, `gps`, or any NORAD catalogue number such as `25544`. Run any script with `--help` for the full list.

## Why it matters, and how we got here

### A very short history of watching satellites

In the *Principia* (1687), Isaac Newton asked you to imagine firing a cannonball horizontally from a very tall mountain. Fire it faster and it lands farther away; fire it fast enough and the ground curves away beneath it as fast as it falls, and it never lands. That is an orbit. Johannes Kepler had already worked out the shape of those paths (ellipses, in laws published between 1609 and 1619), and the maths in this project is still Newton and Kepler plus a few corrections.

Tracking real satellites began on **4 October 1957**, when the Soviet Union launched Sputnik 1. It carried a simple radio transmitter, and amateur radio operators all over the world could hear the beeps as it passed. Within months the first American satellite, Explorer 1 (1958), had discovered the Van Allen radiation belts, and governments realised that if you were going to put things in orbit, you needed to keep track of every one of them. Today the US Space Force's space surveillance network follows tens of thousands of objects and publishes their orbits, and CelesTrak (run by Dr. T.S. Kelso since 1985) makes that data easy to download.

The format those orbits are published in, the **two-line element set** or TLE, goes back to the days of 80-column punched cards. It is only 138 characters, but together with a matching mathematical model called **SGP4** (developed in the 1970s and described in the famous "Spacetrack Report No. 3" of 1980) it can predict where a satellite will be. The catch is that a TLE is a snapshot. Air drag, and the Sun and Moon's gravity, keep changing the real orbit, so a TLE is only good for a few days to a couple of weeks. This project downloads fresh data for that reason, and warns you when it's old.

### The three satellites, and why they look so different on a map

| | Altitude | Time per orbit | Tilt of the orbit (inclination) |
|---|---|---|---|
| **ISS** | about 420 km | 93 minutes | 51.6° |
| **NOAA** polar orbiters | about 800 km | 101 minutes | about 98.7° |
| **GPS** | 20,200 km | 11 h 58 min | about 55° |

<p align="center">
  <img src="docs/images/orbit_scale.png" width="520" alt="The ISS, NOAA, GPS and geostationary orbits drawn to scale around the Earth">
</p>

**The International Space Station** has been orbiting since its first module launched in November 1998 and has been continuously crewed since November 2000. Its 51.6° inclination was chosen so that Russian rockets launching from Kazakhstan (at 46° north) can reach it, while it is still within reach of American launches from Florida. Because the orbit is tilted 51.6° to the equator, the ISS never goes farther north or south than 51.6°. That limit is why its ground track, in the picture above, is a wave that peaks at that latitude. Its low altitude (it's just skimming the top of the atmosphere) means it slowly loses height from drag and has to be pushed back up periodically, which is one reason its orbit data needs to be fresh.

**NOAA's polar-orbiting weather satellites** (the TIROS series began with the first weather satellite in 1960, and satellites carrying the NOAA name followed from 1970) fly in a nearly polar orbit that takes them almost over the poles, so that as the Earth turns beneath them they eventually see everywhere. They use a clever trick called a **sun-synchronous orbit**. The Earth isn't a perfect sphere; it bulges at the equator, and that bulge slowly twists the orbit's plane around the Earth. At an inclination of about 98.7° and an altitude of about 800 km, that twist works out to almost exactly 0.9856° per day, the same rate at which the Earth goes around the Sun. The satellite therefore passes over each place at the same local solar time every day, so the lighting on its pictures is always the same. The rate of the twist is

```
dΩ/dt = -(3/2) · J₂ · (Rₑ / p)² · n · cos(i)
```

where `J₂` measures the Earth's bulge, `p` is the size of the orbit, `n` is the satellite's average angular speed, and `i` is the inclination. For `i` greater than 90° the cosine is negative, so the twist goes the other way, eastward, which is exactly what's needed. Several of these older satellites, including NOAA 15, have been or are being retired from service, but they are still tracked, so they make good examples.

**GPS satellites** orbit at 20,200 km, a height chosen so that one orbit takes half a sidereal day (11 h 58 min). Each GPS satellite therefore retraces the same ground track every sidereal day (about four minutes earlier each solar day), which makes the constellation very predictable. The first GPS satellite launched in 1978 and the full set of 24 was complete in 1995. The satellites are spread over six orbital planes tilted 55° to the equator so that from almost anywhere on Earth at least four are above the horizon at any time, which is the minimum a receiver needs to work out where it is. Their ground track covers the Earth much more slowly than the lower satellites, as the green line above shows.

### What's actually in a TLE?

This is a real ISS TLE from 2008 (the one used in this project's tests), taken apart:

```
1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927
2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537
```

| Field | Value | Meaning |
|---|---|---|
| Satellite number | `25544` | The satellite's NORAD catalogue number |
| Launch designator | `98067A` | Launched in 1998, the 67th launch of the year, piece A |
| Epoch | `08264.51782528` | When the orbit was measured: 2008, day 264.5178 (20 September, 12:25 UTC) |
| Drag term | `-11606-4` | How strongly drag affects it, used by SGP4 |
| Inclination | `51.6416` | Tilt of the orbit to the equator, in degrees |
| RAAN | `247.4627` | Where the orbit crosses the equator going north, as an angle around the equator from a fixed direction in space |
| Eccentricity | `0006703` | How oval the orbit is (0.0006703: very nearly a circle) |
| Argument of perigee | `130.5360` | Where the lowest point of the orbit is: the angle along the orbit from the equator crossing |
| Mean anomaly | `325.0288` | Where the satellite is on the orbit at the epoch, in degrees |
| Mean motion | `15.72125391` | Orbits per day |
| Last digit of each line | `7` | A check digit. If any character is changed, it no longer adds up |

That last row is what this project uses to check the data. The three TLEs in the first version of this project all failed it, so they had been corrupted or altered somewhere along the way. They were also from June 2024, which would have made any prediction today meaningless. The tools now download fresh orbit data and reject anything that doesn't pass.

### What this is useful for

The same calculations are used on real missions every day:

- **Planning contact with a ground station.** A satellite can only talk to a ground station while it's above the horizon, so operators work out the passes in advance (this is what `ISS.py` does).
- **Earth observation.** Where a satellite flies determines what it can image and how often it comes back. Sun-synchronous orbits like NOAA's exist for exactly that reason.
- **Navigation.** A GPS receiver is working out where the satellites are from the same kind of orbit data.
- **Hobby radio and astronomy.** Amateurs receive weather pictures from NOAA satellites on 137 MHz with a simple antenna, and look up when the ISS is going over for photography or just to wave.
- **Learning orbital mechanics.** Seeing why an inclined orbit makes a wave and why a polar orbit makes a spiral is the quickest way to understand it.

## How accurate is it?

- **It depends on how fresh the orbit data is.** Right after the TLE's epoch, a low-orbit prediction is good to about a kilometre or so. The error grows by a few kilometres each day, mostly because of air drag. After a week or two the position along the track can be off by a noticeable amount. The scripts print a warning once the data is more than 14 days from the time you ask about.
- **Maneuvers aren't predicted.** The ISS regularly fires its engines to climb back up and to dodge debris, and no TLE can know about that in advance.
- **"Passes" means above the horizon.** See the note under `ISS.py`: it doesn't know whether the sky is dark or the satellite is lit.
- **Hills, buildings and trees aren't modelled,** which is why the default minimum elevation is 10° rather than 0°.

## Tests

```bash
python tests/test_orbit_tools.py       # no extra packages needed
# or, with pytest installed:
python -m pytest
```

The 44 tests run in about two seconds. They check the TLE reader against a published TLE (checksum, epoch, orbital elements), the downloading and caching (including a dropped connection), the date-line cutting, the pass finder (with the rise and set times worked out exactly), basic orbital mechanics (the latitude limit, the westward shift per orbit, the orbit's period), and all three scripts including bad input. The tests themselves were checked by deliberately breaking the code in several ways (a wrong checksum rule, a wrong date-line threshold, missing interpolation) and confirming each change made a test fail.

One honest limit: the tests don't use the real Skyfield or Cartopy. `tests/fake_skyfield.py` stands in for Skyfield (a simplified circular-orbit model), so the tests prove that this project's logic is right, but not that Skyfield and Cartopy behave as expected. Those two libraries are used in only a handful of places, in `Orbit_Tools.py`, through their standard documented interfaces.

## About the pictures in this README

They were generated offline by `docs/make_readme_figures.py`, which runs the real scripts but with the stand-in propagator from `tests/fake_skyfield.py` and orbital elements from the original project's TLEs. That means the shapes are right (the ISS limited to ±51.6°, NOAA reaching ±81°, GPS moving slowly) but they are illustrations, not real predictions, and there are no coastlines because those come from Cartopy. Run the scripts yourself for the real thing.

## Project layout

```
Ground_Track.py       global ground tracks of several satellites
ISS.py                passes of a satellite over your location
TLE.py                one satellite's orbit data and ground track
Orbit_Tools.py        shared code (TLEs, Skyfield, pass finding, plots)
requirements.txt
tests/                test suite and the Skyfield stand-in
docs/                 README pictures and the script that makes them
```

## Credits and further reading

- [Skyfield](https://rhodesmill.org/skyfield/), by Brandon Rhodes, does the orbit calculations (using the `sgp4` package)
- [CelesTrak](https://celestrak.org), by Dr. T.S. Kelso, supplies the orbit data
- [Cartopy](https://scitools.org.uk/cartopy/) and [Matplotlib](https://matplotlib.org) draw the maps
- D. Vallado, P. Crawford, R. Hujsak and T.S. Kelso, *Revisiting Spacetrack Report #3* (2006): how SGP4 and TLEs work
- D. Vallado, *Fundamentals of Astrodynamics and Applications*: the standard textbook
- Wikipedia's article on the [two-line element set](https://en.wikipedia.org/wiki/Two-line_element_set) is a good plain-language introduction

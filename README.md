
Photovoltaic System Engineering Final Project

Zayne Khouja

18-883, Autumn 2021

Carnegie Mellon University


## Brief: 
This repository contains a python-based, trip-planning utility for use with an off-grid photovoltaic system installed on a van, RV, truck, or other road-trip vehicle.


## Goal: 
Given a long-distance road trip route plotted in Google Maps, this program computes the average daily drive time along your route to ensure that the battery wired to your PV system is fully charged by the end of the day's drive.


## Motivation: 
Van-life is the hot new lifestyle of the decade, yet actually living out of your vehicle requires a consistent source of power to run appliances such as stoves/coolers, power heaters/fans for regulating temperature, or charge your laptop so you can work remotely.  A roof-mounted solar panel wired to a deep-cycle battery, then run through an inverter, can supply passively generated electricity for all these applications.  Yet, depending on the season and where you are on the road, the output of your PV system may be unreliable or insufficient.  This program attempts to quantify how many hours your solar panels need direct sunlight depending on your road trip location so you can drive assured that your battery will be fully charged when you park for the night.


## Usage: 
This tool is meant as a next step in evaluating your road trip route after you've plotted it in Google Maps.  Once your route is designed, simply copy the google maps URL and use it as input to this program.  Specifically, from within the './code/' directory, run 'python project.py 
\<URL\>' where \<URL\> is the google maps url you planned. Additionally, after the URL, you can include an optional argument to indicate the month you plan to road trip in.  Specify this as 'python project.py \<URL\> \<month\>' where \<month\> can be the full month name (ie. "January") or the abbreviation (ie. "Jan").  If no month is provided, data for July is used to demonstrate the best-case result for a road trip in the middle of the summer.

Once executed, this program prints various logging outputs and generates 3 plots.  First, it parses out the separate waypoint locations along the route and plots them, then it interpolates a straight line between these points and plots that route.  Finally, it maps each coordinate to an irradiance and temperature value for the given location from data collected on renewables.ninja and uses this to resolve the current generated by the solar panel. With this info and some stats on your battery, it computes the expected daily drive time to fully charge the panel along your route and plots this data with separate legs of the trip labeled with their drive time in hours per day.

You can then use this final generated plot to determine the pace of your road trip as you traverse the country and understand roughly how long you must drive to make sure you can power all your appliances once you stop driving for the night.


## Implementation: 
This project was written in python and uses Google's geocoding (and reverse geocoding) APIs to parse a google maps route and convert locations to their coordinates.  It then maps these coordinates to the states they fall within.  Then, using state-level data from renewables.ninja, it maps these states to irradiance and temperature values per month, and uses these values to compute drive time.


## Troubleshooting: 
Running this code requires several python packages.  Specifically, you should install the following:
- 'pip install requests' 
- 'brew install geos'
- 'pip install https://github.com/matplotlib/basemap/archive/master.zip'
- 'pip install geopy'
- 'pip install us'
- 'pip install reverse_geocoder'
  
Additionally, many shells don't like arguments with "!" in them, so if you get an error like "bash: !: event not found" when running the code with a google maps URL, run "set +H" to fix the issue and try again. 

More debugging info can be found in the comment in ./code/project.py.  If problems persist, email zck@andrew.cmu.edu or zayne98@gmail.com


## Future improvements: 
Currently, the accuracy of this project is primarily limited by the APIs available to me.  Specifically, I have to use state-level irradiance data.  With a paid subscription to a service like SoDa (http://www.soda-pro.com/), I could use an API to programmatically request specific irradiance data accurate to exact coordinates, for specific panel tilts, and for specific days and times.  I contacted SoDa requesting free, temporary, student access, but they have yet to respond.

Additionally, I convert a google maps route to a series of waypoints connected by straight lines, which serves as a rough approximation of a route but doesn't incorporate turns and irregularities in a route.  Google Map's more powerful API requires payment as well, but with access to that, the route predictions would be more accurate.

Finally, currently only routes within the continental US are supported due to the limitations of my irradiation data.  This could be extended to any route given access to the SoDa API.

Essentially, with paid access to SoDa's and Google Map's API, the output of this project would be really quite accurate and could serve as quite a powerful trip-planning tool.


## Example: 
In the directory './example/' you can find plots demonstrating this program's output for a cross-country route from Augusta, Maine to San Diego, California.  The google maps route can be found at [this link](https://www.google.com/maps/dir/Augusta,+ME/Albany,+NY/Pittsburgh,+PA/Indianapolis,+IN/Oklahoma+City,+OK/Albuquerque,+NM/Phoenix,+AZ/San+Diego,+CA/@37.1101918,-111.5280353,4z/data=!3m1!4b1!4m50!4m49!1m5!1m1!1s0x4cb200fdafacc49d:0x79a3488d64220b2d!2m2!1d-69.7794897!2d44.3106241!1m5!1m1!1s0x89de0a34cc4ffb4b:0xe1a16312a0e728c4!2m2!1d-73.7562317!2d42.6525793!1m5!1m1!1s0x8834f16f48068503:0x8df915a15aa21b34!2m2!1d-79.9958864!2d40.4406248!1m5!1m1!1s0x886b50ffa7796a03:0xd68e9df640b9ea7c!2m2!1d-86.158068!2d39.768403!1m5!1m1!1s0x87ad8a547ef8d281:0x33a21274d14f3a9d!2m2!1d-97.5164276!2d35.4675602!1m5!1m1!1s0x87220addd309837b:0xc0d3f8ceb8d9f6fd!2m2!1d-106.650422!2d35.0843859!1m5!1m1!1s0x872b12ed50a179cb:0x8c69c7f8354a1bac!2m2!1d-112.0740373!2d33.4483771!1m5!1m1!1s0x80d9530fad921e4b:0xd3a21fdfd15df79!2m2!1d-117.1610838!2d32.715738!3e0).  Here is a breakdown of the maps found in the 'example' directory:
- route.png: Screenshot of the route plotted in google maps
- waypoints.png: Extraction of the waypoints from the google maps route plotted on a basemap of the continental US
- interpolation.png: Shows the linear interpolation between route waypoints used to approximate the driving route
- january/march/july/october/november.png: Shows the route map with approximate charge times for each leg of the trip during different months.  Shows how irradiation levels very seasonally and latitudinally. 


# This file contains code that, when given a google maps route (provided as a
#   URL), will determine the average amount of drive time per day to maintain
#   full charge in your solar system's battery while still progressing along the
#   route at the most liesurly pace possible.

# If testing this locally, follow these steps:
#   1. This project requires a number of packages.  You should have "requests",
#       geos, and basemap installed. If you do not, try the following commands:
#           'pip install requests' 
#           'brew install geos'
#           'pip install https://github.com/matplotlib/basemap/archive/master.zip'
#           'pip install geopy'
#           'pip install us'
#           'pip install reverse_geocoder'
#       Note: These instructions probably only work on Linux/Macs.  
#           I'm not sure about Windows
#   2. Plot a route on google maps and copy the URL
#       An example route from Pittsburgh to Washington DC is:
#       https://www.google.com/maps/dir/Pittsburgh,+PA/Washington,+DC,+DC/@39.6690482,-79.6376097,8z/data=!3m1!4b1!4m14!4m13!1m5!1m1!1s0x8834f16f48068503:0x8df915a15aa21b34!2m2!1d-79.9958864!2d40.4406248!1m5!1m1!1s0x89b7c6de5af6e45b:0xc2524522d4885d2a!2m2!1d-77.0368707!2d38.9071923!3e0
#   3. From the command line, navigate to the folder containing this file.
#   4. Run 'python project.py <URL>' replacing URL with your google maps route
#       URL. No '<>' or quotes needed.  Just the URL as you copied it.
#       Example: python project.py https://www.google.com/maps/dir/Pittsburgh,+PA/Washington,+DC,+DC/@39.6690482,-79.6376097,8z/data=!3m1!4b1!4m14!4m13!1m5!1m1!1s0x8834f16f48068503:0x8df915a15aa21b34!2m2!1d-79.9958864!2d40.4406248!1m5!1m1!1s0x89b7c6de5af6e45b:0xc2524522d4885d2a!2m2!1d-77.0368707!2d38.9071923!3e0
#       NOTE: You must have a Wifi connection to run this.  It pings google.com
#
#   This will output... TODO

# NOTE: If you're running this from a bash shell and recieve an error 
#     similar to "bash: !: event not found", you should run 'set +H' 
#     from the command line.  Bash doesn't like all the pesky "!" in the
#     metadata of the URL.
#     Alternatively, just manually delete everything after "/data=" in 
#     the URL to get rid of all the "!"
#
# If you have any other troubles running this, email me (Zayne Khouja)
#   at zck@andrew.cmu.edu

from os import stat
import requests
import sys
import math
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from geopy.geocoders import GoogleV3
import pandas as pd
import reverse_geocoder
import us
import numpy as np


# I created a developer account for google's cloud services.  I will be using
#   googles Geocoding API to convert named locations to their coordinates.
# I believe anyone should be able to use this key as I have not restricted it
# This key will expire on February 21, 2022
GOOGLE_API_KEY = "AIzaSyA7CzmKB2lQpQl-RYVY9cjRriYcjFJJCg0"

# Defines coordinate degree granularity. One coordinate degree is roughly 
#   69 miles (111 km).  The below definition defines how far apart to place 
#   additional coordinates along a route between two locations.  A lower number 
#   will yield more granular data, but will take longer to compute. 
# The unit of this value is lat/long degrees.  A value of 0.25 defines a
#   coordinate roughly every 17 miles (28 km).
INTERPOLATION_GRANULARITY = 0.25

# Defines for plotting coordinates on the map
BIG_DOT = 5
SMALL_DOT = 2

# Standardizing the date and time to use for resolving irradiance
# Using the Date/Time with the highest irradiance (from Homework 1)
# As an approximation for the lower bound of drive time
# This represents July 1st
MONTH = "7"
DAY = "1"
IRRADIANCE_PATH = "./irradiance.csv"



# Converts a given address (formatted as per a google URL's input) to its
#   respective longitude and lattitude coordinate.
def get_coordinates(addr):
    lat, lng = None, None
    base_url = "https://maps.googleapis.com/maps/api/geocode/json"
    endpoint = f"{base_url}?address={addr}&key={GOOGLE_API_KEY}"

    r = requests.get(endpoint)

    if r.status_code not in range(200, 299):
        print("\n*** Google API Response Error!")
        print("*** Check URL formatting, or try again later (API may be down)\n\n")
        quit()

    results = r.json()['results'][0]
    lat = results['geometry']['location']['lat']
    lng = results['geometry']['location']['lng']
   
    return lat, lng



# Given a Google maps route URL as input, resolves a list of coordinates for  
#   the locations in the route
def parse_url(url):
    # Content between "https://www.google.com/maps/dir/" and "/@" represents
    #   the locations of interest.  Rest is metadata that we don't care about
    url = url.split("https://www.google.com/maps/dir/")[1]
    url = url.split("/@")[0]

    locations = url.split("/")
    print("Parsed", str(len(locations)), "Locations:", str(locations), "\n")

    coordinates = [get_coordinates(loc) for loc in locations]
    print("Resolved",str(len(coordinates)), "coordinates as:", str(coordinates), "\n")
    return coordinates



# Groups adjacent coordinates into pairs to define legs of the route
def build_pairs(coords):
    print("Grouping", str(len(coords)), "locations into", str(len(coords) - 1), "pair(s)\n")
    pairs = []
    for i in range(len(coords)-1):
        pairs.append((coords[i], coords[i+1]))
    assert(len(pairs) == (len(coords) - 1))
    return pairs



# Given a list of coordinates, resolves a straight line between each adjacent 
#   pair of points, then defines additional coordinates every 
#   INTERPOLATION_GRANULARITY degrees along the line. (constant defined above)
# This provides a degree of granularity in the irradiance calculations along
#   the route.
def interpolate(coord_pairs):
    interpolated_coords = []
    for ((x1, y1), (x2, y2)) in coord_pairs:
        new_cords = [(x1, y1)]
        dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        slope = (y2 - y1) / (x2 - x1)
        num_coords = int((dist // INTERPOLATION_GRANULARITY) - 1)

        # Will be roughly the same as INTERPOLATION_GRANULARITY, but not exact
        #   because we want equally spaced points
        coord_dist = dist / (num_coords + 1)

        incremental_dist = coord_dist
        for i in range(num_coords):
            dist_ratio = incremental_dist / dist
            new_x = ((1 - dist_ratio) * x1) + (dist_ratio * x2)
            new_y = ((1 - dist_ratio) * y1) + (dist_ratio * y2)
            new_cords.append((new_x, new_y))
            incremental_dist = incremental_dist + coord_dist

        new_cords.append((x2, y2))
        interpolated_coords.append(new_cords)
    return interpolated_coords


# Given a list of (long, lat) coordinates that fall within the continental US, 
#   plot them on a map!
# Also allows plotting of a list of lists of coordinates (so I can call this
#   with more flexible input types)
# If plotting many points, use small dot_size (2).  Use 5 for fewer points
def plot_coordinates(coords, dot_size):
    # If given a list of lists, flatten it to a single list
    if (any(isinstance(coord, list) for coord in coords)):
        coords = [item for sublist in coords for item in sublist]

    m = Basemap(llcrnrlon=-121,llcrnrlat=20,urcrnrlon=-62,urcrnrlat=51,
        projection='lcc',lat_1=32,lat_2=45,lon_0=-95)
    plt.figure(num=1,figsize=(9,7))
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawmapboundary(fill_color='paleturquoise')
    m.fillcontinents(color='peachpuff',lake_color='paleturquoise')

    for (lat, long) in coords:
        (x, y) = m(long, lat)
        m.plot(x,y,marker='o',color='Red',markersize=dot_size)

    plt.show()


# Uses Google's reverse geocoding API to obtain the state a coordinate lands in
# This sometimes returns google's "Plus Codes" for some reason, which don't 
#   carry the state info annoyingly.  Another solution is in the function below
# def reverse_geocode(coord):
#     (lat, long) = coord
#     geolocator = GoogleV3(api_key=GOOGLE_API_KEY)
#     location = geolocator.reverse(str(lat) + ", " + str(long))
#     if location:
#         print(location.raw)
#         loc = str(location[0])
#         # print(loc)
#         loc = loc.split(", USA")[0]
#         loc = loc.split(", ")[-1]
#         return loc[:2]
#     else:
#         print("\n*** Reverse Geocoding API Response Error!")
#         print("*** Check coordinates, or try again later (API may be down)")
#         print("*** Coordinate must land within a USA State\n\n")
#         quit()


# Uses the reverse_geocode library and the US state library to determine the
#   state a coordinate lies in and then convert state name to state abbreviation
# This solution is sadly slower than Google's API, but it is more reliable
def reverse_geocode(coord):
    coord_data = reverse_geocoder.search(coord)
    state_name = us.states.lookup(coord_data[0]["admin1"])
    if (state_name is None):
        return "DC"
    else:
        return state_name.abbr


# Unpacks a list of lists of coordinates and maps each to the US state it lands in
def coords_to_states(coord_lists):
    total_state_list = []
    for coords in coord_lists:
        single_state_list = []
        for coord in coords:
            single_state_list.append(reverse_geocode(coord))
        total_state_list.append(single_state_list)
    print("Resolved states for coordinates as:", str(total_state_list), "\n")
    return total_state_list


# Using a single date and time right now for simplicity purposes, but should
#   be simple enough to build a dynamic time/date calculator for the 
#   progression of a road trip
# I select the maximum irradiance for each state (regardless of date), to 
#   get a lower bound on drive time
# Returns the temperature for that day as well
def state_to_irrad_temp(state_lists):
    df = pd.read_csv(IRRADIANCE_PATH, header = 0, dtype={1: float, 2: float, 
        3: float, 4: 'str'})
    df['time'] =  pd.to_datetime(df['time'])

    total_list = []
    for states in state_lists:
        single_list = []
        for state in states:
            state_df = df[df["state"] == state]
            state_df = state_df.reset_index(drop=True)
            max_entry = state_df.iloc[state_df["irradiance_surface"].idxmax()]
            irrad = max_entry["irradiance_surface"]
            temp = max_entry["temperature"]
            single_list.append((irrad,temp))
        total_list.append(single_list)
    print("Resolved irrad and temp values for states as:", str(total_list), "\n")
    return total_list



# Constants for Max Power Point calculations:
V_OC_0 = 45.9 #V
I_SC_0 = 9.07 #A
GAMMA = -0.32 #%/C
LAMBDA = -0.41 #%/C
DELTA = 0.058
NOCT = 44 #C
N_CONST = 1.1
K_CONST = 1.38 * math.pow(10, -23) #J/K
E_CONST = 1.602 * math.pow(10,-19) #A
R_S = 0.001 #Ohm
R_SH = 1000 #Ohm

# ===============
# Next several functions involve computing the Max Power Point
# ===============

def get_T_cell(G, T_amb):
    return T_amb + (G * (NOCT - 20) / 800) 

def get_V_oc(G, T_cell):
    term1 = V_OC_0 / (1 + (DELTA * np.log(1000 / G)))
    term2 = 1 + (GAMMA * (T_cell - 25))
    return abs(term1 * term2)

def get_I_sc(G, T_cell):
    term1 = I_SC_0 * (G / 1000)
    term2 = 1 + (LAMBDA * (T_cell - 25))
    return abs(term1 * term2)

def get_FF_0(V_oc, T_amb):
    V_T = N_CONST * K_CONST * T_amb / E_CONST
    num = np.log(abs((V_oc / V_T) + 0.72))
    denom = 1 + (V_oc / V_T)
    return 1 - (num / denom)

def get_R_CH(V_oc, I_sc):
    return V_oc / I_sc

def get_r_S(R_CH):
    return R_S / R_CH

def get_r_SH(R_CH):
    return R_SH / R_CH

def get_FF(FF_0, V_oc, I_sc):
    R_CH = get_R_CH(V_oc, I_sc)
    r_S = get_r_S(R_CH)
    r_SH = get_r_SH(R_CH)
    return FF_0 * (1 - r_S) * (1 - (1 / r_SH))


def get_MPP_Isc(G_T_pair):
    (g, T_amb) = G_T_pair
    T_cell = get_T_cell(g, T_amb)

    V_oc = get_V_oc(g, T_cell)
    I_sc = get_I_sc(g, T_cell)

    ff_0 = get_FF_0(V_oc, T_amb)
    ff = get_FF(ff_0, V_oc, I_sc)

    mpp = ff * V_oc * I_sc
    return (mpp, I_sc)

# ===============
# End MPP-related functions
# ===============


# Given the Irradiance and Temperature for the coordinates, resolves the max
#   power point!
def irrad_temp_to_MPP_Isc(irrad_temp_lists):
    total_list = []
    for irrad_temp_list in irrad_temp_lists:
        single_list = []
        for G_T_amb_pair in irrad_temp_list:
            mpp_isc = get_MPP_Isc(G_T_amb_pair)
            single_list.append(mpp_isc)
        total_list.append(single_list)
    print("Resolved MPP and I_sc as:", str(total_list), "\n")
    return total_list


# Constants for a deep cycle battery
# I referenced a Duralast Deep Cycle Marine Battery for these numbers:
# https://www.autozone.com/batteries-starting-and-charging/rv-battery/p/duralast-29dp-dl-group-29-deep-cycle-marine-and-rv-battery/95824_0_0
AMP_HOUR = 65

# Computes the expected charge time to fill a deep cycle battery
def get_charge_time(I_sc):
    return AMP_HOUR / I_sc

# Given the Isc and MPP, maps them to expeced charge times
def MPP_Isc_to_charge_time(mpp_isc_lists):
    total_list = []
    for mpp_isc_list in mpp_isc_lists:
        single_list = []
        for (mpp, I_sc) in mpp_isc_list:
            charge_time = get_charge_time(I_sc)
            single_list.append(charge_time)
        total_list.append(single_list)
    print("Resolved charge times as:", str(total_list), "\n")
    return total_list

# TODO: Figure out how to attach charge times to the map to understand distance
# TODO: Understand why my numbers seem so very wrong :(


# Main driver function to organize computations
def driver(url, do_plot):
    coordinates = parse_url(url)
    if do_plot:
        plot_coordinates(coordinates, BIG_DOT)
    coordinate_pairs = build_pairs(coordinates)
    granular_coordinates = interpolate(coordinate_pairs)
    if do_plot:
        plot_coordinates(granular_coordinates, SMALL_DOT)
    state_map = coords_to_states(granular_coordinates)
    irrad_temp_map = state_to_irrad_temp(state_map)
    MPP_Isc_map = irrad_temp_to_MPP_Isc(irrad_temp_map)
    charge_times_map = MPP_Isc_to_charge_time(MPP_Isc_map)


# Code to parse command-line arguments 
# Pass a random second argument if you want to skip the plotting part
if len(sys.argv) > 1:
    print("=====")
    do_plot = True
    if (len(sys.argv) > 2):
        do_plot = False
    driver(sys.argv[1], do_plot)
    print("Execution complete.")
    print("=====")
else:
    print("=====\nNo arguments provided.")
    print("Please supply the URL from a google maps route.")
    print("Usage: python project.py <URL>\n=====")

    
# Short route for testing
"""
https://www.google.com/maps/dir/Pittsburgh,+PA/Washington,+DC,+DC/@39.6690482,-79.6376097,8z/data=!3m1!4b1!4m14!4m13!1m5!1m1!1s0x8834f16f48068503:0x8df915a15aa21b34!2m2!1d-79.9958864!2d40.4406248!1m5!1m1!1s0x89b7c6de5af6e45b:0xc2524522d4885d2a!2m2!1d-77.0368707!2d38.9071923!3e0
"""

# Long route for testing (part of my winter break return route)
"""
https://www.google.com/maps/dir/Pittsburgh,+PA/Hot+Springs+National+Park/Magazine+Mountain,+Arkansas/Petrified+Forest+National+Park/Mesa+Verde+National+Park,+Mesa+Verde,+CO/Capitol+Reef+National+Park,+Utah/Great+Basin+National+Park,+Nevada/Death+Valley,+CA/San+Jose,+CA/@36.1182533,-119.0157453,4z/data=!3m1!4b1!4m56!4m55!1m5!1m1!1s0x8834f16f48068503:0x8df915a15aa21b34!2m2!1d-79.9958864!2d40.4406248!1m5!1m1!1s0x87cd2a4f3e14f595:0x582d6639762948dd!2m2!1d-93.0423545!2d34.5216915!1m5!1m1!1s0x87cc7df3a3deba49:0xa00510e07fa0fe9e!2m2!1d-93.6449152!2d35.167312!1m5!1m1!1s0x872f9b3e3ba72d9b:0xf026f23d8c8c2ecd!2m2!1d-109.78204!2d35.065931!1m5!1m1!1s0x873960bf2ed7711f:0x79f695a21bf61863!2m2!1d-108.4618335!2d37.2308729!1m5!1m1!1s0x874a00ff07e7a253:0xde3bec53484fff07!2m2!1d-111.1354983!2d38.0877312!1m5!1m1!1s0x80b15c25d7cc0d03:0x3cd4750fafbebd31!2m2!1d-114.2633787!2d38.9299798!1m5!1m1!1s0x80c739a21e8fffb1:0x1c897383d723dd25!2m2!1d-116.9325408!2d36.5322649!1m5!1m1!1s0x808fcae48af93ff5:0xb99d8c0aca9f717b!2m2!1d-121.8863286!2d37.3382082!3e0
"""

# East Coast vertical route
"""
https://www.google.com/maps/dir/Bangor,+ME+04401/Boston,+MA/Philadelphia,+PA/Cape+Charles,+VA/Wilmington,+NC/Charleston,+SC/Jacksonville,+FL/Miami,+FL/@34.9822916,-84.2372354,5z/data=!3m1!4b1!4m50!4m49!1m5!1m1!1s0x4cae4b46101129bd:0x4d0918b0a7af7677!2m2!1d-68.7712257!2d44.8016128!1m5!1m1!1s0x89e3652d0d3d311b:0x787cbf240162e8a0!2m2!1d-71.0588801!2d42.3600825!1m5!1m1!1s0x89c6b7d8d4b54beb:0x89f514d88c3e58c1!2m2!1d-75.1652215!2d39.9525839!1m5!1m1!1s0x89ba5c609acedfa3:0x52caf608b4c59f27!2m2!1d-76.0174336!2d37.267916!1m5!1m1!1s0x89a9f5a20debaed5:0x5e66493884093032!2m2!1d-77.8868117!2d34.2103894!1m5!1m1!1s0x88fe7a42dca82477:0x35faf7e0aee1ec6b!2m2!1d-79.9310512!2d32.7764749!1m5!1m1!1s0x88e5b716f1ceafeb:0xc4cd7d3896fcc7e2!2m2!1d-81.655651!2d30.3321838!1m5!1m1!1s0x88d9b0a20ec8c111:0xff96f271ddad4f65!2m2!1d-80.1917902!2d25.7616798!3e0
"""
# This is where the final function for optimization will live, as well as the initial conditions 
import warnings
import math
import numpy as np
import csv
import random

# Plotting Imports
import shapely
import cartopy.crs as ccrs
import cartopy.geodesic
import matplotlib.pyplot as plt

# Brahe Imports
import brahe as bh
import brahe.data_models as bdm
import brahe.access.access as ba
import src.functions as functions

from shapely.geometry import Polygon, MultiPolygon
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point

from scipy.optimize import minimize, Bounds, differential_evolution 
import pyswarms as ps 

def diff_ev(num_sats, deleted): 
    # Initialize heights and altitudes for satellites 
    n = num_sats
    alt_sats_orig = [650 for _ in range(n)]
    inclination_orig = np.arange(0, 180, 180/n).tolist()
    sats_initial = []
    sats_initial.extend(alt_sats_orig)
    sats_initial.extend(inclination_orig)

    coords = functions.calculate_lat_lon(sats_initial)
    original_constellation = functions.compute_area(coords, alt_sats=alt_sats_orig, inclination = inclination_orig) 

    # print("Original Area: ", original_constellation.area)
    # print("Original Satellites: ", sats_initial)
    # print("Original Coordinates", coords)
    # functions.write_polygon(original_constellation, "original_constellation")
    with open('original_coords.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(coords)

    # lose a satellite 
    alt_sats_remaining = [alt_sats_orig[i] for i in range(len(alt_sats_orig)) if i not in deleted]
    inclination_remaining = [inclination_orig[i] for i in range(len(inclination_orig)) if i not in deleted]

    sats_remaining = []
    sats_remaining.extend(alt_sats_remaining) 
    sats_remaining.extend(inclination_remaining)

    coords_remaining = functions.calculate_lat_lon(sats_remaining)
    remaining_constellation = functions.compute_area(coords_remaining, alt_sats=alt_sats_remaining, inclination = inclination_remaining)


    def obj_func(input_vec):
        coords = functions.calculate_lat_lon(input_vec)
        alt_sats = input_vec[:n-1]
        inclination = input_vec[n-1:]
        constellation = functions.compute_area(coords, alt_sats, inclination)
        area_overlap = Polygon()
        # for sat in constellation: 
        #     
        area_overlap = area_overlap.union(constellation)
        # for x in original_area: 
        area_overlap = area_overlap.intersection(original_constellation)
        
        if area_overlap.area == 0:
            return -1e-9
        else:
            return -area_overlap.area

    # initial guess 
    x0 = sats_remaining

    n = n - len(deleted)
    # Define lower and upper bounds for alt_sats and inclination
    lower_bounds = [500.0]*(n-1) + [0.00000001]*(len(x0)-(n-1))
    upper_bounds = [650.0]*(n-1) + [180.0]*(len(x0)-(n-1))
    bounds = Bounds(lower_bounds, upper_bounds)

    # Diff_Ev + LBFGS-B

    lbfgs_sats = differential_evolution(obj_func, bounds, tol=0.01, maxiter=100, disp = True)


    coords_lbfgs = functions.calculate_lat_lon(lbfgs_sats.x)
    alt_sats_lbfgs = lbfgs_sats.x[:n-1]
    inclination_lbfgs = lbfgs_sats.x[n-1:]
    lbfgs_constellation = functions.compute_area(coords_lbfgs, alt_sats=alt_sats_lbfgs, inclination = inclination_lbfgs) 

    # print("LBFGS-B Satellites", lbfgs_sats.x)
    # print("LBFGS Area: ", lbfgs_constellation.area)
    print("Percent Difference from Original: ", 100*(original_constellation.area - lbfgs_constellation.area)/original_constellation.area)
    # print("LBFGS-B Coordinates", coords_lbfgs)
    # functions.write_polygon(lbfgs_constellation, "lbfgs_constellation")

    lbfgs_percentdiff = 100*(original_constellation.area - lbfgs_constellation.area)/original_constellation.area

    return lbfgs_percentdiff


# number of original satellites
n = 20
# lose this many satellites
missing = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17] 

with open('avg_diffev_percent_diffs.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Num_Missing', 'DiffEv_LBFGS-B'])

    for i in range(len(missing)-1, -1, -1):
        print(i)
        missing_copy = missing[:i]
        print(missing_copy)
        num_missing = len(missing_copy)
        print(num_missing)
        percent_diff = diff_ev(n, missing_copy)
        writer.writerow([num_missing, percent_diff])
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

def optimize(num_sats, deleted): 
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

    # lose a satellite 
    alt_sats_remaining = [alt_sats_orig[i] for i in range(len(alt_sats_orig)) if i not in deleted]
    inclination_remaining = [inclination_orig[i] for i in range(len(inclination_orig)) if i not in deleted]

    sats_remaining = []
    sats_remaining.extend(alt_sats_remaining) 
    sats_remaining.extend(inclination_remaining)

    coords_remaining = functions.calculate_lat_lon(sats_remaining)
    remaining_constellation = functions.compute_area(coords_remaining, alt_sats=alt_sats_remaining, inclination = inclination_remaining)

    # print("Remaining Area: ", remaining_constellation.area)
    # print("Percent Difference: ", 100*(original_constellation.area - remaining_constellation.area)/original_constellation.area)
    remaining_percentdiff = 100*(original_constellation.area - remaining_constellation.area)/original_constellation.area
    # optimize the area_coverage objective function 
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
        
        return -area_overlap.area

    # initial guess 
    x0 = sats_remaining

    n = n - len(deleted)
    # Define lower and upper bounds for alt_sats and inclination
    lower_bounds = [500.0]*(n-1) + [1e-8]*(len(x0)-(n-1))
    upper_bounds = [650.0]*(n-1) + [180.0]*(len(x0)-(n-1))
    bounds = Bounds(lower_bounds, upper_bounds)

    # Differential Evolution + LBFGS-B

    # for method in ['L-BFGS-B']: 
    # for i in range(5):
    #     lbfgs_sats = differential_evolution(obj_func, bounds, tol=0.01, maxiter=100, disp = True)
    #     # lbfgs_sats = minimize(obj_func, x0, bounds=bounds, method=method,
    #     #         options={'xatol': 1e-8, 'disp': True})
    #     if lbfgs_sats.success:
    #         break

    # lbfgs_sats = differential_evolution(obj_func, bounds, tol=0.01, maxiter=100, disp = True)

    # coords_lbfgs = functions.calculate_lat_lon(lbfgs_sats.x)
    # alt_sats_lbfgs = lbfgs_sats.x[:n-1]
    # inclination_lbfgs = lbfgs_sats.x[n-1:]
    # lbfgs_constellation = functions.compute_area(coords_lbfgs, alt_sats=alt_sats_lbfgs, inclination = inclination_lbfgs) 
    # lbfgs_percentdiff = 100*(original_constellation.area - lbfgs_constellation.area)/original_constellation.area

    # print("LBFGS-B Satellites", lbfgs_sats.x)
    # print("LBFGS Area: ", lbfgs_constellation.area)
    # print("Diff Ev Percent Difference from Original: ", 100*(original_constellation.area - lbfgs_constellation.area)/original_constellation.area)
    # # print("LBFGS-B Coordinates", coords_lbfgs)
    # functions.write_polygon(lbfgs_constellation, "lbfgs_constellation")

    # particle swarm optimization
    # Set up hyperparameters
    options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

    # Create bounds
    bounds = [(lower, upper) for lower, upper in zip(lower_bounds, upper_bounds)]
    max_bound = np.array(upper_bounds)
    min_bound = np.array(lower_bounds)
    bounds = (min_bound, max_bound)

    # Initialize swarm
    optimizer = ps.single.GlobalBestPSO(n_particles=1, dimensions=len(x0), options=options, bounds=bounds)

    # Perform optimization
    cost, pso_sats = optimizer.optimize(obj_func, iters=1000)
    coords_pso = functions.calculate_lat_lon(pso_sats)
    alt_sats_pso = pso_sats[:n-1]
    inclination_pso = pso_sats[n-1:]
    pso_constellation = functions.compute_area(coords_pso, alt_sats=alt_sats_pso, inclination = inclination_pso)

    pso_percentdiff = 100*(original_constellation.area - pso_constellation.area)/original_constellation.area
    # print("PSO Satellites: ", pso_sats)
    # print("PSO Area: ", pso_constellation.area)
    print("PSO Percent Difference from Original: ", 100*(original_constellation.area - pso_constellation.area)/original_constellation.area)
    # # print("PSO Coordinates", coords_pso)
    # functions.write_polygon(pso_constellation, "pso_constellation")

    random_sats = []
    for i in range(n-1):
        random_sats.append(random.uniform(500, 650))
    for i in range(n-1):
        random_sats.append(random.uniform(0, 180))

    coords_random = functions.calculate_lat_lon(random_sats)
    alt_sats_random = random_sats[:n-1]
    inclination_random = random_sats[n-1:]
    random_constellation = functions.compute_area(coords_random, alt_sats=alt_sats_random, inclination = inclination_random)

    random_percentdiff = 100*(original_constellation.area - random_constellation.area)/original_constellation.area
    # print("Random Satellites: ", random_sats)
    # print("Random Area: ", random_constellation.area)
    print("Random Percent Difference from Original: ", 100*(original_constellation.area - random_constellation.area)/original_constellation.area)
    # # print("Random Coordinates", coords_random)
    # functions.write_polygon(random_constellation, "random_constellation")

    lbfgs_percentdiff = 0
    return remaining_percentdiff, lbfgs_percentdiff, pso_percentdiff, random_percentdiff
    

# number of original satellites
n = 20
# lose this many satellites
missing = [0] 
num_missing = len(missing)
remaining_percentdiff = []
lbfgs_percentdiff = []
pso_percentdiff = []
random_percentdiff = []

for i in range(1, 31):
    r_percentdiff, l_percentdiff, p_percentdiff, rand_percentdiff = optimize(n, missing)
    remaining_percentdiff.append(r_percentdiff)
    lbfgs_percentdiff.append(l_percentdiff)
    pso_percentdiff.append(p_percentdiff)
    random_percentdiff.append(rand_percentdiff)

print("Average Percent Difference No Optimization",sum(remaining_percentdiff)/len(remaining_percentdiff))
# print("Average Percent Difference DiffEv + LBFGS-B",sum(lbfgs_percentdiff)/len(lbfgs_percentdiff))
print("Average Percent Difference PSO",sum(pso_percentdiff)/len(pso_percentdiff))
print("Average Percent Difference Random",sum(random_percentdiff)/len(random_percentdiff))

with open(f'avg_percent_diffs_{num_missing}.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Remaining'] + [sum(remaining_percentdiff)/len(remaining_percentdiff)])
    # writer.writerow(['LBFGS'] + [sum(lbfgs_percentdiff)/len(lbfgs_percentdiff)])
    writer.writerow(['PSO'] + [sum(pso_percentdiff)/len(pso_percentdiff)])
    writer.writerow(['Random'] + [sum(random_percentdiff)/len(random_percentdiff)])
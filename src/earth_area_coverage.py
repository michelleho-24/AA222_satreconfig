# from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import numpy as np
import math
import shapely
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point
import json
from shapely.geometry import shape
import shapely.geometry as sgeom
from shapely.ops import unary_union

# import functions 
#import brahe as bh
import cartopy.geodesic
from shapely.geometry import Polygon, MultiPolygon
import warnings
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D

def read_csv_file(file_path):
    lat_lon = []
    with open(file_path, 'r') as csvfile: 
        csvreader = csv.reader(csvfile)
        headers_skipped = False
        for row in csvreader: 
            if not headers_skipped: 
                headers_skipped = True
                continue
            lat_lon.append(row)
    return lat_lon

R_EARTH = 6278*1e3 # radius of earth in meters
def compute_earth_interior_angle(ele, alt):
    ele = ele * math.pi / 180.0
    rho = math.asin(R_EARTH/(R_EARTH + alt * 1e3))
    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele
    return lam #returns in radians 

ele = 10.0

def calculate_coverage_area(lat, lon, alt,ele): 
    angle = compute_earth_interior_angle(ele,alt)
    radius = angle*R_EARTH
    coverage_area = Point(lon, lat).buffer(radius)
    return coverage_area

# Input latitude and longitude matrix (dummy for now)
'''file_path = 'satellite_coordinates.csv'
lat_lon_ar = read_csv_file(file_path)
lat_lon = np.array(lat_lon_ar)
#lat_lon = np.array([[1,2,3,4], [5,6,7,8],[9,10,11,12], [13, 14, 15, 16]])

# Assume an input of latitude and longitude input as a matrix where columnes correspond to latitude and longitude for each satellite 
rows, columns = lat_lon.shape
num_points  = rows
num_sat = columns/2
heights = np.array([500,500,500])

constellation_coverage = MultiPolygon()
original_coverage = MultiPolygon()

# While running through every time step
constellation_coverage = []
for i in range(0,int(num_sat)-1):
    x = i*2
    satellite_coverage = [] # Initialize satellite coverage
    lat = lat_lon[:,x] # extract latitude 
    n = x+1
    lon = lat_lon[:,n] # extract longitude 
    alt = heights[i] # extract height
    for t in range(0,num_points): # iterate through all of the data points 
        coverage = calculate_coverage_area(lat[t], lon[t],alt, ele)
        satellite_coverage.append(coverage)

    satellite_overall = Polygon()
    for area in satellite_coverage:
        satellite_overall = satellite_overall.union(area)
    constellation_coverage.append(satellite_overall)

constellation = Polygon()
for sat in constellation_coverage: 
    constellation = constellation.union(sat)'''

def write_polygon(polygon, filename):
    geojson ={
        "type":"Feature",
        "properties":{},
        "geometry": polygon.__geo_interface__
    }
    with open(filename, "w") as f: 
        json.dump(geojson,f)

def read_polygon(filename):
    with open(filename, "r") as f: 
        geojson = json.load(f)

    polygon = shape(geojson['geometry'])
    return polygon

# area_covered = constellation.area
# write_polygon(constellation, "test.geojson")

# const = read_polygon("test.geojson")

def plot_coverage(coverage, filename): 
    fig, ax = plt.subplots(figsize=(10,8))
    ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    ax.stock_img()
    # for polygon in coverage.geoms: 
    ax.plot(*coverage.exterior.xy, color = 'red', transform=ccrs.PlateCarree())
    ax.set_global()
    ax.coastlines()
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Coverage Area')
    plt.grid(True)
    # plt.show()
    plt.savefig(filename, bbox_inches = 'tight', pad_inches = 0.1)
    plt.close()

#plot_coverage(const, "test_coverage.png")

def plot_eddy(sat_coords):
    # og_coords = [550, 550, 550, 550, 550, 550, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0]
    # LGBFS_coords = [553.69411486,550,548.15294257,553.69411486,550,30,60,90,120,150]
    # PSO_coords = [605.33456691,624.50644854,593.98916237,591.4515896 ,625.9709508,60.53297263 ,25.3114063 ,120.00667778 ,26.98097552 ,91.31959186]
    # random_coords = [633, 508, 584, 640, 566, 179, 77, 110, 93, 8]

    elevation_min = 10.0
    # Create the figure
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.stock_img()
    c = 'b' # Set the plot color
    # plot the optimized sats 
    for i in range(len(sat_coords)/2):
        lam = compute_earth_interior_angle(ele=elevation_min, alt=sat_coords[i])

        # Get the satellite position from above
        x_ecef = functions.tle.state_ecef(t) # Get the ECEF state one day later
        x_geod = bh.sECEFtoGEOD(x_ecef[0:3], use_degrees=True) # Need the array slice to get only the position
        lon, lat = x_geod[0], x_geod[1]

        # Plot Groundstation Location
        ax.plot(lon, lat, color=c, marker='o', markersize=3, transform=ccrs.Geodetic())

        # Get a bunch of points in a circle space at the the right angle offset from the sub-satellite point to draw the view cone
        circle_points = cartopy.geodesic.Geodesic().circle(lon=lon, lat=lat, radius=lam*bh.R_EARTH, n_samples=100, endpoint=False)
        geom = shapely.geometry.Polygon(circle_points)
        ax.add_geometries((geom,), crs=ccrs.Geodetic(), facecolor=c, alpha=0.5, edgecolor='none', linewidth=0)

        plt.show()

'''def generate_circle_points(lon, lat, radius_meters, n_samples): 
    geod = pyproj.Geod(ellps='WGS84')  # Define a geodesic with WGS84 ellipsoid
    azimuths = [360 * i / n_samples for i in range(n_samples)]  # Generate azimuths
    lon_array, lat_array, _ = geod.fwd(lon, lat, azimuths, [radius_meters] * n_samples)
    return lon_array, lat_array
'''
def plot_eddy_lat_lon(lat_lon, lat_lon_optimize, alt, alt_optimize, filename, c, b):
#def plot_eddy_lat_lon(lat_lon, alt, filename, c, b):
    elevation_min = 10.0
    fig = plt.figure(figsize=(10,5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.stock_img()
    rows, columns = lat_lon.shape
    R_EARTH = 6378*1e3
    rows2, columns2 = lat_lon_optimize.shape

    all_geoms_original = []
    all_geoms_optimized = []
    for i in range(0,int(columns)//2-1):
        a = alt[i]
        lam = compute_earth_interior_angle(elevation_min, a)

        for j in range(rows):
            index = 2*i
            longitude = float(lat_lon[j,index])
            latitude = float(lat_lon[j, index+1])
            #ax.plot(longitude, latitude, color = b, marker='o', markersize=3, transform=ccrs.Geodetic())
            radius = lam*R_EARTH
            circle_points = cartopy.geodesic.Geodesic().circle(longitude, latitude, radius, 100, endpoint=False)
            circle_points = circle_points.tolist()
            circle_points.append(circle_points[0])
            geom = shapely.geometry.Polygon(circle_points)
            all_geoms_original.append(geom)
            ax.add_geometries([geom], crs=ccrs.Geodetic(), color = b, alpha=0.05, edgecolor='none', linewidth=0, label = "Original")
    
    for i in range(0,int(columns2)//2-1):
        a = alt_optimize[i]
        lam = compute_earth_interior_angle(elevation_min, a)

        for j in range(rows2):
            index = 2*i
            longitude = float(lat_lon_optimize[j,index])
            latitude = float(lat_lon_optimize[j, index+1])
            #ax.plot(longitude, latitude, color = c, marker='o', markersize=3, transform=ccrs.Geodetic())
            radius = lam*R_EARTH
            circle_points = cartopy.geodesic.Geodesic().circle(longitude, latitude, radius, 100, endpoint=False)
            circle_points = circle_points.tolist()
            circle_points.append(circle_points[0])
            geom = shapely.geometry.Polygon(circle_points)
            if not geom.is_valid:
                geom = geom.buffer(0)
            all_geoms_optimized.append(geom)
            ax.add_geometries([geom], crs=ccrs.Geodetic(), color = c, alpha=0.05, edgecolor='none', linewidth=0, label = "Particle Swarm")
    
    #merged_original = unary_union([geom for geom in all_geoms_original if geom.is_valid])
    #merged_optimize = unary_union([geom for geom in all_geoms_optimized if geom.is_valid])
    #ax.add_geometries([merged_original], crs=ccrs.Geodetic(),color = b, alpha = 0.05, edgecolor = 'none')
    #ax.add_geometries([merged_optimize], crs=ccrs.Geodetic(), color = c, alpha=0.05, edgecolor = 'none')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title("Ground Coverage")
    legend_elements = [
        Line2D([0],[0],color = b, lw=2, label = 'Original Coverage'),
        Line2D([0],[0],color = c, lw=2, label = 'Random')]
    ax.legend(handles=legend_elements, loc='upper left')
    
    plt.savefig(filename, dpi = 1000)
    

APPLIED_FILTER_WARNINGS = False

def filter_cartopy_warnings():
    global APPLIED_FILTER_WARNINGS

    if not APPLIED_FILTER_WARNINGS:
        warnings.filterwarnings("ignore", message="Approximating coordinate system")
        APPLIED_FILTER_WARNINGS = True

filter_cartopy_warnings()

original_alt = [650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650,650]
LBFGS_alt = [650.0,650.0,650.0,650.0,650.0,650.0,650.0,650.0]
PSO_alt = [549.0450144888142,591.6218815400181,646.2505174251835,642.7817016433726,584.980256172985,608.209004108612,524.1578103806704,578.1331252812461]
random_alt = [610.5183086432758,595.2437592504358,508.54816400526215,554.8771172503688,627.5118269396709,634.908011254997,576.8588255799604,590.7650583485722]
diff_alt = [517.1017477608356,608.987416630821,596.1562925678791,633.7061670946465,636.2625085230242,566.1753211927348,650.0,558.2418921964143]
filepath_diff = "diffev_coordinates.csv"
filepath_original = "original_coords.csv"
filepath_LBFGS = 'lbfgs_coordinates.csv'
filepath_pso = 'pso_coordinates.csv'
filepath_random = 'random_coordinates.csv'


lat_lon_original = np.array(read_csv_file(filepath_original))
lat_lon_LBFGS = np.array(read_csv_file(filepath_LBFGS))
lat_lon_PSO = np.array(read_csv_file(filepath_pso))
lat_lon_random = np.array(read_csv_file(filepath_random))
lat_lon_diff = np.array(read_csv_file(filepath_diff))

#plot_eddy_lat_lon(lat_lon_original, original_alt, "plot_original.png", 'b', 'r')
#plot_eddy_lat_lon(lat_lon_original, lat_lon_PSO, original_alt, PSO_alt,"plot_pso.png", 'b', 'r')
# plot_eddy_lat_lon(lat_lon_original, lat_lon_random, original_alt, random_alt,"plot_random.png", 'b', 'r')
plot_eddy_lat_lon(lat_lon_original, lat_lon_diff, original_alt, diff_alt,"plot_diff.png", 'b', 'r')


 #buffered_geoms = [geom.buffer(0) for geom in all_geoms]
    #combined_geom = unary_union(all_geoms)
    #globe_poly = shapely.geometry.Polygon([(-180,-90), (180,-90),(180, 90),(-180,90),(-180,-90)])
    #buffered_geom = combined_geom.buffer(0)
    #inner_area = globe_poly.difference(combined_geom)
    #ax.add_geometries([inner_area], crs = ccrs.Geodetic(),color=c, alpha = 0.5, edgecolor = 'none', linewidth = 0)
'''combined_circle_points = np.concatenate(all_circle_points)
    combined_circle_geom = sgeom.MultiPoint(combined_circle_points).convex_hull
    path = Path(combined_circle_geom.exterior.coords)
    patch = PathPatch(path, facecolor=c, edgecolor = 'none', alpha = 0.5)
    ax.add_patch(patch)'''
    # plt.savefig(filename)
#plot_eddy([550, 550, 550, 550, 550, 550, 0.0, 30.0, 60.0, 90.0, 120.0, 150.0]
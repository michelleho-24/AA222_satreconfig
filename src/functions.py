
import warnings
import math
import numpy as np
import shapely
import cartopy.geodesic
import brahe as bh
import brahe.data_models as bdm
import brahe.access.access as ba
from shapely.geometry import Polygon, MultiPolygon
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import csv
from matplotlib.patches import Circle
from shapely.geometry import Point
import json
from shapely.geometry import shape

# Here we can download the latest Earth orientation data and load it.

# Uncomment this line ONCE the data has been downloaded. Recomment it once it has been downloaded.
# bh.utils.download_iers_bulletin_ab()
# Load the latest Earth Orientation Data
# bh.EOP.load('iau2000A_finals_ab.txt')

def filter_cartopy_warnings():
    global APPLIED_FILTER_WARNINGS

    if not APPLIED_FILTER_WARNINGS:
        warnings.filterwarnings("ignore", message="Approximating coordinate system")
        APPLIED_FILTER_WARNINGS = True
        
def sECEFtoECIpatch(epc, x):


    """Transforms an Earth fixed state into an Inertial state

    The transformation is accomplished using the IAU 2006/2000A, CIO-based
    theory using classical angles. The method as described in section 5.5 of
    the SOFA C transformation cookbook.

    Args:
        epc (Epoch): Epoch of transformation
        x (np.ndarray): Earth-fixed state (position, velocity) [*m*; *m/s*]

    Returns:
        x_ecef (np.ndarray): Inertial state (position, velocity)
    """

    # Ensure input is array-like
    x = np.asarray(x)

    # Set state variable size
    dim_x = len(x)
    x_eci = np.zeros((dim_x,))

    # Extract State Components
    r_ecef = x[0:3]

    if dim_x >= 6:
        v_ecef = x[3:6]

    # Compute Sequential Transformation Matrices
    rot = bh.earth_rotation(epc)

    # Create Earth's Angular Rotation Vector
    omega_vec = np.array([0, 0, bh.constants.OMEGA_EARTH]) # Neglect LOD effect

    # Calculate ECEF State
    x_eci[0:3] = ( rot ).T @ r_ecef
    # x_eci[0:3] = (pm @ rot @ bpn).T @ r_ecef

    if dim_x >= 6:
        x_eci[3:6] = (rot ).T @ (v_ecef + bh.utils.fcross(omega_vec, r_ecef))

    return x_eci

def get_tle(epc0, alt, ecc, inc, raan, argp, M, ndt2=0.0, nddt6=0.0, bstar=0.0, norad_id=99999):
    '''Get a TLE object from the given orbital elements

    Args:
    - epc0 (Epoch): Epoch of the orbital elements / state
    - alt (float): Altitude of the orbit [km]
    - ecc (float): Eccentricity of the orbit
    - inc (float): Inclination of the orbit [deg]
    - raan (float): Right Ascension of the Ascending Node [deg]
    - argp (float): Argument of Perigee [deg]
    - M (float): Mean Anomaly [deg]

    Returns:
    - tle (TLE): TLE object for the given orbital elements
    '''
    alt *= 1e3 # Convert to meters

    # Get semi-major axis
    sma = bh.R_EARTH + alt

    # Get mean motion
    n = bh.mean_motion(sma)/(2*np.pi)*86400

    tle_string = bh.tle_string_from_elements(epc0, np.array([n, ecc, inc, raan, argp, M, ndt2, nddt6, bstar]), norad_id)
    tle = bh.TLE(*tle_string)
    return tle

def current_epoch(dt):
    epc = bh.Epoch(2024, 5, 20, 0, 0, 0) + dt
    return epc

def get_tle_mult_sats_initial(height_incl):
    n = len(height_incl) / 2
    num = int(n)
    epc0 = current_epoch(0)
    ecc = 0.001
    raan = 45
    argp = 90
    M = 45
    norad_id = 99999
    
    tle_list = []
    for i in range(num):
        tle_list.append(get_tle(epc0, height_incl[i], ecc, height_incl[num + i], raan, argp, M, norad_id=norad_id))
    
    return tle_list

def get_coords_mult_sats(height_incl, dt):
    n = len(height_incl) / 2
    num = int(n)
    # epc = current_epoch(dt)
    # epc0 = current_epoch(0)
    epc0 = bh.Epoch(2022, 5, 20, 0, 0, 0)
    ecc = 0.001
    raan = 45
    argp = 90
    M = 45
    norad_id = 99999
    t = epc0 + dt
    
    # tle_list = []
    coord_list = []
    for i in range(num):
        tle_current_sat = get_tle(epc0, height_incl[i], ecc, height_incl[num + i], raan, argp, M, norad_id=norad_id)
        # tle_list.append(tle_current_sat)
        x_ecef = tle_current_sat.state_ecef(t)
        x_geod = bh.sECEFtoGEOD(x_ecef[0:3], use_degrees=True) 
        coord_list.append(x_geod[0:2])
    return coord_list

def get_alt_mult_sats(height_incl, dt):
    n = len(height_incl) / 2
    num = int(n)
    # epc = current_epoch(dt)
    # epc0 = current_epoch(0)
    epc0 = bh.Epoch(2022, 5, 20, 0, 0, 0)
    ecc = 0.001
    raan = 45
    argp = 90
    M = 45
    norad_id = 99999
    t = epc0 + dt
    
    # tle_list = []
    alt_list = []
    for i in range(num):
        tle_current_sat = get_tle(epc0, height_incl[i], ecc, height_incl[num + i], raan, argp, M, norad_id=norad_id)
        # tle_list.append(tle_current_sat)
        x_ecef = tle_current_sat.state_ecef(t)
        x_geod = bh.sECEFtoGEOD(x_ecef[0:3], use_degrees=True) 
        alt_list.append(x_geod[2])
    return alt_list

def compute_earth_interior_angle(ele, alt):
    # elevation in degrees and altitude in meters 
    ele = ele * math.pi / 180.0
    # rho = math.asin(bh.R_EARTH/(bh.R_EARTH + (alt*1000)))
    # Calculate the argument for math.asin
    arg = bh.R_EARTH/(bh.R_EARTH + (alt*1000))

    # Check if the argument is within the valid range for math.asin
    if arg < -1:
        arg = -1
    elif arg > 1:
        arg = 1

    # Now it's safe to call math.asin
    rho = math.asin(arg)
    eta = math.asin(math.cos(ele)*math.sin(rho))
    lam = math.pi/2.0 - eta - ele
    return lam #returns in radians 

def calculate_coverage_area(lat, lon, alt, ele): 
    angle = compute_earth_interior_angle(ele,alt)
    radius = angle*bh.R_EARTH
    coverage_area = Point(lon, lat).buffer(radius)
    return coverage_area

def compute_area(coords, alt_sats, inclination):
    alt = np.copy(alt_sats)
    # input is an array 
    lat_lon = np.array(coords)
    # print((lat_lon).size)
    # print(lat_lon.shape)
    rows, columns = lat_lon.shape
    num_points = rows
    # num_sat = columns/2
    num_sat = len(alt_sats)
    # print(coords)
    # print(num_sat)
    # print(num_points)
    ele = 10.0

    # constellation_coverage = MultiPolygon()

    # While running through every time step
    constellation_coverage = []
    for i in range(0,int(num_sat)-1):
        x = i*2
        satellite_coverage = [] # Initialize satellite coverage
        # print(lat_lon.size)
        lon = lat_lon[:,x] # extract latitude 
        n = x+1
        lat = lat_lon[:,n] # extract longitude 
        altitude = alt[i] # extract height
        for t in range(0,num_points): # iterate through all of the data points 
            coverage = calculate_coverage_area(lat[t], lon[t],altitude,ele)
            satellite_coverage.append(coverage)

        satellite_overall = Polygon()
        for area in satellite_coverage:
            satellite_overall = satellite_overall.union(area)
        constellation_coverage.append(satellite_overall)

    constellation = Polygon()
    for sat in constellation_coverage: 
        constellation = constellation.union(sat)
        
    #area_covered = constellation.area
    return constellation

def calculate_lat_lon(sats_initial):
        coords_initial = get_coords_mult_sats(sats_initial, 0) # find initial coordinates of each satellite 
        alt_initial = get_alt_mult_sats(sats_initial, 0)       # find initial altitude of each satellite 
        num_sats = int(len(sats_initial) / 2) # define number of satellites in problem 
        coord_info = {}                       # define emtpy dictionary to carry through coordinate information 
        alt_info = {}                         # define emtpy dictionary to carry through altitude information 
        for i in range(num_sats):
            coord_info[f'coords_sat_{i}'] = []
            alt_info[f'alt_sat_{i}'] = []

        end_time = 60 * 60      # propagate by 60 minutes 
        dt = 60                 # extract coordinate data every 60 seconds 

        height_incl = sats_initial
        coords = []             # empty list for coordinate information 
        alt = []                # empty list for altitude information 
        for time in range(0, end_time + 60, dt):
            coords_current = get_coords_mult_sats(height_incl, time)  # get coord info for current time step
            alt_current = get_alt_mult_sats(height_incl, time)        # get alt info for current time step
            
            coords_row = []   # define empty array for csv writing 

            for i in range(num_sats):
                coord_info[f'coords_sat_{i}'].append(coords_current[i])
                alt_info[f'alt_sat_{i}'].append(alt_current[i])
                coords_row.extend(coords_current[i])            # "flatten" lat/long for each satellite into the same row  
            
            coords.append(coords_row)  
        return coords

def plot_coverage(coverage, filename): 
        # coverage is a polygon 
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

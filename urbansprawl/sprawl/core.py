###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import pandas as pd
import geopandas as gpd
import numpy as np
from osmnx import log
from shapely.geometry import Point

from ..osm.core import get_route_graph, get_processed_osm_data
from .landusemix import compute_grid_landusemix
from .accessibility import compute_grid_accessibility
from .dispersion import compute_grid_dispersion

def get_indices_grid(df_osm_built, df_osm_building_parts, df_osm_pois, step=100):
	"""
	Creates an input geodataframe with points sampled in a regular grid

	Parameters
	----------
	df_osm_built : geopandas.GeoDataFrame
		OSM processed buildings
	df_osm_building_parts : geopandas.GeoDataFrame
		OSM processed building parts
	df_osm_pois : geopandas.GeoDataFrame
		OSM processed points of interest
	step : int
		step to sample the regular grid in meters

	Returns
	----------
	geopandas.GeoDataFrame
		regular grid
	"""
	# Get bounding box
	west, south, east, north = pd.concat( [ df_osm_built, df_osm_building_parts, df_osm_pois ], sort=False ).total_bounds
	return get_indices_grid_from_bbox([west, south, east, north],
                                          step,
                                          df_osm_built.crs)


def get_indices_grid_from_bbox(bounding_box, step=100, crs={"init":"epsg:4326"}):
	"""
	Creates an input geodataframe with points sampled in a regular grid

	Parameters
	----------
        bounding_box : list
            Geographical coordinates in which one has to build the grid
	step : int
	    Step to sample the regular grid in meters

	Returns
	----------
	geopandas.GeoDataFrame
		regular grid
	"""
	# Get bounding box
	west, south, east, north = bounding_box
	# Create indices
	df_indices = gpd.GeoDataFrame( [ Point(i,j) for i in np.arange(west, east, step) for j in np.arange(south, north, step) ], columns=["geometry"] )
	# Set projection
	df_indices.crs = crs
	return df_indices


def process_spatial_indices(city_ref=None, region_args={"polygon":None, "place":None, "which_result":1, "point":None, "address":None, "distance":None, "north":None, "south":None, "east":None, "west":None},
			grid_step = 100,
			process_osm_args = {"retrieve_graph":True, "default_height":3, "meters_per_level":3, "associate_landuses_m2":True, "minimum_m2_building_area":9, "date":None},
			dispersion_args = {'radius_search': 750, 'use_median': False, 'K_nearest': 50},
			landusemix_args = {'walkable_distance': 600, 'compute_activity_types_kde': True, 'weighted_kde': True, 'pois_weight': 9, 'log_weighted': True},
			accessibility_args = {'fixed_distance': True, 'fixed_activities': False, 'max_edge_length': 200, 'max_node_distance': 250,
				'fixed_distance_max_travel_distance': 2000, 'fixed_distance_max_num_activities': 250, 'fixed_activities_min_number': 20},
			indices_computation = {"dispersion":True, "landusemix":True, "accessibility":True} ):
	"""
	Process sprawling indices for an input region of interest
	1) OSM data is retrieved and processed.
		If the city name has already been processed, locally stored data will be loaded
	2) A regular grid is created where indices will be calculated
	3) Sprawling indices are calculated and returned

	Parameters
	----------
	city_ref : str
		Name of input city / region
	grid_step : int
		step to sample the regular grid in meters
	region_args : dict
		contains the information to retrieve the region of interest as the following:
			polygon : shapely Polygon or MultiPolygon
				geographic shape to fetch the land use footprints within
			place : string or dict
				query string or structured query dict to geocode/download
			which_result : int
				result number to retrieve from geocode/download when using query string
			point : tuple
				the (lat, lon) central point around which to construct the region
			address : string
				the address to geocode and use as the central point around which to construct the region
			distance : int
				retain only those nodes within this many meters of the center of the region
			north : float
				northern latitude of bounding box
			south : float
				southern latitude of bounding box
			east : float
				eastern longitude of bounding box
			west : float
				western longitude of bounding box
	process_osm_args : dict
		additional arguments to drive the OSM data extraction process:
			retrieve_graph : boolean
				that determines if the street network for input city has to be retrieved and stored
			default_height : float
				height of buildings under missing data
			meters_per_level : float
				buildings number of levels assumed under missing data
			associate_landuses_m2 : boolean
				compute the total square meter for each land use
			minimum_m2_building_area : float
				minimum area to be considered a building (otherwise filtered)
			date : datetime.datetime
				query the database at a certain timestamp
	dispersion_args : dict
		arguments to drive the dispersion indices calculation
			radius_search: int
				circle radius to consider the dispersion calculation at a local point
			use_median : bool
				denotes whether the median or mean should be used to calculate the indices
			K_nearest : int
				number of neighboring buildings to consider in evaluation
	landusemix_args : dict
		arguments to drive the land use mix indices calculation
			walkable_distance : int
				the bandwidth assumption for Kernel Density Estimation calculations (meters)
			compute_activity_types_kde : bool
				determines if the densities for each activity type should be computed
			weighted_kde : bool
				use Weighted Kernel Density Estimation or classic version
			pois_weight : int
				Points of interest weight equivalence with buildings (squared meter)
			log_weighted : bool
				apply natural logarithmic function to surface weights
	accessibility_args : dict
		arguments to drive the accessibility indices calculation
			fixed_distance : bool
				denotes the cumulative opportunities access to activity land uses given a fixed maximum distance to travel
			fixed_activities : bool
				represents the distance needed to travel in order to reach a certain number of activity land uses
			max_edge_length: int
				maximum length, in meters, to tolerate an edge in a graph (otherwise, divide edge)
			max_node_distance: int
				maximum distance tolerated from input point to closest graph node in order to calculate accessibility values
			fixed_distance_max_travel_distance: int
				(fixed distance) maximum distance tolerated (cut&branch) when searching for the activities
			fixed_distance_max_num_activities: int
				(fixed distance) cut iteration if the number of activities exceeds a threshold
			fixed_activities_min_number: int
				(fixed activities) minimum number of activities required
	indices_computation : dict
		determines what sprawling indices should be computed

	Returns
	----------
	gpd.GeoDataFrame
		returns the regular grid with the indicated sprawling indices
	"""
	try:
		# Process OSM data
		df_osm_built, df_osm_building_parts, df_osm_pois = get_processed_osm_data(city_ref=city_ref, region_args=region_args, kwargs=process_osm_args)
		# Get route graph
		G = get_route_graph(city_ref)

		if (not ( indices_computation.get("accessibility") or indices_computation.get("landusemix") or indices_computation.get("dispersion") ) ):
			log("Not computing any spatial indices")
			return None

		# Get indices grid
		df_indices = get_indices_grid(df_osm_built, df_osm_building_parts, df_osm_pois, grid_step)

		# Compute sprawling indices
		if (indices_computation.get("accessibility")):
			compute_grid_accessibility(df_indices, G, df_osm_built, df_osm_pois, accessibility_args)
		if (indices_computation.get("landusemix")):
			compute_grid_landusemix(df_indices, df_osm_built, df_osm_pois, landusemix_args)
		if (indices_computation.get("dispersion")):
			compute_grid_dispersion(df_indices, df_osm_built, dispersion_args)

		return df_indices

	except Exception as e:
		log("Could not compute the spatial indices. An exception occurred: " + str(e))
		return None

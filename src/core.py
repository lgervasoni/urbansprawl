###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import os.path
import time
import osmnx as ox
import pandas as pd
import geopandas as gpd
import numpy as np

from osmnx.utils import log

from .parameters import get_dataframes_filenames, files

from .osm.osm_overpass import create_landuse_gdf, create_pois_gdf, create_building_parts_gdf, create_buildings_gdf_from_input, retrieve_route_graph
from .osm.osm_tags import columns_osm_tag, height_tags, building_parts_to_filter
from .osm.osm_data import compute_landuse_inference, classify_tag, classify_activity_category
from .osm.osm_surface import compute_landuses_m2
from .osm.osm_utils import load_geodataframe, store_geodataframe, associate_structures, sanity_check_height_tags

from .sprawl.landusemix import compute_grid_landusemix
from .sprawl.accessibility import compute_grid_accessibility
from .sprawl.dispersion import compute_grid_dispersion

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
	from shapely.geometry import Point
	# Get bounding box
	west, south, east, north = pd.concat( [ df_osm_built, df_osm_building_parts, df_osm_pois ] ).total_bounds
	# Create indices
	df_indices = gpd.GeoDataFrame( [ Point(i,j) for i in np.arange(west, east, step) for j in np.arange(south, north, step) ], columns=["geometry"] )
	# Set projection
	df_indices.crs = df_osm_built.crs
	return df_indices


def get_route_graph(city_ref, date="", polygon=None, north=None, south=None, east=None, west=None, force_crs=None):
	""" 
	Wrapper to retrieve city's street network
	Loads the data if stored locally
	Otherwise, it retrieves the graph from OpenStreetMap using the osmnx package
	Input polygon or bounding box coordinates determine the region of interest

	Parameters
	----------
	city_ref : string
		name of the city
	polygon : shapely.Polygon
		polygon shape of input city
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	force_crs : dict
		graph will be projected to input crs

	Returns
	----------
	networkx.multidigraph
		projected graph
	"""
	return retrieve_route_graph(city_ref, date, polygon, north, south, east, west, force_crs)


def get_processed_osm_data(city_ref_file=None, region_args={"polygon":None, "place":None, "which_result":1, "point":None, "address":None, "distance":None, "north":None, "south":None, "east":None, "west":None},
			kwargs={"retrieve_graph":True, "default_height":3, "meters_per_level":3, "associate_landuses_m2":True, "mixed_building_first_floor_activity":True, "minimum_m2_building_area":9, "date":None}):
	"""
	Retrieves buildings, building parts, and Points of Interest associated with a residential/activity land use from OpenStreetMap data for input city
	If a name for input city is given, the data will be loaded (if it was previously stored)
	If no stored files exist, it will query and process the data and store it under the city name
	Queries data for input region (polygon, place, point/address and distance around, or bounding box coordinates)
	Additional arguments will drive the overall process

	Parameters
	----------
	city_ref_file : str
		Name of input city / region
	region_args : dict
		contains the information to retrieve the region of interest as the following:
			polygon : shapely Polygon or MultiPolygon
				geographic shape to fetch the landuse footprints within
			place : string or dict
				query string or structured query dict to geocode/download
			which_result : int
				result number to retrieve from geocode/download when using query string 
			point : tuple
				the (lat, lon) central point around which to construct the graph
			address : string
				the address to geocode and use as the central point around which to construct the graph
			distance : int
				retain only those nodes within this many meters of the center of the graph
			north : float
				northern latitude of bounding box
			south : float
				southern latitude of bounding box
			east : float
				eastern longitude of bounding box
			west : float
				western longitude of bounding box
	kwargs : dict
		additional arguments to drive the process:
			retrieve_graph : boolean
				that determines if the street network for input city has to be retrieved and stored
			default_height : float
				height of buildings under missing data
			meters_per_level : float
				buildings number of levels assumed under missing data
			associate_landuses_m2 : boolean
				compute the total square meter for each land use
			mixed_building_first_floor_activity : Boolean
				if True: Associates building's first floor to activity uses and the rest to residential uses
				if False: Associates half of the building's area to each land use (Activity and Residential)
			minimum_m2_building_area : float
				minimum area to be considered a building (otherwise filtered)
			date : datetime.datetime
				query the database at a certain timestamp

	Returns
	----------
	[ gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame ]
		returns the output geo dataframe containing all buildings, building parts, and points associated to a residential or activity land usage
	
	"""
	log("OSM data requested for city: " + str(city_ref_file) )

	if (city_ref_file):
		geo_poly_file, geo_poly_parts_file, geo_point_file = get_dataframes_filenames(city_ref_file)

	##########################
	### Stored file ?
	##########################
	if ( (city_ref_file) and ( os.path.isfile(geo_poly_file) ) ): # File exists
		log("Found stored files for city "+city_ref_file)
		# Load local GeoDataFrames
		return load_geodataframe(geo_poly_file), load_geodataframe(geo_poly_parts_file), load_geodataframe(geo_point_file)

	# Get keyword arguments for input region of interest
	polygon, place, which_result, point, address, distance, north, south, east, west = region_args.get("polygon"), region_args.get("place"), region_args.get("which_result"), region_args.get("point"), region_args.get("address"), region_args.get("distance"), region_args.get("north"), region_args.get("south"), region_args.get("east"), region_args.get("west")

	### Valid input?
	if not( any( [not (polygon is None), place, point, address, north, south, east, west] ) ):
		log("Error: Must provide at least one type of input")
		return None, None, None

	if ( kwargs.get("date") ): # Non-null date
		date_ = kwargs.get("date").strftime("%Y-%m-%dT%H:%M:%SZ")
		log("Requesting OSM database at timestamp: "+date_)
		# e.g.: [date:"2004-05-06T00:00:00Z"]
		date_query = '[date:"'+date_+'"]'
	else:
		date_query = ""

	##########################
	### Overpass query: Buildings
	##########################
	# Query and update bounding box / polygon
	df_osm_built, polygon, north, south, east, west = create_buildings_gdf_from_input(date=date_query, polygon=polygon, place=place, which_result=which_result, point=point, address=address, distance=distance, north=north, south=south, east=east, west=west)
	df_osm_built["osm_id"] = df_osm_built.index
	df_osm_built.reset_index(drop=True, inplace=True)
	##########################
	### Overpass query: Land use polygons. Aid to perform buildings land use inference
	##########################
	df_osm_lu = create_landuse_gdf(date=date_query, polygon=polygon, north=north, south=south, east=east, west=west)
	df_osm_lu["osm_id"] = df_osm_lu.index
	# Drop useless columns
	columns_of_interest = ["osm_id", "geometry", "landuse"]
	df_osm_lu.drop( [ col for col in list( df_osm_lu.columns ) if not col in columns_of_interest ], axis=1, inplace=True )
	df_osm_lu.reset_index(drop=True, inplace=True)
	##########################
	### Overpass query: POIs
	##########################
	df_osm_pois = create_pois_gdf(date=date_query, polygon=polygon, north=north, south=south, east=east, west=west)
	df_osm_pois["osm_id"] = df_osm_pois.index
	df_osm_pois.reset_index(drop=True, inplace=True)
	##########
	### Overpass query: Building parts. Allow to calculate the real amount of M^2 for each building
	##########
	df_osm_building_parts = create_building_parts_gdf(date=date_query, polygon=polygon, north=north, south=south, east=east, west=west)	
	# Filter rows not needed (roof, etc) and not building already exists in buildings extract
	if ("building" in df_osm_building_parts.columns):
		df_osm_building_parts = df_osm_building_parts[ (~ df_osm_building_parts["building:part"].isin(building_parts_to_filter) ) & (~ df_osm_building_parts["building:part"].isnull() ) & (df_osm_building_parts["building"].isnull()) ]
	else:
		df_osm_building_parts = df_osm_building_parts[ (~ df_osm_building_parts["building:part"].isin(building_parts_to_filter) ) & (~ df_osm_building_parts["building:part"].isnull() ) ]
	df_osm_building_parts["osm_id"] = df_osm_building_parts.index
	df_osm_building_parts.reset_index(drop=True, inplace=True)
	
	log("Done: OSM data request")

	####################################################
	### Sanity check of height tags
	####################################################
	sanity_check_height_tags(df_osm_built)
	sanity_check_height_tags(df_osm_building_parts)

	def remove_nan_dict(x): # Remove entries with nan values
		return { k:v for k, v in x.items() if pd.notnull(v) }

	df_osm_built['height_tags'] = df_osm_built[ [ c for c in height_tags if c in df_osm_built.columns ] ].apply(lambda x: remove_nan_dict(x.to_dict() ), axis=1)
	df_osm_building_parts['height_tags'] = df_osm_building_parts[ [ c for c in height_tags if c in df_osm_building_parts.columns ] ].apply(lambda x: remove_nan_dict(x.to_dict() ), axis=1)

	###########
	### Remove columns which do not provide valuable information
	###########	
	columns_of_interest = columns_osm_tag + ["osm_id", "geometry", "height_tags"]
	df_osm_built.drop( [ col for col in list( df_osm_built.columns ) if not col in columns_of_interest ], axis=1, inplace=True )
	df_osm_building_parts.drop( [ col for col in list( df_osm_building_parts.columns ) if not col in columns_of_interest ], axis=1, inplace=True)

	columns_of_interest = columns_osm_tag + ["osm_id", "geometry"]
	df_osm_pois.drop( [ col for col in list( df_osm_pois.columns ) if not col in columns_of_interest ], axis=1, inplace=True )	

	log('Done: Sanity check + drop unnecesary columns')	

	###########
	### Classification
	###########
	df_osm_built['classification'], df_osm_built['key_value'] = list( zip(*df_osm_built.apply( classify_tag, axis=1) ) )
	df_osm_pois['classification'], df_osm_pois['key_value'] = list( zip(*df_osm_pois.apply( classify_tag, axis=1) ) )
	df_osm_building_parts['classification'], df_osm_building_parts['key_value'] = list( zip(*df_osm_building_parts.apply( classify_tag, axis=1) ) )

	# Remove unnecesary buildings
	df_osm_built.drop( df_osm_built[ df_osm_built.classification.isnull() ].index, inplace=True )
	df_osm_built.reset_index(inplace=True, drop=True)
	# Remove unnecesary POIs
	df_osm_pois.drop( df_osm_pois[ df_osm_pois.classification.isin(["infer","other"]) | df_osm_pois.classification.isnull() ].index, inplace=True )
	df_osm_pois.reset_index(inplace=True, drop=True)
	# Building parts will acquire its containing building land use if it is not available
	df_osm_building_parts.loc[ df_osm_building_parts.classification.isin(["infer","other"]), "classification" ] = None

	log('Done: OSM tags classification')
	
	###########
	### Remove already used tags
	###########
	df_osm_built.drop( [ c for c in columns_osm_tag if c in df_osm_built.columns ], axis=1, inplace=True )
	df_osm_pois.drop( [ c for c in columns_osm_tag if c in df_osm_pois.columns ], axis=1, inplace=True )
	df_osm_building_parts.drop( [ c for c in columns_osm_tag if c in df_osm_building_parts.columns ], axis=1, inplace=True)

	###########
	### Project, drop small buildings and reset indices
	###########
	### Project to UTM coordinates within the same zone
	df_osm_built = ox.project_gdf(df_osm_built)
	df_osm_lu = ox.project_gdf(df_osm_lu, to_crs=df_osm_built.crs)
	df_osm_pois = ox.project_gdf(df_osm_pois, to_crs=df_osm_built.crs)
	df_osm_building_parts = ox.project_gdf(df_osm_building_parts, to_crs=df_osm_built.crs)

	# Drop buildings with an area lower than a threshold
	df_osm_built.drop( df_osm_built[ df_osm_built.geometry.area < kwargs["minimum_m2_building_area"] ].index, inplace=True )

	log('Done: Geometries reprojection')

	####################################################
	### Infer buildings land use (under uncertainty)
	####################################################
	compute_landuse_inference(df_osm_built, df_osm_lu)
	# Free space
	del df_osm_lu

	assert( len( df_osm_built[df_osm_built.key_value =={"inferred":"other"} ] ) == 0 )
	assert( len( df_osm_built[df_osm_built.classification.isnull()] ) == 0 )
	assert( len( df_osm_pois[df_osm_pois.classification.isnull()] ) == 0 )

	log('Done: Land use deduction')

	####################################################
	### Associate for each building, its containing building parts and Points of interest
	####################################################
	associate_structures(df_osm_built, df_osm_building_parts, operation='contains', column='containing_parts')
	associate_structures(df_osm_built, df_osm_pois, operation='intersects', column='containing_poi')

	# Classify activity types
	df_osm_built['activity_category'] = df_osm_built.apply(lambda x: classify_activity_category(x.key_value), axis=1)
	df_osm_pois['activity_category'] = df_osm_pois.apply(lambda x: classify_activity_category(x.key_value), axis=1)
	df_osm_building_parts['activity_category'] = df_osm_building_parts.apply(lambda x: classify_activity_category(x.key_value), axis=1)

	log('Done: Building parts association + activity categorization')

	####################################################
	### Associate effective number of levels, and measure the surface dedicated to each land use per building
	####################################################
	if (kwargs["associate_landuses_m2"]):
		default_height = kwargs["default_height"]
		meters_per_level = kwargs["meters_per_level"]
		mixed_building_first_floor_activity = kwargs["mixed_building_first_floor_activity"]
		compute_landuses_m2(df_osm_built, df_osm_building_parts, df_osm_pois, default_height=default_height, meters_per_level=meters_per_level, mixed_building_first_floor_activity=mixed_building_first_floor_activity)

	df_osm_built.loc[ df_osm_built.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
	df_osm_pois.loc[ df_osm_pois.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
	df_osm_building_parts.loc[ df_osm_building_parts.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan

	if (kwargs["associate_landuses_m2"]):
		# Set the composed classification given, for each building, its containing Points of Interest and building parts classification
		df_osm_built.loc[ df_osm_built.apply(lambda x: x.landuses_m2["activity"]>0 and x.landuses_m2["residential"]>0, axis=1 ), "classification" ] = "mixed"

	log('Done: Land uses surface association')

	##########################
	### Overpass query: Street network graph
	##########################
	if (kwargs["retrieve_graph"]): # Save graph for input city shape
		get_route_graph(city_ref_file, date=date_query, polygon=polygon, north=north, south=south, east=east, west=west, force_crs=df_osm_built.crs)

	log('Done: Street network graph retrieval')

	##########################
	### Store file ?
	##########################
	if ( city_ref_file ): # File exists
		# Save GeoDataFrames
		store_geodataframe(df_osm_built, geo_poly_file)
		store_geodataframe(df_osm_building_parts, geo_poly_parts_file)
		store_geodataframe(df_osm_pois, geo_point_file)
		log("Storing file for city: "+city_ref_file)

	return df_osm_built, df_osm_building_parts, df_osm_pois


def process_spatial_indices(city_ref_file=None, region_args={"polygon":None, "place":None, "which_result":1, "point":None, "address":None, "distance":None, "north":None, "south":None, "east":None, "west":None},
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
	city_ref_file : str
		Name of input city / region
	grid_step : int
		step to sample the regular grid in meters
	region_args : dict
		contains the information to retrieve the region of interest as the following:
			polygon : shapely Polygon or MultiPolygon
				geographic shape to fetch the landuse footprints within
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
				(fixed distance) cut iteration if the number of activites exceeds a threshold
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
		df_osm_built, df_osm_building_parts, df_osm_pois = get_processed_osm_data(city_ref_file=city_ref_file, region_args=region_args, kwargs=process_osm_args)
		# Get route graph
		G = get_route_graph(city_ref_file)

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
		log("Could not compute the spatial indices. An exception occured: " + str(e))
		return None

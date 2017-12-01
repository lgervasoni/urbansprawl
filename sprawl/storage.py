###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

from os.path import dirname, join, exists, isdir
from os import makedirs, remove, listdir
from six import string_types
import pandas as pd
import osmnx as ox
import time

import logging as lg
from .utils import _create_mesh_grid, _create_grid, log
from .osm_shapefile_parse import process_shapefiles

from .dispersion import compute_grid_dispersion, _apply_polygon_closest_distance_neighbor
from .accessibility import compute_grid_accessibility
from .landusemix import compute_grid_landusemix, kde_activity_residential, kde_activity_types

from .parameters import storage_folder

##########################
# Avoid PerformanceWarning when storing Polygon type to HDF5
import warnings
from pandas.io.pytables import PerformanceWarning
warnings.simplefilter(action = "ignore", category = PerformanceWarning)
##########################

##############################################################
### Parameters
##############################################################
# CONSTANTS FOR THE STORAGE-RELATED DEFINITIONS
FILE_EXTENSIONS = '.h5'

# Over write HDF5 data?
DATA_HDF5_OVERWITE = False

### KEYS
# osm_data: Data for building polygons; bbox: Bounding box
KEY_OSM_DF = ['osm_data','bbox']
# Kernel density estimations for 'activity' and 'residential' land usages
KEY_KDES = ['kde_activity_','kde_residential_'] # + step

# Type classification for activities land usages
from .osm_tags import ACTIVITIES_CLASSIFICATION
KEY_KDES_ACTIVITIES = [ 'kde_classification_activity_'+activ.replace('/','_') for activ in ACTIVITIES_CLASSIFICATION ]

# Indices (to index: KEY_INDEX + step)
KEY_DISPERSION = 'grid_dispersion_' # + step
KEY_ACCESSIBILITY = 'grid_accessibility_' # + step
KEY_LANDUSEMIX = 'grid_landusemix_' # + step
# Mesh grid keys (to index: KEY_MESH_GRID + step)
KEY_MESH_GRID = 'xy_grid_' # + step

##############################################################
##############################################################

### CUSTOM EXCEPTION
class NoMethodProvidedException(Exception):
	pass
##############################################################

# LOCAL STORAGE UTILS

def _generate_file_path(city_ref):
	"""
	Generates a file path for local storage data of `city_ref`

	Parameters
	----------
	city_ref : string
		name of the city
	
	Returns
	----------
	string
		path to storage file
	"""
	return join(storage_folder, city_ref + FILE_EXTENSIONS)


def _get_local_data(city_ref, hdfs_key):
	"""
	Get local data stored as `hdfs_key` in city_ref's file

	Parameters
	----------
	city_ref : string
		name of the city
	hdfs_key : string
		key name for data
	
	Returns
	----------
	pandas.DataFrame
		loaded data for input key
	"""
	data = None
	try:
		# Read mode
		storage = pd.HDFStore(_generate_file_path(city_ref), 'r')
		if (hdfs_key in storage): # Read `hdfs_key` objects
			data = pd.read_hdf(storage,hdfs_key)
			log( "Found `"+hdfs_key+"` for `"+city_ref+"` stored locally" )
		storage.close()
		return data
	except (IOError, KeyError) as e: # Local data does not exist for `city_ref`
		return None

def _update_local_data(city_ref, result, hdfs_key, over_write = False):
	"""
	Store data for input data frame under name `hdfs_key`

	Parameters
	----------
	city_ref : string
		name of the city
	result : pandas.DataFrame
		data
	hdfs_key : string
		key name for data
	over_write : boolean
		optionally over-write existing data
	
	Returns
	----------

	"""
	if (not isdir(storage_folder)): makedirs( storage_folder ) # Create cache folder
	if (result is None): return # Empty data?
	try:
		# Write mode: append
		storage = pd.HDFStore(_generate_file_path(city_ref), 'a')
		# Write data?
		if ( (hdfs_key not in storage) or (over_write) ):
			result.to_hdf(_generate_file_path(city_ref),hdfs_key)
			log( "Data `"+hdfs_key+"` for `"+city_ref+"` has been stored" )
		storage.close()
	except Exception as e:
		log( "Error while saving data. Local data not updated. "+str(e) , level=lg.WARNING)

# QUERY UTILS

def load_data(city_ref, hdfs_keys, extra_method=None, extra_args=None, store_data=True):
	"""
	Retrieve the data, either from locally stored data or calculating it through `extra_method`

	Parameters
	----------
	city_ref : string
		name of the city
	hdfs_keys : list of string
		keys of data to be retrieved/stored
	extra_method : function
		method to calculate data if needed
	extra_args : array
		arguments for extra_method
	store_data : boolean
		optionally store data
	
	Returns
	----------
	list of pandas.DataFrame
		data related to input keys
	"""
	# If stand-alone key required
	if ( type(hdfs_keys) != list ): hdfs_keys = [hdfs_keys] 
	# Get the data for each corresponding key
	results = [ _get_local_data(city_ref,key) for key in hdfs_keys ] 
	if ( any( [ element is None for element in results ] ) ): # Data is not stored locally. Method to calculate it is provided?
		if (extra_method == None): raise NoMethodProvidedException() # Method not provided
		results = extra_method(*extra_args) # Compute required results
		# If single object is computed
		if (len(hdfs_keys) == 1): results = [ results ] 
		# Store data
		if (store_data): [ _update_local_data(city_ref,result,key,DATA_HDF5_OVERWITE) for key,result in zip(hdfs_keys,results) ]
	# Return results
	if (len(results) == 1 ): return results[0] # Single object
	else: return results # List of objects


# PUBLIC METHODS

# Related to stored HDFS

def remove_hdfs_matching_keys(city_ref, matching_key):
	"""
	Remove all stored keys which match the regular expression "`matching_key` in key" for the HDFS corresponding to `city_ref`

	Parameters
	----------
	city_ref : string
		name of the city
	matching_key : string
		keys to remove matching input key
	
	Returns
	----------

	"""
	to_remove = [ key for key in get_stored_keys(city_ref) if matching_key in key ]	
	storage = pd.HDFStore(_generate_file_path(city_ref))
	for key in to_remove:
		storage.remove(key)
	storage.close()
	log( "Removed stored keys: "+str(to_remove) )

def get_stored_keys(city_ref):
	"""
	Returns stored keys in HDFS file for `city_ref`

	Parameters
	----------
	city_ref : string
		name of the city
	
	Returns
	----------
	list of string
		list with found keys
	"""
	try:
		storage = pd.HDFStore(_generate_file_path(city_ref), 'r')
		stored_keys = [ key[1:] for key in storage.keys() ] # Remove absolute path "/"
		storage.close()
		return stored_keys
	except Exception as e:
		log( "HDF Error while retrieving stored keys."+str(e), level=lg.ERROR )

def get_stored_cities():
	"""
	Get the name of the cities for which stored data exists

	Parameters
	----------
	
	Returns
	----------
	list of string
		list with the name of the cities found
	"""
	stored_cities = []
	try:
		for file in listdir(storage_folder):
			if file.endswith(FILE_EXTENSIONS):
				stored_cities.append( file[:-len(FILE_EXTENSIONS)] )
	except:
		pass
	return stored_cities

# Load data

###############################
## Interface to call
###############################

def get_osm_data(city_ref, osm_df_generator = process_shapefiles, extra_args = [], get_bbox = False, store_data = True):
	"""
	Retrieves the OpenStreetMap data related to activity and residential land usages for given `city_ref`.
	Loads the data if stored locally.

	Parameters
	----------
	city_ref : string
		name of the city
	osm_df_generator : function
		method to retrieve OSM data
	extra_args : array
		additional arguments for the function
	get_bbox : boolean
		Optionally return the pertinent bounding box
	store_data : boolean
		optionally store the data
	
	Returns
	----------
	pandas.DataFrame, array
		data frame with following columns: ['osm_id', 'geometry', 'classification', 'key_value', 'activity_type', 'closest_d'], and optionally its bounding box
			osm_id: Corresponds to its OpenStreetMap identifier
			geometry: shapely.geometry
			classification: Land usage type {"activity", "residential", "mixed"}
			key_value: OpenStreetMap pair(s) of key(s):value(s) which define(s) its land usage classification
			activity_type: Activity classification
			closest_d: [Optional column] Closest distance of buildings (polygons) to their nearest neighbor in meters
	"""
	start = time.time()

	df, bbox = load_data(city_ref, KEY_OSM_DF, osm_df_generator, extra_args, store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on OSM data retrieval: " + str(end-start) )

	if (get_bbox): return [ df, bbox ]
	else: return df

###############################

def get_mesh_grid(city_ref, step = 100, mesh_grid_generator = _create_mesh_grid, store_data = True):
	""" 
	Retrieve the mesh grid loading polygons relative to `city_ref` and using step parameter
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	mesh_grid_generator : function
		method to calculate mesh grid
	store_data : boolean
		optionally store the data
	
	Returns
	----------
	pandas.Panel
		meshgrid
	"""
	return load_data(city_ref, KEY_MESH_GRID+str(step), mesh_grid_generator, [ get_osm_data(city_ref) , step ], store_data )

def get_route_graph(city_ref, bbox = None, store_data = True, over_write = False):
	""" 
	Retrieves osmnx graph for given `city_ref`
	Loads the data if stored locally. Otherwise, it queries to OpenStreetMap using osmnx package
	If no bounding box is provided, the stored `city_ref` bounding box will be retrieved

	Parameters
	----------
	city_ref : string
		name of the city
	bbox : dict
		bounding box with nort, south, east and west coordinates
	store_data : boolean
		optionally store the data
	over_write : boolean
		optionally over-write stored graph

	Returns
	----------
	networkx multidigraph
		projected graph
	"""
	try:
		if (over_write): raise Exception # Force exception
		G = ox.load_graphml(city_ref+'_network.graphml')
		log( "Found graph for `"+city_ref+"` stored locally" )
	except:
		try:
			if (bbox is None): bbox = get_osm_data(city_ref, get_bbox=True)[1]
			# Create graph network from bounding box
			G = ox.graph_from_bbox(bbox.north, bbox.south, bbox.east, bbox.west, network_type='drive_service')
			# Project graph
			G = ox.project_graph(G)
			if (store_data): # save street network as GraphML file
				ox.save_graphml(G, filename=city_ref+'_network.graphml')
				log( "Data graph for `"+city_ref+"` has been stored" )
		except Exception as e:
			log( "Osmnx graph could not be retrieved."+str(e), level=lg.ERROR )
			return None
	return G

def get_kde_activities_residential(city_ref, step = 100, kde_generator = kde_activity_residential, store_data = True):
	""" 
	Retrieves densities for Activity and Residential land usages of given `city_ref`
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	kde_generator : function
		method to calculate densities
	extra_args : array
		additional arguments for the function to calculate indices
	store_data : boolean
		optionally store the data

	Returns
	----------
	pandas.DataFrame, pandas.DataFrame
		activity and residential calculated densities
	"""
	df = get_osm_data(city_ref)
	XY = get_mesh_grid( city_ref, step = step )
	# Key stored: KEY_LANDUSEMIX + step + RADIUS_SEARCH
	KEYS_ = [ key+str(step) for key in KEY_KDES ]

	start = time.time()

	kde_activity, kde_residential = load_data(city_ref, KEYS_, kde_generator, [ df, XY ], store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on Kernel density estimations: " + str(end-start) )

	return kde_activity, kde_residential


def get_kde_activities_types(city_ref, step = 100, kde_generator = kde_activity_types, store_data = True):
	""" 
	Retrieves densities for the different activity land usages classification of given `city_ref`
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	kde_generator : function
		method to calculate densities
	extra_args : array
		additional arguments for the function to calculate indices
	store_data : boolean
		optionally store the data

	Returns
	----------
	array of pandas.DataFrame, array of string
		returns an array of the rescaled densities for each activity classification, and an array of the activities classification
	"""
	df = get_osm_data(city_ref)
	XY = get_mesh_grid( city_ref, step = step )
	# Key stored: KEY_LANDUSEMIX + step + RADIUS_SEARCH
	KEYS_ = [ key+str(step) for key in KEY_KDES_ACTIVITIES ]

	start = time.time()

	kde_activities = load_data(city_ref, KEYS_, kde_generator, [ df, XY ], store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on activities Kernel density estimations: " + str(end-start) )

	return kde_activities, ACTIVITIES_CLASSIFICATION


def get_dispersion_indicator(city_ref, step = 100, dispersion_generator = compute_grid_dispersion, kw_args={"radius_search":750, "use_median":False}, extra_args = [], store_data = True):
	""" 
	Retrieves dispersion indices for given `city_ref`
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	dispersion_generator : function
		method to calculate indices
	extra_args : array
		additional arguments for the function to calculate indices
	store_data : boolean
		optionally store the data

	Returns
	----------
	pandas.DataFrame
		grid with calculated indices
	"""
	df = get_osm_data(city_ref)

	if ( not 'closest_d' in df.columns ): # Column for closest distance exists?
		# Compute for each building the distance to its closest building
		_apply_polygon_closest_distance_neighbor(df)
		# Remove stored data related to data frame
		remove_hdfs_matching_keys(city_ref, KEY_OSM_DF[0])
		# Store (Over-write) data including closest distance computation
		_update_local_data(city_ref, df, KEY_OSM_DF[0])
	
	XY = get_mesh_grid( city_ref, step = step )

	# Key stored: KEY_DISPERSION + step
	KEY_ = KEY_DISPERSION+str(step)

	start = time.time()

	# Calculate dispersion grid
	dispersion_grid = load_data(city_ref, KEY_, dispersion_generator, [ XY, df, kw_args, extra_args ], store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on Dispersion indices: " + str(end-start) )

	return dispersion_grid

def get_accessibility_indicator(city_ref, step = 100, accessibility_generator = compute_grid_accessibility, kw_args={}, extra_args = [], store_data = True):
	""" 
	Retrieves accessibility indices for given `city_ref`
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	accessibility_generator : function
		method to calculate indices
	extra_args : array
		additional arguments for the function to calculate indices
	store_data : boolean
		optionally store the data

	Returns
	----------
	pandas.DataFrame
		grid with calculated indices
	"""
	df, bbox = get_osm_data(city_ref, get_bbox=True)
	G = get_route_graph(city_ref, bbox)
	XY = get_mesh_grid( city_ref, step = step )

	# Key stored: KEY_ACCESSIBILITY + step
	KEY_ = KEY_ACCESSIBILITY+str(step)

	start = time.time()

	# Calculate accessibility grid
	accessibility_grid = load_data(city_ref, KEY_, accessibility_generator, [ XY, G, df, kw_args, extra_args ] , store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on Accessibility indices: " + str(end-start) )

	return accessibility_grid

def get_landusemix_indicator(city_ref, step = 100, landusemix_generator = compute_grid_landusemix, kw_args={}, extra_args = [], store_data = True):
	""" 
	Retrieves land use mix indices for given `city_ref`
	Loads the data if stored locally

	Parameters
	----------
	city_ref : string
		name of the city
	step : int
		grid step in meters
	landusemix_generator : function
		method to calculate indices
	extra_args : array
		additional arguments for the function to calculate indices
	store_data : boolean
		optionally store the data

	Returns
	----------
	pandas.DataFrame
		grid with calculated indices
	"""
	XY = get_mesh_grid( city_ref, step = step )
	kde_activities, kde_residential = get_kde_activities_residential(city_ref, step=step)
	
	# Key stored: KEY_LANDUSEMIX + step
	KEY_ = KEY_LANDUSEMIX+str(step)

	start = time.time()

	# Calculate land use mix grid
	landusemix_grid = load_data(city_ref, KEY_, landusemix_generator, [ XY, kde_activities, kde_residential, kw_args, extra_args ], store_data )

	end = time.time()
	if ((end - start) > 5): # Computation, not loading
		log( "Time on Land use mix indices: " + str(end-start) )

	return landusemix_grid
###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import osmnx as ox
import pandas as pd
import geopandas as gpd
import numpy as np
from scipy import spatial

from osmnx.utils import log

from .osm_tags import key_classification, landuse_classification, activity_classification

####################################################################################
# Under uncertainty: Residential assumption?
RESIDENTIAL_ASSUMPTION_UNCERTAINTY = True
####################################################################################

############################################
### Tag land use classification
############################################

def aggregate_classification(classification_list):
	""" 
	Aggregate into a unique classification given an input list of classifications

	Parameters
	----------
	classification_list : list
		list with the input land use classifications
	
	Returns
	----------
	string
		returns the aggregated classification
	"""
	if ("other" in classification_list): # other tag -> Non-interesting building
		classification = None
	elif ( ("activity" in classification_list) and ("residential" in classification_list) ): # Mixed
		classification = "mixed"
	elif ( "mixed" in classification_list ): # Mixed
		classification = "mixed"
	elif ("activity" in classification_list): # Activity
		classification = "activity"
	elif ("residential" in classification_list): # Residential
		classification = "residential"
	elif ("infer" in  classification_list): # To infer
		classification = "infer"
	else: # No valuable classification
		classification = None

	return classification	

def classify_tag(tags, return_key_value=True):
	""" 
	Classify the land use of input OSM tag in `activity`, `residential`, `mixed`, None, or `infer` (to infer later)

	Parameters
	----------
	tags : dict
		OpenStreetMap tags
	
	Returns
	----------
	string, dict
		returns the classification, and a dict relating `key`:`value` defining its classification
	"""
	# key_value: Dictionary of osm key : osm value
	classification, key_value = [], {}

	for key, value in key_classification.items():
		# Get the corresponding key tag (without its land use)
		key_tag = key.replace("activity_","").replace("residential_","").replace("other_","").replace("infer_","")

		if tags.get(key_tag) in value:
			# First part of key defines the land use
			new_classification = key.split("_")[0]
			# Add the new classification
			classification.append( new_classification )
			# Associate the key-value
			key_value[key_tag] = tags.get(key_tag)

	classification = aggregate_classification(classification)

	if (return_key_value):
		return classification, key_value
	else:
		return classification

############################################
### Land use inference
############################################

def classify_landuse_inference(land_use):
	""" 
	Classify input land use into a defined category: `other`, `activity`, `residential`, or None

	Parameters
	----------
	land_use : string
		input land use tag
	
	Returns
	----------
	string
		returns the land use classification
	"""
	for key, value in landuse_classification.items():
		# key: Classification ; value: keys contained in the classifcation
		if (land_use in value):
			return key
	# Uncertain case
	if (RESIDENTIAL_ASSUMPTION_UNCERTAINTY): # Undefined use. Assumption: Residential
		return "residential"
	else:
		return None # No tag

def compute_landuse_inference(df_buildings, df_landuse):
	""" 
	Compute land use inference for building polygons with no information
	The inference is done using polygons with defined land use
	A building polygon's land use is inferred by means of adopting the land use of the smallest encompassing polygon with defined land use

	Parameters
	----------
	df_buildings : geopandas.GeoDataFrame
		input buildings
	df_landuse : geopandas.GeoDataFrame
		land use polygons to aid inference procedure

	Returns
	----------
	
	"""
	# Get those indices which need to be inferred, and keep geometry column only
	df_buildings_to_infer = df_buildings.loc[ df_buildings['classification'] == 'infer', ["geometry"] ]
	# Add land use polygon's area
	df_landuse['area'] = df_landuse.apply(lambda x: x.geometry.area, axis=1)

	# Get geometries to infer within Land use polygons matching
	sjoin = gpd.sjoin(df_buildings_to_infer, df_landuse, op='within')

	# Add index column to sort values
	sjoin['index'] = sjoin.index
	# Sort values by index, then by area
	sjoin.sort_values(by=['index','area'], inplace=True)
	# Drop duplicates. Keep first (minimum computing area)
	sjoin.drop_duplicates(subset=['index'], keep='first', inplace=True)

	##### Set key:value and classification
	# Set default value: inferred:None
	df_buildings.loc[ df_buildings_to_infer.index, "key_value" ] = df_buildings.loc[ df_buildings_to_infer.index].apply(lambda x: {"inferred":None} , axis=1)
	# Set land use for those buildings within a defined land use polygon
	df_buildings.loc[ sjoin.index, "key_value" ] = sjoin.apply(lambda x: {'inferred':x.landuse}, axis=1)

	# Set classification
	df_buildings.loc[ df_buildings_to_infer.index, "classification" ] = df_buildings.loc[ df_buildings_to_infer.index, "key_value" ].apply(lambda x: classify_landuse_inference(x.get("inferred")) )

	# Remove useless rows
	df_buildings.drop( df_buildings[ df_buildings.classification.isin([None,"other"]) ].index, inplace=True)
	df_buildings.reset_index(inplace=True,drop=True)
	assert( len( df_buildings[df_buildings.classification.isnull()] ) == 0 )

############################################
### Activity type classification
############################################

def value_activity_category(x):
	""" 
	Classify the activity of input activity value

	Parameters
	----------
	x : string
		activity value
	
	Returns
	----------
	string
		returns the activity classification
	"""
	for key, value in activity_classification.items():
		if x in value:
			return key
	return None

def key_value_activity_category(key, value):
	""" 
	Classify the activity of input pair key:value

	Parameters
	----------
	key : string
		key dict
	value : string
		value dict
	
	Returns
	----------
	string
		returns the activity classification
	"""
	# Note that some values repeat for different keys (e.g. shop=fuel and amenity=fuel), but they do not belong to the same activity classification
	return {
		'shop': 'shop',
		'leisure': 'leisure/amenity',
		'amenity': 'leisure/amenity',
		'man_made' : 'commercial/industrial',
		'industrial' : 'commercial/industrial',
		'landuse' : value_activity_category(value),
		'inferred' : value_activity_category(value), # Inferred cases adopted land use values
		'building' : value_activity_category(value),
		'building:use' : value_activity_category(value),
		'building:part' : value_activity_category(value)
	}.get(key, None)

def classify_activity_category(key_values):
	""" 
	Classify input activity category into `commercial/industrial`, `leisure/amenity`, or `shop`

	Parameters
	----------
	key_values : dict
		contain pairs of key:value relating to its usage
	
	Returns
	----------
	string
		returns the activity classification
	"""
	####################
	### Categories: commercial/industrial, leisure/amenity, shop
	####################
	categories = set( [ key_value_activity_category(key,value) for key,value in key_values.items() ] )
	categories.discard(None)
	return list(categories)

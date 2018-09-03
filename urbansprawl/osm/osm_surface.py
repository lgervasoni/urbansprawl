###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import osmnx as ox
import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

from osmnx.utils import log

from .osm_tags import height_tags, activity_classification
from .osm_data import aggregate_classification

############################################
### Land uses surface association
############################################

def get_composed_classification(building, df_pois):
	""" 
	Retrieve the composed classification given the building's containing Points of Interest

	Parameters
	----------
	building : geopandas.GeoSeries
		input building
	df_pois : geopandas.GeoDataFrame
		Points of Interest contained in the building

	Returns
	----------
	geopandas.GeoSeries
		returns a composed classification building
	"""
	# POIs aggregated classification
	pois_classification = aggregate_classification( df_pois.classification.values )
	# Composed building-POIs classification
	composed_classification = aggregate_classification( [building.classification, pois_classification] )
	# Composed activity categories
	try:
		composed_activity_category = list( set( [element for list_ in df_pois.activity_category for element in list_] + building.activity_category ) )
	except: # df_pois.activity_category.isnull().all() Returns True
		composed_activity_category = building.activity_category
	# Create a Series for a new row with composed classification
	composed_row = pd.Series( [building.geometry,composed_classification,composed_activity_category,building.building_levels], index=["geometry","classification","activity_category","building_levels"])
	return composed_row

def sum_landuses(x, landuses_m2, default_classification = None, mixed_building_first_floor_activity=True):
	""" 
	Associate to each land use its correspondent surface use for input building
	Mixed-uses building Option 1:
		First floor: Activity use
		Rest: residential use
	Mixed-uses building Option 2:
		Half used for Activity uses, the other half Residential use


	Parameters
	----------
	x : geopandas.GeoSeries
		input building
	landuses_m2 : dict
		squared meter surface associated to each land use
	default_classification : pandas.Series
		main building land use classification and included activity types
	mixed_building_first_floor_activity : Boolean
		if True: Associates building's first floor to activity uses and the rest to residential uses
		if False: Associates half of the building's area to each land use (Activity and Residential)

	Returns
	----------

	"""
	# Empty
	if ( not x.get("geometry")): return
	# Mixed building assumption: First level for activity uses, the rest residential use	
	if (x["classification"] is "activity"): # Sum activity use
		landuses_m2["activity"] += x["geometry"].area * x["building_levels"]
		# Sum activity category m2
		area_per_activity_category = x["geometry"].area * x["building_levels"] / len( x["activity_category"] )
		for activity_type in x["activity_category"]:
			landuses_m2[activity_type] += area_per_activity_category		
	elif (x["classification"] is "mixed"): # Sum activity and residential use

		if (x["building_levels"] > 1) and (mixed_building_first_floor_activity): # More than one level
			# First floor
			landuses_m2["residential"] += x["geometry"].area * ( x["building_levels"] - 1 )
			# Rest of the building
			landuses_m2["activity"] += x["geometry"].area
			area_per_activity_category = x["geometry"].area
		
		else: # One level building
			landuses_m2["residential"] += x["geometry"].area * x["building_levels"] / 2.
			landuses_m2["activity"] += x["geometry"].area * x["building_levels"] / 2.
			area_per_activity_category = ( x["geometry"].area * x["building_levels"] / 2. ) / len( x["activity_category"] )
		
		# Sum activity category m2		
		for activity_type in x["activity_category"]:
			landuses_m2[activity_type] += area_per_activity_category			
	elif (x["classification"] == "residential"): # Sum residential use
		landuses_m2["residential"] += x["geometry"].area * x["building_levels"]
	else: 
		# Row does not contain a classification, use given default classification creating a new dict
		dict_x = {"classification":default_classification.classification, "geometry":x.geometry, "building_levels":x.building_levels, "activity_category":default_classification.activity_category}
		# Recursive call
		sum_landuses(dict_x, landuses_m2, mixed_building_first_floor_activity=mixed_building_first_floor_activity)


def calculate_landuse_m2(building, mixed_building_first_floor_activity=True):
	""" 
	Calculate the total squared meters associated to residential and activity uses for input building
	In addition, surface usage for each activity types is performed

	Parameters
	----------
	building : geopandas.GeoSeries
		input building
	mixed_building_first_floor_activity : Boolean
		if True: Associates building's first floor to activity uses and the rest to residential uses
		if False: Associates half of the building's area to each land use (Activity and Residential)

	Returns
	----------
	dict
		contains the total associated surface to each land use key
	"""
	# Initialize
	landuse_m2 = {}
	landuse_m2["activity"] = 0
	landuse_m2["residential"] = 0
	for activity_type in list( activity_classification.keys() ):
		landuse_m2[activity_type] = 0
	
	# Get the composed classification from input building + containing POIs
	building_composed_classification = get_composed_classification(building, building.pois_full_parts)
	
	def no_min_level_geometry(building_parts):
		"""
		Returns building parts with no min. level associated
		"""
		def no_min_level_tag(x): # Buildings starts from a specific num level?
			if ( x.get("building:min_level") or x.get("min_level") or x.get("building:min_height") or x.get("min_height") ):
				return True
			else:
				return False
		# Get the geometries of the contained buildings with no height/level tags available
		geometries = building_parts.loc[ building_parts.height_tags.apply(lambda x: no_min_level_tag(x) ) ].geometry
		
		# Create the union of those geometries
		no_min_level_geom = Polygon()
		for shape in geometries.values:
			no_min_level_geom = no_min_level_geom.union(shape)

		# Return the final shape
		return no_min_level_geom
	
	# Remove from the main building geometry, those building parts geometries that do not contain a minimum level/height: Avoid duplicating first level surface
	building_composed_classification.geometry = building_composed_classification.geometry.difference( no_min_level_geometry(building.full_parts) )

	# Sum land uses for main building
	sum_landuses(building_composed_classification, landuse_m2, mixed_building_first_floor_activity=mixed_building_first_floor_activity)
	
	# Sum land uses for building parts. If no classification given, use the building's land use
	building.full_parts.apply(lambda x: sum_landuses(x, landuse_m2, building_composed_classification[["classification","activity_category"]], mixed_building_first_floor_activity=mixed_building_first_floor_activity), axis=1)
	
	return landuse_m2

def associate_levels(df_osm, default_height, meters_per_level):
	""" 
	Calculate the effectiver number of levels for each input building
	Under missing tag data, default values are used
	A column ['building_levels'] is added to the data frame

	Parameters
	----------
	df_osm : geopandas.GeoDataFrame
		input data frame
	default_height : float
		default building height in meters
	meters_per_level : float
		default meters per level

	Returns
	----------

	"""
	def levels_from_height(height, meters_per_level):
		"""
		Returns estimated number of levels given input height (meters)
		"""
		levels = abs( round( height / meters_per_level ) )
		if (levels >= 1):
			return levels
		else: # By default: 1 level
			assert( height > 0 )
			return 1

	def associate_level(x):
		"""
		Associates the number of levels to input building given its height tags information
		Returns the absolute value in order to consider the cases of underground levels
		"""
		# Buildings starts from a specific num level?
		if ( x.get("building:min_level") ): # building:min_level
			min_level = x["building:min_level"]
		elif ( x.get("min_level") ): # min_level
			min_level = x["min_level"]
		######################### Height based
		elif ( x.get("building:min_height") ): # min_level
			min_level = levels_from_height( x["building:min_height"], meters_per_level )
		elif ( x.get("min_height") ): # min_level
			min_level = levels_from_height( x["min_height"], meters_per_level )
		else:
			min_level = 0

		######################### Levels based
		if ( x.get("building:levels") ): # Number of building:levels given
			number_levels = abs( x["building:levels"] - min_level )
		elif ( x.get("levels") ): # Number of levels given
			number_levels = abs( x["levels"] - min_level )
		######################### Height based
		elif ( x.get("building:height") ): # building:height given
			number_levels = abs( levels_from_height( x["building:height"], meters_per_level ) - min_level )
		elif ( x.get("height") ): # height given
			number_levels = abs( levels_from_height( x["height"], meters_per_level ) - min_level )
		else: # No information given
			number_levels = levels_from_height(default_height, meters_per_level)

		assert( number_levels >= 0 )
		if (number_levels == 0): # By default at least 1 level 
			number_levels = 1
		return number_levels

	df_osm["building_levels"] = df_osm.height_tags.apply( lambda x: associate_level(x) )

def classification_sanity_check(building):
	"""
	Performs a sanity check in order to achieve coherence between the building's classification and the amount of M^2 associated to each land use
	Example: A building's classification could be 'residential', but contains its building parts (occupying 100% of the area, then all land uses M^2 associated to this land use) contain an acitivty use
	This would impose a problem of coherence between the classification and its surface land use

	Parameters
	----------
	building : geopandas.GeoSeries
		one row denoting the building's information

	Returns
	----------
	string
		returns the coherent classification
	"""
	if ( building.landuses_m2["residential"] > 0 ):
		if ( building.landuses_m2["activity"] > 0 ): # Mixed use
			return "mixed" 
		else: # Residential use
			return "residential"
	else: # Activity use
		return "activity"

def compute_landuses_m2(df_osm_built, df_osm_building_parts, df_osm_pois, default_height=6, meters_per_level=3, mixed_building_first_floor_activity=True):
	"""
	Determine the effective number of levels per building or building parts
	Calculate the amount of squared meters associated to residential and activity uses per building
	In addition, surface usage for each activity types is performed

	Parameters
	----------
	df_osm_built : geopandas.GeoDataFrame
		OSM Buildings
	df_osm_building_parts : geopandas.GeoDataFrame
		OSM building parts
	df_osm_pois : geopandas.GeoDataFrame
		OSM Points of interest
	default_height : float
		default building height in meters
	meters_per_level : float
		default meters per level
	mixed_building_first_floor_activity : Boolean
		if True: Associates building's first floor to activity uses and the rest to residential uses
		if False: Associates half of the building's area to each land use (Activity and Residential)

	Returns
	----------

	"""
	# Associate the number of levels to each building / building part
	associate_levels(df_osm_built, default_height=default_height, meters_per_level=meters_per_level)
	associate_levels(df_osm_building_parts, default_height=default_height, meters_per_level=meters_per_level)

	##################
	# Calculate for each building, the M^2 associated to each land usage considering building parts (area calculated given UTM coordinates projection assumption)
	##################
	
	# Associate the complete data frame of containing building parts
	col_interest = ["geometry","activity_category", "classification", "key_value", "height_tags", "building_levels"]
	df_osm_built["full_parts"] = df_osm_built.containing_parts.apply(lambda x: df_osm_building_parts.loc[x, col_interest ] if isinstance(x, list) else df_osm_building_parts.loc[ [], col_interest ] )
	
	# Associate the complete POIs contained in buildings
	col_interest = ["geometry","activity_category", "classification", "key_value"]
	df_osm_built["pois_full_parts"] = df_osm_built.containing_poi.apply(lambda x: df_osm_pois.loc[x, col_interest] if isinstance(x, list) else df_osm_pois.loc[ [], col_interest] )

	# Calculate m2's for each land use, plus for each activity category
	df_osm_built["landuses_m2"] = df_osm_built.apply(lambda x: calculate_landuse_m2(x, mixed_building_first_floor_activity=mixed_building_first_floor_activity), axis=1 )

	# Drop added full parts
	df_osm_built.drop( ["full_parts"], axis=1, inplace=True )
	df_osm_built.drop( ["pois_full_parts"], axis=1, inplace=True )

	# Sanity check: For each building land use classfication, its M^2 associated to these land uses must be greater than 1
	df_osm_built["classification"] = df_osm_built.apply(lambda x: classification_sanity_check(x), axis=1 )

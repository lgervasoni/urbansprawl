###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import geopandas as gpd
import pandas as pd
import numpy as np
import osmnx as ox
import os.path
import time

from osmnx.utils import log

from .utils import get_aggregated_squares, get_population_df_filled_empty_squares

from .data_extract import get_extract_population_data
from ..osm.osm_surface import get_composed_classification
# Filenames
from ..parameters import get_population_urban_features_filename, get_population_training_validating_filename
# Sprawl indices
from ..sprawl.dispersion import compute_grid_dispersion
from ..sprawl.landusemix import compute_grid_landusemix
from ..sprawl.accessibility import compute_grid_accessibility


def compute_urban_features(df_osm_built, df_osm_pois, x_square):
	"""
	Computes a set of urban features for an input square, given as data the buildings and Points of Interest

	Parameters
	----------
	df_osm_built : geopandas.GeoDataFrame
		input buildings
	df_osm_pois : geopandas.GeoDataFrame
		input points of interest
	x_square : geopandas.GeoSeries
		geometry square where urban features will be calculated

	Returns
	----------
	geopandas.GeoSeries
		geometry with updated urban features
	"""
	try: # Building fully-contained or partially-contained in X's square geometry
		df_osm_built = gpd.sjoin( df_osm_built, gpd.GeoDataFrame([x_square], crs=df_osm_built.crs), op='intersects')
	except: # Empty intersection
		return x_square
	# Empty
	if ( df_osm_built.empty ): return x_square

	### Apply percentage of building presence within square: 1 if fully contained, 0.5 if half the building contained, ...
	df_osm_built['building_ratio'] = df_osm_built.apply( lambda x: x.geometry.intersection(x_square.geometry).area / x.geometry.area, axis=1)
	
	### Total M^2
	x_square['m2_total_residential'] = df_osm_built.apply(lambda x: x.landuses_m2['residential'] * x.building_ratio, axis=1).sum()
	x_square['m2_total_activity'] = df_osm_built.apply(lambda x: x.landuses_m2['activity'] * x.building_ratio, axis=1).sum()

	### M^2 footprint 
	df_osm_built_selection = df_osm_built[ df_osm_built.composed_classification.isin(['residential']) ]
	x_square['m2_footprint_residential'] = 0 if df_osm_built_selection.empty else df_osm_built_selection.apply(lambda x: x.geometry.area * x.building_ratio, axis=1).sum()
	df_osm_built_selection = df_osm_built[ df_osm_built.composed_classification.isin(['activity']) ]
	x_square['m2_footprint_activity'] = 0 if df_osm_built_selection.empty else df_osm_built_selection.apply(lambda x: x.geometry.area * x.building_ratio, axis=1).sum()
	df_osm_built_selection = df_osm_built[ df_osm_built.composed_classification.isin(['mixed']) ]
	x_square['m2_footprint_mixed'] = 0 if df_osm_built_selection.empty else df_osm_built_selection.apply(lambda x: x.geometry.area * x.building_ratio, axis=1).sum()

	### Containing building types
	count_ratio_by_classification = df_osm_built.groupby('composed_classification').building_ratio.sum()
	x_square['num_built_activity'] = count_ratio_by_classification.get('activity', 0)
	x_square['num_built_residential'] = count_ratio_by_classification.get('residential',0)
	x_square['num_built_mixed'] = count_ratio_by_classification.get('mixed',0)

	### Containing POIs
	df_osm_pois_selection = df_osm_pois[ df_osm_pois.classification.isin(["activity","mixed"]) ]
	x_square['num_activity_pois'] = len( gpd.sjoin( df_osm_pois_selection, gpd.GeoDataFrame([x_square], crs=df_osm_pois.crs), op='intersects' ) )
	
	### Total number of levels
	x_square['num_levels'] = df_osm_built.apply(lambda x: x.building_levels * x.building_ratio , axis=1).sum()
	x_square['num_buildings'] = df_osm_built.building_ratio.sum()

	### Percentage of built-up area
	x_square['built_up_relation'] = df_osm_built.apply(lambda x: x.geometry.area * x.building_ratio, axis=1).sum() / x_square.geometry.area 
		
	return x_square

def compute_full_urban_features(city_ref, df_osm_built=None, df_osm_pois=None, graph=None, df_insee=None, data_source=None, kwargs={"max_dispersion":15}):
	"""
	Computes a set of urban features for each square where population count data exists

	Parameters
	----------
	city_ref : string
		city reference name
	df_osm_built : geopandas.GeoDataFrame
		input buildings
	df_osm_pois : geopandas.GeoDataFrame
		input points of interest
	graph : 
	x_square : geopandas.GeoSeries
		geometry square where urban features will be calculated

	Returns
	----------
	geopandas.GeoSeries
		geometry with updated urban features
	"""

	# Population extract exists?
	if ( os.path.exists( get_population_urban_features_filename(city_ref, data_source) ) ):
		log("Urban features from population gridded data exist for input city: "+city_ref)
		# Read from GeoJSON (default projection coordinates)
		df_insee_urban_features_4326 = gpd.read_file( get_population_urban_features_filename(city_ref, data_source) )
		# Project to UTM coordinates
		return ox.project_gdf(df_insee_urban_features_4326)

	# Required arguments
	assert( not df_osm_built is None )
	assert( not df_osm_pois is None )
	assert( not graph is None )
	assert( not df_insee is None )

	# Copy data frame in order to modify it
	#df_insee_urban_features = df_insee.copy()
	# Data frame + creation of empty squares with 0 count population
	df_insee_urban_features = get_population_df_filled_empty_squares(df_insee)
	##################
	### Sprawling indices
	##################
	df_insee_urban_features['geometry_squares'] = df_insee_urban_features.geometry
	df_insee_urban_features['geometry'] = df_insee_urban_features.geometry.centroid
	
	'''
	compute_grid_accessibility(df_insee_urban_features, graph, df_osm_built, df_osm_pois)
	'''
	
	# Compute land uses mix + densities estimation
	compute_grid_landusemix(df_insee_urban_features, df_osm_built, df_osm_pois)
	# Dispersion indices
	compute_grid_dispersion(df_insee_urban_features, df_osm_built)

	if (kwargs.get("max_dispersion")): # Set max bounds for dispersion values
		df_insee_urban_features.loc[ df_insee_urban_features.dispersion > kwargs.get("max_dispersion"), "dispersion" ] = kwargs.get("max_dispersion")
	
	# Set back original geometries
	df_insee_urban_features['geometry'] = df_insee_urban_features.geometry_squares
	df_insee_urban_features.drop('geometry_squares', axis=1, inplace=True)
	
	##################
	### Additional urban features
	##################	
	log("Composed classification (Buildings + parts + POIs) calculation")
	start = time.time()

	# Calculate the composed classification
	df_osm_built.containing_poi = df_osm_built.containing_poi.apply(lambda x: x if isinstance(x,list) else [])
	df_osm_built.activity_category = df_osm_built.activity_category.apply(lambda x: x if isinstance(x,list) else [])

	df_osm_built['composed_classification'] = df_osm_built.apply(lambda x: get_composed_classification( x, df_osm_pois.loc[x.containing_poi] ).classification, axis=1 )
	df_osm_built.loc[ df_osm_built.containing_poi.apply(lambda x: len(x)==0 ), "containing_poi" ] = np.nan
	df_osm_built.loc[ df_osm_built.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan

	end = time.time()
	log("Composed classification (Buildings + parts + POIs) calculation time: "+str(end-start) )

	# Compute the urban features for each square
	log("Urban features calculation")
	start = time.time()

	df_insee_urban_features = df_insee_urban_features.apply(lambda x: compute_urban_features(df_osm_built, df_osm_pois, x), axis=1)
	# FillNA, set CRS
	df_insee_urban_features.fillna(0, inplace=True)
	df_insee_urban_features.crs = df_insee.crs

	end = time.time()
	log("Urban features calculation time: "+str(end-start) )

	# Save to GeoJSON file (no projection conserved, then use EPSG 4326)
	ox.project_gdf(df_insee_urban_features, to_latlong=True).to_file( get_population_urban_features_filename(city_ref, data_source), driver='GeoJSON' )
	
	return df_insee_urban_features


def get_training_testing_data(city_ref, df_insee_urban_features=None):
	"""
	Returns the Y and X arrays for training/testing population downscaling estimates.

	Y contains vectors with the correspondent population densities
	X contains vectors with normalised urban features
	X_columns columns referring to X values
	Numpy arrays are stored locally

	Parameters
	----------
	cities_selection : string
		list of cities to select
	cities_skip : string
		list of cities to skip (retrieve the rest)

	Returns
	----------
	np.array, np.array, np.array
		Y vector, X vector, X column names vector
	"""
	# Population extract exists?
	if ( os.path.exists( get_population_training_validating_filename(city_ref) ) ):
		log("Urban population training+validation data/features exist for input city: " + city_ref)
		# Read from Numpy.Arrays
		data = np.load( get_population_training_validating_filename(city_ref) )
		# Project to UTM coordinates
		return data["Y"], data["X"], data["X_columns"]

	log("Urban population training+validation data/features calculation for city: " + city_ref)

	# Select columns to normalise
	columns_to_normalise = [col for col in df_insee_urban_features.columns if "num_" in col or "m2_" in col or "dispersion" in col or "accessibility" in col]
	# Normalise selected columns
	df_insee_urban_features.loc[:,columns_to_normalise] = df_insee_urban_features.loc[:,columns_to_normalise].apply(lambda x: x / x.max() , axis=0)

	# By default, idINSPIRE for created squares (0 population count) is 0: Change for 'CRS' string: Coherent with squares aggregation procedure (string matching)
	df_insee_urban_features.loc[ df_insee_urban_features.idINSPIRE == 0, "idINSPIRE" ] = "CRS"

	# Aggregate 5x5 squares: Get all possible aggregations (step of 200 meters = length of individual square)
	aggregated_df_insee_urban_features = get_aggregated_squares(ox.project_gdf(df_insee_urban_features, to_crs="+init=epsg:3035"), step=200., conserve_squares_info=True)

	# X values: Vector <x1,x2, ... , xn> with normalised urban features
	X_values = []
	# Y values: Vector <y1, y2, ... , ym> with normalised population densities. m=25
	Y_values = []

	# For each <Indices> combination, create a X and Y vector
	for idx in aggregated_df_insee_urban_features.indices:
		# Extract the urban features in the given 'indices' order (Fill to 0 for non-existent squares)
		square_info = df_insee_urban_features.reindex( idx ).fillna(0)
		# Y input (Ground truth): Population densities
		population_densities = (square_info["pop_count"] / square_info["pop_count"].sum() ).values

		if (all (pd.isna(population_densities)) ): # If sum of population count is 0, remove (NaN values)
			continue

		# X input: Normalised urban features
		urban_features = square_info[ [col for col in square_info.columns if col not in ['idINSPIRE','geometry','pop_count']] ].values

		# Append X, Y
		X_values.append(urban_features)
		Y_values.append(population_densities)

	# Get the columns order referenced in each X vector
	X_values_columns = df_insee_urban_features[  [col for col in square_info.columns if col not in ['idINSPIRE','geometry','pop_count']]  ].columns
	X_values_columns = np.array(X_values_columns)

	# To Numpy Array
	X_values = np.array(X_values)
	Y_values = np.array(Y_values)

	# Save to file
	np.savez( get_population_training_validating_filename(city_ref), X=X_values, Y=Y_values, X_columns=X_values_columns)
	
	log("Urban population training+tesing data/features done")

	return Y_values, X_values, X_values_columns

def get_Y_X_features_population_data(cities_selection=None, cities_skip=None):
	"""
	Returns the Y and X arrays for training/testing population downscaling estimates.
	It gathers either a selection of cities or all stored cities but a selected list to skip

	Y contains vectors with the correspondent population densities
	X contains vectors with normalised urban features
	X_columns columns referring to X values
	Numpy arrays are previously stored

	Parameters
	----------
	cities_selection : string
		list of cities to select
	cities_skip : string
		list of cities to skip (retrieve the rest)

	Returns
	----------
	np.array, np.array, np.array
		Y vector, X vector, X column names vector
	"""
	arr_X, arr_Y = [], []
	
	# Get the complete training-testig dataset
	for Y_X_data_city in os.listdir("data/training"):
		# Only if it contains a valid extension
		if ('.npz' not in Y_X_data_city): continue
		
		# Get city's name
		city_ref = Y_X_data_city.replace('_X_Y.npz','')
		
		# Only retrieve data from cities_selection (if ever given)
		if ( (cities_selection is not None) and (city_ref not in cities_selection) ): 
			log('Skipping city:', city_ref)
			continue
			
		# Skip cities data from from cities_skip (if ever given)
		if ( (cities_skip is not None) and (city_ref in cities_skip) ): 
			log('Skipping city:', city_ref)
			continue
		
		log('Retrieving data:', city_ref)
		
		# Get stored data
		city_Y, city_X, city_X_cols = get_training_testing_data(city_ref)
		# Append values
		arr_Y.append(city_Y)
		arr_X.append(city_X)
		
	# Assumption: All generated testing-training data contain the same X columns
	return np.concatenate(arr_Y), np.concatenate(arr_X), city_X_cols
###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
###################################################################################################

import geopandas as gpd
import pandas as pd
import numpy as np
import osmnx as ox
import os.path

from osmnx.utils import log

from .data_extract import get_extract_population_data
from ..osm.osm_surface import get_composed_classification
# Filenames
from ..parameters import get_population_urban_features_filename
# Sprawl indices
from ..sprawl.dispersion import compute_grid_dispersion
from ..sprawl.landusemix import compute_grid_landusemix
from ..sprawl.accessibility import compute_grid_accessibility


def compute_urban_features(df_osm_built, df_osm_pois, x_square):
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
	#count_by_classification = df_join.groupby('classification').size()
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

def compute_full_urban_features(df_osm_built=None, df_osm_pois=None, graph=None, df_insee=None, city_ref=None, data_source=None):
	""" Compute the urban features for each INSEE square
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


	df_insee_urban_features = df_insee.copy()
	##################
	### Sprawling indices
	##################
	df_insee_urban_features['geometry_squares'] = df_insee_urban_features.geometry
	df_insee_urban_features['geometry'] = df_insee_urban_features.geometry.centroid
	
	'''
	compute_grid_accessibility(df_insee_urban_features, graph, df_osm_built, df_osm_pois)
	'''
	compute_grid_dispersion(df_insee_urban_features, df_osm_built)
	# Compute land uses densities estimation
	compute_grid_landusemix(df_insee_urban_features, df_osm_built, df_osm_pois)
	
	# Set back original geometries
	df_insee_urban_features['geometry'] = df_insee_urban_features.geometry_squares
	df_insee_urban_features.drop('geometry_squares', axis=1, inplace=True)
	
	##################
	### Additional urban features
	##################	
	# Calculate the composed classification
	df_osm_built.containing_poi = df_osm_built.containing_poi.apply(lambda x: x if isinstance(x,list) else [])
	df_osm_built.activity_category = df_osm_built.activity_category.apply(lambda x: x if isinstance(x,list) else [])

	df_osm_built['composed_classification'] = df_osm_built.apply(lambda x: get_composed_classification( x, df_osm_pois.loc[x.containing_poi] ).classification, axis=1 )
	df_osm_built.loc[ df_osm_built.containing_poi.apply(lambda x: len(x)==0 ), "containing_poi" ] = np.nan
	df_osm_built.loc[ df_osm_built.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
	# Compute the urban features for each square
	df_insee_urban_features = df_insee_urban_features.apply(lambda x: compute_urban_features(df_osm_built, df_osm_pois, x), axis=1)
	# FillNA, set CRS
	df_insee_urban_features.fillna(0, inplace=True)
	df_insee_urban_features.crs = df_insee.crs

	# Save to GeoJSON file (no projection conserved, then use EPSG 4326)
	ox.project_gdf(df_insee_urban_features, to_latlong=True).to_file( get_population_urban_features_filename(city_ref, data_source), driver='GeoJSON' )
	
	return df_insee_urban_features
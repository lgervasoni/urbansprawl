###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

from shapely.geometry import Polygon, GeometryCollection
import geopandas as gpd
import pandas as pd
import os
import numpy as np
import osmnx as ox

from osmnx.utils import log
from ..parameters import get_population_extract_filename

DATA_SOURCES = ['insee','gpw']

##############################
### I/O for population data
##############################

def get_df_extract(df_data, poly_gdf, operation = "within"):
	"""
	Indexes input geo-data frame within an input region of interest
	If the region of interest is given as a polygon, its bounding box is indexed

	Parameters
	----------
	df_data : geopandas.GeoDataFrame
		input data frame to index
	poly_gdf : geopandas.GeoDataFrame
		geodataframe containing the region of interest in form of polygon
	operation : string
		the desired spatial join operation: 'within' or 'intersects'

	Returns
	----------
	geopandas.GeoDataFrame
		returns the population data frame indexed within the region of interest
	"""
	# Project to same system coordinates
	poly_gdf = ox.project_gdf(poly_gdf, to_crs=df_data.crs)
	# Spatial join
	df_extract = gpd.sjoin(df_data, poly_gdf, op=operation)
	# Keep original columns
	df_extract = df_extract[ df_data.columns ]
	return df_extract

def get_population_df(pop_shapefile, pop_data_file, data_source, to_crs, poly_gdf):
	"""
	Read the population shapefile from input filename/s
	Index the data within the bounding box
	Project to desired CRS

	Parameters
	----------
	pop_shapefile : string
		population count shapefile
	pop_data_file : string
		population data additional file (required for INSEE format)
	data_source : string
		desired population data source
	to_crs : dict
		desired coordinate reference system
	poly_gdf : geopandas.GeoDataFrame
		geodataframe containing the region of interest in form of polygon

	Returns
	----------
	geopandas.GeoDataFrame
		returns the indexed and projected population data frame
	"""
	#######################################
	### Load GPW/INSEE population data
	#######################################
	# Read population data
	df_pop = gpd.read_file(pop_shapefile)
		
	### Extract region of interest (epsg 4326)
	# Filter geometries not contained in bounding box
	df_pop = get_df_extract(df_pop, poly_gdf)

	if (data_source is 'insee'):
		#######################################
		### Additional step for INSEE data
		#######################################	
		# Read dbf files
		data_pop = gpd.read_file(pop_data_file)
		# Get columns of interest
		data_pop = data_pop[["idINSPIRE","ind_c"]]
		df_pop = df_pop[["geometry","idINSPIRE"]]
		# Inner join to obtain population count data associated to each geometry
		df_pop = pd.merge(df_pop, data_pop, how='inner', on='idINSPIRE')
	
	# Rename population count column
	df_pop.rename(columns={"ind_c":"pop_count", "DN":"pop_count"}, inplace=True)

	return ox.project_gdf(df_pop, to_crs=to_crs)

def get_extract_population_data(city_ref, data_source, pop_shapefile=None, pop_data_file=None, to_crs={'init': 'epsg:4326'}, df_osm_built=None):
	"""
	Get data population extract of desired data source for input city
	The population data frame is projected to the desired coordiante reference system
	Stores the extracted shapefile
	Returns the stored population data for input 'data source' and 'city reference' if it was previously stored

	Parameters
	----------
	city_ref : string
		name of input city
	data_source : string
		desired population data source
	pop_shapefile : string
		population count shapefile
	pop_data_file : string
		population data additional file (required for INSEE format)
	to_crs : dict
		desired coordinate reference system
	df_osm_built : geopandas.GeoDataFrame
		buildings for input region of interest

	Returns
	----------
	geopandas.GeoDataFrame
		returns the extracted population data
	"""
	# Input data source type given?
	assert( data_source in DATA_SOURCES )

	# Population extract exists?
	if ( os.path.exists( get_population_extract_filename(city_ref, data_source) ) ):
		log("Population extract exists for input city: "+city_ref)
		return gpd.read_file( get_population_extract_filename(city_ref, data_source) )

	# Input shape given?
	assert( not ( (np.all(df_osm_built is None) ) and (polygon is None) ) )
	# Input population shapefile given?
	assert( not pop_shapefile is None )
	# All input files given?
	assert( not ( (data_source == 'insee') and (pop_data_file is None) ) )

	# Get buildings convex hull
	polygon = GeometryCollection( df_osm_built.geometry.values.tolist() ).convex_hull
	# Convert to geo-dataframe with defined CRS
	poly_gdf = gpd.GeoDataFrame([polygon], columns=["geometry"], crs=df_osm_built.crs)
	
	# Compute extract
	df_pop = get_population_df(pop_shapefile, pop_data_file, data_source, to_crs, poly_gdf)
	
	# Save to shapefile
	df_pop.to_file( get_population_extract_filename(city_ref, data_source), driver='ESRI Shapefile' )
	return df_pop	


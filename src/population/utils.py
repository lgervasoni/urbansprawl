###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
###################################################################################################

import geopandas as gpd
import pandas as pd
import osmnx as ox
import numpy as np

from shapely.geometry import Polygon
from shapely.geometry import Point

def get_aggregated_squares(df_insee, step=1000.):
	"""
	Aggregates input population data in squares of 5x5
	Assumption: Input squares 200m by 200m
	INSEE data contains column 'idINSPIRE' which denotes in epsg:3035 the northing/easting coordinates of the south-west box endpoint
	Output: Aggregated squares of 1km by 1km

	Parameters
	----------
	df_insee : geopandas.GeoDataFrame
		INSEE population data
	step : float
		sampling step (default of 1 kilometer)

	Returns
	----------
	geopandas.GeoDataFrame
		returns the aggregated population data
	"""
	def get_northing_easting(x): # Extract northing and easting coordinates
		north, east = x.idINSPIRE.split("N")[1].split("E")
		x["north"] = int(north)
		x["east"] = int(east)
		return x

	def index_square(x, df_insee, offset_index):
		squares = df_insee.cx[ x.geometry.x - offset_index: x.geometry.x + offset_index , x.geometry.y - offset_index: x.geometry.y + offset_index]
		aggregated_polygon = Polygon()
		for geom in squares.geometry:
			aggregated_polygon = aggregated_polygon.union(geom)
		x["square_geometry"] = aggregated_polygon
		x["pop_count"] = squares.pop_count.sum()
		return x

	# Get northing and easting coordinates
	coordinates = df_insee.apply(lambda x: get_northing_easting(x), axis=1 )[["north","east"]]

	# North, east coordinates denote the south-west box endpoint: +500 meters to obtain the centroid of the 5x5 squared-box
	north_min, north_max = coordinates.north.min() + 500., coordinates.north.max() + 500.
	east_min, east_max = coordinates.east.min() + 500., coordinates.east.max() + 500.

	# Create mesh grid: One point for each square's centroid: Each square has an extent of 1km by 1km
	xv, yv = np.meshgrid( np.arange(east_min, east_max, step), np.arange(north_min, north_max, step) )
	points = [ Point(x,y) for x,y in zip( xv.ravel(), yv.ravel() ) ]
	df_squares = gpd.GeoDataFrame( points, columns=[ "geometry" ], crs="+init=epsg:3035" )
	# Project
	df_squares = ox.project_gdf(df_squares, to_crs = df_insee.crs)

	# Index, for each square centroid, +- 400 meters to achieve squares of 5 by 5
	df_squares = df_squares.apply(lambda x: index_square( x, df_insee, offset_index=400 ), axis=1 )
	# Update geometry
	df_squares['geometry'] = df_squares.square_geometry
	df_squares.drop( 'square_geometry', axis=1, inplace=True )
	# Remove empty squares
	df_squares.drop(df_squares[ df_squares.geometry.area == 0 ].index, axis=0, inplace=True)
	df_squares.reset_index(drop=True, inplace=True)
	# Set final CRS
	df_squares.crs = df_insee.crs
	return df_squares


def population_downscaling_validation(df_osm_built, df_insee):
	"""
	Validates the poplation downscaling estimation by means of aggregating the sum of buildings estimated population lying within each population square
	Allows to compare the real population count with the estimated population lying within each square
	Updates new column 'pop_estimation' for each square in the population data frame

	Parameters
	----------
	df_osm_built : geopandas.GeoDataFrame
		input buildings with computed population count
	df_insee : geopandas.GeoDataFrame
		INSEE population data

	Returns
	----------

	"""
	df_osm_built['geom'] = df_osm_built.geometry
	df_osm_built_residential = df_osm_built[ df_osm_built.apply(lambda x: x.landuses_m2['residential'] > 0, axis = 1) ]
	df_insee.crs = df_osm_built_residential.crs

	# Intersecting gridded population - buildings
	sjoin = gpd.sjoin( df_insee, df_osm_built_residential, op='intersects')
	# Calculate area within square (percentage of building with the square)
	sjoin['pop_estimation'] = sjoin.apply(lambda x: x.population * (x.geom.intersection(x.geometry).area / x.geom.area), axis=1 )
	
	# Initialize
	df_insee['pop_estimation'] = np.nan
	sum_pop_per_square = sjoin.groupby(sjoin.index)['pop_estimation'].sum()
	
	df_insee.loc[ sum_pop_per_square.index, "pop_estimation" ] = sum_pop_per_square.values
	# Drop unnecessary column
	df_osm_built.drop('geom', axis=1, inplace=True)
	# Set to 0 nan values
	df_insee.loc[ df_insee.pop_estimation.isnull(), "pop_estimation" ] = 0

	# Compute absolute and relative error
	df_insee["absolute_error"] = df_insee.apply(lambda x: abs(x.pop_count - x.pop_estimation), axis=1)
	df_insee["relative_error"] = df_insee.apply(lambda x: abs(x.pop_count - x.pop_estimation) / x.pop_count, axis=1)
###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

from scipy import spatial
import numpy as np
import pandas as pd
import time

from osmnx.utils import log

##############################################################
### Dispersion indices methods
##############################################################

def closest_building_distance_median( point_ref, tree, df_closest_d, radius_search ):
	""" 
	Dispersion metric at point_ref
	Computes the median of the closest distance to another building for each building within a radius search
	Uses the input KDTree to accelerate calculations

	Parameters
	----------
	point_ref : shapely.Point
		calculate indice at input point
	tree : scipy.spatial.KDTree
		KDTree of buildings centroid
	df : pandas.DataFrame 
		data frame of buildings with closest distance calculation
	radius_search : float
		circle radius to consider the dispersion calculation at a local point

	Returns
	----------
	float
		value of dispersion at input point
	"""
	# Query buildings within radius search
	indices = tree.query_ball_point( point_ref, radius_search )
	# No dispersion value
	if (len(indices) == 0): return np.NaN
	# Calculate median of closest distance values. If no information is available, NaN is set
	return df_closest_d.loc[ indices ].median()

def closest_building_distance_average( point_ref, tree, df_closest_d, radius_search ):
	""" 
	Dispersion metric at point_ref
	Computes the mean of the closest distance to another building for each building within a radius search
	Uses the input KDTree to accelerate calculations

	Parameters
	----------
	point_ref : shapely.Point
		calculate indice at input point
	tree : scipy.spatial.KDTree
		KDTree of buildings centroid
	df : pandas.DataFrame 
		data frame of buildings with closest distance calculation
	radius_search : int
		circle radius to consider the dispersion calculation at a local point

	Returns
	----------
	float
		value of dispersion at input point
	"""
	# Query buildings within radius search
	indices = tree.query_ball_point( point_ref, radius_search )
	# No dispersion value
	if (len(indices) == 0): return np.NaN
	# Calculate mean of closest distance values. If no information is available, NaN is set
	return df_closest_d.loc[ indices ].mean()


##############################################################
### Dispersion indices calculation
##############################################################

def compute_grid_dispersion(df_indices, df_osm_built, kwargs={"radius_search":750, "use_median":True, "K_nearest":50} ):
	""" 
	Creates grid and calculates dispersion indices.

	Parameters
	----------
	df_indices : geopandas.GeoDataFrame
		data frame containing the (x,y) reference points to calculate indices
	df_osm_built : geopandas.GeoDataFrame
		data frame containing the building's geometries
	kw_args: dict
		additional keyword arguments for the indices calculation
			radius_search: int
				circle radius to consider the dispersion calculation at a local point
			use_median : bool
				denotes whether the median or mean should be used to calculate the indices
			K_nearest : int
				number of neighboring buildings to consider in evaluation

	Returns
	----------
	geopandas.GeoDataFrame
		data frame with the added column for dispersion indices
	"""
	log("Dispersion calculation")
	start = time.time()

	# Get radius search: circle radius to consider the dispersion calculation at a local point
	radius_search = kwargs["radius_search"]
	# Use the median or mean computation ?
	use_median = kwargs["use_median"]

	# Assign dispersion calculation method
	if (kwargs["use_median"]):
		_calculate_dispersion = closest_building_distance_median
	else:
		_calculate_dispersion = closest_building_distance_average
	
	# Calculate the closest distance for each building within K_nearest centroid buildings
	_apply_polygon_closest_distance_neighbor(df_osm_built, K_nearest = kwargs["K_nearest"] )
	
	# For dispersion calculation approximation, create KDTree with buildings centroid
	coords_data = [ point.coords[0] for point in df_osm_built.loc[ df_osm_built.closest_d.notnull() ].geometry.apply(lambda x: x.centroid) ]
	# Create KDTree
	tree = spatial.KDTree( coords_data )
	
	# Compute dispersion indices
	index_column = "dispersion"
	df_indices[index_column] = df_indices.geometry.apply(lambda x: _calculate_dispersion(x, tree, df_osm_built.closest_d, radius_search ) )
	
	# Remove added column
	df_osm_built.drop('closest_d', axis=1, inplace=True)

	end = time.time()
	log("Dispersion calculation time: "+str(end-start))
	

def _apply_polygon_closest_distance_neighbor(df_osm_built, K_nearest = 50):
	""" 
	Computes for each polygon, the distance to the (approximated) nearest neighboring polygon
	Approximation is done using distance between centroids to K nearest neighboring polygons, then evaluating the real polygon distance
	A column `closest_d` is added in the data frame

	Parameters
	----------
	df_osm_built: geopandas.GeoDataFrame
		data frame containing the building's geometries
	K_nearest: int
		number of neighboring polygons to evaluate

	Returns
	----------

	"""
	def get_closest_indices(tree, x, K_nearest):
		# Query the closest buidings considering their centroid
		return tree.query( x.centroid.coords[0] , k=K_nearest+1)[1][1:]
	def compute_closest_distance(x, buildings):
		# Minimum distance of all distances between reference building 'x' and the other buildings
		return (buildings.apply(lambda b: x.distance(b) ) ).min()

	# Use all elements to get the exact closest neighbor?
	if ( (K_nearest == -1) or (K_nearest >= len(df_osm_built)) ): K_nearest = len(df_osm_built)-1

	# Get separate list for coordinates
	coords_data = [ geom.centroid.coords[0] for geom in df_osm_built.geometry ]
	# Create KD Tree using polygon's centroid
	tree = spatial.KDTree( coords_data )

	# Get the closest buildings indices
	df_osm_built['closest_buildings'] = df_osm_built.geometry.apply(lambda x: get_closest_indices(tree, x, K_nearest) )
	# Compute the minimum real distance for the closest buildings
	df_osm_built['closest_d'] = df_osm_built.apply(lambda x: compute_closest_distance(x.geometry,df_osm_built.geometry.loc[x.closest_buildings]), axis=1 )
	# Drop unnecessary column
	df_osm_built.drop('closest_buildings', axis=1, inplace=True)
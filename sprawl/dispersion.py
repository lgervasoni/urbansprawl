###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

from scipy import spatial
import numpy as np
import pandas as pd
import time

from .utils import _create_grid, log


##############################################################
### Parameters
##############################################################
# Mean or median computation?
DEF_Use_median = False

# Local radius search in meters
DEF_Radius_search = 750

# Number of neighbors to consider for given input point. -1 indicates to perform a full brute-force search (no approximation done)
K_NEAREST_NEIGHBORS_DISPERSION = 100

##############################################################
### Dispersion indices methods
##############################################################

def closest_building_distance( point_ref, tree, df, radius_search, use_median ):
	""" 
	Dispersion metric at point_ref
	Computes the median/mean of the closest distance to another building for each building within a radius search
	Uses the input KDTree to accelerate calculations (approximation)

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

	# Find indices within radius of distance
	idx_within_radius = tree.query_ball_point( point_ref, radius_search )
	# Calculate mean of closest distance values. If no information is available, NaN is set
	if (use_median):
		return df.closest_d.iloc[idx_within_radius].median()
	else:
		return df.closest_d.iloc[idx_within_radius].mean()


#### Assign land use mix method
_calculate_dispersion = closest_building_distance

##############################################################
### Dispersion indices calculation
##############################################################

def compute_grid_dispersion(XX_YY, df, kwargs={}, extra_args=[]):
	""" 
	Creates grid and calculates dispersion indices.

	Parameters
	----------
	XX_YY : pandas.Panel 
		meshgrid with (x,y) reference points to calculate indices
	df : pandas.DataFrame
		input data frame. Must contain [x_utm_centroid,y_utm_centroid] columns
	kw_args: dict
		additional keyword arguments for the indices calculation
	extra_args : array
		additional arguments for the indices calculation

	Returns
	----------
	pandas.DataFrame
		grid result of dispersion calculation
	"""
	log("Dispersion calculation")
	start = time.time()

	# Create grid
	grid = _create_grid( XX_YY )

	# Get radius search: circle radius  to consider the dispersion calculation at a local point
	radius_search = kwargs.get("radius_search",DEF_Radius_search)
	# Use the median or mean computation ?
	use_median = kwargs.get("use_median",DEF_Use_median)

	# Get geometries with computed closest distance
	df_closest_d = df.loc[df.closest_d.notnull()]

	# For dispersion calculation approximation, create KDTree with buildings centroid
	coords_data = [ point.coords[0] for point in df_closest_d.geometry.apply(lambda x: x.centroid) ]
	# Create KDTree
	tree = spatial.KDTree( coords_data )

	# Iterate
	rows, cols = grid.shape[0] , grid.shape[1]

	# Percentage of calculation
	percentage_100 = [round(p) for p in (np.linspace(0,1,num=21) * rows) ]
	log("Rows, columns: "+str(rows)+" "+str(cols))
	log("Percentage completed:")

	for i in range(rows): #Rows
		percentage_achieved = float(i) /rows
		if (i in percentage_100):
			log( str(round(percentage_achieved,2)) )

		for j in range(cols): #Cols
			# Point reference (UTM coordinates)
			point_ref = XX_YY["x",i,j], XX_YY["y",i,j] # xx[i,j], yy[i,j]
			# Calculate value
			grid[i,j] = _calculate_dispersion( point_ref, tree, df_closest_d, radius_search, use_median )

	end = time.time()
	log("Dispersion calculation time: "+str(end-start))

	return pd.DataFrame(grid)


def _apply_polygon_closest_distance_neighbor(df, K_nearest = None):
	""" 
	Computes for each polygon, the distance to the (approximated) nearest neighboring polygon
	Approximation is done using distance between centroids to K nearest neighboring polygons, then evaluating the real polygon distance
	A column `closest_d` is added in the data frame

	Parameters
	----------
	df: pandas.DataFrame
		must contain the following columns: [x_utm_centroid,y_utm_centroid,polygon]
	K_nearest: int
		number of neighboring polygons to evaluate

	Returns
	----------

	"""
	# Initialize values
	df.loc[:,'closest_d'] = np.nan

	# Get buildings: Geometries which do not correspond to Points (Polygon, MultiPolygon, ...)
	# Assumption: Any filtered polygon with activity or residential use is a building (not always tagged this way)
	df_buildings = df[[x.type != "Point" for x in df['geometry'] ]]

	# If not given as a parameter, use utils.K_NEAREST_NEIGHBORS
	if (K_nearest == None): K_nearest = K_NEAREST_NEIGHBORS_DISPERSION
	# Use all elements to get the exact closest neighbor? OR, K_nearest value higher than the number of elements to analyse?
	if ( (K_nearest == -1) or (K_nearest >= len(df_buildings)) ): K_nearest = len(df_buildings)-1


	# Get separate list for coordinates
	#xy = list(zip(*[ geom.centroid.coords[0] for geom in df_buildings.geometry ] ))
	coords_data = [ geom.centroid.coords[0] for geom in df_buildings.geometry ]
	# Create KD Tree using polygon's centroid
	tree = spatial.KDTree( coords_data )

	# For each polygon: Query the {distance,polygon_index} to the nearest polygons according to their centroid
	query_arr = [ tree.query(i, k=K_nearest+1) for i in coords_data ]
	# Get the indices of the nearest (approximation) K polygon's
	query_indices = [ query[1][1:] for query in query_arr ]

	def closest_distance(poly_ref, indices, df_buildings):
		# Minimum distance of all distances between Reference Polygon and Indices Polygons
		return min( [ poly_ref.distance( df_buildings.geometry.iloc[idx] ) for idx in indices ] )
	# Foreach {Poly,KnearestNeighborPolygons}: Evaluate the (min) real distance relative to the approximated nearest polygons
	# Only calculated for building geometries
	df.loc[ [x.type != "Point" for x in df['geometry'] ] , 'closest_d'] = [ closest_distance(poly_ref,indices,df_buildings) for poly_ref,indices in zip(df_buildings.geometry,query_indices) ]
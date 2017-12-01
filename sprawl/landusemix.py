###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

import math
import numpy as np
import pandas as pd
import time

#from seaborn.distributions import _statsmodels_bivariate_kde
#from sklearn.neighbors.kde import KernelDensity
import statsmodels.nonparametric.api as smnp

from .utils import _create_grid, log


##############################################################
### Parameters
##############################################################
# Distance in meters to use as bandwidth. Denominated a walkable distance
WALKABLE_DISTANCE = 400

# Rescale PDF to [0;1]. Sum is greater or equal to 1. If False, the PDF is normalised
RESCALE_PDF = True


##############################################################
### Utils
##############################################################

def rescale_kde(kde_pdf):
	""" 
	Rescales the values of the probability density function
	Elements are divided by its maximum element
	Values are rescaled to [0,1]
	Note that the sum is greater or equal to 1

	Parameters
	----------
	kde_pdf : numpy.matrix
		evaluated probability density function of KDE

	Returns
	----------
	pandas.DataFrame
		rescaled densities
	"""
	# Convert to pandas Data Frame
	df_kde_pdf = pd.DataFrame(kde_pdf)
	max_element = df_kde_pdf.values.max()
	df_kde_pdf = df_kde_pdf.applymap(lambda x: x/max_element)
	return df_kde_pdf

def normalise_kde(kde_pdf):
	""" 
	Normalises the values of the probability density function
	Sum of elements is 1

	Parameters
	----------
	kde_pdf : numpy.matrix
		evaluated probability density function of KDE

	Returns
	----------
	pandas.DataFrame
		normalised densities
	"""
	# Convert to pandas Data Frame
	df_kde_pdf = pd.DataFrame(kde_pdf)
	sum_elements = df_kde_pdf.values.sum()
	df_kde_pdf = df_kde_pdf.applymap(lambda x: x/sum_elements)
	return df_kde_pdf

if (RESCALE_PDF):
	transform_kde = rescale_kde
else:
	transform_kde = normalise_kde

##############################################################
### Land use mix indices methods
##############################################################

def metric_phi_entropy(x,y):
	""" 
	Shannon's entropy metric
	Based on article "Comparing measures of urban land use mix, 2013"

	Parameters
	----------
	x : float
		probability related to land use X
	y : float
		probability related to land use Y

	Returns
	----------
	float
		entropy value
	"""
	if (x <= 0 or y <= 0): return np.nan
	# Sum = 1
	x_,y_ = x/(x+y), y/(x+y)
	phi_value = - ( ( x_*math.log(x_) ) + ( y_*math.log(y_) ) ) / math.log(2)
	return phi_value

#### Assign land use mix method
_land_use_mix = metric_phi_entropy

##############################################################
### Land use mix indices calculation
##############################################################

def compute_grid_landusemix(XX_YY, kde_activities, kde_residential, kw_args={}, extra_args=[]):
	""" 
	Calculate land use mix indices on input grid

	Parameters
	----------
	XX_YY : pandas.Panel
		meshgrid with (x,y) reference points to calculate indices
	kde_activities : pandas.DataFrame
		Activity land use densities
	kde_residential : pandas.DataFrame
		Residential land use densities
	kw_args: dict
		additional keyword arguments for the indices calculation
	extra_args : array
		additional arguments for the indices calculation

	Returns
	----------
	pandas.DataFrame
		land use mix indices
	"""
	log("Land use mix calculation")
	start = time.time()

	# Create grid
	grid = _create_grid( XX_YY )

	# KDE's: Activity and Residential uses as numpy.matrix
	np_kde_activities = np.matrix(kde_activities)
	np_kde_residential = np.matrix(kde_residential)

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
			grid[i,j] = _land_use_mix(np_kde_activities[i,j] , np_kde_residential[i,j])
	
	end = time.time()
	log("Land use mix calculation time: "+str(end-start))

	return pd.DataFrame(grid)
	
####

def kde_activity_residential(df, XX_YY):
	""" 
	Compute Kernel Density Estimations for activity and residential land uses

	Parameters
	----------
	df : pandas.DataFrame
		data with column `classification` of land usage
	XX_YY : pandas.Panel
		meshgrid to evaluate the probability density function

	Returns
	----------
	numpy.matrix, numpy.matrix
		returns rescaled activity and residential densities
	"""
	# Get Activities and Residential land uses
	df_activities = df.loc[ df.classification.isin(["activity","mixed"]) ]
	df_residential = df.loc[ df.classification.isin(["residential","mixed"]) ]
	# Compute KDE's
	return _kde(df_activities, XX_YY), _kde(df_residential,XX_YY)

def kde_activity_types(df, XX_YY):
	""" 
	Compute Kernel Density Estimations for a certain classification of activity land uses

	Parameters
	----------
	df : pandas.DataFrame
		data with column `activity_type` containing the activity classification
	XX_YY : pandas.Panel
		meshgrid to evaluate the probability density function

	Returns
	----------
	array of numpy.matrix
		returns an array of the rescaled densities for each activity classification
	"""	
	# Import activities classification
	from .osm_tags import ACTIVITIES_CLASSIFICATION
	densities_activities = []

	for activity_type in ACTIVITIES_CLASSIFICATION: # For each classification
		# Get data
		df_activity_type = df.loc[ df.activity_type.isin([activity_type]) ]
		# Estimate densities
		densities_activities.append( _kde(df_activity_type,XX_YY) )

	# Return estimated densities
	return densities_activities

def _kde(df, XX_YY):
	"""
	Evaluate the probability density function using Kernel Density Estimation of input geo-localized data
	KDE's bandwidth related to walkable-distances

	Parameters
	----------
	df : pandas.DataFrame
		input data with column [geometry] containing shapely geometries
	XX_YY : pandas.Panel
		meshgrid to evaluate the probability density function

	Returns
	----------
	pandas.DataFrame
		
	"""
	# Get coordinates
	_x, _y = list( zip(*[ geom.centroid.coords[0] for geom in df.geometry ] ) )

	# WALKABLE_DISTANCE meters bandwidth
	bandwidth = WALKABLE_DISTANCE
	# Kernel Density Estimation: Continuous estimation for X and Y dimensions
	kde = smnp.KDEMultivariate( [_x, _y] , "cc", [ bandwidth, bandwidth ])

	# Get mesh grid x and y values to evaluate the PDF
	x = XX_YY.x.values.ravel()
	y = XX_YY.y.values.ravel()
	# Evaluate Probability Density Function
	kde_pdf = kde.pdf( [x,y] ).reshape( XX_YY.x.shape )
	
	return transform_kde(kde_pdf)

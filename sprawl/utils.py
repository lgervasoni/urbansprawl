###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

import numpy as np
import pandas as pd
import utm
import scipy.spatial.distance as distance 
import geopandas as gpd
import math
import logging as lg
import os
import sys
import unicodedata

from . import parameters

##############################################################
### Parameters
##############################################################

# Numpy meshgrid indexing type
INDEXING = 'xy'

##############################################################
### Utils
##############################################################
	
def _create_bounds(axis_values):
	""" 
	Create minimum and maximum bounds given input values

	Parameters
	----------
	axis_values : array
		values

	Returns
	----------
	float, float
		minimum and maximum corresponding values
	"""
	return axis_values.min(), axis_values.max()

def _create_grid_axis(axis_values, step):
	""" 
	Create a range of values using input step

	Parameters
	----------
	axis_values : array
		values
	step: float
		step augmentation value

	Returns
	----------
	numpy array
		range of values
	"""
	_min, _max = _create_bounds(axis_values)
	return np.arange(_min, _max, step)

def _create_mesh_grid(df, step):
	""" 
	Creates a mesh grid covering geometries bounding box

	Parameters
	----------
	df : pandas.DataFrame
		data frame containing 'geometry' column
	step: float
		step augmentation value

	Returns
	----------
	pandas.Panel
		meshgrid containing x and y values
	"""
	# Extract x and y coordinates from geometries centroid
	x,y = list( zip(*[ point.coords[0] for point in df.geometry.apply(lambda x: x.centroid) ] ) )
	# Create Panel
	panel = pd.Panel( np.meshgrid( _create_grid_axis( pd.Series(x), step), _create_grid_axis( pd.Series(y), step), indexing=INDEXING) )
	# Set "x" and "y" items
	panel.items = ["x","y"]
	return panel

def _create_grid(grid_panel):
	""" 
	Creates a matrix grid initialized with zeros
	[row,column] indexing format is assumed

	Parameters
	----------
	grid_panel : pandas.Panel
		panel containing [ x_y_slice, rows, columns ]

	Returns
	----------
	numpy.matrix
		meshgrid containing x and y values
	"""
	# Return initialized grid matrix
	return np.zeros( grid_panel.shape[1:3] )

##############################################################

def get_cities_surfaces(city_refs, shapefile_dirs):
	""" 
	Compute the surface covered by the shapefiles of input city/es

	Parameters
	----------
	city_refs: list or string
		city or cities name/s
	shapefile_dir: string
		corresponding city shapefile/s

	Returns
	----------
	list
		calculated surfaces for corresponding cities
	"""
	# Only one element?
	if ( (type(city_refs) == str) and (type(shapefile_dirs) == str) ):
		city_refs = [city_refs]
		shapefile_dirs = [shapefile_dirs]

	cities_surfaces = []
	# For each input city and corresponding shapefile
	for city_ref, shapefile_dir in zip(city_refs, shapefile_dirs):
		# Read shapefile
		df_shape = gpd.read_file(shapefile_dir)
		# Get bounding box
		lon1, lat1, lon2, lat2 = df_shape.total_bounds

		# Get UTM zone number corresponding to the center of the bounding box
		zone_number = utm.latlon_to_zone_number( (lat1+lat2)/2 , (lon1+lon2)/2 )

		# Get X and Y bound coordinates of bounding box
		x1, y1 = utm.from_latlon(lat1, lon1, force_zone_number=zone_number)[0:2]
		x2, y2 = utm.from_latlon(lat2, lon2, force_zone_number=zone_number)[0:2]

		# Compute surface in squared meters
		surface_ = abs(x1-x2) * abs(y1-y2)

		# Append to list
		cities_surfaces.append(surface_)
		# Log results
		log( "City surface: "+city_ref+"\t"+str(surface_) )
	return cities_surfaces

##############################################################

def valid_shapefile_city(city_ref,shapefile_dir):
	""" 
	Determine if input shapefile can be read by geopandas

	Parameters
	----------
	city_ref: string
		city name
	shapefile_dir: string
		path to shapefile

	Returns
	----------
	boolean
		returns True if the file could be read, False otherwise
	"""
	try:
		gpd.read_file(shapefile_dir)
		log("Ok: "+city_ref)
		return True
	except:
		log("Error: "+city_ref,level=lg.ERROR)
		return False

##############################################################
# Logger
##############################################################

def log(message, level=None, name=None, filename=None):
	"""
	Write a message to the log file and/or print to the the console.
	Parameters
	----------
	message : string
		the content of the message to log
	level : int
		one of the logger.level constants
	name : string
		name of the logger
	filename : string
		name of the log file
	Returns
	-------
	None
	"""
	if level is None:
		level = parameters.log_level
	if name is None:
		name = parameters.log_name
	if filename is None:
		filename = parameters.log_filename

	# if logging to file is turned on
	if parameters.log_file:
		# get the current logger (or create a new one, if none), then log message at requested level
		logger = get_logger(level=level, name=name, filename=filename)
		if level == lg.DEBUG:
			logger.debug(message)
		elif level == lg.INFO:
			logger.info(message)
		elif level == lg.WARNING:
			logger.warning(message)
		elif level == lg.ERROR:
			logger.error(message)

	# if logging to console is turned on, convert message to ascii and print to the console
	if parameters.log_console:
		# capture current stdout, then switch it to the console, print the message, then switch back to what had been the stdout
		# this prevents logging to notebook - instead, it goes to console
		standard_out = sys.stdout
		sys.stdout = sys.__stdout__
		
		# convert message to ascii for console display so it doesn't break windows terminals
		message = unicodedata.normalize('NFKD', make_str(message)).encode('ascii', errors='replace').decode()
		print(message)
		
		sys.stdout = standard_out

def get_logger(level=None, name=None, filename=None):
	"""
	Create a logger or return the current one if already instantiated.
	Parameters
	----------
	level : int
		one of the logger.level constants
	name : string
		name of the logger
	filename : string
		name of the log file
	Returns
	-------
	logger.logger
	"""
	import datetime as dt

	if level is None:
		level = parameters.log_level
	if name is None:
		name = parameters.log_name
	if filename is None:
		filename = parameters.log_filename

	logger = lg.getLogger(name)

	# if a logger with this name is not already set up
	if not getattr(logger, 'handler_set', None):

		# get today's date and construct a log filename
		todays_date = dt.datetime.today().strftime('%Y_%m_%d')
		log_filename = '{}/{}_{}.log'.format(parameters.logs_folder, filename, todays_date)

		# if the logs folder does not already exist, create it
		if not os.path.exists(parameters.logs_folder):
			os.makedirs(parameters.logs_folder)

		# create file handler and log formatter and set them up
		handler = lg.FileHandler(log_filename, encoding='utf-8')
		formatter = lg.Formatter('%(asctime)s %(levelname)s %(name)s %(message)s')
		handler.setFormatter(formatter)
		logger.addHandler(handler)
		logger.setLevel(level)
		logger.handler_set = True

	return logger

def config(data_folder=parameters.storage_folder,
		   logs_folder=parameters.logs_folder,
		   imgs_folder=parameters.images_folder,
		   log_file=parameters.log_file,
		   log_console=parameters.log_console,
		   log_level=parameters.log_level,
		   log_name=parameters.log_name,
		   log_filename=parameters.log_filename):
	"""
	Configure the framework by setting different variable values
	Parameters
	---------
	data_folder : string
		where to save and load data files
	logs_folder : string
		where to write the log files
	imgs_folder : string
		where to save figures
	log_file : bool
		if true, save log output to a log file in logs_folder
	log_console : bool
		if true, print log output to the console
	log_level : int
		one of the logger.level constants
	log_name : string
		name of the logger
	log_filename : string
		name of the file to save the log
	Returns
	-------
	
	"""
	# set each variable to the passed-in parameter value
	parameters.storage_folder = data_folder
	parameters.images_folder = imgs_folder
	parameters.logs_folder = logs_folder
	parameters.log_console = log_console
	parameters.log_file = log_file
	parameters.log_level = log_level
	parameters.log_name = log_name
	parameters.log_filename = log_filename

	# if logging is turned on, log that we are configured
	if parameters.log_file or parameters.log_console:
		log('Log file initialized')


def make_str(value):
	"""
	Convert a passed-in value to unicode if Python 2, or string if Python 3.
	Parameters
	----------
	value : any
		the value to convert to unicode/string
	Returns
	-------
	unicode or string
	"""
	try:
		# for python 2.x compatibility, use unicode
		return unicode(value)
	except:
		# python 3.x has no unicode type, so if error, use str type
		return str(value)

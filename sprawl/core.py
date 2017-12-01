###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

from .parameters import shapefiles_folder, DEF_dispersion_step, DEF_landusemix_step, DEF_accessibility_step
from .utils import log
from .storage import get_osm_data, get_dispersion_indicator, get_landusemix_indicator, get_accessibility_indicator
import logging as lg

def get_input_shapefiles(shapefiles_dir = shapefiles_folder):
	"""
	Get available city shapefiles in local folder

	Parameters
	----------
	shapefiles_dir : string
		path to the folder containing folder corresponding to the shapefiles of each city
	
	Returns
	----------
	array, array, array
		Returns three arrays, respectively for the cities name, their polygon shapefiles path, and their point shapefiles path
	"""
	city_refs = []
	poly_files = []
	point_files = []

	import os
	for root, dirs, files in os.walk( shapefiles_dir ):
		for name in dirs:
			city_refs.append(name)
			log("Found shapefile for city: "+name)
		for name in files:
			if ( ("osm_polygon" in name) and (name.endswith(".shp")) ):
				poly_files.append( os.path.join(root, name) )
			if ( ("osm_point" in name) and (name.endswith(".shp")) ):
				point_files.append( os.path.join(root, name) )
	return city_refs, poly_files, point_files

def process_spatial_indices(city_ref, file_poly=None, file_point=None, K_nearest=50, step={"dispersion":100,"landusemix":100,"accessibility":200}, kw_args={}):
	"""
	Calculates urban sprawl spatial indices for input city
	Only the indices with associated step are processed
	For instance, if "dispersion" key does not appear in step parameter, its calculation will be omitted
	If no polygon and point shapefiles are provided, OSM data must have been stored before
	Results are stored locally

	Parameters
	----------
	city_ref : string
		name of the city to analyse
	file_poly: string
		path to the polygon shapefile
	file_point: string
		path to the point shapefile
	K_nearest: int
		number of nearest polygon neighbors to include in the land use inference of uncertain objects
	step: dict
		step to be used to calculate the different spatial indices
	kw_args: dict
		keyword arguments for spatial indices calculation

	Returns
	----------
	"""	
	log("*** Processing spatial indices for city: "+city_ref)
	log("Step: "+str(step), lg.INFO)
	log("K_nearest: "+str(K_nearest), lg.INFO)

	# Force computation of data frame
	K_nearest = 50
	df, bbox = get_osm_data(city_ref, extra_args = [file_poly, file_point, K_nearest] , get_bbox=True)

	# Computation of indices
	if ( step.get("dispersion",None) ):
		grid_dispersion = get_dispersion_indicator(city_ref, step.get("dispersion",DEF_dispersion_step), kw_args=kw_args )
	if ( step.get("landusemix",None) ):
		grid_landusemix = get_landusemix_indicator(city_ref, step.get("landusemix",DEF_landusemix_step), kw_args=kw_args )
	if ( step.get("accessibility",None) ):
		grid_accessibility = get_accessibility_indicator(city_ref, step.get("accessibility",DEF_accessibility_step), kw_args=kw_args )

	log("*** Finished processing spatial indices for city: "+city_ref)
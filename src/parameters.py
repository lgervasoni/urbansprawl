###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

storage_folder = 'data'
images_folder = 'images'

# Format for load/save the geo-data ['geojson','shp']
geo_format = 'geojson' # 'shp'
geo_driver = 'GeoJSON' # 'ESRI Shapefile'


def get_dataframes_filenames(city_ref_file):
	"""
	Get data frame file names for input city

	Parameters
	----------
	city_ref_file : string
		name of input city

	Returns
	----------
	[ string, string, string ]
		returns filenames for buildings, building parts, and points of interest
	
	"""
	import os
	if not(os.path.isdir(storage_folder)): 
		os.makedirs(storage_folder)
	geo_poly_file = storage_folder+"/"+city_ref_file+"_poly."+geo_format
	geo_poly_parts_file = storage_folder+"/"+city_ref_file+"_poly_parts."+geo_format
	geo_point_file = storage_folder+"/"+city_ref_file+"_poi."+geo_format
	return geo_poly_file, geo_poly_parts_file, geo_point_file

def get_population_extract_filename(city_ref_file, data_source):
	"""
	Get data population extract filename for input city

	Parameters
	----------
	city_ref_file : string
		name of input city
	data_source : string
		desired population data source

	Returns
	----------
	string
		returns the population extract filename
	
	"""
	# Folder exists?
	import os
	if not(os.path.isdir(storage_folder + "/" + data_source)): 
		os.makedirs(storage_folder + "/" + data_source)
	return storage_folder + "/" + data_source + "/" + city_ref_file + "_population.shp"
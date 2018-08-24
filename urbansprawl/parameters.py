###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

storage_folder = 'data'
images_folder = 'images'

# Format for load/save the geo-data ['geojson','shp']
geo_format = 'geojson' # 'shp'
geo_driver = 'GeoJSON' # 'ESRI Shapefile'

### Files referring to population gridded data
files = {}
# Geometries for INSEE population data
files["insee_shapefile"] = "../200m-carreaux-metropole/carr_pop4326.shp"
# dbf file with attributes
files["insee_data_file"] = "../200m-carreaux-metropole/car_m.dbf"
# Worlde-wide gridded population data: Converted to shapefile (EPSG 4326)
files["gpw_population_world"] = "../gpwv4/gpw-v4.shp"


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

def get_population_urban_features_filename(city_ref_file, data_source):
	"""
	Get population urban features extract filename for input city
	Force GeoJSON format: Shapefiles truncate column names

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
	return storage_folder + "/" + data_source + "/" + city_ref_file + "_urban_features." + geo_format

def get_population_training_validating_filename(city_ref_file, data_source="training"):
	"""
	Get population normalised urban features extract and population densities filename for input city
	Stored in Numpy.Arrays

	Parameters
	----------
	city_ref_file : string
		name of input city

	Returns
	----------
	string
		returns the numpy stored/storing filename
	
	"""
	# Folder exists?
	import os
	if not(os.path.isdir(storage_folder + "/" + data_source)): 
		os.makedirs(storage_folder + "/" + data_source)
	return storage_folder + "/" + data_source + "/" + city_ref_file + "_X_Y.npz"
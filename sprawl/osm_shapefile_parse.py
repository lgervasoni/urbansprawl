###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

import geopandas as gpd
import pandas as pd
from scipy import spatial
import utm
import numpy as np
import time

from .utils import log
from .osm_tags import landuse_notResidentialActivity, landuse_activities, landuse_residential
from .osm_tags import amenities_activities, leisure_activies, shop_activities, building_activities, building_residential, building_civic_amenity, building_commercial


####################################################################################
# Under uncertainty: Residential assumption?
RESIDENTIAL_ASSUMPTION_UNCERTAINTY = True

# Tags that define an activity or residential use
activity_uses = {"amenity":amenities_activities, "leisure":leisure_activies, "shop":shop_activities, "building":building_activities}
residential_uses = {"building":building_residential}

# Land use tags that aid to define either an activity, residential or other use
infer_notResidentialActivity = {"landuse":landuse_notResidentialActivity}
infer_activity = {"landuse":landuse_activities}
infer_residential = {"landuse":landuse_residential}
####################################################################################

######################
### Classify activity type
### Categories: commercial/industrial, leisure/amenity, shop

def value_activity_category(x):
	""" 
	Classify the activity of input activity value

	Parameters
	----------
	x : string
		activity value
	
	Returns
	----------
	string
		returns the activity classification
	"""
	if (x in shop_activities):
		return 'shop'
	if (x in amenities_activities):
		return 'leisure/amenity'
	if (x in leisure_activies):
		return 'leisure/amenity'
	if (x in building_civic_amenity):
		return 'leisure/amenity'
	if (x in building_commercial):
		return 'commercial/industrial'
	# By default: commercial/industrial 
	# Assumption: Otherwise it would have been likely to contain a shop, leisure or amenity tag
	return 'commercial/industrial'

def key_value_activity_category(key, value):
	""" 
	Classify the activity of input pair key:value

	Parameters
	----------
	key : string
		key dict
	value : string
		value dict
	
	Returns
	----------
	string
		returns the activity classification
	"""
	return {
		'shop': 'shop',
		'leisure': 'leisure/amenity',
		'amenity': 'leisure/amenity',
		'commercial' : 'commercial/industrial',
		'industrial' : 'commercial/industrial',
		'landuse' : value_activity_category(value),
		'inferred' : value_activity_category(value), # Inferred cases adopted land use values
		'building': value_activity_category(value)
	}.get(key, None)

def classify_activity_category(key_values):
	""" 
	Classify input activity category into `commercial/industrial`, `leisure/amenity`, or `shop`

	Parameters
	----------
	key_values : dict
		contain pairs of key:value relating to its usage
	
	Returns
	----------
	string
		returns the activity classification
	"""
	####################
	### Categories: commercial/industrial, leisure/amenity, shop
	#######
	categories = [ key_value_activity_category(key,value) for key,value in key_values.items() ]
	if (len(categories)==1): # Unique key_value
		return categories[0]
	if (len(categories)==2): # Mixed uses: activityKey, residentialKey
		# In case of mixed uses, find the one related to the activityKey (avoid the default value if other category exists)
		default_value = "commercial/industrial"
		return next( (categ for categ in categories if categ != default_value), default_value)
	assert(False) # Length of categories must be 1 or 2

######################
### Land use inference

def classify_landuse_inference(land_use):
	""" 
	Classify input land use into a defined category: `other`, `activity`, `residential`, or None

	Parameters
	----------
	land_use : string
		input land use tag
	
	Returns
	----------
	string
		returns the land use classification
	"""
	if (land_use in infer_notResidentialActivity["landuse"]): # Neither residential nor activity
		return "other"
	if (land_use in infer_activity["landuse"]): # Activity use
		return "activity"
	if (land_use in infer_residential["landuse"]): # Residential use
		return "residential"
	# Uncertain case
	if (RESIDENTIAL_ASSUMPTION_UNCERTAINTY): # Undefined use. Assumption: Residential
		return "residential"
	else:
		return None # No tag

def classify_tag(tags):
	""" 
	Classify the land use of input OSM tag in `activity`, `residential`, `mixed`, None, or `infer` (to infer later)

	Parameters
	----------
	tags : dict
		OpenStreetMap tags
	
	Returns
	----------
	string, dict
		returns the classification, and a dict relating `key`:`value` defining its classification
	"""
	# Classification: "residential", "activity", "mixed", "infer", or None
	# key_value: Dictionary of osm key : osm value
	classification, key_value = None, {}
	# Check activity tag
	for key in activity_uses.keys():
		if tags.get(key) in activity_uses[key]: # Found activity tag?
			classification = "activity"
			key_value[key] = tags[key]
			break
	# Check residential tag
	for key in residential_uses.keys():
		if tags.get(key) in residential_uses[key]: # Found residential tag?
			if (classification == "activity"): classification = "mixed"
			else: classification = "residential"
			key_value[key] = tags[key]
			break
	# Infer building use?
	if ( (classification is None) and (not tags.get("building") is None)  ): # No classifiaction found; Building tag
		classification = "infer"
	return classification, key_value

def infer_landuse(x_geometry, tree, aid_inference, K_nearest):
	""" 
	Infer the land use for input geometry
	Input geometry's land use is inferred by means of adopting the land use of the smallest encompassing polygon with defined land use

	Parameters
	----------
	x_geometry : shapely.Geometry
		input geometry
	tree : scipy.spatial.KDTree
		KDTree with polygons centroid
	aid_inference : pandas.DataFrame
		polygons with defined land use
	K_nearest : int
		number of nearest polygon neighbors to include in the inference of input geometry

	Returns
	----------
	string
		land use inferred
	"""
	x_point = x_geometry.centroid.coords[0]
	distances, indices = tree.query(x_point, k=K_nearest)
	# Initialize values
	area = float("inf")
	land_usage = None
	# Assumption: Smaller area have more expressiveness than bigger polygon areas
	# Get the land usage of the polygon containing x with defined land usage which minimizes its area
	for element in aid_inference.iloc[indices].itertuples(): # Iterate for closest queried elements
		if ( (element.geometry.contains(x_geometry) ) and (element.area < area) ): # If x is contained; Area small enough?
				area = element.area
				land_usage = element.landuse
	return land_usage

def compute_landuse_inference(df, K_nearest):
	""" 
	Compute land use inference for building polygons with no information
	The inference is done using polygons with defined land use
	A building polygon's land use is inferred by means of adopting the land use of the smallest encompassing polygon with defined land use
	A KDTree with the polygons centroid is used to perform a faster approximation of the smallest encompassing polygon, consindering the nearest `K_nearest` number of polygons
	The higher `K_nearest`, the better the approximation as well as the higher the computational load

	Parameters
	----------
	df : geopandas.GeoDataFrame
		geopandas data
	K_nearest : int
		number of nearest neighbors used in land use inference approximation

	Returns
	----------
	
	"""
	idx_to_infer = df[ df['classification'] == 'infer' ].index
	# Buildings which need to be inferred
	to_infer_geom = df.loc[ idx_to_infer, 'geometry' ]

	# Get polygons with defined land use -> Aid the inference process
	aid_inference = df[ df['landuse'].notnull() ][['landuse','geometry']].copy()
	aid_inference.reset_index(inplace=True, drop=True)
	aid_inference['area'] = aid_inference.geometry.area

	# K_nearest higher than the number of available polygons?
	if (K_nearest >= len(aid_inference)): K_nearest = len(aid_inference)

	# Get separate list for coordinates
	coords_data = [ point.coords[0] for point in aid_inference.centroid ]
	# Create input tree
	tree = spatial.KDTree( coords_data )
	
	# Update key:value with inferred land use for each element to infer
	df.loc[idx_to_infer , "key_value"] = to_infer_geom.apply( lambda x: {"inferred":infer_landuse(x, tree, aid_inference, K_nearest)} )
	# Set the classification use
	df.loc[idx_to_infer , "classification"] = df.loc[idx_to_infer , "key_value"].apply(lambda x: classify_landuse_inference(x.get("inferred")) )

	# Remove useless rows
	df.drop( df[ df.classification.isin([None,"other"]) ].index ,inplace=True )

######################

def project_to_utm(df, force_zone_number = None):
	""" 
	Project input geopandas data frame to UTM coordinates
	Bounding box center is used to define the zone number

	Parameters
	----------
	df : geopandas.GeoDataFrame
		geopandas data
	force_zone_number : int
		optionally force the UTM zone number projection

	Returns
	----------
	
	"""
	import utm
	if (force_zone_number):
		zone_number = force_zone_number
	else:
		# Get latitude and longitude bounds
		lon1, lat1, lon2, lat2 = df.total_bounds
		# Get UTM zone number corresponding to the center of the bounding box
		zone_number = utm.latlon_to_zone_number( (lat1+lat2)/2 , (lon1+lon2)/2 )
	# Define projection
	projection = "+proj=utm +zone="+str(zone_number)+" +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
	df.to_crs(projection, inplace=True)

def get_bounding_box(bbox1, bbox2):
	""" 
	Computes encompasing bounding box of two input bounding box

	Parameters
	----------
	bbox1 : list
		bounding box with longitude, latitude, longitude, latitude coordinates
	bbox2 : list
		bounding box with longitude, latitude, longitude, latitude coordinates

	Returns
	----------
	list
		encompassing bounding box with longitude, latitude, longitude, latitude coordinates
	"""
	lon1, lat1, lon2, lat2 = bbox1
	lon1_, lat1_, lon2_, lat2_ = bbox2
	return ( min(lon1,lon1_),min(lat1,lat1_), max(lon2,lon2_),max(lat2,lat2_) )

def clear_useless_columns(df):
	""" 
	Convert osm_id to long type. Drop useless rows and columns. Reset data frame index

	Parameters
	----------
	df : pandas.DataFrame
		input data frame

	Returns
	----------
	pandas.DataFrame
		returns modified data frame
	"""
	# Set long type to OSM identification
	df.osm_id = df.osm_id.apply(lambda x: int(x))
	
	# Drop rows which does not contain any important classification
	df = df[ df.classification.notnull() ]
	###df.dropna(axis=0, subset=['classification'], inplace=True)
	df = df[ df.key_value != {"inferred":"other"} ]
	###df.drop( df.index[df.key_value == {"inferred":"other"} ] ,inplace=True)
	
	# Drop columns which correspond to original OSM tags
	columns_of_interest = ['geometry','osm_id','classification','key_value']
	columns_drop = [ col for col in list(df.columns) if not col in columns_of_interest ]
	df.drop(columns_drop, axis=1, inplace=True)
	
	# Reset index
	df.reset_index(inplace=True,drop=True)
	return df

def process_shapefiles(polygon_shapefile, point_shapefile, K_nearest = 40):
	""" 
	Process the input polygon and point shapefiles (OSM2PGSQL)
	1) Read the input shapefiles
	2) Project to UTM coordinates
	3) Classify into residential, activity, or mixed land use
	4) Perform land use inference for undefined building polygons
	5) Clear useless rows and columns

	Parameters
	----------
	polygon_shapefile : string
		input polygons shapefile path
	point_shapefile: string
		input points shapefile path
	K_nearest : int
		number of nearest neighbors used in land use inference approximation

	Returns
	----------
	pandas.DataFrame, dict
		returns extracted data frame and its bounding box
	"""
	log("Processing OSM shapefile data")

	# Read the files using geopandas
	df_poly = gpd.read_file(polygon_shapefile)
	df_point = gpd.read_file(point_shapefile)

	# Convert to lowercase
	df_poly.rename(columns=lambda x: x.lower(), inplace=True)
	df_point.rename(columns=lambda x: x.lower(), inplace=True)

	### Remove columns which do not provide valuable information
	from .osm_tags import columns_osm_tag
	columns_of_interest = columns_osm_tag + ["osm_id","geometry"]
	df_poly.drop( [ col for col in list( df_poly.columns ) if not col in columns_of_interest ], axis=1, inplace=True)
	df_point.drop( [ col for col in list( df_point.columns ) if not col in columns_of_interest ], axis=1, inplace=True)

	# Get bounding box in latitude longitude coordinates {lon1,lat2,lon2,lat2}
	lon1,lat1,lon2,lat2 = get_bounding_box(df_poly.total_bounds, df_point.total_bounds)
	zone_number = utm.latlon_to_zone_number( (lat1+lat2)/2 , (lon1+lon2)/2 )

	# Project input geo data frame to UTM coordinates using same zone number
	project_to_utm(df_poly, force_zone_number=zone_number)
	project_to_utm(df_point, force_zone_number=zone_number)

	# Python 2+3 compatibility: Set none values
	def setNoneValues(x): # Set None values to blank '' strings
		if (x == ''): return None
		return x
	df_poly = df_poly.applymap(lambda x: setNoneValues(x) )
	df_point = df_point.applymap(lambda x: setNoneValues(x) )

	# Classify tags
	df_poly['classification'], df_poly['key_value'] = list( zip(*df_poly.apply( classify_tag, axis=1) ) )
	df_point['classification'], df_point['key_value'] = list( zip(*df_point.apply( classify_tag, axis=1) ) )
	
	# Remove points to infer. Inference is made only on polygons
	df_point = df_point[ ~ df_point.classification.isin(["infer"]) ]

	start = time.time()

	# Land use inference: Polygons
	compute_landuse_inference(df_poly, K_nearest=K_nearest)

	end = time.time()
	log("OSM data retrieval: Time on land use inference: "+str(end-start))

	# Clear useless columns and rows
	df_poly = clear_useless_columns(df_poly)
	df_point = clear_useless_columns(df_point)

	# Bounding box Pandas series format
	bbox = pd.Series( { "north":float(lat2), "south":float(lat1), "east":float(lon2), "west":float(lon1) } )

	# Append to same data frame, and convert to pandas Data Frame
	df = df_poly.append(df_point, ignore_index=True)

	# Classify activities
	df['activity_type'] = np.nan
	df.loc[df["classification"].isin(["activity","mixed"]), "activity_type" ] = df.loc[df["classification"].isin(["activity","mixed"]), "key_value" ].apply(classify_activity_category)

	return pd.DataFrame(df), bbox

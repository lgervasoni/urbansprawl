###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import time
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon

from osmnx.core import consolidate_subdivide_geometry, get_polygons_coordinates, overpass_request, get_osm_filter
from osmnx.plot import save_and_show
from osmnx.projection import project_geometry
from osmnx.utils import log, geocode
import logging as lg
import osmnx as ox

#######################################################################
### Buildings
#######################################################################

def create_buildings_gdf_from_input(date="", polygon=None, place=None, which_result=1, point=None, address=None, distance=None, north=None, south=None, east=None, west=None):
	""" 
	Retrieve OSM buildings according to input data
	Queries data for input region (polygon, place, point/address and distance around, or bounding box coordinates)
	Updates the used polygon/bounding box to determine the region of interest	

	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the landuse footprints within
	place : string or dict
		query string or structured query dict to geocode/download
	which_result : int
		result number to retrieve from geocode/download when using query string 
	point : tuple
		the (lat, lon) central point around which to construct the graph
	address : string
		the address to geocode and use as the central point around which to construct the graph
	distance : int
		retain only those nodes within this many meters of the center of the graph
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box

	Returns
	----------
	[ geopandas.GeoDataFrame, shapely.Polygon, float, float, float, float ]
		retrieved buildings, region of interest polygon, and region of interest bounding box
	"""
	##########################
	### Osmnx query: Buildings
	##########################
	if (not polygon is None):  # Polygon
		log("Input type: Polygon")
		# If input geo data frame, extract polygon shape
		if ( type(polygon) is gpd.GeoDataFrame ):
			assert( polygon.shape[0] == 1 )
			polygon = polygon.geometry[0]
		df_osm_built = buildings_from_polygon(date, polygon)
	
	elif ( all( [point,distance] ) ):  # Point + distance
		log("Input type: Point")
		df_osm_built = buildings_from_point(date, point, distance=distance)
		# Get bounding box
		west, south, east, north = df_osm_built.total_bounds
	
	elif ( all( [address,distance] ) ):  # Address
		log("Input type: Address")
		df_osm_built = buildings_from_address(date, address, distance=distance)
		# Get bounding box
		west, south, east, north = df_osm_built.total_bounds

	elif (place):  # Place
		log("Input type: Place")
		if (which_result is None): which_result = 1
		df_osm_built = buildings_from_place(date, place, which_result=which_result)
		# Get encompassing polygon
		poly_gdf = ox.gdf_from_place(place, which_result=which_result)
		polygon = poly_gdf.geometry[0]
	
	elif ( all( [north,south,east,west] ) ): # Bounding box
		log("Input type: Bounding box")
		# Create points in specific order
		p1 = (east,north)
		p2 = (west,north)
		p3 = (west,south)
		p4 = (east,south)	
		polygon = Polygon( [p1,p2,p3,p4] )
		df_osm_built = buildings_from_polygon(date, polygon)	
	else:
		log("Error: Must provide at least one input")
		return
	return df_osm_built, polygon, north, south, east, west

def osm_bldg_download(date="", polygon=None, north=None, south=None, east=None, west=None,
					  timeout=180, memory=None, max_query_area_size=50*1000*50*1000):
	"""
	Download OpenStreetMap building footprint data.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the building footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	max_query_area_size : float
		max area for any part of the geometry, in the units the geometry is in:
		any polygon bigger will get divided up for multiple queries to API
		(default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are
		meters))
	Returns
	-------
	list
		list of response_json dicts
	"""

	# check if we're querying by polygon or by bounding box based on which
	# argument(s) where passed into this function
	by_poly = polygon is not None
	by_bbox = not (north is None or south is None or east is None or west is None)
	if not (by_poly or by_bbox):
		raise ValueError('You must pass a polygon or north, south, east, and west')

	response_jsons = []

	# pass server memory allocation in bytes for the query to the API
	# if None, pass nothing so the server will use its default allocation size
	# otherwise, define the query's maxsize parameter value as whatever the
	# caller passed in
	if memory is None:
		maxsize = ''
	else:
		maxsize = '[maxsize:{}]'.format(memory)

	# define the query to send the API
	if by_bbox:
		# turn bbox into a polygon and project to local UTM
		polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
		geometry_proj, crs_proj = project_geometry(polygon)

		# subdivide it if it exceeds the max area size (in meters), then project
		# back to lat-long
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		log('Requesting building footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
		start_time = time.time()

		# loop through each polygon rectangle in the geometry (there will only
		# be one if original bbox didn't exceed max area size)
		for poly in geometry:
			# represent bbox as south,west,north,east and round lat-longs to 8
			# decimal places (ie, within 1 mm) so URL strings aren't different
			# due to float rounding issues (for consistent caching)
			west, south, east, north = poly.bounds
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};((way["building"]({south:.8f},'
							  '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["building"]'
							  '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')
			query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all building footprints data within bounding box from '
			   'API in {:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(geometry), time.time()-start_time))

	elif by_poly:
		# project to utm, divide polygon up into sub-polygons if area exceeds a
		# max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
		geometry_proj, crs_proj = project_geometry(polygon)
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		polygon_coord_strs = get_polygons_coordinates(geometry)
		log('Requesting building footprints data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
		start_time = time.time()

		# pass each polygon exterior coordinates in the list to the API, one at
		# a time
		for polygon_coord_str in polygon_coord_strs:
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};(way'
							  '(poly:"{polygon}")["building"];(._;>;);relation'
							  '(poly:"{polygon}")["building"];(._;>;));out;')
			query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all building footprints data within polygon from API in '
			   '{:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(polygon_coord_strs), time.time()-start_time))

	return response_jsons


def create_buildings_gdf(date="", polygon=None, north=None, south=None, east=None,
						 west=None, retain_invalid=False):
	"""
	Get building footprint data from OSM then assemble it into a GeoDataFrame.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the building footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	responses = osm_bldg_download(date, polygon, north, south, east, west)

	vertices = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='node':
				vertices[result['id']] = {'lat' : result['lat'],
										  'lon' : result['lon']}

	buildings = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='way':
				nodes = result['nodes']
				try:
					polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])
				except Exception:
					log('Polygon has invalid geometry: {}'.format(nodes))
				building = {'nodes' : nodes,
							'geometry' : polygon}

				if 'tags' in result:
					for tag in result['tags']:
						building[tag] = result['tags'][tag]

				buildings[result['id']] = building

	gdf = gpd.GeoDataFrame(buildings).T
	gdf.crs = {'init':'epsg:4326'}

	if not retain_invalid:
		# drop all invalid geometries
		gdf = gdf[gdf['geometry'].is_valid]

	return gdf


def buildings_from_point(date, point, distance, retain_invalid=False):
	"""
	Get building footprints within some distance north, south, east, and west of
	a lat-long point.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	point : tuple
		a lat-long point
	distance : numeric
		distance in meters
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	bbox = bbox_from_point(point=point, distance=distance)
	north, south, east, west = bbox
	return create_buildings_gdf(date=date, north=north, south=south, east=east, west=west, retain_invalid=retain_invalid)


def buildings_from_address(date, address, distance, retain_invalid=False):
	"""
	Get building footprints within some distance north, south, east, and west of
	an address.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	address : string
		the address to geocode to a lat-long point
	distance : numeric
		distance in meters
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	# geocode the address string to a (lat, lon) point
	point = geocode(query=address)

	# get buildings within distance of this point
	return buildings_from_point(date, point, distance, retain_invalid=retain_invalid)


def buildings_from_polygon(date, polygon, retain_invalid=False):
	"""
	Get building footprints within some polygon.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : Polygon
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	return create_buildings_gdf(date=date, polygon=polygon, retain_invalid=retain_invalid)


def buildings_from_place(date, place, which_result=1, retain_invalid=False):
	"""
	Get building footprints within the boundaries of some place.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	place : string
		the query to geocode to get geojson boundary polygon
	which_result : int
		result number to retrieve from geocode/download when using query string 
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""
	city = ox.gdf_from_place(place, which_result=which_result)
	polygon = city['geometry'].iloc[0]
	return create_buildings_gdf(date=date, polygon=polygon, retain_invalid=retain_invalid)

#######################################################################
### Street network graph
#######################################################################

def retrieve_route_graph(city_ref, date="", polygon=None, north=None, south=None, east=None, west=None, force_crs=None):
	""" 
	Retrieves street network graph for given `city_ref`
	Loads the data if stored locally
	Otherwise, it retrieves the graph from OpenStreetMap using the osmnx package
	Input polygon or bounding box coordinates determine the region of interest

	Parameters
	----------
	city_ref : string
		name of the city
	date : string
		query the database at a certain timestamp
	polygon : shapely.Polygon
		polygon shape of input city
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	force_crs : dict
		graph will be projected to input crs

	Returns
	----------
	networkx.multidigraph
		projected graph
	"""
	try:
		G = ox.load_graphml(city_ref+'_network.graphml')
		log( "Found graph for `"+city_ref+"` stored locally" )
	except:
		try:
			if (not polygon is None):
				G = graph_from_polygon(polygon, network_type='drive_service', date=date)
			elif ( all( [north,south,east,west] ) ):
				G = graph_from_bbox(north, south, east, west, network_type='drive_service', date=date)
			else: # No inputs
				log("Need an input to retrieve graph")
				assert(False)

			# Project graph
			G = ox.project_graph(G, to_crs=force_crs)
			
			# Save street network as GraphML file
			ox.save_graphml(G, filename=city_ref+'_network.graphml')
			log( "Graph for `"+city_ref+"` has been retrieved and stored" )
		except Exception as e:
			log( "Osmnx graph could not be retrieved."+str(e), level=lg.ERROR )
			return None
	return G

def graph_from_polygon(polygon, network_type='all_private', simplify=True,
					   retain_all=False, truncate_by_edge=False, name='unnamed',
					   timeout=180, memory=None, date="",
					   max_query_area_size=50*1000*50*1000,
					   clean_periphery=True, infrastructure='way["highway"]'):
	"""
	Create a networkx graph from OSM data within the spatial boundaries of the
	passed-in shapely polygon.
	Parameters
	----------
	polygon : shapely Polygon or MultiPolygon
		the shape to get network data within. coordinates should be in units of
		latitude-longitude degrees.
	network_type : string
		what type of street network to get
	simplify : bool
		if true, simplify the graph topology
	retain_all : bool
		if True, return the entire graph even if it is not connected
	truncate_by_edge : bool
		if True retain node if it's outside bbox but at least one of node's
		neighbors are within bbox
	name : string
		the name of the graph
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	date : string
		query the database at a certain timestamp
	max_query_area_size : float
		max size for any part of the geometry, in square degrees: any polygon
		bigger will get divided up for multiple queries to API
	clean_periphery : bool
		if True (and simplify=True), buffer 0.5km to get a graph larger than
		requested, then simplify, then truncate it to requested spatial extent
	infrastructure : string
		download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
		infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
	Returns
	-------
	networkx multidigraph
	"""

	# verify that the geometry is valid and is a shapely Polygon/MultiPolygon
	# before proceeding
	if not polygon.is_valid:
		raise ValueError('Shape does not have a valid geometry')
	if not isinstance(polygon, (Polygon, MultiPolygon)):
		raise ValueError('Geometry must be a shapely Polygon or MultiPolygon')

	if clean_periphery and simplify:
		# create a new buffered polygon 0.5km around the desired one
		buffer_dist = 500
		polygon_utm, crs_utm = project_geometry(geometry=polygon)
		polygon_proj_buff = polygon_utm.buffer(buffer_dist)
		polygon_buffered, _ = project_geometry(geometry=polygon_proj_buff, crs=crs_utm, to_latlong=True)

		# get the network data from OSM,  create the buffered graph, then
		# truncate it to the buffered polygon
		response_jsons = osm_net_download(polygon=polygon_buffered, network_type=network_type,
										  timeout=timeout, memory=memory,
										  max_query_area_size=max_query_area_size,
										  infrastructure=infrastructure)
		G_buffered = ox.create_graph(response_jsons, name=name, retain_all=True, network_type=network_type)
		G_buffered = ox.truncate_graph_polygon(G_buffered, polygon_buffered, retain_all=True, truncate_by_edge=truncate_by_edge)

		# simplify the graph topology
		G_buffered = ox.simplify_graph(G_buffered)

		# truncate graph by polygon to return the graph within the polygon that
		# caller wants. don't simplify again - this allows us to retain
		# intersections along the street that may now only connect 2 street
		# segments in the network, but in reality also connect to an
		# intersection just outside the polygon
		G = ox.truncate_graph_polygon(G_buffered, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

		# count how many street segments in buffered graph emanate from each
		# intersection in un-buffered graph, to retain true counts for each
		# intersection, even if some of its neighbors are outside the polygon
		G.graph['streets_per_node'] = ox.count_streets_per_node(G_buffered, nodes=G.nodes())

	else:
		# download a list of API responses for the polygon/multipolygon
		response_jsons = osm_net_download(polygon=polygon, network_type=network_type,
										  timeout=timeout, memory=memory,
										  max_query_area_size=max_query_area_size,
										  infrastructure=infrastructure)

		# create the graph from the downloaded data
		G = ox.create_graph(response_jsons, name=name, retain_all=True, network_type=network_type)

		# truncate the graph to the extent of the polygon
		G = ox.truncate_graph_polygon(G, polygon, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

		# simplify the graph topology as the last step. don't truncate after
		# simplifying or you may have simplified out to an endpoint beyond the
		# truncation distance, in which case you will then strip out your entire
		# edge
		if simplify:
			G = ox.simplify_graph(G)

	log('graph_from_polygon() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
	return G

def graph_from_bbox(north, south, east, west, network_type='all_private',
					simplify=True, retain_all=False, truncate_by_edge=False,
					name='unnamed', timeout=180, memory=None, date="",
					max_query_area_size=50*1000*50*1000, clean_periphery=True,
					infrastructure='way["highway"]'):
	"""
	Create a networkx graph from OSM data within some bounding box.
	Parameters
	----------
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	network_type : string
		what type of street network to get
	simplify : bool
		if true, simplify the graph topology
	retain_all : bool
		if True, return the entire graph even if it is not connected
	truncate_by_edge : bool
		if True retain node if it's outside bbox but at least one of node's
		neighbors are within bbox
	name : string
		the name of the graph
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	date : string
		query the database at a certain timestamp
	max_query_area_size : float
		max size for any part of the geometry, in square degrees: any polygon
		bigger will get divided up for multiple queries to API
	clean_periphery : bool
		if True (and simplify=True), buffer 0.5km to get a graph larger than
		requested, then simplify, then truncate it to requested spatial extent
	infrastructure : string
		download infrastructure of given type (default is streets (ie, 'way["highway"]') but other
		infrastructures may be selected like power grids (ie, 'way["power"~"line"]'))
	Returns
	-------
	networkx multidigraph
	"""

	if clean_periphery and simplify:
		# create a new buffered bbox 0.5km around the desired one
		buffer_dist = 500
		polygon = Polygon([(west, north), (west, south), (east, south), (east, north)])
		polygon_utm, crs_utm = project_geometry(geometry=polygon)
		polygon_proj_buff = polygon_utm.buffer(buffer_dist)
		polygon_buff, _ = project_geometry(geometry=polygon_proj_buff, crs=crs_utm, to_latlong=True)
		west_buffered, south_buffered, east_buffered, north_buffered = polygon_buff.bounds

		# get the network data from OSM then create the graph
		response_jsons = osm_net_download(north=north_buffered, south=south_buffered,
										  east=east_buffered, west=west_buffered,
										  network_type=network_type, timeout=timeout,
										  memory=memory, date=date,
										  max_query_area_size=max_query_area_size,
										  infrastructure=infrastructure)
		G_buffered = ox.create_graph(response_jsons, name=name, retain_all=retain_all, network_type=network_type)
		G = ox.truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=True, truncate_by_edge=truncate_by_edge)

		# simplify the graph topology
		G_buffered = ox.simplify_graph(G_buffered)

		# truncate graph by desired bbox to return the graph within the bbox
		# caller wants
		G = ox.truncate_graph_bbox(G_buffered, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

		# count how many street segments in buffered graph emanate from each
		# intersection in un-buffered graph, to retain true counts for each
		# intersection, even if some of its neighbors are outside the bbox
		G.graph['streets_per_node'] = ox.count_streets_per_node(G_buffered, nodes=G.nodes())

	else:
		# get the network data from OSM
		response_jsons = osm_net_download(north=north, south=south, east=east,
										  west=west, network_type=network_type,
										  timeout=timeout, memory=memory, date=date,
										  max_query_area_size=max_query_area_size,
										  infrastructure=infrastructure)

		# create the graph, then truncate to the bounding box
		G = ox.create_graph(response_jsons, name=name, retain_all=retain_all, network_type=network_type)
		G = ox.truncate_graph_bbox(G, north, south, east, west, retain_all=retain_all, truncate_by_edge=truncate_by_edge)

		# simplify the graph topology as the last step. don't truncate after
		# simplifying or you may have simplified out to an endpoint
		# beyond the truncation distance, in which case you will then strip out
		# your entire edge
		if simplify:
			G = ox.simplify_graph(G)

	log('graph_from_bbox() returning graph with {:,} nodes and {:,} edges'.format(len(list(G.nodes())), len(list(G.edges()))))
	return  G

def osm_net_download(polygon=None, north=None, south=None, east=None, west=None,
					 network_type='all_private', timeout=180, memory=None, date="",
					 max_query_area_size=50*1000*50*1000, infrastructure='way["highway"]'):
	"""
	Download OSM ways and nodes within some bounding box from the Overpass API.
	Parameters
	----------
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the street network within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	network_type : string
		{'walk', 'bike', 'drive', 'drive_service', 'all', 'all_private'} what
		type of street network to get
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	date : string
		query the database at a certain timestamp
	max_query_area_size : float
		max area for any part of the geometry, in the units the geometry is in:
		any polygon bigger will get divided up for multiple queries to API
		(default is 50,000 * 50,000 units [ie, 50km x 50km in area, if units are
		meters])
	infrastructure : string
		download infrastructure of given type. default is streets, ie,
		'way["highway"]') but other infrastructures may be selected like power
		grids, ie, 'way["power"~"line"]'
	Returns
	-------
	response_jsons : list
	"""

	# check if we're querying by polygon or by bounding box based on which
	# argument(s) where passed into this function
	by_poly = polygon is not None
	by_bbox = not (north is None or south is None or east is None or west is None)
	if not (by_poly or by_bbox):
		raise ValueError('You must pass a polygon or north, south, east, and west')

	# create a filter to exclude certain kinds of ways based on the requested
	# network_type
	osm_filter = get_osm_filter(network_type)
	response_jsons = []

	# pass server memory allocation in bytes for the query to the API
	# if None, pass nothing so the server will use its default allocation size
	# otherwise, define the query's maxsize parameter value as whatever the
	# caller passed in
	if memory is None:
		maxsize = ''
	else:
		maxsize = '[maxsize:{}]'.format(memory)

	# define the query to send the API
	# specifying way["highway"] means that all ways returned must have a highway
	# key. the {filters} then remove ways by key/value. the '>' makes it recurse
	# so we get ways and way nodes. maxsize is in bytes.
	if by_bbox:
		# turn bbox into a polygon and project to local UTM
		polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
		geometry_proj, crs_proj = project_geometry(polygon)

		# subdivide it if it exceeds the max area size (in meters), then project
		# back to lat-long
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		log('Requesting network data within bounding box from API in {:,} request(s)'.format(len(geometry)))
		start_time = time.time()

		# loop through each polygon rectangle in the geometry (there will only
		# be one if original bbox didn't exceed max area size)
		for poly in geometry:
			# represent bbox as south,west,north,east and round lat-longs to 8
			# decimal places (ie, within 1 mm) so URL strings aren't different
			# due to float rounding issues (for consistent caching)
			west, south, east, north = poly.bounds
			query_template = date+'[out:json][timeout:{timeout}]{maxsize};({infrastructure}{filters}({south:.8f},{west:.8f},{north:.8f},{east:.8f});>;);out;'
			query_str = query_template.format(north=north, south=south,
											  east=east, west=west,
											  infrastructure=infrastructure,
											  filters=osm_filter,
											  timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		log('Got all network data within bounding box from API in {:,} request(s) and {:,.2f} seconds'.format(len(geometry), time.time()-start_time))

	elif by_poly:
		# project to utm, divide polygon up into sub-polygons if area exceeds a
		# max size (in meters), project back to lat-long, then get a list of
		# polygon(s) exterior coordinates
		geometry_proj, crs_proj = project_geometry(polygon)
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		polygon_coord_strs = get_polygons_coordinates(geometry)
		log('Requesting network data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
		start_time = time.time()

		# pass each polygon exterior coordinates in the list to the API, one at
		# a time
		for polygon_coord_str in polygon_coord_strs:
			query_template = date+'[out:json][timeout:{timeout}]{maxsize};({infrastructure}{filters}(poly:"{polygon}");>;);out;'
			query_str = query_template.format(polygon=polygon_coord_str, infrastructure=infrastructure, filters=osm_filter, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		log('Got all network data within polygon from API in {:,} request(s) and {:,.2f} seconds'.format(len(polygon_coord_strs), time.time()-start_time))

	return response_jsons

#######################################################################
### Land use
#######################################################################

def osm_landuse_download(date="", polygon=None, north=None, south=None, east=None, west=None,
					  timeout=180, memory=None, max_query_area_size=50*1000*50*1000):
	"""
	Download OpenStreetMap landuse footprint data.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the landuse footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	max_query_area_size : float
		max area for any part of the geometry, in the units the geometry is in:
		any polygon bigger will get divided up for multiple queries to API
		(default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are
		meters))
	Returns
	-------
	list
		list of response_json dicts
	"""

	# check if we're querying by polygon or by bounding box based on which
	# argument(s) where passed into this function
	by_poly = polygon is not None
	by_bbox = not (north is None or south is None or east is None or west is None)
	if not (by_poly or by_bbox):
		raise ValueError('You must pass a polygon or north, south, east, and west')

	response_jsons = []

	# pass server memory allocation in bytes for the query to the API
	# if None, pass nothing so the server will use its default allocation size
	# otherwise, define the query's maxsize parameter value as whatever the
	# caller passed in
	if memory is None:
		maxsize = ''
	else:
		maxsize = '[maxsize:{}]'.format(memory)

	# define the query to send the API
	if by_bbox:
		# turn bbox into a polygon and project to local UTM
		polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
		geometry_proj, crs_proj = project_geometry(polygon)

		# subdivide it if it exceeds the max area size (in meters), then project
		# back to lat-long
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		log('Requesting landuse footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
		start_time = time.time()

		# loop through each polygon rectangle in the geometry (there will only
		# be one if original bbox didn't exceed max area size)
		for poly in geometry:
			# represent bbox as south,west,north,east and round lat-longs to 8
			# decimal places (ie, within 1 mm) so URL strings aren't different
			# due to float rounding issues (for consistent caching)
			west, south, east, north = poly.bounds
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};((way["landuse"]({south:.8f},'
							  '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["landuse"]'
							  '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')
			query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all landuse footprints data within bounding box from '
			   'API in {:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(geometry), time.time()-start_time))

	elif by_poly:
		# project to utm, divide polygon up into sub-polygons if area exceeds a
		# max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
		geometry_proj, crs_proj = project_geometry(polygon)
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		polygon_coord_strs = get_polygons_coordinates(geometry)
		log('Requesting landuse footprints data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
		start_time = time.time()

		# pass each polygon exterior coordinates in the list to the API, one at
		# a time
		for polygon_coord_str in polygon_coord_strs:
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};(way'
							  '(poly:"{polygon}")["landuse"];(._;>;);relation'
							  '(poly:"{polygon}")["landuse"];(._;>;));out;')
			query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all landuse footprints data within polygon from API in '
			   '{:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(polygon_coord_strs), time.time()-start_time))

	return response_jsons

def create_landuse_gdf(date="", polygon=None, north=None, south=None, east=None,
						 west=None, retain_invalid=False):
	"""
	Get landuse footprint data from OSM then assemble it into a GeoDataFrame.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the landuse footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	retain_invalid : bool
		if False discard any landuse footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	responses = osm_landuse_download(date, polygon, north, south, east, west)

	vertices = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='node':
				vertices[result['id']] = {'lat' : result['lat'],
										  'lon' : result['lon']}

	landuses = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='way':
				nodes = result['nodes']
				try:
					polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])
				except Exception:
					log('Polygon has invalid geometry: {}'.format(nodes))
				landuse = {'nodes' : nodes,
							'geometry' : polygon}

				if 'tags' in result:
					for tag in result['tags']:
						landuse[tag] = result['tags'][tag]

				landuses[result['id']] = landuse

	gdf = gpd.GeoDataFrame(landuses).T
	gdf.crs = {'init':'epsg:4326'}

	if not retain_invalid:
		# drop all invalid geometries
		gdf = gdf[gdf['geometry'].is_valid]

	return gdf

#######################################################################
### Points of interest
#######################################################################

def osm_pois_download(date="", polygon=None, north=None, south=None, east=None, west=None,
					  timeout=180, memory=None, max_query_area_size=50*1000*50*1000):
	"""
	Download OpenStreetMap POIs footprint data.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the POIs footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	max_query_area_size : float
		max area for any part of the geometry, in the units the geometry is in:
		any polygon bigger will get divided up for multiple queries to API
		(default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are
		meters))
	Returns
	-------
	list
		list of response_json dicts
	"""

	# check if we're querying by polygon or by bounding box based on which
	# argument(s) where passed into this function
	by_poly = polygon is not None
	by_bbox = not (north is None or south is None or east is None or west is None)
	if not (by_poly or by_bbox):
		raise ValueError('You must pass a polygon or north, south, east, and west')

	response_jsons = []

	# pass server memory allocation in bytes for the query to the API
	# if None, pass nothing so the server will use its default allocation size
	# otherwise, define the query's maxsize parameter value as whatever the
	# caller passed in
	if memory is None:
		maxsize = ''
	else:
		maxsize = '[maxsize:{}]'.format(memory)

	# define the query to send the API
	if by_bbox:
		# turn bbox into a polygon and project to local UTM
		polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
		geometry_proj, crs_proj = project_geometry(polygon)

		# subdivide it if it exceeds the max area size (in meters), then project
		# back to lat-long
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		log('Requesting POIs footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
		start_time = time.time()

		# loop through each polygon rectangle in the geometry (there will only
		# be one if original bbox didn't exceed max area size)
		for poly in geometry:
			# represent bbox as south,west,north,east and round lat-longs to 8
			# decimal places (ie, within 1 mm) so URL strings aren't different
			# due to float rounding issues (for consistent caching)
			west, south, east, north = poly.bounds
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};((node["amenity"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f}););(node["leisure"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f}););(node["office"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f}););(node["shop"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f}););(node["sport"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f}););(node["building"]({south:.8f},'
				'{west:.8f},{north:.8f},{east:.8f});));out;')
			query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all POIs footprints data within bounding box from '
			   'API in {:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(geometry), time.time()-start_time))

	elif by_poly:
		# project to utm, divide polygon up into sub-polygons if area exceeds a
		# max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
		geometry_proj, crs_proj = project_geometry(polygon)
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		polygon_coord_strs = get_polygons_coordinates(geometry)
		log('Requesting POIs footprints data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
		start_time = time.time()

		# pass each polygon exterior coordinates in the list to the API, one at
		# a time
		for polygon_coord_str in polygon_coord_strs:
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};('
				'(node["amenity"](poly:"{polygon}"););'
				'(node["leisure"](poly:"{polygon}"););'
				'(node["office"](poly:"{polygon}"););'
				'(node["shop"](poly:"{polygon}"););'
				'(node["sport"](poly:"{polygon}"););'
				'(node["building"](poly:"{polygon}");));out;')
			query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all POIs footprints data within polygon from API in '
			   '{:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(polygon_coord_strs), time.time()-start_time))

	return response_jsons

def create_pois_gdf(date="", polygon=None, north=None, south=None, east=None,
						 west=None, retain_invalid=False):
	"""
	Get POIs footprint data from OSM then assemble it into a GeoDataFrame.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the POIs footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	retain_invalid : bool
		if False discard any POIs footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	responses = osm_pois_download(date, polygon, north, south, east, west)

	vertices = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='node':

				point = Point( result['lon'], result['lat'] )

				POI = {'geometry' : point}

				if 'tags' in result:
					for tag in result['tags']:
						POI[tag] = result['tags'][tag]

				vertices[result['id']] = POI

	gdf = gpd.GeoDataFrame(vertices).T
	gdf.crs = {'init':'epsg:4326'}

	if not retain_invalid:
		try:
			# drop all invalid geometries
			gdf = gdf[gdf['geometry'].is_valid]
		except: # Empty data frame
			# Create a one-row data frame with null information (avoid later Spatial-Join crash)
			if (polygon is not None): # Polygon given
				point = polygon.centroid
			else: # Bounding box
				point = Point( (east+west)/2. , (north+south)/2. )
			data = {"geometry":[point], "osm_id":[0]}
			gdf = gpd.GeoDataFrame(data, crs={'init': 'epsg:4326'})

	return gdf

#######################################################################
### OSM Building parts
#######################################################################

def osm_bldg_part_download(date="", polygon=None, north=None, south=None, east=None, west=None,
					  timeout=180, memory=None, max_query_area_size=50*1000*50*1000):
	"""
	Download OpenStreetMap building parts footprint data.
	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the building footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	timeout : int
		the timeout interval for requests and to pass to API
	memory : int
		server memory allocation size for the query, in bytes. If none, server
		will use its default allocation size
	max_query_area_size : float
		max area for any part of the geometry, in the units the geometry is in:
		any polygon bigger will get divided up for multiple queries to API
		(default is 50,000 * 50,000 units (ie, 50km x 50km in area, if units are
		meters))
	Returns
	-------
	list
		list of response_json dicts
	"""

	# check if we're querying by polygon or by bounding box based on which
	# argument(s) where passed into this function
	by_poly = polygon is not None
	by_bbox = not (north is None or south is None or east is None or west is None)
	if not (by_poly or by_bbox):
		raise ValueError('You must pass a polygon or north, south, east, and west')

	response_jsons = []

	# pass server memory allocation in bytes for the query to the API
	# if None, pass nothing so the server will use its default allocation size
	# otherwise, define the query's maxsize parameter value as whatever the
	# caller passed in
	if memory is None:
		maxsize = ''
	else:
		maxsize = '[maxsize:{}]'.format(memory)

	# define the query to send the API
	if by_bbox:
		# turn bbox into a polygon and project to local UTM
		polygon = Polygon([(west, south), (east, south), (east, north), (west, north)])
		geometry_proj, crs_proj = project_geometry(polygon)

		# subdivide it if it exceeds the max area size (in meters), then project
		# back to lat-long
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		log('Requesting building part footprints data within bounding box from API in {:,} request(s)'.format(len(geometry)))
		start_time = time.time()

		# loop through each polygon rectangle in the geometry (there will only
		# be one if original bbox didn't exceed max area size)
		for poly in geometry:
			# represent bbox as south,west,north,east and round lat-longs to 8
			# decimal places (ie, within 1 mm) so URL strings aren't different
			# due to float rounding issues (for consistent caching)
			west, south, east, north = poly.bounds
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};((way["building:part"]({south:.8f},'
							  '{west:.8f},{north:.8f},{east:.8f});(._;>;););(relation["building:part"]'
							  '({south:.8f},{west:.8f},{north:.8f},{east:.8f});(._;>;);));out;')
			query_str = query_template.format(north=north, south=south, east=east, west=west, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all building part footprints data within bounding box from '
			   'API in {:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(geometry), time.time()-start_time))

	elif by_poly:
		# project to utm, divide polygon up into sub-polygons if area exceeds a
		# max size (in meters), project back to lat-long, then get a list of polygon(s) exterior coordinates
		geometry_proj, crs_proj = project_geometry(polygon)
		geometry_proj_consolidated_subdivided = consolidate_subdivide_geometry(geometry_proj, max_query_area_size=max_query_area_size)
		geometry, _ = project_geometry(geometry_proj_consolidated_subdivided, crs=crs_proj, to_latlong=True)
		polygon_coord_strs = get_polygons_coordinates(geometry)
		log('Requesting building part footprints data within polygon from API in {:,} request(s)'.format(len(polygon_coord_strs)))
		start_time = time.time()

		# pass each polygon exterior coordinates in the list to the API, one at
		# a time
		for polygon_coord_str in polygon_coord_strs:
			query_template = (date+'[out:json][timeout:{timeout}]{maxsize};(way'
							  '(poly:"{polygon}")["building:part"];(._;>;);relation'
							  '(poly:"{polygon}")["building:part"];(._;>;));out;')
			query_str = query_template.format(polygon=polygon_coord_str, timeout=timeout, maxsize=maxsize)
			response_json = overpass_request(data={'data':query_str}, timeout=timeout)
			response_jsons.append(response_json)
		msg = ('Got all building part footprints data within polygon from API in '
			   '{:,} request(s) and {:,.2f} seconds')
		log(msg.format(len(polygon_coord_strs), time.time()-start_time))

	return response_jsons



def create_building_parts_gdf(date="", polygon=None, north=None, south=None, east=None,
						 west=None, retain_invalid=False):
	"""
	Get building footprint data from OSM then assemble it into a GeoDataFrame.
	If no building parts are retrieved, a default (null-data) point located at the centroid of the region of interest is created

	Parameters
	----------
	date : string
		query the database at a certain timestamp
	polygon : shapely Polygon or MultiPolygon
		geographic shape to fetch the building footprints within
	north : float
		northern latitude of bounding box
	south : float
		southern latitude of bounding box
	east : float
		eastern longitude of bounding box
	west : float
		western longitude of bounding box
	retain_invalid : bool
		if False discard any building footprints with an invalid geometry
	Returns
	-------
	GeoDataFrame
	"""

	responses = osm_bldg_part_download(date, polygon, north, south, east, west)

	vertices = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='node':
				vertices[result['id']] = {'lat' : result['lat'],
										  'lon' : result['lon']}

	buildings = {}
	for response in responses:
		for result in response['elements']:
			if 'type' in result and result['type']=='way':
				nodes = result['nodes']
				try:
					polygon = Polygon([(vertices[node]['lon'], vertices[node]['lat']) for node in nodes])
				except Exception:
					log('Polygon has invalid geometry: {}'.format(nodes))
				building = {'nodes' : nodes,
							'geometry' : polygon}

				if 'tags' in result:
					for tag in result['tags']:
						building[tag] = result['tags'][tag]

				buildings[result['id']] = building

	gdf = gpd.GeoDataFrame(buildings).T
	gdf.crs = {'init':'epsg:4326'}

	if not retain_invalid:
		try:
			# drop all invalid geometries
			gdf = gdf[gdf['geometry'].is_valid]
		except: # Empty data frame
			# Create a one-row data frame with null information (avoid later Spatial-Join crash)
			if (polygon is not None): # Polygon given
				point = polygon.centroid
			else: # Bounding box
				point = Point( (east+west)/2. , (north+south)/2. )
			# Data as records
			data = {"geometry":[point], "osm_id":[0], "building:part":["yes"], "height":[""]}
			gdf = gpd.GeoDataFrame(data, crs={'init': 'epsg:4326'})

	return gdf
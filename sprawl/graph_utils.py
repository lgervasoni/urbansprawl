###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

import numpy as np
import pandas as pd
import osmnx as ox
from shapely.geometry import LineString

def graph_bounding_box(G):
	""" 
	Compute the graph bounding box

	Parameters
	----------
	G : networkx multidigraph
		input graph

	Returns
	----------
	list
		values for xmin, xmax, ymin, ymax
	"""
	# Initialize values
	xmin, xmax = G.nodes(data=True)[0][1]["x"], G.nodes(data=True)[0][1]["x"]
	ymin, ymax = G.nodes(data=True)[0][1]["y"], G.nodes(data=True)[0][1]["y"]
	# Iterate
	for node in G.nodes(data=True):
		xmin = min(xmin, node[1]["x"])
		xmax = max(xmax, node[1]["x"])
		ymin = min(ymin, node[1]["y"])
		ymax = max(ymax, node[1]["y"])
	return xmin, xmax, ymin, ymax

def get_nearest_node_utm(G, point, return_dist=False):
	"""
	Return the nearest graph node to some specified point in UTM coordinates
	
	Parameters
	----------
	G : networkx multidigraph
		input graph
	point : tuple
		the (x, y) point for which we will find the nearest node in the graph
	return_dist : bool
		optionally also return the distance between the point and the nearest node
	
	Returns
	-------
	int or tuple
		corresponding node or tuple (int node, float distance)
	""" 
	# dump graph node coordinates into a pandas dataframe indexed by node id with x and y columns
	coords = np.array([[node, data['x'], data['y']] for node, data in G.nodes(data=True)])
	df = pd.DataFrame(coords, columns=['node', 'x', 'y']).set_index('node')
	# Point coordinates
	p_x, p_y = point
	distances = df.apply(lambda x: np.sqrt( ( x.x - p_x)**2 + ( x.y- p_y)**2), axis=1 )
	
	# nearest node's ID is the index label of the minimum distance
	nearest_node = distances.idxmin()
	
	# if caller requested return_dist, return distance between the point and the nearest node as well
	if return_dist:
		return int(nearest_node), distances.loc[nearest_node]
	else:
		return int(nearest_node)

def cut_in_two(line):
	"""
	Cuts input line into two lines of equal length

	Parameters
	----------
	line : shapely.LineString
		input line

	Returns
	----------
	list (LineString, LineString, Point)
		two lines and the middle point cutting input line
	"""
	from shapely.geometry import Point, LineString
	# Get final distance value
	distance = line.length/2
	# Cuts a line in two at a distance from its starting point
	if distance <= 0.0 or distance >= line.length:
		return [LineString(line)]
	coords = list(line.coords)
	for i, p in enumerate(coords):
		pd = line.project(Point(p))
		if pd == distance:
			return [LineString(coords[:i+1]), LineString(coords[i:]), pd]
		if pd > distance:
			cp = line.interpolate(distance)
			return [ LineString(coords[:i] + [(cp.x, cp.y)]), LineString([(cp.x, cp.y)] + coords[i:]), cp]

class NodeCounter:
	"""
	Node negative counter. Utils for node osmid creation. Start on -1 and it auto decrements
	"""
	def __init__(self):
		self._num = 0
	def get_num(self):
		self._num -= 1
		return self._num

def verify_divide_edge(G, u, v, key, data, node_creation_counter, MAX_EDGE_LENGTH):
	"""
	Verify if edge(u,v)[key] length is higher than a certain threshold
	In this case, divide edge(u,v) in two edges of equal length
	Assign negative values to the edges new osm id
	Call recursively to continue dividing each of the lines if necessary

	Parameters
	----------
	G : networkx multidigraph
		input graph
	u : node
		origin node
	v : node
		destination node
	key : int
		(u,v) arc identifier
	data : dict
		arc data
	node_creation_counter : NodeCounter
		node identifier creation
	MAX_EDGE_LENGTH : float
		maximum tolerated edge length

	Returns
	----------

	"""
	# Input: Two communidated nodes (u, v)
	if ( data["length"] <= MAX_EDGE_LENGTH ): # Already satisfy condition?
		return
	
	# Get geometry connecting (u,v)
	if ( data.get("geometry",None) ): # Geometry exists
		geometry = data["geometry"]
	else: # Real geometry is a straight line between the two nodes
		P_U = G.node[u]["x"], G.node[u]["y"]
		P_V = G.node[v]["x"], G.node[v]["y"]
		geometry = LineString( (P_U, P_V) )
	
	# Get geometries for edge(u,middle), edge(middle,v) and node(middle)
	line1, line2, middle_point = cut_in_two(geometry)
	
	# Copy edge(u,v) data to conserve attributes. Modify its length
	data_e1 = data.copy()
	data_e2 = data.copy()
	# Associate correct length
	data_e1["length"] = line1.length
	data_e2["length"] = line2.length
	# Assign geometries
	data_e1["geometry"] = line1
	data_e2["geometry"] = line2
	
	# Create new node: Middle distance of edge
	x,y = list(middle_point.coords)[0]
	# Set a new unique osmid: Negative (coherent with osm2pgsql, created objects contain negative osmid)
	node_osmid = node_creation_counter.get_num()
	node_data = {'osmid':node_osmid, 'x':x, 'y':y}
	
	# Add middle node with its corresponding data
	G.add_node(node_osmid,node_data)
	
	# Add edges (u,middle) and (middle,v)
	G.add_edge(u, node_osmid, attr_dict=data_e1)
	G.add_edge(node_osmid, v, attr_dict=data_e2)
	
	# Remove edge (u,v)
	G.remove_edge(u,v,key=key)
	
	# Recursively verify created edges and divide if necessary. Use last added key to identify the edge
	last_key = len( G[u][node_osmid] ) -1
	verify_divide_edge(G, u, node_osmid, last_key, data_e1, node_creation_counter, MAX_EDGE_LENGTH)
	last_key = len( G[node_osmid][v] ) -1
	verify_divide_edge(G, node_osmid, v, last_key, data_e2, node_creation_counter, MAX_EDGE_LENGTH)
	

def divide_long_edges_graph(G, MAX_EDGE_LENGTH):
	"""
	Divide all edges with a higher length than input threshold by means of dividing the arcs and creating new nodes

	Parameters
	----------
	G : networkx multidigraph
		input graph
	MAX_EDGE_LENGTH : float
		maximum tolerated edge length

	Returns
	----------

	"""
	# Negative osm_id indicate created nodes
	node_creation_counter = NodeCounter()
	# Iterate on edges
	for u, v, key, data in G.edges(data=True, keys=True):
		if ( data["length"] > MAX_EDGE_LENGTH ):
			# Divide the edge (u,v) recursively
			verify_divide_edge(G, u, v, key, data, node_creation_counter, MAX_EDGE_LENGTH)

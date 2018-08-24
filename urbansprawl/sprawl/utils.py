###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import numpy as np
import pandas as pd
import networkx as nx
import math
from shapely.geometry import LineString
from scipy.spatial.distance import cdist


def WeightedKernelDensityEstimation(X, Weights, bandwidth, Y, max_mb_per_chunk = 1000):
    """ 
    Computes a Weighted Kernel Density Estimation

    Parameters
    ----------
    X : array
        input points
    Weights : array
        array of weights associated to points
    bandwidth : float
        bandwidth for kernel density estimation
    Y : array
        points where density estimations will be performed

    Returns
    ----------
    pd.Series
        returns an array of the estimated densities rescaled between [0;1]
    """
    def get_megabytes_pairwise_distances_allocation(X, Y):
    	# Calculate MB needed to allocate pairwise distances
    	return len(X) * len(Y) * 8 * 1e-6
    
    # During this procedure, pairwise euclidean distances are computed between inputs points X and points to estimate Y
    # For this reason, Y is divided in chunks to avoid big memory allocations. At most, X megabytes per chunk are allocated for pairwise distances
    Y_split = np.array_split( Y, math.ceil( get_megabytes_pairwise_distances_allocation(X,Y) / max_mb_per_chunk ) )
    
    """
    ### Step by step
    # Weighed KDE: Sum{ Weight_i * K( (X-Xi) / h) }
    W_norm = np.array( Weights / np.sum(Weights) )
    cdist_values = cdist( Y, X, 'euclidean') / bandwidth
    Ks = np.exp( -.5 * ( cdist_values ) ** 2  )
    PDF = np.sum( Ks * W_norm, axis=1)
    """
    """
    ### Complete version. Memory consuming
    PDF = np.sum( np.exp( -.5 * ( cdist( Y, X, 'euclidean') / bandwidth ) ** 2  ) * ( np.array( Weights / np.sum(Weights) ) ), axis=1)
    """

    ### Divide Y in chunks to avoid big memory allocations
    PDF = np.concatenate( [ np.sum( np.exp( -.5 * ( cdist( Y_i, X, 'euclidean') / bandwidth ) ** 2  ) * ( np.array( Weights / np.sum(Weights) ) ), axis=1) for Y_i in Y_split ] )
    # Rescale
    return pd.Series( PDF / PDF.sum() )


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

def verify_divide_edge(G, u, v, key, data, node_creation_counter, max_edge_length):
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
	max_edge_length : float
		maximum tolerated edge length

	Returns
	----------

	"""
	# Input: Two communidated nodes (u, v)
	if ( data["length"] <= max_edge_length ): # Already satisfy condition?
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
	# Set a new unique osmid: Negative (as in OSM2PGSQL, created objects contain negative osmid)
	node_osmid = node_creation_counter.get_num()
	node_data = {'osmid':node_osmid, 'x':x, 'y':y}
	
	# Add middle node with its corresponding data
	G.add_node(node_osmid)
	nx.set_node_attributes(G, {node_osmid : node_data } )
	
	# Add edges (u,middle) and (middle,v)
	G.add_edge(u, node_osmid)
	nx.set_edge_attributes(G, { (u, node_osmid, 0): data_e1 } )
	G.add_edge(node_osmid, v)
	nx.set_edge_attributes(G, { (node_osmid, v, 0): data_e2 } )
	
	# Remove edge (u,v)
	G.remove_edge(u,v,key=key)
	
	# Recursively verify created edges and divide if necessary. Use last added key to identify the edge
	last_key = len( G[u][node_osmid] ) -1
	verify_divide_edge(G, u, node_osmid, last_key, data_e1, node_creation_counter, max_edge_length)
	last_key = len( G[node_osmid][v] ) -1
	verify_divide_edge(G, node_osmid, v, last_key, data_e2, node_creation_counter, max_edge_length)
	

def divide_long_edges_graph(G, max_edge_length):
	"""
	Divide all edges with a higher length than input threshold by means of dividing the arcs and creating new nodes

	Parameters
	----------
	G : networkx multidigraph
		input graph
	max_edge_length : float
		maximum tolerated edge length

	Returns
	----------

	"""
	# Negative osm_id indicate created nodes
	node_creation_counter = NodeCounter()

	for u, v, key, data in list( G.edges(data=True, keys=True) ):
		if ( data["length"] > max_edge_length ):
			# Divide the edge (u,v) recursively
			verify_divide_edge(G, u, v, key, data, node_creation_counter, max_edge_length)
###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

from scipy import spatial
import numpy as np
import pandas as pd
import networkx as nx
import osmnx as ox
from bisect import bisect
import multiprocessing
import time

from .graph_utils import get_nearest_node_utm, divide_long_edges_graph
from .utils import _create_grid, log

##############################################################
### Parameters
##############################################################
# Maximum length, in meters, to tolerate an edge in a graph (otherwise, divide edge)
MAX_EDGE_LENGTH = 150
# Maximum distance tolerated from input point to closest graph node in order to calculate accessibility values
MAX_NODE_DISTANCE = 250
# Count number of activities in fixed distance?
USE_COUNT_ACTIVITIES_FIXED_DISTANCE = True

#### Fixed distance 
# Maximum distance tolerated (cut&branch) when searching for the activities
MAX_DISTANCE_TO_TRAVEL = 1000
# Cut iteration if the number of activites exceeds a threshold
MAX_NUM_ACTIVITIES_TRAVEL = 250

#### Fixed number of opportunities

# FIXED_NUMBER_ACTIVITIES cost: Minimum number of activities required
ACTIVITIES_MIN_NUMBER = 10

#### Parallelization parameters
# By default, number of cores
NUM_PROCESSES = multiprocessing.cpu_count()
# Number of chunks when dividing data for parallelization
NUM_CHUNKS = NUM_PROCESSES * 2
# Use parallel version?
USE_PARALLEL = True

##############################################################
### Accessibility indices method
##############################################################

def get_count_activities_fixed_distance(G, N0):
	""" 
	Accessibility indices metric 
	Based on counting the number of (activity) opportunities given a fixed maximum distance to travel

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	N0 : int
		initial node id

	Returns
	----------
	int
		returns the number of reached activities
	"""
	# Initialize data structures
	visited_nodes = []
	neighboring_nodes_id = []
	neighboring_nodes_cost = []
	num_activities_travelled = 0

	N_visit = N0

	# Pre-compute the shortest path length from source node N0 to other nodes; using lengths of roads as weight
	shortest_path_length_N0_ = nx.single_source_dijkstra_path_length(G, source=N0, cutoff=MAX_DISTANCE_TO_TRAVEL, weight="length")

	while ( True ):
		# Store visited node
		visited_nodes.append(N_visit)

		# Update travelled activities
		num_activities_travelled += G.node[N_visit]["num_activities"]

		# Reached sufficient number of activities
		if ( num_activities_travelled >= MAX_NUM_ACTIVITIES_TRAVEL ): return MAX_NUM_ACTIVITIES_TRAVEL
		
		# Add to neighboring_nodes the neighbors of visited node
		for N_i in G.neighbors(N_visit):
			if ( (not N_i in neighboring_nodes_id) and (not N_i in visited_nodes) ): # Not stored/visited already
				# Store neighboring nodes, ordered by their distance cost
				cost = shortest_path_length_N0_.get(N_i)

				if (cost): # If path within MAX_DISTANCE_TO_TRAVEL exists, add neighboring node
					idx_to_insert = bisect( neighboring_nodes_cost, cost )
					# Insert in ordered list
					neighboring_nodes_id.insert(idx_to_insert, N_i)
					neighboring_nodes_cost.insert(idx_to_insert, cost)
		
		# If neighborings nodes exist: Continue iteration
		# If no neighboring nodes: Reached maximum distance tolerated, cut the iteration
		if (neighboring_nodes_id): # If not empty
			# Update next node to visit
			N_visit = neighboring_nodes_id.pop(0)
			# Pop cost associated to N_visit
			neighboring_nodes_cost.pop(0)
		else: # Empty neighbors
			return num_activities_travelled
	
	return np.nan

def get_minimum_cost_activities_travel(G, N0):
	""" 
	Accessibility indices metric 
	Based on the minimum radius travel cost to accomplish a certain quantity of activities

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	N0 : int
		initial node id

	Returns
	----------
	float
		returns the computed radius cost length
	"""
	# Initialize data structures
	visited_nodes = []
	neighboring_nodes_id = []
	neighboring_nodes_cost = []
	activities_travelled = 0
	
	N_visit = N0

	while ( not activities_travelled >= ACTIVITIES_MIN_NUMBER ):
		# Store visited node
		visited_nodes.append(N_visit)

		# Update travelled activities
		activities_travelled += G.node[N_visit]["num_activities"]

		# Add to neighboring_nodes the neighbors of visited node
		for N_i in G.neighbors(N_visit):
			if ( (not N_i in neighboring_nodes_id) and (not N_i in visited_nodes) ): # Not stored/visited already
				# Store neighboring nodes, ordered by their distance cost
				cost = nx.shortest_path_length(G,N0,N_i,weight="length")
				idx_to_insert = bisect( neighboring_nodes_cost, cost )
				# Insert in ordered list
				neighboring_nodes_id.insert(idx_to_insert, N_i)
				neighboring_nodes_cost.insert(idx_to_insert, cost)
		
		if (neighboring_nodes_id): # If not empty
			# Update next node to visit
			N_visit = neighboring_nodes_id.pop(0)
			cost_travel = neighboring_nodes_cost.pop(0)

			# Reached maximum distance tolerated. Cut iteration
			if (cost_travel > MAX_DISTANCE_TO_TRAVEL):
				return MAX_DISTANCE_TO_TRAVEL
		else: # Empty neighbors
			return np.nan
	
	# Accomplished. End node: visited_nodes[-1]
	return nx.shortest_path_length(G,N0,visited_nodes[-1],weight="length")

#############

if (USE_COUNT_ACTIVITIES_FIXED_DISTANCE):
	get_accessibility_ = get_count_activities_fixed_distance
else:
	get_accessibility_ = get_minimum_cost_activities_travel

def _calculate_accessibility(G, point_ref ):
	""" 
	Calculate accessibility value at point_ref according to chosen metric
	If no graph node exists nearby input point reference, NaN is set

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	point_ref: shapely.Point
		reference point to calculate accesisibility

	Returns
	----------
	int/float
		resulting accessibility value for input point reference
	"""
	# Find closest node to point_ref
	N0, distance = get_nearest_node_utm(G, point_ref, return_dist=True)
	# Distance to closest node too high?
	if (distance > MAX_NODE_DISTANCE): return np.nan
	return get_accessibility_(G,N0)

##############################################################
### Compute accessibility grid
##############################################################

def compute_grid_accessibility(XX_YY, G, df, kw_args={'max_edge_length':150}, extra_args=[]):
	""" 
	Calculate accessibility values at point_ref

	Parameters
	----------
	XX_YY : pandas.Panel
		meshgrid with (x,y) reference points to calculate indices
	G : networkx multidigraph
		input graph to calculate accessibility
	df : pandas.DataFrame
		data with geo-referenced activity land uses
	kw_args: dict
		additional keyword arguments for the indices calculation
	extra_args : array
		additional arguments for the indices calculation

	Returns
	----------
	int
		number of activities found within a radius distance using the street network
	"""
	log("Accessibility calculation")
	start = time.time()

	# Create grid
	grid = _create_grid( XX_YY )

	# Divide long edges
	divide_long_edges_graph(G, kw_args.get('max_edge_length',MAX_EDGE_LENGTH) )
	# Get rows with defined activity type
	df_activities = df.loc[ df.activity_type.notnull() ]	
	# Associate them to its closest node in the graph
	associate_activities_closest_node(G,df_activities)

	# Iterate
	rows, cols = grid.shape[0] , grid.shape[1]

	log("Number of points in grid:"+str(rows*cols) )

	# Allocate tuple
	grid_tuple_pt_ref = np.ndarray(grid.shape,dtype=tuple)
	# Initialize values
	for i in range(rows): #Rows
		for j in range(cols): #Cols
			# (x,y) tuples
			grid_tuple_pt_ref[i,j] = ( XX_YY["x",i,j], XX_YY["y",i,j] )

	grid_accessibility = perform_parallel_accessibility(G,grid_tuple_pt_ref)

	end = time.time()
	log("Accessibility calculation time: "+str(end-start))

	return pd.DataFrame(grid_accessibility)

##############################################################
### Parallel implementation
##############################################################

def _calculate_accessibility_parallel(G, grid_acc_pt_ref ):
	""" 
	Calculate accessibility values for input chunk

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	grid_acc_pt_ref : numpy.array
		chunk of reference points to calculate accessibility indices

	Returns
	----------
	array
		returns accessibility indices for input chunk
	"""
	#_calculate_accessibility_vectorized = np.vectorize(_calculate_accessibility)
	#return _calculate_accessibility_vectorized( G, grid_acc_pt_ref )
	log("+")
	return [ _calculate_accessibility(G,pt_ref) for pt_ref in grid_acc_pt_ref ]

def perform_parallel_accessibility(G, grid_tuples_pt_ref):
	""" 
	Calculate accessibility indices
	Multi-processing implementation, dividing the grid to process in chunks

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	grid_tuples_pt_ref : numpy.ndarray
		input array with reference point (x,y) tuples

	Returns
	----------
	numpy.matrix
		returns accessibility indices
	"""
	log("Parallel version")

	shape_ = grid_tuples_pt_ref.shape
	# Chunk
	grid_tuple_chunks = np.array_split(grid_tuples_pt_ref.ravel(), NUM_CHUNKS)

	log("Number of elements to process: "+str(shape_[0]*shape_[1]) )
	log("Dividing in "+str(NUM_CHUNKS)+" chunks")
	log("Using "+str(NUM_PROCESSES)+" processes")

	# Initialize parallel 
	# Multi-processing
	from multiprocessing import Pool
	# Multi-threading
	#from multiprocessing.pool import ThreadPool as Pool

	from functools import partial
	import multiprocessing

	pool = Pool(processes=NUM_PROCESSES)
	func = partial(_calculate_accessibility_parallel, G)
	result = pool.map(func, grid_tuple_chunks)
	#result = pool.imap_unordered(func, grid_tuple_chunks)
	pool.close()
	pool.join()

	# Concatenate back again and reshape
	grid_result = np.concatenate( result ).reshape( shape_ )	
	return grid_result

###############
# Utils
###############

def associate_activities_closest_node(G, df_activities ):
	""" 
	Associates the number of existing activities to their closest nodes in the graph

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	df_activities : pandas.DataFrame
		data selection of activity land uses

	Returns
	----------

	"""
	# Initialize number of activity references
	for u, data in G.nodes(data=True):
		data["num_activities"] = 0
	
	# Initialize KDTree of graph nodes
	coords = np.array([[node, data['x'], data['y']] for node, data in G.nodes(data=True)])
	df_nodes = pd.DataFrame(coords, columns=['node', 'x', 'y'])
	# zip coordinates
	data = list(zip( df_nodes["x"].ravel(), df_nodes["y"].ravel() ))
	# Create input tree
	tree = spatial.KDTree( data )
	
	for idx_activity,row in df_activities.iterrows():
		# Get the centroid of the geometry
		P = row.geometry.centroid
		# Get the closest node and its distance
		distance, idx_node = tree.query( (P.x,P.y) )
		node = df_nodes.loc[idx_node, "node"]

		# Append one more activity
		G.node[node]["num_activities"] += 1
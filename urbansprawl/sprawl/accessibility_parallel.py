###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import sys
import numpy as np
import pandas as pd
import networkx as nx
from bisect import bisect
import time

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
	#print(df.loc[124550])
	distances = df.apply(lambda x: np.sqrt( ( x.x - p_x)**2 + ( x.y- p_y)**2), axis=1 )
	
	# nearest node's ID is the index label of the minimum distance
	nearest_node = distances.idxmin()
	
	# if caller requested return_dist, return distance between the point and the nearest node as well
	if return_dist:
		return int(nearest_node), distances.loc[nearest_node]
	else:
		return int(nearest_node)

##############################################
### Accessibility indices calculation
##############################################

def get_count_activities_fixed_distance(G, point_ref, arguments):
	""" 
	Calculate accessibility value at point_ref according to chosen metric
	If no graph node exists nearby input point reference, NaN is set
	Based on counting the number of (activity) opportunities given a fixed maximum distance to travel

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	point_ref: shapely.Point
		reference point to calculate accesisibility

	Returns
	----------
	int
		returns the number of reached activities
	"""
	# Find closest node to point_ref
	N0, distance = get_nearest_node_utm(G, point_ref.coords[0], return_dist=True)
	# Distance to closest node too high?
	if (distance > arguments.max_node_distance): return np.nan

	# Initialize data structures
	visited_nodes = set()
	neighboring_nodes_id = []
	neighboring_nodes_cost = []
	num_activities_travelled = 0
	N_visit = N0

	# Pre-compute the shortest path length from source node N0 to other nodes; using lengths of roads as weight
	shortest_path_length_N0_ = nx.single_source_dijkstra_path_length(G, source=N0, cutoff=arguments.fixed_distance_max_travel_distance, weight="length")

	while ( True ):
		# Store visited node
		visited_nodes.add(N_visit)

		# Update travelled activities
		num_activities_travelled += G.node[N_visit]["num_activities"]

		# Reached sufficient number of activities
		if ( num_activities_travelled >= arguments.fixed_distance_max_num_activities ): return arguments.fixed_distance_max_num_activities
		
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
		
		if (neighboring_nodes_id): # If neighborings nodes exist: Continue iteration
			# Update next node to visit
			N_visit = neighboring_nodes_id.pop(0)
			# Pop cost associated to N_visit
			neighboring_nodes_cost.pop(0)
		else: # If no neighboring nodes: Reached maximum distance tolerated, cut the iteration
			return num_activities_travelled
	
	return np.nan



def get_minimum_cost_activities_travel(G, point_ref, arguments):
	""" 
	Calculate accessibility value at point_ref according to chosen metric
	If no graph node exists nearby input point reference, NaN is set
	Based on the minimum radius travel cost to accomplish a certain quantity of activities

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	point_ref: shapely.Point
		reference point to calculate accesisibility

	Returns
	----------
	float
		returns the computed radius cost length
	"""
	# Find closest node to point_ref
	N0, distance = get_nearest_node_utm(G, point_ref.coords[0], return_dist=True)
	# Distance to closest node too high?
	if (distance > arguments.max_node_distance): return np.nan

	# Initialize data structures
	visited_nodes = []
	neighboring_nodes_id = []
	neighboring_nodes_cost = []
	activities_travelled = 0
	
	N_visit = N0

	while ( not activities_travelled >= arguments.fixed_activities_min_number ):
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
			if (cost_travel > arguments.fixed_activities_max_travel_distance):
				return arguments.fixed_activities_max_travel_distance
		else: # Empty neighbors
			return np.nan
	
	# Accomplished. End node: visited_nodes[-1]
	return nx.shortest_path_length(G,N0,visited_nodes[-1],weight="length")


def main(argv):
	""" 
	Main program to drive the accessibility indices calculation

	Parameters
	----------
	argv : array
		arguments to drive the calculation

	Returns
	----------

	"""
	start = time.time()

	# Load graph
	G = nx.read_gpickle(argv[1])

	# Load indices points
	indices = pd.read_pickle( argv[2] )

	# Load indices calculation arguments
	arguments = pd.read_pickle( argv[3] )

	if ( ( len(argv) > 4 ) and (argv[4] == "memory_test" ) ): # Test memory used for current subprocess
		import os
		import psutil
		process = psutil.Process(os.getpid())
		Allocated_process_MB = process.memory_info().rss / 1000 / 1000
		Free_system_MB = psutil.virtual_memory().available / 1000 / 1000
		#print( MB_total )
		#print( Free_system_MB )
		max_processes = int( Free_system_MB / Allocated_process_MB )
		print(max_processes)
		return


	if (arguments.fixed_activities):
		_calculate_accessibility = get_minimum_cost_activities_travel
	elif (arguments.fixed_distance):
		_calculate_accessibility = get_count_activities_fixed_distance
	else:
		assert(False)

	# Calculate accessibility
	indices["accessibility"] = indices.geometry.apply(lambda x: _calculate_accessibility(G, x, arguments) )

	# Store results
	indices.to_pickle( argv[2].replace('points','indices') )

	end = time.time()
	print( "Time:",str(end-start) )


if __name__ == "__main__":
	main(sys.argv)
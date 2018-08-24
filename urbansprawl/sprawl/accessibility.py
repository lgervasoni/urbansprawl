###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

from scipy import spatial
import numpy as np
import pandas as pd
import geopandas as gpd
import networkx as nx
import time
import os
import subprocess
import shutil

from multiprocessing import cpu_count

from osmnx.utils import log
from .utils import divide_long_edges_graph

##############################################################
### Compute accessibility grid
##############################################################

def compute_grid_accessibility(df_indices, G, df_osm_built, df_osm_pois, kw_args={'fixed_distance':True,'fixed_activities':False,'max_edge_length':200,'max_node_distance':250,
				'fixed_distance_max_travel_distance':2000, 'fixed_distance_max_num_activities':250, 'fixed_activities_min_number': 20, 'fixed_activities_max_travel_distance':5000} ):
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
			fixed_distance : bool
				denotes the cumulative opportunities access to activity land uses given a fixed maximum distance to travel
			fixed_activities : bool
				represents the distance needed to travel in order to reach a certain number of activity land uses
			max_edge_length: int
				maximum length, in meters, to tolerate an edge in a graph (otherwise, divide edge)
			max_node_distance: int
				maximum distance tolerated from input point to closest graph node in order to calculate accessibility values
			fixed_distance_max_travel_distance: int
				(fixed distance) maximum distance tolerated (cut&branch) when searching for the activities
			fixed_distance_max_num_activities: int
				(fixed distance) cut iteration if the number of activites exceeds a threshold
			fixed_activities_min_number: int
				(fixed activities) minimum number of activities required
			fixed_activities_max_travel_distance : int
				(fixed activities) maximum distance tolerated (cut&branch) when searching for the activities


	Returns
	----------
	int
		number of activities found within a radius distance using the street network
	"""
	log("Accessibility calculation")
	start = time.time()

	# Assert that only one option is set
	assert( kw_args["fixed_distance"] ^ kw_args["fixed_activities"] )

	# Arguments to pandas.Series
	kw_arguments = pd.Series(kw_args)

	##############
	### Prepare input data for indices calculation in parallel call
	##############
	# Temporary folder to pickle data
	if ( not os.path.exists("temp") ): os.makedirs("temp")
	# Number of CPU cores on your system
	num_cores = cpu_count() 
	# Prepare input data: As many chunks of data as cores
	prepare_data(G, df_osm_built, df_osm_pois, df_indices, num_cores, kw_arguments )

	#This command could have multiple commands separated by a new line \n
	parallel_code = os.path.realpath(__file__).replace(".py","_parallel.py")
	command_call = "python " + parallel_code + " temp/graph.gpickle temp/points_NUM_CHUNK.pkl temp/arguments.pkl"
	
	##############
	### Verify amount of memory used per subprocess
	##############
	p = subprocess.Popen(command_call.replace("NUM_CHUNK",str(0)) + " memory_test", stdout=subprocess.PIPE, shell=True)
	output, err = p.communicate()
	p.wait()

	# Max number of subprocess allocations given its memory consumption
	numbers = [ numb for numb in str(output) if numb in ["0","1","2","3","4","5","6","7","8","9"] ]
	max_processes = int( ''.join(numbers) )
	log("Max processes to allocate (consdering memory availability):"+str(max_processes) )
	log("Numbero of available cores:"+str(num_cores) )

	##############
	### Set chunks to run in parallel: If more core than allowed processes, divide chunks to run at most X processes
	##############
	if (num_cores > max_processes):	# Run parallel-chunks at a splitted pace, to avoid memory swap		
		chunks_run = np.array_split( list( range(num_cores) ), max_processes )
	else: # Run all chunks in parallel
		chunks_run = [ list( range(num_cores) ) ]


	# Parallel implementation
	for chunk in chunks_run: # Run full chunk
		Ps_i = []
		for i in chunk: # Run each index
			p = subprocess.Popen(command_call.replace("NUM_CHUNK",str(i)), stdout=subprocess.PIPE, shell=True)
			Ps_i.append( p )
		
		# Get the output
		Output_errs = [ p.communicate() for p in Ps_i ]
		
		# This makes the wait possible
		Ps_status = [ p.wait() for p in Ps_i ]

		# Output for chunk
		for output, err in Output_errs:
			log ( str(output) )

	# Associate data by getting the chunk results concatenated
	index_column = "accessibility"
	df_indices[index_column] = pd.concat( [ pd.read_pickle("temp/indices_NUM_CHUNK.pkl".replace("NUM_CHUNK",str(i)) ) for i in range(num_cores) ], ignore_index=True ).accessibility

	# Delete temporary folder
	shutil.rmtree('temp')

	end = time.time()
	log("Accessibility calculation time: "+str(end-start))

###############
# Utils
###############

def prepare_data(G, df_osm_built, df_osm_pois, df_indices, num_processes, kw_arguments):
	""" 
	Pickles data to a temporary folder in order to achieve parallel accessibiltiy calculation
	A new subprocess will be created in order to minimize memory requirements

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	df_osm_built : geopandas.GeoDataFrame
		buildings data
	df_osm_pois : geopandas.GeoDataFrame
		buildings data
	df_indices : geopandas.GeoDataFrame
		data frame where indices will be calculated
	num_processes : int
		number of data chunks to create
	kw_arguments : pandas.Series
		additional keyword arguments

	Returns
	----------

	"""
	# Divide long edges
	divide_long_edges_graph(G, kw_arguments.max_edge_length )
	log("Graph long edges shortened")

	# Get activities
	df_built_activ = df_osm_built[ df_osm_built.classification.isin(["activity","mixed"]) ]
	df_pois_activ = df_osm_pois[ df_osm_pois.classification.isin(["activity","mixed"]) ]
	
	# Associate them to its closest node in the graph
	associate_activities_closest_node(G, df_built_activ, df_pois_activ )
	log("Activities associated to graph nodes")

	# Nodes dict
	for n, data in G.nodes.data(data=True):
		# Remove useless keys
		keys_ = list( data.keys() )
		[ data.pop(k) for k in keys_ if k not in ["x","y","num_activities"] ]
	
	# Edges dict
	for u, v, data in G.edges.data(data=True, keys=False):
		# Remove useless keys
		keys_ = list( data.keys() )
		[ data.pop(k) for k in keys_ if k not in ["length","key"] ]

	try:
		G.graph.pop("streets_per_node")
	except:
		pass
	# Pickle graph
	nx.write_gpickle(G, "temp/graph.gpickle")

	# Prepare input indices points
	data_split = np.array_split(df_indices, num_processes)
	for i in range(num_processes):
		data_split[i].to_pickle("temp/points_"+str(i)+".pkl")
	# Pickle arguments
	kw_arguments.to_pickle("temp/arguments.pkl")

def associate_activities_closest_node(G, df_activities_built, df_activities_pois ):
	""" 
	Associates the number of existing activities to their closest nodes in the graph

	Parameters
	----------
	G : networkx multidigraph
		input graph to calculate accessibility
	df_activities_built : pandas.DataFrame
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

	def associate_to_node(tree, point, G):
		distance, idx_node = tree.query( (point.x,point.y) )
		G.node[ df_nodes.loc[ idx_node, "node"] ]["num_activities"] += 1
	
	# Associate each activity to its closest node
	df_activities_built.apply(lambda x: associate_to_node(tree, x.geometry.centroid, G) , axis=1)
	df_activities_pois.apply(lambda x: associate_to_node(tree, x.geometry.centroid, G) , axis=1)
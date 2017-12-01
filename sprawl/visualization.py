###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
import collections
import numpy as np
import pandas as pd
import osmnx as ox
import os
import math
import logging as lg
from .utils import log
from .parameters import images_folder
from .storage import get_osm_data, get_stored_keys, load_data, get_kde_activities_residential, KEY_MESH_GRID, KEY_DISPERSION, KEY_ACCESSIBILITY, KEY_LANDUSEMIX, KEY_KDES, KEY_KDES_ACTIVITIES

##############################################################
### Parameters
##############################################################

# Number of colors in contourf call
num_zlevels = 15

# Save to files
SAVE_IMAGE = True
FILE_EXTENSION = 'png'
# When saving to file, DPI of image
DPI = 100

# Remove X and Y tick values during plot
REMOVE_XY_TICK_LABELS = True
# Add title to image
SAVE_WITH_TITLE = True

### Default parameters

# Plots
DEF_plot_landuses = True
DEF_plot_buildings = True
DEF_plot_histogram = True
DEF_plot_dispersion = True
DEF_plot_accessibility = True
DEF_plot_landusemix = True
DEF_plot_kdes = True
DEF_plot_kde_activity_types = True
DEF_plot_graph = True

# Divide long edges for graph
DEF_divide_long_edges = True

# Maximum indices values
DEF_accessibility_vminmax = None
DEF_dispersion_vminmax = None

# Maximum histogram value to analyse in bins
DEF_max_histogram_x = 100

# Bubble plot
DEF_desired_num_bubbles=5000
DEF_Bubble_Size_scaling=150

# Graph node size
DEF_Node_size = 4

# Figure size
DEF_figsize = (12,8)

# Make maximum values indices plot relative to cities maximum values?
PLOT_RELATIVE_MAX_VALUES = False
##############################################################

##############################################################
### Utils
##############################################################

def _generate_file_path(city_ref):
	"""
	Generates a file path for local storage data of `city_ref`

	Parameters
	----------
	city_ref : string
		name of the city
	
	Returns
	----------
	string
		path to storage file
	"""
	return os.path.join(images_folder, city_ref + '.' + FILE_EXTENSION)

### Creation of a nonlinear color mapping
class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""

    name = 'nlcmap'

    def __init__(self, cmap, levels):
        self.cmap = cmap
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels/ self.levels.max()
        self.levmax = self.levels.max()
        self.levmin = self.levels.min()
        self._y = np.linspace(self.levmin, self.levmax, len(self.levels))

    def __call__(self, xi, alpha=1.0, **kw):
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi/self.levmax, alpha)


def get_step_bubbles(rows, cols, DESIRED_NUM_BUBBLES):
	# Initialize step
	step = 1
	number_bubbles = math.ceil( float(rows) / float(step) ) * math.ceil( float(cols) / float(step) )
	DIFF_REAL_BUBBLES = abs(DESIRED_NUM_BUBBLES - number_bubbles)

	while (True):
		# Test number of bubbles with following step
		step_ = step + 1
		number_bubbles_ = math.ceil( float(rows) / float(step_) ) * math.ceil( float(cols) / float(step_) )
		DIFF_REAL_BUBBLES_ = abs(DESIRED_NUM_BUBBLES - number_bubbles_)

		# Are we closer to the desired number of bubbles?
		if (DIFF_REAL_BUBBLES_ > DIFF_REAL_BUBBLES): break
		# Update for following iteration
		step, number_bubbles, DIFF_REAL_BUBBLES = step_, number_bubbles_, DIFF_REAL_BUBBLES_
	
	return step

def bubble_group(matrix, DESIRED_NUM_BUBBLES):
	rows, cols = matrix.shape
	step = get_step_bubbles(rows, cols, DESIRED_NUM_BUBBLES)
	
	# Get rows and columns for grouped bubble matrix
	rows_ = int( math.ceil( float(rows) / float(step) ) )
	cols_ = int( math.ceil( float(cols) / float(step) ) )
	bubble_grouped = np.zeros( (rows_,cols_) , dtype=float)
	
	for i in range(0,rows,step):
		for j in range(0,cols,step):
			bubble_grouped[ int(i/step), int(j/step)] = matrix[i:i+step, j:j+step].mean()
	return bubble_grouped

def rescale_grid(size_grid):
	size_grid = size_grid / size_grid.max()

##############################################################
### Generic contour grid plot
##############################################################

def plot_contour_grid(grid_xy, grid_result, figsize, title = None, Z_levels=None, tick_label_last_element_plus=False):
	""" Contour plot for grids
	http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.contourf
	"""	
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	# Get array for X and Y values (xx and yy correspond to meshgrid)
	X,Y = grid_xy["x",0,:] , grid_xy["y",:,0]
	fig, ax = plt.subplots(figsize=figsize)

	# Ax
	if (REMOVE_XY_TICK_LABELS):
		ax.set_xticklabels([])
		ax.set_yticklabels([])
	
	# Color maps
	cmap_lin = cm.jet
	cmap_nonlin = nlcmap(cmap_lin, Z_levels)

	# Contour f
	cp = ax.contourf( X, Y, grid_result, Z_levels, cmap=cmap_nonlin, antialiased=True)
	
	# Color bar
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	cbar = fig.colorbar( cp, cax=cax )

	# Add a '+' for last element in colorbar
	if (tick_label_last_element_plus):
		last_element = cbar.ax.get_yticklabels()[-1].get_text()
		yticklabels = cbar.ax.get_yticklabels()
		yticklabels[-1].set_text( str(last_element) + "+" )
		cbar.ax.set_yticklabels(yticklabels)

	# Additional
	if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title)
	if (SAVE_IMAGE): plt.savefig( _generate_file_path(title) , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
	plt.show()



##############################################################
### Bubble plot
##############################################################

def plot_bubble(XY, color_grid, size_grid, figsize, DESIRED_NUM_BUBBLES=3500, Bubble_Size_scaling=100, title=None):
	""" Bubble plot for grid with values within [0,1]
	"""
	from mpl_toolkits.axes_grid1 import make_axes_locatable

	X, Y = bubble_group( np.matrix(XY.x), DESIRED_NUM_BUBBLES ), bubble_group( np.matrix(XY.y), DESIRED_NUM_BUBBLES )
	COLOR_GRID = bubble_group( np.matrix(color_grid), DESIRED_NUM_BUBBLES)
	SIZE_GRID = bubble_group( np.matrix(size_grid), DESIRED_NUM_BUBBLES)
	# Re scale between 0 and 1
	rescale_grid(SIZE_GRID)
	###

	fig, ax = plt.subplots(figsize=figsize)
	# Axes
	ax.set_ylim(Y.min(), Y.max())
	ax.set_xlim(X.min(), X.max())
	if (REMOVE_XY_TICK_LABELS):
		ax.set_xticklabels([])
		ax.set_yticklabels([])
	# Scatter plot
	cmap = plt.get_cmap('jet')
	sc = plt.scatter(X,Y, c=COLOR_GRID, s=SIZE_GRID*Bubble_Size_scaling, cmap=cmap , antialiased=True, edgecolors='black', linewidths=0.5)

	# Color bar
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	plt.colorbar(sc, cax=cax)

	# Additional
	if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title)
	if (SAVE_IMAGE): plt.savefig( _generate_file_path(title) , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
	plt.show()

##############################################################
### Histogram analysis
##############################################################

def plot_distance_histogram(df, figsize, title = None, max_x = None):
	""" Plot the histogram of the computed closest distances
	"""
	import math
	df_poly = df.iloc[ df.closest_d.dropna().index ]
	if (max_x is None): max_x = int( df_poly.closest_d.max() / 10 )
	bins_ = 20
	histo, bins = np.histogram(df_poly.closest_d,bins=bins_, range=(0,max_x))
	
	# Log histogram and bins
	log( "Histogram; Size: " + str(histo) + "\t" + str( len(histo) ) )
	log( "Bins; Size:" + str( bins ) + "\t" + str( len(bins) ) )

	fig, ax = plt.subplots(figsize=figsize)
	if (REMOVE_XY_TICK_LABELS):
		ax.set_xticklabels([])
	#ax.set_yticklabels([])
	ax.set_xlabel('Distance (meters) to closest neighbor')
	ax.set_ylabel('Frequency')

	ax.bar( bins[:-1], histo, tick_label=bins[:-1], width=4 )

	if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title+": Histogram of nearest neighboring distance")
	if (SAVE_IMAGE): plt.savefig( _generate_file_path(title) , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
	plt.show()	

##############################################################
### Plot activity and residential geometry uses 
##############################################################

def plot_landuse_scatter(df, figsize, title=None, alphas=[0.5,0.1],little_margin=0.01):

	df_activities = df.loc[ df.classification.isin(["activity","mixed"]) ]
	df_residential = df.loc[ df.classification.isin(["residential","mixed"]) ]

	x_act, y_act = list( zip(*[ geom.centroid.coords[0] for geom in df_activities.geometry ]) )
	x_res, y_res = list( zip(*[ geom.centroid.coords[0] for geom in df_residential.geometry ]) )

	colors=['b','r']

	miny,maxy = min( min(y_act), min(y_res) ), max( max(y_act), max(y_res) )
	minx,maxx = min( min(x_act), min(x_res) ), max( max(x_act), max(x_res) )

	log( "Activity uses: "+str( len(df_activities) ) )
	log( "Residential uses: "+str( len(df_residential) ) )
	alphas = [ min( 1., 1000./len(df_activities) ), min( 1., 1000./len(df_residential) ) ]

	for x,y,color,alpha,category in zip([x_act,x_res],[y_act,y_res],colors,alphas,["Activity points","Residential points"]):
		fig, ax = plt.subplots(figsize=figsize)

		if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title+": "+category)
		if (REMOVE_XY_TICK_LABELS):
			ax.set_xticklabels([])
			ax.set_yticklabels([])
		ax.set_ylim( miny - little_margin, maxy + little_margin )
		ax.set_xlim( minx - little_margin, maxx + little_margin )
		#ax.set_aspect('equal')
		if ((alpha < 0) or (alpha > 1)): alpha = 0.1
		plt.scatter(x, y, alpha=alpha, color=color)
		if (SAVE_IMAGE): plt.savefig( _generate_file_path(title+": "+category) , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
		plt.show()	

##############################################################
### Related to polygon buildings
##############################################################

def plot_polygons_scatter(df, figsize, title=None):
	""" Scatter plot using polygon's centroid
	Bubble size proportional to nearest neighbor distance
	"""
	# Bubble size multiplication factor (Scatter plot context)
	BUBBLE_MULTIPLICATION_FACTOR = 1
	# Bubble size offset (sum; Scatter plot context)
	BUBBLE_OFFSET = 1

	# Get only geometries with closest distance computed
	df_poly = df.iloc[ df.closest_d.dropna().index ]
	# Get separate list for coordinates
	x,y = list( zip(*[ point.centroid.coords[0] for point in df_poly.geometry ] ) )

	fig, ax = plt.subplots(figsize=figsize)
	# Scatter plot: Bubble size according to distance? 
	ax.scatter( x , y, s= df_poly.closest_d * BUBBLE_MULTIPLICATION_FACTOR + BUBBLE_OFFSET, alpha=0.55 )
	ax.set_xlim( float(min(x) ), float(max(x) ) )
	ax.set_ylim( float(min(y) ), float(max(y) ) )
	if (REMOVE_XY_TICK_LABELS):
		ax.set_xticklabels([])
		ax.set_yticklabels([])

	if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title+": Bubble buildings nearest neighboring distance")
	if (SAVE_IMAGE): plt.savefig( _generate_file_path(title+": Dispersion scatter plot") , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
	plt.show()	
	
def plot_polygons(df, figsize, title=None):
	""" Plot input polygons
	"""
	from matplotlib.collections import PatchCollection
	from descartes import PolygonPatch
	from shapely.geometry import Polygon, MultiPolygon, shape
	# Polygons
	df_polys = df[[x.type != "Point" for x in df['geometry'] ]]
	mp = MultiPolygon(  list( df_polys.geometry.ravel() )  )
	num_colours = 50
	cmap = plt.cm.get_cmap("hsv", num_colours+1)
	# Bounds
	minx, miny, maxx, maxy = mp.bounds
	fig, ax = plt.subplots(figsize=figsize)
	ax.set_xlim(minx , maxx )
	ax.set_ylim(miny , maxy )
	# Patches
	patches = []
	for idx, p in enumerate(mp):		
		patches.append(PolygonPatch(p, fc=cmap(idx%num_colours), ec='#555555', alpha=1., zorder=1))
	ax.add_collection(PatchCollection(patches, match_original=True))
	if (REMOVE_XY_TICK_LABELS):
		ax.set_xticklabels([])
		ax.set_yticklabels([])
	if ((not title is None) and (SAVE_WITH_TITLE) ): ax.set_title(title+": Buildings")
	if (SAVE_IMAGE): plt.savefig( _generate_file_path(title+": Buildings") , format=FILE_EXTENSION, bbox_inches='tight', dpi=DPI)
	plt.show()

##############################################################
### Graph
##############################################################

def plot_graph_length_color(G, fig_size, title=None):
	# show the simplified network with edges colored by edge length
	if (False): # Edge color by length
		ec = ox.get_edge_colors_by_attr(G, attr='length')
	else: # Edge color yellow
		ec = 'y'
	fig_width, fig_height = fig_size
	fig, ax = ox.plot_graph(G, node_color='b', node_edgecolor='k', node_size=DEF_Node_size, node_zorder=3, edge_color=ec, 
		edge_linewidth=3, fig_height=fig_height,fig_width=fig_width, dpi=DPI,save=SAVE_IMAGE,filename=title)
	plt.show()

def plot_graph_city( city_ref, fig_size, divide_long_edges=True ):
	from sprawl.storage import get_route_graph
	from .graph_utils import divide_long_edges_graph
	from .accessibility import MAX_EDGE_LENGTH
	# Get graph
	G = get_route_graph(city_ref)
	
	if (divide_long_edges): # Divide long edges?
		# Divide long edges
		divide_long_edges_graph(G,MAX_EDGE_LENGTH)

	# Plot G
	plot_graph_length_color(G, fig_size, title=city_ref)

##############################################################
### Generic call to plot stored key grids
##############################################################

def plot_grids( city_ref, KEYS_, fig_size, vminmax=None ):
	# For each calculated grid (different step)
	keys_grid = [ key for key in get_stored_keys(city_ref) if KEY_MESH_GRID in key ]
	for key_grid in keys_grid:
		# Load grid
		grid_xy = load_data(city_ref, key_grid)
		# Get step number used
		step = key_grid.split("_")[-1]
		
		desired_keys = []
		# Can be one or more keys, depending on type
		if (type(KEYS_) != list): KEYS_ = [ KEYS_ ]
		
		for key_ in KEYS_: # For each key we can find within grid step
			desired_keys = desired_keys + [ key for key in get_stored_keys(city_ref) if (key_+step in key) ]

		# Empty
		if (len(desired_keys) == 0): continue

		for desired_key in desired_keys:
			# Default: False
			tick_label_last_element_plus = False

			# Load grid accessibility
			grid = load_data(city_ref, desired_key)

			
			if (vminmax): # Set max value if minmax values not set
				Z_levels = np.linspace(start= vminmax[0], stop= vminmax[1], num= num_zlevels).tolist()
			else:
				if ("accessibility" in desired_key): # Pre define Z values
					ACCESSIBILITY_REGULAR_COLOR_MAX_VAL = 250

					# Z levels
					Z_levels = [0,5,10,15,20,25,30,40,50, 75, 100, 125, 150, 200, ACCESSIBILITY_REGULAR_COLOR_MAX_VAL]
					# Set last element '+'
					tick_label_last_element_plus = True
				elif ("dispersion" in desired_key): # Pre define Z values
					DISPERSION_REGULAR_COLOR_MAX_VAL = 35
					def Upper_limit(x,MAX_D):
						if ( np.isnan(x) ): return x
						return min(x,MAX_D)
					# Define upper limit values for grid
					grid = grid.applymap( lambda x: Upper_limit(x,DISPERSION_REGULAR_COLOR_MAX_VAL ) )
					
					# Z levels
					#Z_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, DISPERSION_REGULAR_COLOR_MAX_VAL]
					#Z_levels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, DISPERSION_REGULAR_COLOR_MAX_VAL]
					Z_levels = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, DISPERSION_REGULAR_COLOR_MAX_VAL]
					# Set last element '+'
					tick_label_last_element_plus = True
				else: # Default: [0,1] values
					Z_levels = np.linspace(start= 0, stop=np.matrix(grid).max(), num= num_zlevels).tolist()


			# Title: Cit_ref: KEY_OF_ANALYSIS
			text_association = {"accessibility":"Accessibility","landusemix":"Land use mix","dispersion":"Dispersion","activity":"Activity densities","residential":"Residential densities"}
			
			# Special case for activity types classification
			def activity_classification_tag(key):
				if ('shop' in key):
					return 'Activities (Shop)'
				if ('leisure' in key):
					return 'Activities (Leisure and amenities)'
				if ('commercial' in key):
					return 'Activities (Commercial and industrial)'

			title = city_ref+": "+ text_association.get( desired_key.split('_')[1], activity_classification_tag(desired_key) )

			plot_contour_grid(grid_xy, grid, fig_size, title=title, Z_levels=Z_levels, tick_label_last_element_plus=tick_label_last_element_plus)

			if ("landusemix" in desired_key):
				# grid: land use mix
				# size_grid: Obtain sum of KDE's land usages
				kde_act, kde_res = get_kde_activities_residential(city_ref, step)
				size_grid = kde_act + kde_res
				plot_bubble(grid_xy, grid, size_grid, fig_size, DESIRED_NUM_BUBBLES=DEF_desired_num_bubbles, Bubble_Size_scaling=DEF_Bubble_Size_scaling, title=city_ref+": Bubble land use mix", )

##############################################################
### Maximum indices values relative to ...
##############################################################

def find_max_values(cities, KEYS_):
	MAX_VALUE = 0

	for city_ref in cities:
		# For each calculated grid (different step)
		keys_grid = [ key for key in get_stored_keys(city_ref) if KEY_MESH_GRID in key ]
		for key_grid in keys_grid:
			# Load grid
			grid_xy = load_data(city_ref, key_grid)
			# Get step number used
			step = key_grid.split("_")[-1]
			
			desired_keys = []
			# Can be one or more keys, depending on type
			if (type(KEYS_) != list): KEYS_ = [ KEYS_ ]
			
			for key_ in KEYS_: # For each key we can find within grid step
				desired_keys = desired_keys + [ key for key in get_stored_keys(city_ref) if (key_+step in key) ]

			# Empty
			if (len(desired_keys) == 0): continue

			for desired_key in desired_keys:
				# Load grid
				grid = load_data(city_ref, desired_key)
				# Get maximum value
				MAX_VALUE = max( MAX_VALUE , np.nanmax( grid.as_matrix() ) )
	return MAX_VALUE

##############################################################
### City visualization
##############################################################

def visualize_city(city_ref, kwargs = { 'figsize':(12,8), 'plot_buildings':True, 'plot_histogram':True, 'plot_kdes':True, 'plot_kde_activity_types':True, 'plot_graph':True, 'plot_accessibility':True, 
		'plot_landusemix':True, 'plot_dispersion':True, 'plot_landuses':True, 'divide_long_edges':True, 'dispersion_vminmax':None, 'accessibility_vminmax':None, 'landusemix_vminmax':(0,1),
		'max_histogram_x':100
		} ):
	""" Complete visualization for input data city
	"""	
	log( "Visualization: "+city_ref )

	if (SAVE_IMAGE): # Save images?
		if (not(os.path.exists(images_folder))):
			os.makedirs(images_folder)

	if (PLOT_RELATIVE_MAX_VALUES):
		if (kwargs.get("dispersion_vminmax",DEF_dispersion_vminmax) is None):
			kwargs["dispersion_vminmax"] = ( 0, find_max_values([city_ref], KEY_DISPERSION) )
		if (kwargs.get("accessibility_vminmax",DEF_accessibility_vminmax) is None):
			kwargs["accessibility_vminmax"] = ( 0 , find_max_values([city_ref], KEY_ACCESSIBILITY) )
		if (kwargs.get("landusemix_vminmax",None) ):
			kwargs["landusemix_vminmax"] = (0,1)

	df = get_osm_data(city_ref)	
	
	if (kwargs.get("plot_landuses",DEF_plot_landuses)):
		plot_landuse_scatter(df, kwargs.get("figsize",DEF_figsize), title=city_ref)
	
	### Grid land use mix
	if (kwargs.get("plot_kdes",DEF_plot_kdes)):
		plot_grids( city_ref, KEY_KDES, kwargs.get("figsize",DEF_figsize) )

	if (kwargs.get("plot_kde_activity_types",DEF_plot_kde_activity_types)):
		plot_grids( city_ref, KEY_KDES_ACTIVITIES, kwargs.get("figsize",DEF_figsize) )

	if (kwargs.get("plot_landusemix",DEF_plot_landusemix)):
		plot_grids( city_ref, KEY_LANDUSEMIX, kwargs.get("figsize",DEF_figsize) )

	### Buildings
	if (kwargs.get("plot_buildings",DEF_plot_buildings)):
		plot_polygons(df, kwargs.get("figsize",DEF_figsize), title=city_ref)

	### Grid dispersion
	if (kwargs.get("plot_dispersion",DEF_plot_dispersion)):
		plot_grids( city_ref, KEY_DISPERSION, kwargs.get("figsize",DEF_figsize), vminmax=kwargs.get("dispersion_vminmax",DEF_dispersion_vminmax) )
		if ('closest_d' in df.columns):
			log( "Mean nearest neighbor: " + str( df.closest_d.mean() ) )
			log( "Median nearest neighbor: " + str( df.closest_d.median() ) )
			plot_polygons_scatter(df, kwargs.get("figsize",DEF_figsize), title=city_ref)

	### Closest distances histogram plot
	if (kwargs.get("plot_histogram",DEF_plot_histogram)):
		if ('closest_d' in df.columns):
			plot_distance_histogram(df, kwargs.get("figsize",DEF_figsize), max_x = kwargs.get("max_histogram_x",DEF_max_histogram_x), title=city_ref)

	### Grid accessibility
	if (kwargs.get("plot_accessibility",DEF_plot_accessibility)):
		plot_grids( city_ref, KEY_ACCESSIBILITY, kwargs.get("figsize",DEF_figsize), vminmax=kwargs.get("accessibility_vminmax",DEF_accessibility_vminmax) )
	
	if (kwargs.get("plot_graph",DEF_plot_graph)):
		plot_graph_city( city_ref, kwargs.get("figsize",DEF_figsize), divide_long_edges= kwargs.get("divide_long_edges",DEF_divide_long_edges) )


	


###############

def visualize_cities(cities_ref, kwargs = { 'figsize':(12,8), 'plot_buildings':True, 'plot_histogram':True, 'plot_kdes':True, 'plot_kde_activity_types':True, 'plot_graph':True, 'plot_accessibility':True, 
		'plot_landusemix':True, 'plot_dispersion':True, 'plot_landuses':True, 'divide_long_edges':True, 'dispersion_vminmax':None, 'accessibility_vminmax':None, 'landusemix_vminmax':(0,1),
		'max_histogram_x':100
		} ):
	""" Complete visualization for input cities name
	Assumption: Results are already stored
	"""

	if (PLOT_RELATIVE_MAX_VALUES):
		if (kwargs.get("dispersion_vminmax",DEF_dispersion_vminmax) is None):
			kwargs["dispersion_vminmax"] = ( 0, find_max_values([city_ref], KEY_DISPERSION) )
		if (kwargs.get("accessibility_vminmax",DEF_accessibility_vminmax) is None):
			kwargs["accessibility_vminmax"] = ( 0 , find_max_values([city_ref], KEY_ACCESSIBILITY) )
		if (kwargs.get("landusemix_vminmax",None) ):
			kwargs["landusemix_vminmax"] = (0,1)

	# For each city: Complete visualization
	for city_ref in cities_ref:
		visualize_city(city_ref, kwargs )
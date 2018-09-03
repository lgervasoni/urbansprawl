###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import geopandas as gpd


def proportional_population_downscaling(df_osm_built, df_insee):
	"""
	Performs a proportional population downscaling considerig the surface dedicated to residential land use
	Associates the estimated population to each building in column 'population'

	Parameters
	----------
	df_osm_built : geopandas.GeoDataFrame
		input buildings with computed residential surface
	df_insee : geopandas.GeoDataFrame
		INSEE population data

	Returns
	----------

	"""
	df_osm_built['geom'] = df_osm_built.geometry
	df_osm_built_residential = df_osm_built[ df_osm_built.apply(lambda x: x.landuses_m2['residential'] > 0, axis = 1) ]
	
	# Same projection ?
	assert( df_insee.crs.get('datum') == df_osm_built_residential.crs.get('datum') )
	assert( df_insee.crs.get('proj') == df_osm_built_residential.crs.get('proj') )
	assert( df_insee.crs.get('zone') == df_osm_built_residential.crs.get('zone') )
	assert( df_insee.crs.get('units') == df_osm_built_residential.crs.get('units') )
	# Loading/saving using geopandas loses the 'ellps' key
	df_insee.crs = df_osm_built_residential.crs

	# Intersecting gridded population - buildings
	sjoin = gpd.sjoin( df_insee, df_osm_built_residential, op='intersects')
	# Calculate area within square (percentage of building with the square)
	sjoin['residential_m2_within'] = sjoin.apply(lambda x: x.landuses_m2['residential'] * (x.geom.intersection(x.geometry).area / x.geom.area), axis=1 )
	# Initialize
	df_insee['residential_m2_within'] = 0
	# Sum residential area within square
	sum_m2_per_square = sjoin.groupby(sjoin.index)['residential_m2_within'].sum()
	# Assign total residential area within each square
	df_insee.loc[ sum_m2_per_square.index, "residential_m2_within" ] = sum_m2_per_square.values
	# Get number of M^2 / person
	df_insee[ "m2_per_person" ] = df_insee.apply(lambda x: x.residential_m2_within / x.pop_count, axis=1)

	def population_building(x, df_insee):
		# Sum of: For each square: M2 of building within square / M2 per person
		return ( x.get('m2',[]) / df_insee.loc[ x.get('idx',[]) ].m2_per_person ).sum()
	# Index: Buildings , Values: idx:Indices of gridded square population, m2: M2 within that square
	buildings_square_m2_association = sjoin.groupby('index_right').apply(lambda x: {'idx':list(x.index), 'm2':list(x.residential_m2_within)} )
	# Associate
	df_osm_built.loc[ buildings_square_m2_association.index, "population" ] = buildings_square_m2_association.apply(lambda x: population_building(x,df_insee) )
	# Drop unnecessary column
	df_osm_built.drop('geom', axis=1, inplace=True)


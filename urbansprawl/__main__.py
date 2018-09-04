###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

import sys

# Usage: "python urbansprawl [method] [arguments]"

def main():
	print(sys.argv)

	if (sys.argv[1] == "osm_data"):
		print("Retrieve OSM data")
		# TODO CALL

	if (sys.argv[1] == "indices"):
		print("Retrieve sprawl indices")
		# TODO CALL

	if (sys.argv[1] == "population_downscaling"):
		print("Perform population downscaling")
		# TODO CALL

if __name__== "__main__":
	main()
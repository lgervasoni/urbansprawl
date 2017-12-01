###################################################################################################
# Repository: https://gitlab.inria.fr/gervason/urbansprawl
###################################################################################################

##############################################################
# This file contain parameters associated to the framework
##############################################################

import logging as lg

##############################################################
### Utils
##############################################################
# Logger parameters
logs_folder = 'data'
log_console = False
log_file = False
log_level = lg.INFO
log_name = 'urban'
log_filename = 'urban'

storage_folder = 'data'
images_folder = 'images'
shapefiles_folder = 'input_data'

##############################################################
### Spatial indices
##############################################################

############
### Default step to calculate the different spatial indices
DEF_dispersion_step = 100
DEF_landusemix_step = 100
DEF_accessibility_step = 200
############
###################################################################################################
# Repository: https://github.com/lgervasoni/urbansprawl
# MIT License
###################################################################################################

from .data_extract import get_extract_population_data
from .downscaling import proportional_population_downscaling
from .urban_features import compute_full_urban_features, get_training_testing_data, get_Y_X_features_population_data
from .utils import get_aggregated_squares, population_downscaling_validation
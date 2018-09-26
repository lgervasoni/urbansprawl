"""urbansprawl package
"""

# OpenStreetMap data
from .osm.core import get_route_graph, get_processed_osm_data

# Spatial urban sprawl indices
from .sprawl.core import compute_grid_landusemix, compute_grid_accessibility, compute_grid_dispersion
from .sprawl.core import get_indices_grid, process_spatial_indices

# Disaggrated population estimates
from .population.core import get_extract_population_data, compute_full_urban_features, get_training_testing_data, get_Y_X_features_population_data
from .population.core import get_aggregated_squares, proportional_population_downscaling, population_downscaling_validation


__version__ = '1.1'

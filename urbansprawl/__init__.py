"""urbansprawl package
"""

# OpenStreetMap data
from .osm.core import get_route_graph, get_processed_osm_data

# Spatial urban sprawl indices
from .sprawl.core import compute_grid_landusemix, compute_grid_accessibility, compute_grid_dispersion
from .sprawl.core import get_indices_grid, process_spatial_indices

__version__ = '1.1'

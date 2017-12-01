# Urbansprawl

The urbansprawl project provides an open framework to assess the urban sprawl phenomenon.
It uses OpenStreetMap data to calculate its sprawling indices, divided in Accessibility, Land use mix, and Dispersion.


**For more details, refer to:**

* Gervasoni Luciano, Bosch Martí, Fenet Serge, and Sturm Peter. 2016. "[A framework for evaluating urban land use mix from crowd-sourcing data](https://hal.inria.fr/hal-01396792)." 2nd International Workshop on Big Data for Sustainable Development.

* Gervasoni Luciano, Bosch Martí, Fenet Serge, and Sturm Peter. 2017. "[LUM_OSM: une plateforme pour l'évaluation de la mixité urbaine à partir de données participatives](https://hal.inria.fr/hal-01548341)." GAST Workshop, Conférence Extraction et Gestion de Connaissances (EGC 2017).

* Gervasoni Luciano, Bosch Martí, Fenet Serge, and Sturm Peter. 2017. "[Calculating spatial urban sprawl indices using open data](https://hal.inria.fr/hal-01535469)." 15th International Conference on Computers in Urban Planning and Urban Management.

NOTE: The implementation for the previous publications can be found at version 1.0

## Dependencies

urbansprawl works with Python 2+3.

- Python dependencies:
```sh
osmnx scikit-learn
```

* Using anaconda:
```sh
conda update -c conda-forge --all
conda install -c conda-forge osmnx scikit-learn
```

## Example

OpenStreetMap data is retrieved using the Overpass API.

Results are depicted for the city of **Lyon, France**:

- Locations of residential and activity land uses are retrieved

* Activity uses:

<img src="examples/images/Lyon_activities.png" width="550" height="500">

* Residential uses:

<img src="examples/images/Lyon_residential.png" width="550" height="500">

- Densities for each land use are estimated:

<img src="examples/images/Lyon_densities.png" width="1200" height="500">

* Activity uses can be further classified:

<img src="examples/images/Lyon_activities_densities.png" width="1200" height="500">

- Street network:

<img src="examples/images/Lyon_graph.png" width="550" height="500">

**Sprawling indices:**

- Land use mix indices:

<img src="examples/images/Lyon_Landusemix.png" width="550" height="500">

- Accessibility indices:

<img src="examples/images/Lyon_Accessibility.png" width="550" height="500">

- Dispersion indices:

<img src="examples/images/Lyon_Dispersion.png" width="550" height="500">

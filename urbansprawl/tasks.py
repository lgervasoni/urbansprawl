"""This module aims at recovering OpenStreetMap data through Overpass API

To accomplish a task, the following command must be run on the terminal:

```
python -m luigi --local-scheduler --module urbansprawl.tasks <Task> <params>
```

with `Task` one of the class defined below, and `params` the corresponding
parameters.

This computation is done locally (because of the ̀--local-scheduler` option). It
can be done on a server, by first launching an instance of the luigi daemon :
```
luigid
̀``

and then by running the previous command without the `--local-scheduler`
option. The task dependency graph and some miscellaneous information about the
tasks are visible at `localhost:8082` URL address.

"""

from datetime import date
import geopandas as gpd
import luigi
import numpy as np
import os
import osmnx
import pandas as pd

from urbansprawl.osm.osm_overpass import (create_buildings_gdf,
                                          create_building_parts_gdf,
                                          create_pois_gdf,
                                          create_landuse_gdf,
                                          retrieve_route_graph)
from urbansprawl.osm.osm_utils import (sanity_check_height_tags,
                                       associate_structures)
from urbansprawl.osm.osm_data import (classify_tag,
                                      classify_activity_category,
                                      compute_landuse_inference)
from urbansprawl.osm.osm_surface import compute_landuses_m2

# Columns of interest corresponding to OSM keys
OSM_TAG_COLUMNS = [ "amenity", "landuse", "leisure", "shop", "man_made",
                    "building", "building:use", "building:part" ]
COLUMNS_OF_INTEREST = OSM_TAG_COLUMNS + ["osm_id", "geometry", "height_tags"]
COLUMNS_OF_INTEREST_POIS = OSM_TAG_COLUMNS + ["osm_id", "geometry"]
COLUMNS_OF_INTEREST_LANDUSES = ["osm_id", "geometry", "landuse"]
HEIGHT_TAGS = [ "min_height", "height", "min_level", "levels",
                "building:min_height", "building:height", "building:min_level",
                "building:levels", "building:levels:underground" ]
MINIMUM_M2_BUILDING_AREA = 9.0

def define_filename(description, city, date, datapath, geoformat):
    """Build a distinctive filename regarding a given `description`, `city`,
    `date` (ISO-formatted), ̀datapath` and a `geoformat` for the file extension

    Parameters
    ----------
    description : str
        Describe the file content in one word
    city : str
        City of interest, used for the queries to Overpass API
    date : str
        Date of the Overpass query, in ISO format
    datapath : str
        Path of the file on the file system
    geoformat : str
        File extension, *i.e.* GeoJSON

    Returns
    -------
    str
        Full path name on the file system
    """
    os.makedirs(datapath, exist_ok=True)
    filename = "{}-{}.{}".format(description, date, geoformat)
    return os.path.join(datapath, city, filename)

def set_list_as_str(l):
    """Small utility function to transform list in string

    Parameters
    ----------
    l : list
        Input list

    Returns
    -------
    str
        Stringified version of the input list, with items separated with a comma
    """
    if type(l) == list:
        return ','.join(str(e) for e in l)

def clean_list_in_geodataframe_column(gdf, column):
    """Stringify items of `column` within ̀gdf`, in order to allow its
    serialization

    Parameters
    ----------
    gdf : GeoDataFrame
        Input data structure
    column : str
        Column to modify

    Returns
    -------
    GeoDataFrame
        Modified input structure, with a fixed `column` (contains stringified items)
    """
    if column in gdf.columns:
        gdf[column] = gdf[column].apply(lambda x: set_list_as_str(x))
    return gdf


class GetBoundingBox(luigi.Task):
    """Extract the bounding box around a given `city`

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetBoundingBox
    --city valence-drome
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    (default: `./data`)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")

    def output(self):
        """Indicates the task result destination onto the file system
        """
        path = os.path.join(self.datapath, self.city)
        os.makedirs(path, exist_ok=True)
        return luigi.LocalTarget(os.path.join(path, "bounding_box.geojson"))

    def run(self):
        """Main operations of the Luigi task
        """
        city_gdf = osmnx.gdf_from_place(self.city, which_result=1)
        city_gdf.to_file(self.output().path, driver="GeoJSON")


class GetBuildings(luigi.Task):
    """Give a raw version of OpenStreetMap buildings through an Overpass API
    query

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetBuildings
    --city valence-drome --date-query 2017-01-01T1200
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())

    def requires(self):
        """Gives the task(s) that are needed to accomplish the current one. It
        refers implicitely to the project dependency graph.
        """
        return GetBoundingBox(self.city, self.datapath)

    def output(self):
        output_path = define_filename("buildings",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        city_gdf = gpd.read_file(self.input().path)
        north, south, east, west = city_gdf.loc[0, ["bbox_north", "bbox_south",
                                                    "bbox_east", "bbox_west"]]
        date = "[date:'" + str(self.date_query) + "']"
        buildings = create_buildings_gdf(date=date,
                                         north=north, south=south,
                                         east=east, west=west)
        buildings.drop(["nodes"], axis=1, inplace=True)
        buildings.to_file(self.output().path, driver="GeoJSON")

class GetBuildingParts(luigi.Task):
    """Give a raw version of OpenStreetMap building parts through an Overpass API
    query

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetBuildingParts
    --city valence-drome --date-query 2017-01-01T1200
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())

    def requires(self):
        return GetBoundingBox(self.city, self.datapath)

    def output(self):
        output_path = define_filename("building-part",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        city_gdf = gpd.read_file(self.input().path)
        north, south, east, west = city_gdf.loc[0, ["bbox_north", "bbox_south",
                                                    "bbox_east", "bbox_west"]]
        date = "[date:'" + str(self.date_query) + "']"
        building_parts = create_building_parts_gdf(date=date,
                                                   north=north, south=south,
                                                   east=east, west=west)
        columns_to_drop = [col for col in list(building_parts.columns)
                           if not col in COLUMNS_OF_INTEREST]
        building_parts.drop(["nodes"], axis=1, inplace=True)
        building_parts.to_file(self.output().path, driver="GeoJSON")


class GetPOIs(luigi.Task):
    """Give a raw version of OpenStreetMap Points of Interest (POIs) through an
    Overpass API query

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetPOIs
    --city valence-drome --date-query 2017-01-01T1200
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())

    def requires(self):
        return GetBoundingBox(self.city, self.datapath)

    def output(self):
        output_path = define_filename("pois",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        city_gdf = gpd.read_file(self.input().path)
        north, south, east, west = city_gdf.loc[0, ["bbox_north", "bbox_south",
                                                    "bbox_east", "bbox_west"]]
        date = "[date:'" + str(self.date_query) + "']"
        pois = create_pois_gdf(date=date,
                               north=north, south=south,
                               east=east, west=west)
        columns_to_drop = [col for col in list(pois.columns)
                           if not col in COLUMNS_OF_INTEREST_POIS]
        pois.drop(columns_to_drop, axis=1, inplace=True)
        pois["osm_id"] = pois.index
        pois.reset_index(drop=True, inplace=True)
        pois.to_file(self.output().path, driver="GeoJSON")


class CreateLandUse(luigi.Task):
    """Query the OpenStreetMap land uses with Overpass API

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks CreateLandUse
    --city valence-drome
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())

    def requires(self):
        return GetBoundingBox(self.city, self.datapath)

    def output(self):
        output_path = define_filename("land-use",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        city_gdf = gpd.read_file(self.input().path)
        north, south, east, west = city_gdf.loc[0, ["bbox_north", "bbox_south",
                                                    "bbox_east", "bbox_west"]]
        date = "[date:'" + str(self.date_query) + "']"
        landuses = create_landuse_gdf(date=date,
                                       north=north, south=south,
                                       east=east, west=west)
        landuses = landuses[["landuse", "geometry"]]
        landuses["osm_id"] = landuses.index
        columns_to_drop = [col for col in list(landuses.columns)
                           if not col in COLUMNS_OF_INTEREST_LANDUSES]
        landuses.drop(columns_to_drop, axis=1, inplace=True)
        landuses.reset_index(drop=True, inplace=True)
        landuses.to_file(self.output().path, driver="GeoJSON")

class SanityCheck(luigi.Task):
    """Check buildings and building parts GeoDataFrames, especially their
    height tags

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks SanityCheck
    --city valence-drome --table buildings
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    table : str
        Structure to check, either `buildings` or `building-parts`
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    table = luigi.Parameter(default="buildings")

    def requires(self):
        if self.table == "buildings":
            return GetBuildings(self.city, self.datapath,
                                self.geoformat, self.date_query)
        elif self.table == "building-parts":
            return GetBuildingParts(self.city, self.datapath,
                                    self.geoformat, self.date_query)
        else:
            raise ValueError(("Please provide a valid table name (either "
                              "'buildings' or 'building-parts')."))

    def output(self):
        output_path = define_filename("checked-" + self.table,
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        gdf = gpd.read_file(self.input().path)
        sanity_check_height_tags(gdf)
        def remove_nan_dict(x):
            """Remove entries with nan values
            """
            return {k:v for k, v in x.items() if pd.notnull(v)}
        gdf['height_tags'] = gdf[[c for c in HEIGHT_TAGS
                                  if c in gdf.columns]].apply(lambda x:
                                                              remove_nan_dict(x.to_dict() ), axis=1)
        columns_to_drop = [col for col in list(gdf.columns)
                           if not col in COLUMNS_OF_INTEREST]
        gdf.drop(columns_to_drop, axis=1, inplace=True)
        gdf["osm_id"] = gdf.index
        gdf.reset_index(drop=True, inplace=True)
        gdf.to_file(self.output().path, driver="GeoJSON")


class GetClassifiedInfo(luigi.Task):
    """Classify each building, building part or POI record as "residential",
    "activity" or "mixed" according to the associated tags

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetClassifiedInfo
    --city valence-drome --table buildings
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    table : str
        Structure of interest, either `buildings`, `building-parts` or `pois`
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    table = luigi.Parameter(default="buildings")

    def requires(self):
        if self.table == "buildings":
            return SanityCheck(self.city, self.datapath,
                               self.geoformat, self.date_query, self.table)
        elif self.table == "building-parts":
            return SanityCheck(self.city, self.datapath,
                               self.geoformat, self.date_query, self.table)
        elif self.table == "pois":
            return GetPOIs(self.city, self.datapath,
                           self.geoformat, self.date_query)
        else:
            raise ValueError(("Please provide a valid table name (either "
                              "'buildings', 'building-parts', 'pois')."))

    def output(self):
        output_path = define_filename("classified-" + self.table,
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        gdf = gpd.read_file(self.input().path)
        gdf['classification'], gdf['key_value'] = list( zip(*gdf.apply(classify_tag, axis=1)) )
        if self.table == "buildings":
            gdf.drop(gdf[gdf.classification.isnull()].index, inplace=True)
            gdf.reset_index(inplace=True, drop=True)
        elif self.table == "building-parts":
	    # Building parts will acquire its containing building land use
            # if it is not available
            gdf.loc[gdf.classification.isin(["infer", "other"]),
                    "classification"] = None
        elif self.table == "pois":
            gdf.drop(gdf[gdf.classification.isin(["infer", "other"]) | gdf.classification.isnull()].index, inplace=True)
            gdf.reset_index(inplace=True, drop=True)
        else:
            raise ValueError(("Please provide a valid table name (either "
                              "'buildings', 'building-parts', 'pois')."))
        # Drop tag-related columns
        columns_to_drop = [col for col in OSM_TAG_COLUMNS if col in gdf.columns]
        gdf.drop(columns_to_drop, axis=1, inplace=True)
        gdf.to_file(self.output().path, driver="GeoJSON")


class SetupProjection(luigi.Task):
    """Fix the GeoDataFrames projections, so as to ensure that every
    GeoDataFrames has the same projection

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks SetupProjection
    --city valence-drome --table buildings-parts
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    table : str
        Structure of interest, either `buildings`, `building-parts` or `pois`
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    table = luigi.Parameter(default="buildings")
    srid = luigi.Parameter(default=4326)

    def requires(self):
        if self.table in ["buildings", "building-parts", "pois"]:
            return GetClassifiedInfo(self.city, self.datapath,
                                     self.geoformat, self.date_query,
                                     self.table)
        elif self.table == "land-uses":
            return CreateLandUse(self.city, self.datapath,
                                 self.geoformat, self.date_query)
        else:
            raise ValueError(("Please provide a valid table name (either "
                              "'buildings', 'building-parts', "
                              "'pois' or 'land-uses')."))

    def output(self):
        output_path = define_filename("reprojected-" + self.table,
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        gdf = gpd.read_file(self.input().path)
	### Project to UTM coordinates within the same zone
        if self.table == "buildings":
            gdf = osmnx.project_gdf(gdf)
            gdf.drop(gdf[gdf.geometry.area < MINIMUM_M2_BUILDING_AREA].index,
                     inplace=True)
        else:
            gdf = osmnx.project_gdf(gdf,
                                    to_crs={'init': "epsg:{}".format(self.srid)})
        gdf.to_file(self.output().path, driver="GeoJSON")


class InferLandUse(luigi.Task):
    """Infer land use of each OpenStreetMap building thanks to land use
    information

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks InferLandUse
    --city valence-drome --date-query 2017-01-01T1200 --srid 4326
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    srid = luigi.Parameter(default=4326)

    def requires(self):
        return {"buildings": SetupProjection(self.city, self.datapath,
                                             self.geoformat, self.date_query,
                                             "buildings", self.srid),
                "land-uses": SetupProjection(self.city, self.datapath,
                                             self.geoformat, self.date_query,
                                             "land-uses", self.srid)}

    def output(self):
        output_path = define_filename("infered-buildings",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        buildings = gpd.read_file(self.input()["buildings"].path)
        land_uses = gpd.read_file(self.input()["land-uses"].path)
        compute_landuse_inference(buildings, land_uses)
        assert(len(buildings[buildings.key_value=={"inferred":"other"} ]) == 0)
        assert(len(buildings[buildings.classification.isnull()]) == 0)
        buildings.to_file(self.output().path, driver="GeoJSON")


class AssociateStructures(luigi.Task):
    """Associate OpenStreetMap buildings with respectively building parts and
    POIs : enrich the building GeoDataFrame with building parts and POIs data

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks
    AssociateStructures --city valence-drome
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    srid = luigi.Parameter(default=4326)

    def requires(self):
        return {"buildings": InferLandUse(self.city, self.datapath,
                                          self.geoformat, self.date_query,
                                          self.srid),
                "building-parts": SetupProjection(self.city, self.datapath,
                                                  self.geoformat,
                                                  self.date_query,
                                                  "building-parts", self.srid),
                "pois": SetupProjection(self.city, self.datapath,
                                        self.geoformat, self.date_query,
                                        "pois", self.srid)}

    def output(self):
        output_path = define_filename("associated-buildings",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        buildings = gpd.read_file(self.input()["buildings"].path)
        building_parts = gpd.read_file(self.input()["building-parts"].path)
        pois = gpd.read_file(self.input()["pois"].path)
        associate_structures(buildings, building_parts,
                             operation='contains', column='containing_parts')
        associate_structures(buildings, pois,
                             operation='intersects', column='containing_poi')
        buildings.to_file(self.output().path, driver="GeoJSON")


class ComputeLandUse(luigi.Task):
    """Compute land use per building type (residential, activity or mixed)

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks ComputeLandUse
    --city valence-drome --default-heights 6 --meters-per-level 2
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    default_height : int
        Default building height, in meters (default: 3 meters)
    meters_per_level : int
        Default height per level, in meter (default: 3 meters)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    srid = luigi.Parameter(default=4326)
    default_height = luigi.Parameter(3)
    meters_per_level = luigi.Parameter(3)

    def requires(self):
        return {"buildings": AssociateStructures(self.city, self.datapath,
                                                 self.geoformat,
                                                 self.date_query, self.srid),
                "building-parts": SetupProjection(self.city, self.datapath,
                                                  self.geoformat,
                                                  self.date_query,
                                                  "building-parts", self.srid),
                "pois": SetupProjection(self.city, self.datapath,
                                        self.geoformat, self.date_query,
                                        "pois", self.srid)}

    def output(self):
        output_path = define_filename("buildings-with-computed-land-use",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        buildings = gpd.read_file(self.input()["buildings"].path)
        building_parts = gpd.read_file(self.input()["building-parts"].path)
        pois = gpd.read_file(self.input()["pois"].path)
        buildings['activity_category'] = buildings.apply(lambda x: classify_activity_category(x.key_value), axis=1)
        building_parts['activity_category'] = building_parts.apply(lambda x: classify_activity_category(x.key_value), axis=1)
        pois['activity_category'] = pois.apply(lambda x: classify_activity_category(x.key_value), axis=1)
        compute_landuses_m2(buildings,
                            building_parts,
                            pois,
                            default_height=self.default_height,
                            meters_per_level=self.meters_per_level,
                            mixed_building_first_floor_activity=True)
        buildings.loc[buildings.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
        building_parts.loc[building_parts.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
        pois.loc[pois.activity_category.apply(lambda x: len(x)==0 ), "activity_category" ] = np.nan
        # Set the composed classification given, for each building, its containing Points of Interest and building parts classification
        buildings.loc[buildings.apply(lambda x: x.landuses_m2["activity"]>0 and x.landuses_m2["residential"]>0, axis=1 ), "classification" ] = "mixed"
        buildings.to_file(self.output().path, driver="GeoJSON")


class GetRouteGraph(luigi.Task):
    """Retrieve routing graph for the given city, through its encompassing
    bounding box

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks GetRouteGraph
    --city valence-drome --srid 4326
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    srid = luigi.Parameter(default=4326)

    def requires(self):
        return GetBoundingBox(self.city, self.datapath)

    def output(self):
        output_path = define_filename("route-graph",
                                      self.city,
                                      self.date_query.isoformat(),
                                      self.datapath,
                                      self.geoformat)
        return luigi.LocalTarget(output_path)

    def run(self):
        city_gdf = gpd.read_file(self.input().path)
        north, south, east, west = city_gdf.loc[0, ["bbox_north", "bbox_south",
                                                    "bbox_east", "bbox_west"]]
        date = "[date:'" + str(self.date_query) + "']"
        retrieve_route_graph(self.city, date=date,
                             north=north, south=south,
                             east=east, west=west,
                             force_crs={'init': "epsg:{}".format(self.srid)})
        # `retrieve_route_graph` does not return any result
        # in order to consider the task as finished,
        # we may save the bounding box as an output
        city_gdf.to_file(self.output().path, driver="GeoJSON")


class MasterTask(luigi.Task):
    """Generic task that launches every final task

    Example:
    ```
    python -m luigi --local-scheduler --module urbansprawl.tasks MasterTask
    --city valence-drome --date-query 2017-01-01T1200 --srid 4326
    --default-height 3 --meters-per-level 3
    ```

    Attributes
    ----------
    city : str
        City of interest
    datapath : str
        Indicates the folder where the task result has to be serialized
    geoformat : str
        Output file extension (by default: `GeoJSON`)
    date_query : str
        Date to which the OpenStreetMap data must be recovered (format:
    AAAA-MM-DDThhmm)
    srid : int
        Geographical projection (default 4326, *i.e.* WGS84)
    default_height : int
        Default building height, in meters (default: 3 meters)
    meters_per_level : int
        Default height per level, in meter (default: 3 meters)
    """
    city = luigi.Parameter()
    datapath = luigi.Parameter("./data")
    geoformat = luigi.Parameter("geojson")
    date_query = luigi.DateMinuteParameter(default=date.today())
    srid = luigi.Parameter(default=4326)
    default_height = luigi.Parameter(default=3)
    meters_per_level = luigi.Parameter(default=3)

    def requires(self):
        yield ComputeLandUse(self.city, self.datapath,
                             self.geoformat, self.date_query, self.srid,
                             self.default_height, self.meters_per_level)
        yield GetRouteGraph(self.city, self.datapath,
                            self.geoformat, self.date_query, self.srid)

    def complete(self):
        return False

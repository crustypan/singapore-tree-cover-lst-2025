"""Measure distances from 500 m cell centers to mapped Singapore water."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rasterio.warp import transform_geom
from shapely import distance, force_2d, make_valid, points, union_all
from shapely.geometry import MultiLineString, Polygon, shape
from shapely.strtree import STRtree

from screen_landsat import GRID_CRS


MAP = Path("data/snapshot/grass_tree_substitution/sla_national_map_polygon.geojson.gz")
MAP_SHA256 = "80791c1e271c378f26ab8c618eb019848c89d0c77306adbaf559f744c5c9048b"
PRIMARY_AREA_KM2 = 0.25
ALTERNATE_AREA_KM2 = 0.10


def polygons(geometry):
    if geometry.is_empty:
        return
    if geometry.geom_type == "Polygon":
        yield geometry
    elif geometry.geom_type in ("MultiPolygon", "GeometryCollection"):
        for part in geometry.geoms:
            yield from polygons(part)


def water_geometry(land_parts, hydro):
    """Assemble shoreline and inland water from projected polygon layers."""
    coast_polygons = list(polygons(union_all(land_parts)))
    shoreline = MultiLineString([p.exterior for p in coast_polygons])
    # Include interior holes in the comparison footprint. The hydrographic
    # layer also delineates named reservoirs that reach beyond the coast layer.
    footprint = union_all([Polygon(p.exterior) for p in coast_polygons])
    names = [name for name, _ in hydro]
    inland = [g if "RESERVOIR" in name else g.intersection(footprint) for name, g in hydro]
    components = [
        p for p in polygons(union_all(inland)) if p.area >= ALTERNATE_AREA_KM2 * 1e6
    ]
    components.sort(key=lambda p: (-p.area, *p.bounds))
    if shoreline.is_empty or not any(p.area >= PRIMARY_AREA_KM2 * 1e6 for p in components):
        raise ValueError("The source must provide a coastline and qualifying inland water")
    index = STRtree(inland)
    inventory = []
    for component in components:
        hits = index.query(component, predicate="intersects")
        ranked = sorted(
            ((inland[i].intersection(component).area, names[i]) for i in hits if names[i]),
            reverse=True,
        )
        area_km2 = component.area / 1e6
        inventory.append({
            "name_largest_named_overlap": ranked[0][1] if ranked and ranked[0][0] > 0 else "",
            "area_km2": area_km2,
            "primary": area_km2 >= PRIMARY_AREA_KM2,
            "easting_m": component.representative_point().x,
            "northing_m": component.representative_point().y,
        })
    return shoreline, components, pd.DataFrame(inventory)


def read_water(path):
    with gzip.open(path, "rb") as f:
        source = f.read()
    if hashlib.sha256(source).hexdigest() != MAP_SHA256:
        raise ValueError("The National Map Polygon input differs from the archived SLA download")
    coast, hydro = [], []
    for feature in json.loads(source)["features"]:
        properties = feature["properties"]
        name = properties["NAME"] or ""
        layer = properties["FOLDERPATH"]
        if feature["geometry"] is None or layer not in ("Layers/Coastal_Outlines", "Layers/Hydrographic"):
            continue
        if "MALAYSIA" in name:
            continue
        projected = transform_geom("EPSG:4326", GRID_CRS, feature["geometry"])
        geometry = union_all(list(polygons(make_valid(force_2d(shape(projected))))))
        if layer == "Layers/Coastal_Outlines":
            coast.append(geometry)
        else:
            hydro.append((name, geometry))
    return water_geometry(coast, hydro)


def distances(grid, coastline, components, inventory):
    centers = points(grid.easting_m, grid.northing_m)
    return grid.assign(
        coast_km=np.asarray(distance(centers, coastline)) / 1000,
        inland_km=np.asarray(distance(centers, union_all(
            [p for p, primary in zip(components, inventory.primary) if primary]
        ))) / 1000,
        inland_010_km=np.asarray(distance(centers, union_all(components))) / 1000,
    )

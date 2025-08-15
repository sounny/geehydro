# GEE Hydro

This repository contains resources for the Google Earth Engine Hydro project. Our goal is to create surface hydrology tools that run directly in the cloud using the Earth Engine JavaScript API. The initial focus is generating Flow Direction Rasters (FDR) from digital elevation data using the D8 algorithm.

## Contents

- `index.html` – project overview and reference links.
- `fdr` – example script implementing the D8 flow direction algorithm.
- `FlowDirectionLegend.png` – legend image describing direction codes.
- `poster_dem_inspector.pdf` – research poster describing early work.
- `flowacc.js` – compute flow accumulation within a drawn AOI.
- `watershed.js` – delineate an upstream watershed from a pour point.
- `streams.js` – derive a stream network from flow accumulation thresholds.
- `basins.js` – generate drainage basins by propagating stream labels upstream.
- `demfill.js` – fill DEM depressions to enforce drainage.

## Running the example

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Copy the contents of the `fdr` file into a new script.
3. Run the script to visualize the computed flow direction raster alongside the HydroSHEDS reference product.

The script demonstrates how to compute flow direction from the HydroSHEDS DEM and can serve as a starting point for additional hydrology tools such as stream and watershed delineation.

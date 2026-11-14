# Optional independent Earth Engine data check

The public-mirror audit can resolve the annual-coverage decision when even the most generous treatment of missing dates yields less than 80% eligible land. The following export can be used to measure coverage on the additional dates and independently compare the inputs where the two sources overlap.

1. Open [Earth Engine Code Editor](https://code.earthengine.google.com/) in a browser signed in to an enabled Earth Engine project.
2. Paste the contents of [export_earth_engine.js](../scripts/export_earth_engine.js) and click **Run**. Inspect the printed dates and native projection.
3. In **Tasks**, start `singapore_mod11a1_061_2025_native_daily`. Download the completed GeoTIFF from the `SG_LST_2025` Google Drive folder. Keep it in `data/raw/`.
4. Run the coverage audit with the file path and a new results directory.

```bash
.venv/bin/python -B scripts/audit_coverage.py \
  --earth-engine data/raw/singapore_mod11a1_061_2025_native_daily.tif \
  --output results/modis_earth_engine
```

The exporter includes the UTC calendar date, raw temperature DN and raw QA flag for every available granule. The audit verifies date order, alignment with the independently read sinusoidal grid, and the overlapping QA-passing values. It recomputes all monthly and quarterly counts from scratch. The older 23-band UTM aggregate does not have the raw daily layers needed by this verifier.

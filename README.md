# EDM code workspace

This folder is a cleaned project-oriented version of the scripts that were previously mixed under `code/`.

## Layout

- `R/`: shared path/bootstrap helpers.
- `projects/jp_flu_subtypes/`: Japan influenza subtype UIC, surrogate testing, and S-map analysis.
- `projects/okinawa_flu_weather/`: Okinawa subprefecture influenza-weather UIC/S-map analysis and figures.
- `projects/okinawa_id_weather/`: Okinawa infectious disease and weather UIC/S-map analysis.
- `projects/okinawa_infectious_diseases/`: subprefecture infectious-disease S-map summaries and heatmaps.
- `projects/archived_development/original_tree/`: original in-house scripts kept for traceability.
- `references/sample_code/`: external reference/sample repositories copied from the old `code/sample_code`.

## Running scripts

Run R from inside `edm_code` or any subdirectory below it. The bootstrap code locates `edm_code`, then sets the working directory to the parent workspace so existing relative paths such as `data/...` and `result/...` keep working.

If the data/result root moves, set:

```r
Sys.setenv(MICROBIOME_DYNAMICS_ROOT = "/path/to/microbiome dynamics")
```

## Review notes

- Absolute `setwd("/Users/...")` calls were removed from the main revised scripts.
- The newest integrated scripts in `code/EDM_infectious_RProject` were used as the canonical starting point where available.
- Older date-stamped experiments are kept in `legacy/` or the archived original tree rather than mixed with current scripts.
- External sample code is preserved under `references/` and should not be treated as project source code.

# Okinawa infectious diseases and weather

Canonical scripts:

- `scripts/01_uic_smap_okinawa_id_subpref.R`
- `scripts/02_uic_smap_okinawa_id_weather_parallel.R`
- `scripts/02_uic_smap_okinawa_id_weather_serial.R`

Supporting script:

- `scripts/function_UIC_smap.R`: converted from a free-standing code fragment into `run_uic_smap_pair()`.

Earlier integrated code is retained in `legacy/`. The interaction-network notebook is in `notebooks/`.

Main data/output paths are workspace-level `data/oknw/` and `result/oknw_ID/`.

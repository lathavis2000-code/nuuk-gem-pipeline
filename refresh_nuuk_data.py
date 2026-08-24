"""
Nuuk GEM -> Earth Engine monthly refresh pipeline.
-----------------------------------------------------------------
Fetches the latest real CTD data from the GEM MarineBasis API, reprocesses
it into the same monthly Surface/Mid/Deep aggregates the Nuuk Diagnostic
Tool already uses, and pushes it as a live Earth Engine Table Asset - so
the GEE tool can query current data instead of the hardcoded snapshot
baked into the script at build time.

Run manually:
    GEM_API_KEY=xxx EE_ASSET_ID=users/you/nuuk_monthly python refresh_nuuk_data.py

Run via GitHub Actions: see .github/workflows/refresh_nuuk_data.yml
(reads GEM_API_KEY and EE_SERVICE_ACCOUNT_KEY from repo secrets).

REAL LIMITS, disclosed directly:
- ee.batch.Export.table.toAsset() FAILS if an asset already exists at the
  target ID - this script deletes the old asset first (if present) before
  re-exporting, which is the standard pattern for a "refresh in place"
  table. There is a brief window where the asset does not exist.
- Export tasks are asynchronous. This script polls until COMPLETED/FAILED
  so a scheduled CI run can tell whether the refresh actually succeeded,
  not just that it was submitted.
- This script does NOT validate the fetched data beyond what pandas
  parsing already enforces (numeric coercion, date parsing). If GEM ever
  changes their CSV's column names or units, this will need updating -
  it will fail loudly (KeyError) rather than silently ingest garbage.
"""

import os
import sys
import time
import requests
import pandas as pd
import ee

# ---------------------------------------------------------------
# CONFIG - all from environment variables (GitHub Actions secrets,
# or your own shell) so no credentials ever live in this file.
# ---------------------------------------------------------------
GEM_API_KEY = os.environ.get("GEM_API_KEY")
GEM_URL = "https://api.g-e-m.dk/api/dataset/10.17897/KMEK-TK21/csv"
EE_ASSET_ID = os.environ.get("EE_ASSET_ID", "users/YOUR_USERNAME/nuuk_monthly_data")
EE_SERVICE_ACCOUNT_EMAIL = os.environ.get("EE_SERVICE_ACCOUNT_EMAIL")
EE_SERVICE_ACCOUNT_KEY_PATH = os.environ.get("EE_SERVICE_ACCOUNT_KEY_PATH", "ee_key.json")

VARS_MAP = {
    "turbidity": "Turbidity (FTU)",
    "fluorescence": "Fluorescence (app. \u00b5g/l)",
    "salinity": "Salinity",
    "temperature": "Temperature (C)",
    "density": "Density, Sigma-theta (kg m-3)",
}
DEPTH_BANDS = {
    "surface": lambda p: p < 20,
    "mid": lambda p: (p >= 20) & (p <= 300),
    "deep": lambda p: p > 300,
}


def fetch_gem_csv(dest_path="nuuk_ctd_latest.csv"):
    """Real, live pull from the GEM MarineBasis API - same endpoint your
    Streamlit app already used for its 'Fetch latest real data' button."""
    if not GEM_API_KEY:
        raise SystemExit("GEM_API_KEY environment variable is not set - cannot authenticate to GEM.")
    print(f"[fetch] Requesting {GEM_URL} ...")
    headers = {"X-API-KEY": GEM_API_KEY}
    r = requests.get(GEM_URL, headers=headers, timeout=120)
    r.raise_for_status()
    with open(dest_path, "wb") as f:
        f.write(r.content)
    print(f"[fetch] Saved {len(r.content):,} bytes to {dest_path}")
    return dest_path


def process(csv_path):
    """Identical monthly-aggregation logic to the Streamlit app and the
    GEE tool's original data extraction - Surface <20 dbar, Mid 20-300
    dbar, Deep >300 dbar, monthly means per variable."""
    print(f"[process] Reading {csv_path} ...")
    df = pd.read_csv(csv_path, sep="\t")
    df["datetime"] = pd.to_datetime(
        df["Date (YYYY-MM-DD)"] + " " + df["Time (HH:MM:SS)"], errors="coerce"
    )
    df = df.dropna(subset=["datetime"]).set_index("datetime").sort_index()
    df["Fluorescence (app. \u00b5g/l)"] = df["Fluorescence (app. \u00b5g/l)"].clip(lower=0)

    if df.empty:
        raise SystemExit("[process] No valid rows after parsing - GEM's CSV format may have changed.")

    print(f"[process] {len(df):,} real rows, {df.index.min()} to {df.index.max()}")

    all_months = pd.period_range(df.index.min(), df.index.max(), freq="M").astype(str)
    out = pd.DataFrame({"month": all_months})
    out = out.set_index("month")

    for band_name, cond in DEPTH_BANDS.items():
        sub = df[cond(df["Pressure (avr. Db)"])]
        for var_name, col in VARS_MAP.items():
            s = sub[col].resample("MS").mean()
            s.index = s.index.to_period("M").astype(str)
            out[f"{var_name}_{band_name}"] = s

    out = out.reset_index()
    n_real_values = out.drop(columns=["month"]).notna().sum().sum()
    print(f"[process] Built {len(out)} monthly rows, {int(n_real_values)} real (non-null) values total.")
    return out


def upload_to_ee(df, asset_id):
    """Deletes any existing asset at this ID, then exports the fresh
    monthly table as a real Earth Engine Table Asset. Small collection
    (a few hundred rows) - built and exported directly, no GCS staging
    step needed."""
    print(f"[ee] Authenticating with service account {EE_SERVICE_ACCOUNT_EMAIL} ...")
    if not EE_SERVICE_ACCOUNT_EMAIL:
        raise SystemExit("EE_SERVICE_ACCOUNT_EMAIL environment variable is not set.")
    credentials = ee.ServiceAccountCredentials(EE_SERVICE_ACCOUNT_EMAIL, EE_SERVICE_ACCOUNT_KEY_PATH)
    ee.Initialize(credentials)

    print(f"[ee] Checking for an existing asset at {asset_id} ...")
    try:
        ee.data.getAsset(asset_id)
        print(f"[ee] Existing asset found - deleting it (export fails on a pre-existing assetId).")
        ee.data.deleteAsset(asset_id)
    except ee.EEException:
        print(f"[ee] No existing asset - this is the first run, or it was already cleared.")

    features = []
    for _, row in df.iterrows():
        props = row.dropna().to_dict()  # real values only - no fabricated nulls sent as 0
        features.append(ee.Feature(None, props))
    fc = ee.FeatureCollection(features)
    print(f"[ee] Built a FeatureCollection with {len(features)} real monthly rows.")

    task = ee.batch.Export.table.toAsset(
        collection=fc,
        description="nuuk_monthly_refresh",
        assetId=asset_id,
    )
    task.start()
    print(f"[ee] Export task started: {task.id}. Polling for completion ...")

    while task.active():
        time.sleep(10)
        status = task.status()
        print(f"[ee] ... state={status.get('state')}")

    final_status = task.status()
    if final_status.get("state") != "COMPLETED":
        raise SystemExit(f"[ee] Export FAILED: {final_status}")
    print(f"[ee] Export COMPLETED. Live asset ready at: {asset_id}")


if __name__ == "__main__":
    csv_path = fetch_gem_csv()
    processed = process(csv_path)
    upload_to_ee(processed, EE_ASSET_ID)
    print("=== Monthly refresh complete. ===")

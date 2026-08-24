# Nuuk GEM -> Earth Engine monthly refresh pipeline

Three pieces, wired together. All three need the SAME asset ID to agree.

## 1. Get a real Earth Engine service account (one-time setup)

The GitHub Actions job can't do an interactive `ee.Authenticate()` login -
it needs a service account:

1. In Google Cloud Console (same project as your Earth Engine account),
   create a service account.
2. Grant it access to Earth Engine (add it as a "Writer" on your GEE
   Cloud Project, or register it directly at
   https://code.earthengine.google.com/register - "Service account").
3. Create a JSON key for that service account and download it.

## 2. Set GitHub repo secrets

In your repo -> Settings -> Secrets and variables -> Actions, add:

| Secret name                | Value                                              |
|-----------------------------|-----------------------------------------------------|
| `GEM_API_KEY`               | your real GEM MarineBasis API key                   |
| `EE_SERVICE_ACCOUNT_EMAIL`  | the service account's email (ends in `.gserviceaccount.com`) |
| `EE_SERVICE_ACCOUNT_KEY_JSON` | the ENTIRE contents of the downloaded JSON key file |
| `EE_ASSET_ID`               | e.g. `users/YOUR_USERNAME/nuuk_monthly_data`        |

## 3. Update the GEE script to match

In `nuuk_diagnostic_tool.js`, find this line near the top of the
`OPTIONAL LIVE DATA OVERRIDE` section:

```javascript
var LIVE_ASSET_ID = 'users/YOUR_USERNAME/nuuk_monthly_data';
```

Change it to the EXACT same value you put in the `EE_ASSET_ID` secret.
These two must match, or the tool will just keep using the hardcoded
snapshot (harmlessly - it won't error, it just won't find anything live).

## 4. Test it once, manually, before trusting the schedule

Push this repo to GitHub, then go to the Actions tab -> "Refresh Nuuk GEM
data -> Earth Engine asset" -> "Run workflow" (the `workflow_dispatch`
trigger lets you fire it on demand). Watch the log. If it fails, the
error will tell you which piece is misconfigured (GEM auth, EE auth, or
the asset export itself).

Once that manual run succeeds, open the GEE tool and check the sidebar's
"DATA SOURCES" section - it should say "LIVE asset merged (N months...)"
instead of "hardcoded snapshot."

## 5. Let it run monthly

After the manual test passes, the `cron: '0 3 1 * *'` schedule in the
workflow file takes over automatically - runs at 03:00 UTC on the 1st of
every month, matching GEM's own monthly CTD cast cadence. No further
action needed.

## What this does NOT do

- Does not touch the Coupling Engine's satellite-source path - that's
  still live MODIS-Aqua fetches, unrelated to this pipeline.
- Does not validate the fetched data beyond what pandas parsing already
  enforces. If GEM changes their CSV format, the script fails loudly
  (visible in the Action's log) rather than silently ingesting garbage.
- Does not delete or modify the hardcoded 2005-2024 snapshot already in
  the GEE script - it only adds to / refreshes on top of it in memory,
  every time the tool loads.

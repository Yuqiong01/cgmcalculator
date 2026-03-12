# cgmcalculator Changelog

All notable changes to this package will be documented in this file.

---

## v1.0.2 (2026-03-12)

### New features
- Added `famm_step_mins` argument to `agpanalyze()` to allow flexible response-grid resolution for functional regression.
- The response matrix used for `refund::pffr()` can now be constructed at user-defined temporal resolutions (e.g., 60, 30, or 15 minutes) instead of the previous fixed hourly grid.

### Improvements
- The default behavior of `pffr_yind` now follows the actual response-grid centers determined by `famm_step_mins`.
- Updated documentation for `agpanalyze()`, including parameter descriptions, return values, and examples.

### Notes
- Internal object names such as `MIMS_hour_mat` are retained for backward compatibility, although the response grid is no longer restricted to hourly sampling.
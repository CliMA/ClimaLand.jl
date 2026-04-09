"""
Quick diagnostic plots of the four CalLMIP NetCDF submission files.

Usage (from any node — pure Python, no Julia needed):
    python3 experiments/callmip_uq_dk_sor/plot_callmip_netcdf.py

Outputs: experiments/callmip_uq_dk_sor/output_evaluation/callmip_nc_*.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc
from datetime import datetime, timedelta

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NC_DIR     = os.path.join(SCRIPT_DIR, "callmip_output")
OUT_DIR    = os.path.join(SCRIPT_DIR, "output_evaluation")
FLUX_NC    = os.path.join(SCRIPT_DIR, "..", "..",
                          "DK_Sor",
                          "DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Unit conversions ──────────────────────────────────────────────────────────
kgC_s_to_gC_d = 1e3 * 86400   # kg C m⁻² s⁻¹ → gC m⁻² d⁻¹

# ── Load the four NetCDF files ────────────────────────────────────────────────
FILES = {
    "Cal_Prior":     "ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Prior.nc",
    "Cal_Posterior": "ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Cal_Posterior.nc",
    "Val_Prior":     "ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Temporal_Prior.nc",
    "Val_Posterior": "ClimaLand.CalLMIP1.0_Expt1_DK-Sor_Val_Temporal_Posterior.nc",
}

data = {}
dates = {}
for label, fname in FILES.items():
    path = os.path.join(NC_DIR, fname)
    ds = nc.Dataset(path, "r")
    # time axis
    t_raw = ds.variables["time"][:]
    t_units = ds.variables["time"].units   # e.g. "days since 2003-01-02"
    t0_str = t_units.split("since ")[-1].strip()
    try:
        t0 = datetime.strptime(t0_str, "%Y-%m-%d")
    except ValueError:
        t0 = datetime.strptime(t0_str[:10], "%Y-%m-%d")
    dates[label] = [t0 + timedelta(days=float(d)) for d in t_raw]

    d = {}
    for v in ds.variables:
        if v == "time":
            continue
        arr = np.array(ds.variables[v][:], dtype=float)
        arr[arr > 1e30] = np.nan
        arr[arr < -1e30] = np.nan
        d[v] = arr
    data[label] = d
    # Print variable inventory
    print(f"\n=== {label} ({fname}) ===")
    print(f"  n_days = {len(dates[label])}")
    for v, arr in d.items():
        finite = arr[np.isfinite(arr)]
        if len(finite) == 0:
            print(f"  {v:30s}  ALL NaN")
        else:
            print(f"  {v:30s}  min={finite.min():.4g}  max={finite.max():.4g}  "
                  f"mean={finite.mean():.4g}  n_nan={np.sum(np.isnan(arr))}")
    ds.close()

# ── Load FLUXNET observations ─────────────────────────────────────────────────
obs_nee = obs_qle = obs_qh = obs_dates = None
if os.path.isfile(FLUX_NC):
    ds_flux = nc.Dataset(FLUX_NC, "r")
    t_raw   = ds_flux.variables["time"][:]
    t_units = ds_flux.variables["time"].units
    t0_str = t_units.split("since ")[-1].strip()[:10]
    t0 = datetime.strptime(t0_str, "%Y-%m-%d")
    # Units may be seconds, hours, or days since epoch
    if "second" in t_units:
        obs_dates = [t0 + timedelta(seconds=float(d)) for d in t_raw]
    elif "hour" in t_units:
        obs_dates = [t0 + timedelta(hours=float(d)) for d in t_raw]
    else:  # days
        obs_dates = [t0 + timedelta(days=float(d)) for d in t_raw]

    def _load(key):
        v = ds_flux.variables[key][:]
        v = np.ma.filled(v, np.nan).astype(float)
        v[np.abs(v) > 1e10] = np.nan
        return v

    obs_nee = _load("NEE_daily")        # gC m⁻² d⁻¹
    obs_qle = _load("Qle_daily")        # W m⁻²
    obs_qh  = _load("Qh_daily")         # W m⁻²
    ds_flux.close()
    print(f"\nFLUXNET loaded: {len(obs_dates)} days")

# ── Helper: combine cal + val arrays onto a single timeline ───────────────────
def combine(var, unit_scale=1.0):
    """Return (all_dates, prior, posterior) arrays over full sim period."""
    d_all, p_all, q_all = [], [], []
    for period in ("Cal", "Val"):
        prior_key = f"{period}_Prior"
        post_key  = f"{period}_Posterior"
        d_all.extend(dates[prior_key])
        prior_arr = data[prior_key].get(var, np.full(len(dates[prior_key]), np.nan))
        post_arr  = data[post_key].get(var, np.full(len(dates[post_key]),  np.nan))
        p_all.extend(prior_arr * unit_scale)
        q_all.extend(post_arr  * unit_scale)
    return np.array(d_all), np.array(p_all), np.array(q_all)

def combine_unc(var, pct, unit_scale=1.0):
    """Return percentile band over full period (posterior only)."""
    key = f"{var}_{pct}"
    out = []
    for period in ("Cal", "Val"):
        pk = f"{period}_Posterior"
        arr = data[pk].get(key, np.full(len(dates[pk]), np.nan))
        out.extend(arr * unit_scale)
    return np.array(out)

def obs_in_range(all_dates):
    """Subset obs to the model date range."""
    if obs_dates is None:
        return None, None, None
    d0, d1 = min(all_dates), max(all_dates)
    mask = [(d >= d0) and (d <= d1) for d in obs_dates]
    mask = np.array(mask)
    return (np.array(obs_dates)[mask],
            obs_nee[mask],
            obs_qle[mask],
            obs_qh[mask])

# ── Figure 1: Full time series — NEE, Qle, Qh ────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(14, 11), sharex=True)
fig.suptitle("ClimaLand DK-Sor — CalLMIP NetCDF submission: full time series",
             fontsize=12, fontweight="bold")

VARS = [
    ("NEE",  kgC_s_to_gC_d,  "NEE (gC m⁻² d⁻¹)"),
    ("Qle",  1.0,             "Qle (W m⁻²)"),
    ("Qh",   1.0,             "Qh (W m⁻²)"),
]

all_dates, prior_nee, post_nee = combine("NEE", kgC_s_to_gC_d)
# cal/val split date
cal_end = max(dates["Cal_Prior"])

for ax, (var, scale, ylabel) in zip(axes, VARS):
    all_dates_v, prior_v, post_v = combine(var, scale)

    # Uncertainty bands (posterior)
    p05 = combine_unc(var, "p05", scale)
    p95 = combine_unc(var, "p95", scale)
    p50 = combine_unc(var, "p50", scale)

    # Check if bands are non-trivial (not all NaN)
    has_bands = not np.all(np.isnan(p05))

    if has_bands:
        ax.fill_between(all_dates_v, p05, p95, alpha=0.25, color="steelblue",
                        label="CES 90% band")
        ax.plot(all_dates_v, p50, color="steelblue", lw=0.8, alpha=0.7,
                linestyle="--", label="CES median")

    ax.plot(all_dates_v, prior_v, color="grey",   lw=1.0, label="Prior")
    ax.plot(all_dates_v, post_v,  color="darkorange", lw=1.2, label="Posterior (EKI)")

    # FLUXNET obs
    obs_info = obs_in_range(all_dates_v)
    if obs_info[0] is not None:
        od, on, oq, oh = obs_info
        obs_arr = on if var == "NEE" else (oq if var == "Qle" else oh)
        ax.scatter(od, obs_arr, s=3, c="black", alpha=0.4, zorder=5, label="FLUXNET")

    # Cal/Val divider
    ax.axvline(cal_end, color="black", linestyle=":", lw=1.2)
    ax.text(cal_end, ax.get_ylim()[1] if ax.get_ylim()[1] != 1.0 else 1,
            " Val→", fontsize=8, va="top", ha="left")

    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(loc="upper left", fontsize=8, ncol=4)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(True, alpha=0.3)

axes[-1].set_xlabel("Date")
plt.tight_layout()
out_path = os.path.join(OUT_DIR, "callmip_nc_timeseries.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {out_path}")

# ── Figure 2: Seasonal cycle ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle("ClimaLand DK-Sor — CalLMIP NetCDF: mean annual cycle",
             fontsize=12, fontweight="bold")

for ax, (var, scale, ylabel) in zip(axes, VARS):
    all_dates_v, prior_v, post_v = combine(var, scale)
    p05 = combine_unc(var, "p05", scale)
    p95 = combine_unc(var, "p95", scale)
    p50 = combine_unc(var, "p50", scale)
    months = np.array([d.month for d in all_dates_v])

    prior_mc  = [np.nanmean(prior_v[months == m]) for m in range(1, 13)]
    post_mc   = [np.nanmean(post_v[months == m])  for m in range(1, 13)]
    p05_mc    = [np.nanmean(p05[months == m])      for m in range(1, 13)]
    p95_mc    = [np.nanmean(p95[months == m])      for m in range(1, 13)]
    p50_mc    = [np.nanmean(p50[months == m])      for m in range(1, 13)]

    mo = np.arange(1, 13)
    has_bands = not np.all(np.isnan(p05))
    if has_bands:
        ax.fill_between(mo, p05_mc, p95_mc, alpha=0.25, color="steelblue",
                        label="CES 90% band")
        ax.plot(mo, p50_mc, color="steelblue", lw=1.0, linestyle="--",
                label="CES median")

    ax.plot(mo, prior_mc, color="grey",       lw=1.5, linestyle="--", label="Prior")
    ax.plot(mo, post_mc,  color="darkorange", lw=2.0, label="Posterior (EKI)")

    # FLUXNET obs seasonal cycle
    if obs_dates is not None:
        obs_arr = (obs_nee if var=="NEE" else obs_qle if var=="Qle" else obs_qh)
        obs_months = np.array([d.month for d in obs_dates])
        obs_mc = [np.nanmean(obs_arr[obs_months == m]) for m in range(1, 13)]
        ax.plot(mo, obs_mc, color="black", lw=1.5, linestyle="-.", label="FLUXNET")

    ax.set_xticks(mo)
    ax.set_xticklabels(["J","F","M","A","M","J","J","A","S","O","N","D"])
    ax.set_ylabel(ylabel)
    ax.set_title(var)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="black", lw=0.5)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "callmip_nc_seasonal.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

# ── Figure 3: All ALMA variables overview (posterior only) ────────────────────
ALMA_GRID = [
    ("GPP",           kgC_s_to_gC_d, "GPP (gC m⁻² d⁻¹)"),
    ("Reco",          kgC_s_to_gC_d, "Reco (gC m⁻² d⁻¹)"),
    ("TVeg",          86400*1e3,      "TVeg (g m⁻² d⁻¹)"),
    ("ESoil",         86400*1e3,      "ESoil (g m⁻² d⁻¹)"),
    ("Qg",            1.0,            "Qg (W m⁻²)"),
    ("AvgSurfT",      1.0,            "AvgSurfT (K)"),
    ("SoilMoist",     1.0,            "SoilMoist (kg m⁻²)"),
    ("LAI",           1.0,            "LAI (m² m⁻²)"),
    ("TotAbovBioMass",1.0,            "TotAbovBioMass (kg C m⁻²)"),
    ("TotSoilCarb",   1.0,            "TotSoilCarb (kg C m⁻²)"),
]

fig, axes = plt.subplots(5, 2, figsize=(14, 18), sharex=False)
fig.suptitle("ClimaLand DK-Sor — CalLMIP NetCDF: ALMA variables (prior vs posterior)",
             fontsize=12, fontweight="bold")
axes_flat = axes.ravel()

for ax, (var, scale, ylabel) in zip(axes_flat, ALMA_GRID):
    all_dates_v, prior_v, post_v = combine(var, scale)
    ax.plot(all_dates_v, prior_v, color="grey",       lw=0.8, alpha=0.8, label="Prior")
    ax.plot(all_dates_v, post_v,  color="darkorange", lw=1.0, label="Posterior (EKI)")
    ax.axvline(cal_end, color="black", linestyle=":", lw=1.0)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(var, fontsize=10)
    ax.legend(fontsize=8)
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(True, alpha=0.3)
    # Flag all-NaN variables
    if np.all(np.isnan(post_v)):
        ax.text(0.5, 0.5, "ALL NaN", transform=ax.transAxes,
                ha="center", va="center", color="red", fontsize=14, fontweight="bold")

# Hide unused subplot
for ax in axes_flat[len(ALMA_GRID):]:
    ax.set_visible(False)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "callmip_nc_alma_vars.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

# ── Figure 4: CES uncertainty bands zoom — show how flat they are ─────────────
fig, axes = plt.subplots(3, 1, figsize=(14, 10))
fig.suptitle("ClimaLand DK-Sor — CalLMIP Posterior NetCDF: CES uncertainty bands detail",
             fontsize=12, fontweight="bold")

for ax, (var, scale, ylabel) in zip(axes, VARS):
    all_dates_v, prior_v, post_v = combine(var, scale)
    p05 = combine_unc(var, "p05", scale)
    p50 = combine_unc(var, "p50", scale)
    p95 = combine_unc(var, "p95", scale)

    ax.fill_between(all_dates_v, p05, p95, alpha=0.3, color="steelblue",
                    label="CES 90% band (p05–p95)")
    ax.plot(all_dates_v, p50, color="steelblue",   lw=1.0, linestyle="--", label="CES median (p50)")
    ax.plot(all_dates_v, post_v, color="darkorange", lw=1.5, label="EKI optimal")

    if obs_dates is not None:
        od, on, oq, oh = obs_in_range(all_dates_v)
        obs_arr = on if var == "NEE" else (oq if var == "Qle" else oh)
        ax.scatter(od, obs_arr, s=3, c="black", alpha=0.35, zorder=5, label="FLUXNET")

    ax.axvline(cal_end, color="black", linestyle=":", lw=1.0)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(fontsize=8, ncol=4)
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

axes[-1].set_xlabel("Date")
plt.tight_layout()
out_path = os.path.join(OUT_DIR, "callmip_nc_ces_bands.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")

print("\n=== All plots saved to output_evaluation/ ===")

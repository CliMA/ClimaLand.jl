"""Fig 5: prior model vs obs monthly climatology, all 21 calibration sites.
Three PNGs of 7 sites each (rows) x NEE/LE/H (cols). Run after PRIORS_DONE."""
import numpy as np, datetime as dt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob, os

BLUE, ORANGE = "#2a78d6", "#eb6834"
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "text.color": INK, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "font.size": 8.5, "axes.titlesize": 9,
    "axes.spines.top": False, "axes.spines.right": False})
D = "/home/renatob/ClimaLand.jl/experiments/callmip_phase1b/output_deliverables"
FLUXDIR = "/net/sampo/data1/renatob/callmip_artifacts/callmip_phase1_20260818/Data/Phase1b"
FIG = "/home/renatob/ClimaLand.jl/experiments/callmip_phase1b/figures"

def site_clim(site):
    nc = Dataset(glob.glob(f"{D}/ClimaLand.v1_Phase1b_Scen1_{site}_Cal_Prior.nc")[0])
    y0 = int(nc["time"].units.split("since ")[1][:4]) + 1
    d0 = dt.date(y0 - 1, 12, 31)
    mm = np.array([(d0 + dt.timedelta(days=int(x))).month for x in np.array(nc["time"][:]).astype(int)])
    def mc(var, scale=1.0, sign=1.0):
        x = np.ma.filled(nc[var][:], np.nan).flatten(); x[x >= 1e37] = np.nan
        x = sign * x * scale
        return np.array([np.nanmean(x[mm == m]) for m in range(1, 13)])
    model = {"NEE": mc("nep", 1e3 * 86400, -1.0), "LE": mc("hfls"), "H": mc("hfss")}
    nc.close()
    fo = Dataset(glob.glob(f"{FLUXDIR}/{site}_daily_aggregated_*_Flux.nc")[0])
    fy0 = int(str(fo["time"].units).split("since ")[1][:4])
    fm = np.array([(dt.date(fy0, 1, 1) + dt.timedelta(seconds=float(s))).month for s in fo["time"][:].flatten()])
    def oc(var):
        x = np.ma.filled(fo[var][:], np.nan).flatten(); x[x >= 1e19] = np.nan
        return np.array([np.nanmean(x[fm == m]) if np.isfinite(x[fm == m]).sum() >= 5 else np.nan for m in range(1, 13)])
    obs = {"NEE": oc("NEE_daily"), "LE": oc("Qle_daily"), "H": oc("Qh_daily")}
    fo.close()
    return model, obs

SITES = sorted(os.path.basename(f).split("_")[3] for f in glob.glob(f"{D}/ClimaLand.v1_Phase1b_Scen1_*_Cal_Prior.nc"))
print(f"{len(SITES)} sites with prior files")
months = np.arange(1, 13)
groups = [SITES[i:i + 7] for i in range(0, len(SITES), 7)]
for gi, group in enumerate(groups):
    fig, axes = plt.subplots(len(group), 3, figsize=(10, 1.75 * len(group) + 1.2))
    axes = np.atleast_2d(axes)
    for r, site in enumerate(group):
        model, obs = site_clim(site)
        for c, k in enumerate(("NEE", "LE", "H")):
            ax = axes[r, c]
            ax.plot(months, model[k], color=BLUE, lw=1.8, marker="o", ms=3, mec=SURF, mew=0.5)
            ax.plot(months, obs[k], color=ORANGE, lw=0, marker="o", ms=4.5, mec=SURF, mew=0.7)
            if k == "NEE": ax.axhline(0, color=INK2, lw=0.7)
            if r == 0: ax.set_title(("NEE (gC m$^{-2}$ d$^{-1}$)", "LE (W m$^{-2}$)", "H (W m$^{-2}$)")[c], loc="left")
            if c == 0: ax.set_ylabel(site, fontsize=9, fontweight="bold")
            ax.set_xticks([1, 7]); ax.set_xticklabels(["Jan", "Jul"] if r == len(group) - 1 else ["", ""])
    fig.suptitle(f"Prior model (blue line) vs tower obs (orange dots) — monthly climatology, sites {group[0]}…{group[-1]}\nuncalibrated default parameters; obs months with <5 valid days omitted",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{FIG}/fig5{'abc'[gi]}_prior_vs_obs_all_sites.png", dpi=140)
    plt.close(fig)
    print(f"fig5{'abc'[gi]} written ({len(group)} sites)")

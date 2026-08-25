"""Fig 5: model vs obs monthly climatology, all 21 calibration sites.
Three PNGs of 7 sites each (rows) x NEE/LE/H (cols).
Each panel is annotated with R^2 (squared Pearson r) and RMSE computed on the
MONTHLY TIME SERIES (every month with >=5 valid obs days is one point; model
monthly means are all-day) -- the same metric as the prior-vs-posterior skill
gate, so the numbers are directly comparable across calibration stages.
Usage: python3 make_fig5_all_sites.py [Prior|Posterior]"""
import sys
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

STAGE = sys.argv[1] if len(sys.argv) > 1 else "Prior"

def monthly_metrics(model_daily, model_ym, obs_daily, obs_ym, obs_ok):
    """R^2 and RMSE on paired monthly means over all (year, month) cells."""
    cells = sorted(set(map(tuple, model_ym)) & set(map(tuple, obs_ym[obs_ok])))
    m, o = [], []
    for ym in cells:
        msel = (model_ym[:, 0] == ym[0]) & (model_ym[:, 1] == ym[1])
        osel = (obs_ym[:, 0] == ym[0]) & (obs_ym[:, 1] == ym[1]) & obs_ok
        if osel.sum() >= 5 and np.isfinite(model_daily[msel]).any():
            m.append(np.nanmean(model_daily[msel]))
            o.append(np.nanmean(obs_daily[osel]))
    m, o = np.array(m), np.array(o)
    if len(m) < 3:
        return np.nan, np.nan, 0
    r = np.corrcoef(m, o)[0, 1]
    return r * r, float(np.sqrt(np.mean((m - o) ** 2))), len(m)

def site_clim(site):
    nc = Dataset(glob.glob(f"{D}/ClimaLand.v1_Phase1b_Scen1_{site}_Cal_{STAGE}.nc")[0])
    y0 = int(nc["time"].units.split("since ")[1][:4]) + 1
    d0 = dt.date(y0 - 1, 12, 31)
    mm = np.array([(d0 + dt.timedelta(days=int(x))).month for x in np.array(nc["time"][:]).astype(int)])
    dts = np.array([d0 + dt.timedelta(days=int(x)) for x in np.array(nc["time"][:]).astype(int)])
    model_ym = np.array([[d.year, d.month] for d in dts])
    def mdaily(var, scale=1.0, sign=1.0):
        x = np.ma.filled(nc[var][:], np.nan).flatten(); x[x >= 1e37] = np.nan
        return sign * x * scale
    md = {"NEE": mdaily("nep", 1e3 * 86400, -1.0), "LE": mdaily("hfls"), "H": mdaily("hfss")}
    model = {k: np.array([np.nanmean(v[mm == m]) for m in range(1, 13)]) for k, v in md.items()}
    nc.close()
    fo = Dataset(glob.glob(f"{FLUXDIR}/{site}_daily_aggregated_*_Flux.nc")[0])
    fy0 = int(str(fo["time"].units).split("since ")[1][:4])
    fdts = np.array([dt.date(fy0, 1, 1) + dt.timedelta(seconds=float(s)) for s in fo["time"][:].flatten()])
    fm = np.array([d.month for d in fdts])
    obs_ym = np.array([[d.year, d.month] for d in fdts])
    def odaily(var):
        x = np.ma.filled(fo[var][:], np.nan).flatten(); x[x >= 1e19] = np.nan
        return x
    od = {"NEE": odaily("NEE_daily"), "LE": odaily("Qle_daily"), "H": odaily("Qh_daily")}
    obs = {k: np.array([np.nanmean(v[fm == m]) if np.isfinite(v[fm == m]).sum() >= 5 else np.nan for m in range(1, 13)]) for k, v in od.items()}
    fo.close()
    stats = {k: monthly_metrics(md[k], model_ym, od[k], obs_ym, np.isfinite(od[k])) for k in md}
    return model, obs, stats

SITES = sorted(os.path.basename(f).split("_")[3] for f in glob.glob(f"{D}/ClimaLand.v1_Phase1b_Scen1_*_Cal_{STAGE}.nc"))
print(f"{len(SITES)} sites with prior files")
months = np.arange(1, 13)
all_stats = {}
groups = [SITES[i:i + 7] for i in range(0, len(SITES), 7)]
for gi, group in enumerate(groups):
    fig, axes = plt.subplots(len(group), 3, figsize=(10, 1.75 * len(group) + 1.2))
    axes = np.atleast_2d(axes)
    for r, site in enumerate(group):
        model, obs, stats = site_clim(site)
        all_stats[site] = stats
        for c, k in enumerate(("NEE", "LE", "H")):
            ax = axes[r, c]
            ax.plot(months, model[k], color=BLUE, lw=1.8, marker="o", ms=3, mec=SURF, mew=0.5)
            ax.plot(months, obs[k], color=ORANGE, lw=0, marker="o", ms=4.5, mec=SURF, mew=0.7)
            if k == "NEE": ax.axhline(0, color=INK2, lw=0.7)
            r2, rmse, nm = stats[k]
            txt = f"R$^2$={r2:.2f} RMSE={rmse:.2g} (n={nm})" if np.isfinite(r2) else "n<3 months"
            ax.text(0.985, 0.97, txt, transform=ax.transAxes, ha="right", va="top",
                    fontsize=6.5, color=INK2)
            if r == 0: ax.set_title(("NEE (gC m$^{-2}$ d$^{-1}$)", "LE (W m$^{-2}$)", "H (W m$^{-2}$)")[c], loc="left")
            if c == 0: ax.set_ylabel(site, fontsize=9, fontweight="bold")
            ax.set_xticks([1, 7]); ax.set_xticklabels(["Jan", "Jul"] if r == len(group) - 1 else ["", ""])
    fig.suptitle(f"{STAGE} model (blue line) vs tower obs (orange dots) — monthly climatology, sites {group[0]}…{group[-1]}\nR$^2$/RMSE computed on monthly time series (n = valid months); obs months with <5 valid days omitted",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{FIG}/fig5{'abc'[gi]}_{STAGE.lower()}_vs_obs_all_sites.png", dpi=140)
    plt.close(fig)
    print(f"fig5{'abc'[gi]} written ({len(group)} sites)")

with open(f"{FIG}/skill_{STAGE.lower()}.csv", "w") as io:
    io.write("site,flux,r2,rmse,n_months\n")
    for site in SITES:
        for k in ("NEE", "LE", "H"):
            r2, rmse, nm = all_stats[site][k]
            io.write(f"{site},{k},{r2:.4f},{rmse:.4f},{nm}\n")
print(f"skill_{STAGE.lower()}.csv written")

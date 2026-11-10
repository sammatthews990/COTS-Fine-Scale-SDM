"""Combine dated ADRIA Figures 1/4 and reflow the other figures for portrait pages.

Run through script 07. Bands and observation aggregation follow the source
plot_calibration_report.py, with COTS/coral aligned by reef and shared legends.
All raw members are checked against the frozen source metadata before plotting.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd
import numpy as np


def checked(path, expected):
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Frozen member changed: {path}")
    return path


def reflow_supplements(a, meta, tr, obs):
    """Replot source summaries at portrait width with readable eight-point type."""
    colors = ("#0d366b", "#5b2a86")
    grid = "#e1e0d9"
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8})
    extra = []
    for key, run_key in (("trajectories_extra.csv", "run"), ("cull/trajectories_extra.csv", "cull_run")):
        path = a.adria / meta[run_key].replace("\\", "/") / "trajectories_extra.csv"
        extra.append(pd.read_csv(checked(path, meta["inputs"][key])))
    target_path = a.adria / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
    targets = pd.read_csv(checked(target_path, meta["inputs"]["targets_reef_year.csv"]))
    all_tr = [pd.concat([t, e]) for t, e in zip(tr, extra)]

    def style(ax):
        ax.grid(color=grid, linewidth=.6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)

    def save(fig, stem, bottom=.15, top=.96):
        fig.tight_layout(rect=(0, bottom, 1, top))
        for ext in ("png", "pdf"):
            fig.savefig(a.output / f"{stem}.{ext}", dpi=300, facecolor="white")
        plt.close(fig)

    def calibration_lines(ax, reef, scale=1, recent=False):
        for frame, color, ls in zip(all_tr, colors, ("-", "--")):
            d = frame[frame.reef_name == reef].pivot(index="year", columns="treatment", values="adults_ha") * scale
            d = d.loc[2015 if recent else 1992:2026]
            base = d[[n for n in d if "_sx" not in n and "_sh" not in n]]
            lo = .05 if recent else 0
            ax.fill_between(d.index, d.min(axis=1).clip(lower=lo), d.max(axis=1), color=color, alpha=.10, linewidth=0)
            ax.fill_between(base.index, base.min(axis=1).clip(lower=lo), base.max(axis=1), color=color, alpha=.22, linewidth=0)
            ax.plot(base.index, base.mean(axis=1).clip(lower=lo), color=color, linestyle=ls, linewidth=1.2)
        style(ax)

    cal_legend = [Line2D([], [], color=colors[0], label="No recorded culling"),
                  Line2D([], [], color=colors[1], linestyle="--", label="Recorded culling"),
                  Patch(color="gray", alpha=.22, label="Five-seed range"),
                  Patch(color="gray", alpha=.10, label="Parameter range")]
    reefs = list(dict.fromkeys(tr[0].reef_name)) + ["Martin Reef"]
    fig, axes = plt.subplots(3, 2, figsize=(6.26, 7.4))
    for ax, reef in zip(axes.flat, reefs):
        calibration_lines(ax, reef, recent=True)
        for src, marker in (("SALAD", "o"), ("cull_converted", "s")):
            d = targets[(targets.target_reef == reef) & (targets.target_source == src) & (targets.target_density_ha > 0)]
            ax.scatter(d.report_year, d.target_density_ha, marker=marker, color="#eb6834", s=16, zorder=5)
        ax.set(title=reef, xlim=(2014.5, 2026.5), ylim=(.1, 300), yscale="log", xlabel="Year", ylabel="Adult COTS ha$^{-1}$")
        ax.set_xticks([2015, 2020, 2025])
    axes.flat[-1].set_visible(False)
    fig.legend(handles=cal_legend + [Line2D([], [], color="#eb6834", marker="o", linestyle="none", label="SALAD diver density"),
               Line2D([], [], color="#eb6834", marker="s", linestyle="none", label="Cull-derived density")],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    save(fig, "fig23_metacommunity_recent_density", bottom=.13)

    others = {"Martin Reef": "Martin Reef", "Linnet Reef": "Linnet Reef", "Yonge Reef": "Yonge Reef",
              "Carter Reef": "Carter Reef", "No Name Reef": "No Name Reef", "Ribbon No 10 Reef": None}
    fig, axes = plt.subplots(3, 2, figsize=(6.26, 7.2))
    for ax, (reef, oname) in zip(axes.flat, others.items()):
        calibration_lines(ax, reef, scale=.015)
        if oname:
            o = obs[obs.reef_name == oname].groupby("year").cotsptow.mean()
            ax.scatter(o.index, o, s=12, color="#1f1e1d", zorder=5)
        ax.axhline(.22, color="gray", linestyle=":", linewidth=.8)
        ax.set(title=reef, xlim=(1991.5, 2026.5), xlabel="Year", ylabel="COTS per manta tow")
        ax.set_ylim(bottom=0)
        ax.set_xticks([1995, 2005, 2015, 2025])
    fig.legend(handles=cal_legend + [Line2D([], [], color="#1f1e1d", marker="o", linestyle="none", label="LTMP manta"),
               Line2D([], [], color="gray", linestyle=":", label="0.22 COTS per tow")],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    save(fig, "fig24_metacommunity_additional_reefs", bottom=.13)

    pm = json.loads(a.projections.read_text())
    frames, culls = [], []
    for relative, sha in pm["inputs"].items():
        p = checked(a.adria / "sandbox/calibration/runs" / relative, sha)
        d = pd.read_csv(p)
        if p.name == "culls.csv":
            culls.append(d[(d.scope == "domain") & d.year.between(2027, 2050)])
        else:
            frames.append(d[d.reef_name.isin(pm["reefs"]) & d.year.between(2015, 2050)])
    pattern = re.compile(r"_pj(\d{3})g(\d{2})W_cullreef_(fc[NFT]\d{0,3})_xc(_rsAJ)?_aprmak_")
    traits = []
    for arm in pd.concat(frames).treatment.unique():
        m = pattern.search(arm)
        if m and not m.group(4):
            traits.append({"treatment": arm, "ssp": "ssp"+m.group(1), "cull": m.group(3)})
    traits = pd.DataFrame(traits)
    assert len(traits) == pm["arms_used"]
    df = pd.concat(frames).merge(traits, on="treatment")
    cull = pd.concat(culls).merge(traits, on="treatment")
    cull["removed"] = cull.removed_pre_adults + cull.removed_post_adults
    removals = cull.groupby("treatment").removed.sum()
    strategies = {"fcN": ("No further culling", "#2a78d6"), "fcF100": ("Recent effort ×1", "#eb6834"),
                  "fcF300": ("Recent effort ×3", "#1baf7a"), "fcT100": ("Targeted effort ×1", "#eda100")}
    for ssp in ("ssp245", "ssp370"):
        fig, axes = plt.subplots(5, 2, figsize=(6.26, 8.4))
        for row, reef in enumerate(pm["reefs"]):
            for col, var in enumerate(("adults_ha", "coral_cover")):
                ax = axes[row, col]
                for code, (label, color) in strategies.items():
                    g = df[(df.ssp == ssp) & (df.cull == code) & (df.reef_name == reef)].groupby("year")[var]
                    ax.plot(g.median().index, g.median(), color=color, linewidth=1.1)
                    if code in ("fcN", "fcF300"):
                        ax.fill_between(g.median().index, g.quantile(.25), g.quantile(.75), color=color, alpha=.15, linewidth=0)
                ax.axvline(2026.5, color="gray", linestyle=":", linewidth=.8)
                if col == 0:
                    ax.axhline(.22/.015, color="gray", linewidth=.7)
                ax.set(title=reef.replace(" Reef", ""), xlim=(2015, 2050))
                ax.set_xticks([2020, 2030, 2040, 2050])
                style(ax)
        axes[0, 0].set_ylabel("Adults ha$^{-1}$")
        axes[0, 1].set_ylabel("Coral cover (fraction)")
        for ax in axes[-1]:
            ax.set_xlabel("Year")
        fig.legend(handles=[Line2D([], [], color=c, label=l) for l, c in strategies.values()],
                   loc="lower center", ncol=2, frameon=False, fontsize=8)
        save(fig, "fig25_metacommunity_projection_ssp245" if ssp == "ssp245" else "fig26_metacommunity_projection_ssp370", bottom=.08, top=.99)

    rows = []
    for (ssp, code, arm, reef), d in df.groupby(["ssp", "cull", "treatment", "reef_name"]):
        future = d[d.year.between(2027, 2050)]
        rows.append({"ssp": ssp, "cull": code, "arm": arm, "reef": reef,
                     "years": (future.adults_ha >= .22/.015).sum(),
                     "peak": d[d.year.between(2027, 2031)].adults_ha.max(),
                     "coral": future.coral_cover.mean(), "removed": removals[arm]/1000})
    rows = pd.DataFrame(rows)
    keys = ["ssp", "cull", "reef"]
    actual = rows.groupby(keys)[["years", "peak", "coral"]].median().sort_index()
    expected = pd.read_csv(a.projections.parent / "culling_effects_2027_2050.csv").rename(columns={"reef_name": "reef"})
    expected = expected.set_index(keys).sort_index()
    assert actual.index.equals(expected.index)
    assert np.allclose(actual.to_numpy(), expected[["years_in_outbreak", "peak_2027_31_ha", "coral_mean"]].to_numpy(), rtol=1e-12, atol=1e-12)
    print("Reflowed projection summaries exactly reproduce the frozen source table.")
    means = rows.groupby(["ssp", "cull", "arm"])[["years", "peak", "coral", "removed"]].mean().reset_index()
    order = ["fcN", "fcF100", "fcF150", "fcF200", "fcF300", "fcT100", "fcT150", "fcT200", "fcT300"]
    labels = ["None", "Recent ×1", "Recent ×1.5", "Recent ×2", "Recent ×3", "Targeted ×1", "Targeted ×1.5", "Targeted ×2", "Targeted ×3"]
    fig, axes = plt.subplots(2, 2, figsize=(6.26, 5.7))
    for ax, var, title in zip(axes.flat, ("years", "peak", "coral", "removed"),
                             ("Years in irruption, 2027–50", "Peak adults ha$^{-1}$, 2027–31", "Mean coral cover, 2027–50", "Domain removals (thousands)")):
        for offset, ssp, color in ((-.12, "ssp245", "#2a78d6"), (.12, "ssp370", "#eb6834")):
            g = means[means.ssp == ssp].groupby("cull")[var]
            m, lo, hi = [f().reindex(order) for f in (g.median, lambda: g.quantile(.1), lambda: g.quantile(.9))]
            ax.errorbar([i+offset for i in range(9)], m, yerr=[m-lo, hi-m], color=color, fmt="o", ms=3, linewidth=.8)
        ax.set_xticks(range(9), labels, rotation=55, ha="right")
        ax.set_title(title, loc="left")
        style(ax)
    fig.legend(handles=[Line2D([], [], color=c, marker="o", label=s) for c, s in (("#2a78d6", "SSP2-4.5"), ("#eb6834", "SSP3-7.0"))],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    save(fig, "fig27_metacommunity_projection_effects", bottom=.08)


def main():
    parser = argparse.ArgumentParser()
    for key in ("calibration", "projections", "adria", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    a = parser.parse_args()
    meta = json.loads(a.calibration.read_text())
    for mpath in (a.calibration, a.projections):
        m = json.loads(mpath.read_text())
        for name, sha in m["outputs"].items():
            checked(mpath.parent / name, sha)
    paths = {n: a.adria / "sandbox/data" / n for n in
             ("reef_cots.csv", "reef_manta.csv", "reef_photo_transect.csv")}
    paths.update({"trajectories.csv": a.adria / meta["run"].replace("\\", "/") / "trajectories.csv",
                  "cull/trajectories.csv": a.adria / meta["cull_run"].replace("\\", "/") / "trajectories.csv"})
    for n, p in paths.items():
        checked(p, meta["inputs"][n])
    obs, manta, photo = [pd.read_csv(paths[n], low_memory=False) for n in
                         ("reef_cots.csv", "reef_manta.csv", "reef_photo_transect.csv")]
    photo = photo[(photo.variable == "HARD CORAL") & (photo.purpose == "GROUP_LEVEL") & (photo.project_code == "LTMP")]
    tr = [pd.read_csv(paths[n]) for n in ("trajectories.csv", "cull/trajectories.csv")]
    reefs = {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef",
             "North Direction Reef": "North Direction Island", "Eyrie Reef": "Eyrie Reef"}
    # Shared report blue/orange; purple distinguishes the culling diagnostic.
    blue, purple, orange, grid = "#0d366b", "#5b2a86", "#eb6834", "#e1e0d9"
    plt.rcParams.update({"font.family": "Arial", "font.size": 9.3, "axes.titlesize": 10.3})
    fig, axes = plt.subplots(4, 2, figsize=(7.25, 8.0), sharex=True)
    for row, (reef, oname) in enumerate(reefs.items()):
        for col, variable in enumerate(("adults_ha", "coral_cover")):
            ax = axes[row, col]
            for t, color, ls in zip(tr, (blue, purple), ("-", "--")):
                pivot = t[t.reef_name == reef].pivot(index="year", columns="treatment", values=variable).loc[1992:2026]
                if col == 0:
                    pivot *= .015
                base = pivot[[n for n in pivot if "_sx" not in n and "_sh" not in n]]
                assert base.shape[1] == 5 and pivot.shape[1] == 65
                ax.fill_between(pivot.index, pivot.min(axis=1), pivot.max(axis=1), color=color, alpha=.10, linewidth=0)
                ax.fill_between(base.index, base.min(axis=1), base.max(axis=1), color=color, alpha=.22, linewidth=0)
                ax.plot(base.index, base.mean(axis=1), color=color, linestyle=ls, linewidth=1.2)
            if col == 0:
                o = obs[obs.reef_name == oname].groupby("year").cotsptow.mean()
                ax.scatter(o.index, o, color="#1f1e1d", s=10, zorder=5)
                ax.axhline(.22, color="#898781", linestyle=":", linewidth=.8)
                ax.set_ylim(bottom=0)
            else:
                m = manta[manta.reef_name == oname].sort_values("report_year")
                ax.errorbar(m.report_year, m["median"], yerr=[m["median"]-m.lower, m.upper-m["median"]],
                            fmt="o", ms=2.6, color="#1f1e1d", ecolor="#898781", elinewidth=.7, zorder=5)
                ph = photo[photo.domain_name == oname].groupby("report_year")["median"].mean()
                ax.scatter(ph.index, ph, color=orange, marker="D", s=12, zorder=6)
                ax.set_ylim(0, .8)
            ax.set_title(f"{chr(65 + row*2 + col)}  {reef.replace(' Reef', '')}", loc="left", pad=3)
            ax.set_xlim(1991.5, 2026.5)
            ax.set_xticks([1995, 2005, 2015, 2025])
            ax.grid(color=grid, linewidth=.6)
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(labelsize=9.3)
    fig.text(.28, .98, "COTS per manta tow", ha="center", weight="bold", fontsize=10.3)
    fig.text(.77, .98, "Hard coral cover (fraction)", ha="center", weight="bold", fontsize=10.3)
    for ax in axes[-1]:
        ax.set_xlabel("Year")
    handles = [Line2D([], [], color=blue, label="No recorded culling"),
               Line2D([], [], color=purple, linestyle="--", label="Recorded culling"),
               Patch(color="#777777", alpha=.22, label="Five-seed range"),
               Patch(color="#777777", alpha=.10, label="Parameter range"),
               Line2D([], [], color="#1f1e1d", marker="o", linestyle="none", markersize=3, label="LTMP manta"),
               Line2D([], [], color=orange, marker="D", linestyle="none", markersize=3, label="LTMP 9 m photo")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .006), ncol=3, frameon=False, fontsize=9.3)
    fig.subplots_adjust(left=.075, right=.985, top=.94, bottom=.13, hspace=.34, wspace=.23)
    for ext in ("png", "pdf"):
        fig.savefig(a.output / f"fig22_metacommunity_cots_coral.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    reflow_supplements(a, meta, tr, obs)
    print("Built combined calibration figure and five portrait supplementary figures (PNG/PDF).")


if __name__ == "__main__":
    main()

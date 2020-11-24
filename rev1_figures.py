from matplotlib.patches import Patch
from neuro_helper.abstract.map import HierarchyName, Space
from neuro_helper.dataframe import calc_percentage_change
from neuro_helper.generic import combine_topo_map
from neuro_helper.hcp.fmri.generic import task_order
from neuro_helper.map import WangTemplateMap
from neuro_helper.plot import *
from scipy.io import loadmat
from scipy.stats import stats
from statsmodels.stats.multitest import multipletests

import hcp_mf as mf
import hcp_lzc as lzc
import seaborn as sns
import matplotlib.pyplot as plt
import ptitprince as pt
import numpy as np
import pandas as pd
from config import *

font_scale = 1.1
sns.set(font_scale=font_scale, style="whitegrid")
LH_colors = paired_colors()
rest_task_colors = triple_colors()
task_colors = triple_colors()[1:]
tpt = cole_tpt


def map_regions():
    import cifti
    import os
    from neuro_helper.assets.manager import get_full_path, AssetCategory
    img, (lbl_axis, brain) = cifti.read(cole_tpt.file_full_path)
    regions = lbl_axis.label.item()
    nets = cole_tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER)
    img[img > 360] = 0

    lh_colors = paired_colors(True)
    lh_out = {0: ("???", (1.0, 1.0, 1.0, 0.0))}
    net_out = {0: ("???", (1.0, 1.0, 1.0, 0.0))}
    for index in range(1, 361):
        name, c = regions[index]
        net, lh = name.split("_")
        net_parts = net.split("-")
        net, rgn = ("".join(net_parts[:2]), net_parts[2]) if len(net_parts) == 3 else net_parts
        is_l = net in nets["L"]
        lh_out[index] = name, lh_colors[0 if is_l else 1]
        net_out[index] = name, cole_tpt.net_colors[cole_tpt.net_order.index(net)] + [1.0]

    for lbl, out in zip(["lh", "net"], [lh_out, net_out]):
        # noinspection PyTypeChecker
        cifti.write(f"figures/cole.{lbl}.dlabel.nii", img, (cifti.Label([out]), brain))
        os.system(f"wb_command -cifti-separate figures/cole.{lbl}.dlabel.nii COLUMN "
                  f"-label CORTEX_LEFT figures/cole.{lbl}.L.label.gii "
                  f"-label CORTEX_RIGHT figures/cole.{lbl}.R.label.gii")

        os.system(f"wb_command -label-to-border "
                  f"{get_full_path(AssetCategory.HCP1200, 'anat.midthickness.59k.L.surf.gii')} "
                  f"figures/cole.{lbl}.L.label.gii figures/cole.{lbl}.L.border")
        os.system(f"wb_command -label-to-border "
                  f"{get_full_path(AssetCategory.HCP1200, 'anat.midthickness.59k.R.surf.gii')} "
                  f"figures/cole.{lbl}.R.label.gii figures/cole.{lbl}.R.border")


def lzc_rest_net():
    df = lzc.gen_long_data(tpt) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    pt.RainCloud(data=df, x="network", y="metric", order=tpt.net_order, ax=ax, offset=0.1,
                 pointplot=True, palette=tpt.net_colors, scale="width")
    ax.set(xlabel="", ylabel=f"LZC")
    ax.set_xticklabels(tpt.net_labels(), rotation=90)
    print(savefig(fig, "rev1.lzc.rest.network", low=False))


def lzc_rest_lh():
    df = lzc.gen_long_data(tpt) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="metric", order=["L", "H"], ax=ax, offset=0.1,
                 pointplot=True, palette=LH_colors, orient="h")
    ax.set(ylabel="", xlabel=f"LZC", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels)
    print(savefig(fig, "rev1.lzc.rest.lh", low=False))


def lzc_rest_net_3():
    df = lzc.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    pt.RainCloud(data=df, x="network", y="metric", order=tpt.net_order, ax=ax, offset=0.1,
                 pointplot=True, palette=tpt.net_colors, scale="width")
    ax.set(xlabel="", ylabel=f"LZC")
    ax.set_xticklabels(tpt.net_labels(True), rotation=90)
    print(savefig(fig, "rev1.lzc.rest.network.3", low=False))


def lzc_rest_lh_3():
    df = lzc.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(cole_tpt_32.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="metric", order=["L", "H"], ax=ax, offset=0.1,
                 pointplot=True, palette=LH_colors, orient="h")
    ax.set(ylabel="", xlabel=f"LZC", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels)
    print(savefig(fig, "rev1.lzc.rest.lh.3", low=False))


def lzc_task_net():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, False), task_colors)]
    df = lzc.gen_long_data(tpt) \
        .and_filter(NOTtask="REST") \
        .groupby(["task", "region", "network"]).mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    pt.RainCloud(data=df, x="network", y="metric", hue="task", order=tpt.net_order, alpha=.65, dodge=True,
                 move=.1, bw=.2, width_viol=.7, hue_order=task_order(tpt.space, False), ax=ax, offset=0.1,
                 pointplot=True, palette=task_colors, scale="width")
    ax.set(xlabel="", ylabel=f"LZC")
    ax.set_xticklabels(tpt.net_labels(True), rotation=90)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.118, -0.035, 0.79, 1))
    print(savefig(fig, "rev1.lzc.task.network", low=False, extra_artists=(lgn,)))


def lzc_task_lh():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, False), task_colors)]
    df = lzc.gen_long_data(tpt) \
        .and_filter(NOTtask="REST") \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="metric", hue="task", order=["L", "H"], palette=task_colors, move=.2,
                 hue_order=task_order(tpt.space, False), ax=ax, offset=0.1, alpha=.65, dodge=True, bw=.2, width_viol=.7,
                 pointplot=True, orient="h")
    ax.set(ylabel="", xlabel=f"LZC", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.108, 0.05, 0.81, 1))
    print(savefig(fig, "rev1.lzc.task.lh", low=False, extra_artists=(lgn,)))


def lzc_pchange_net():
    shared_subjs = lzc.find_shared_subjects(tpt, task_order(tpt.space, True))
    df = lzc.gen_long_data(tpt).and_filter(subject=shared_subjs) \
        .groupby(["task", "subject", "region", "network"]).mean().reset_index() \
        .groupby(["subject", "network", "region"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_3", 1) \
        .groupby(["task", "region", "network"]).mean().reset_index()

    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, False), task_colors)]

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    sns.barplot(data=df, x="network", y="pchange", hue="task", order=tpt.net_order,
                hue_order=task_order(tpt.space, False), ax=ax, palette=task_colors)
    # pt.RainCloud(data=df, x="network", y="pchange", hue="task", order=tpt.net_order, alpha=.65, dodge=True,
    #              move=.1, bw=.2, width_viol=.7, hue_order=task_order(tpt.space, False), ax=ax, offset=0.1,
    #              pointplot=True, palette=task_colors, scale="width")
    ax.set(xlabel="", ylabel="", ylim=[-25, 25])
    ax.set_xticklabels(tpt.net_labels(True), rotation=90)
    txt0 = ax.text(-0.09, 0.5, "LZC change from REST (%)", rotation=90, transform=ax.transAxes, ha='center',
                   va='center')
    txt1 = ax.text(-0.06, 0.75, "Rest > Task", rotation=90, transform=ax.transAxes, ha='center', va='center')
    txt2 = ax.text(-0.06, 0.25, "Rest < Task", rotation=90, transform=ax.transAxes, ha='center', va='center')
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.118, -0.035, 0.79, 1))
    print(savefig(fig, "rev1.lzc.pchange.network.bar", low=False, extra_artists=(lgn, txt0, txt1, txt2)))


def lzc_pchange_lh():
    shared_subjs = lzc.find_shared_subjects(tpt, task_order(tpt.space, True))
    df = lzc.gen_long_data(tpt).and_filter(subject=shared_subjs) \
        .groupby(["task", "subject", "region", "network"]).mean().reset_index() \
        .groupby(["subject", "network", "region"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_3", 1) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, False), task_colors)]

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="pchange", hue="task", order=["L", "H"], palette=task_colors, move=.2,
                 hue_order=task_order(tpt.space, False), ax=ax, offset=0.1, alpha=.65, dodge=True, bw=.2, width_viol=.7,
                 pointplot=True, orient="h")
    ax.set(ylabel="", xlabel="", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels, xlim=[-20, 20])
    txt0 = ax.text(0.5, -0.3, "LZC change from REST (%)", transform=ax.transAxes, ha='center', va='center')
    txt1 = ax.text(0.75, -0.2, "Rest > Task", transform=ax.transAxes, ha='center', va='center')
    txt2 = ax.text(0.25, -0.2, "Rest < Task", transform=ax.transAxes, ha='center', va='center')
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.108, 0.05, 0.81, 1))
    print(savefig(fig, "rev1.lzc.pchange.lh", low=False, extra_artists=(lgn, txt0, txt1, txt2)))


def lzc_wang():
    wang_tpt = WangTemplateMap(Space.K59)
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(wang_tpt.space, True), rest_task_colors)]
    legend_handles_pchange = [Patch(facecolor=color, edgecolor=color, label=task)
                              for task, color in zip(task_order(wang_tpt.space, False), task_colors)]

    df = lzc.gen_long_data(wang_tpt, show_warning=True) \
        .groupby(["task", "subject", "region"]).mean().reset_index() \
        .groupby(["subject", "region"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_2", 1)

    fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex='col')
    ax = axs[0]
    sns.barplot(data=df, x="region", y="metric", hue="task",
                hue_order=task_order(wang_tpt.space, True), order=wang_tpt.net_order, ax=ax, palette=rest_task_colors)
    ax.set(xlabel="", ylabel="Mean LZC \u00B1 %95 CI", ylim=[0.6, None])
    ax.set_xticklabels(wang_tpt.net_labels(True), rotation=90)
    ax.legend(handles=legend_handles)

    ax = axs[1]
    sns.barplot(data=df, x="region", y="pchange", hue="task",
                hue_order=task_order(wang_tpt.space, True), order=wang_tpt.net_order, ax=ax, palette=task_colors)
    ax.set(xlabel="", ylabel="", ylim=[-15, 15])
    ax.set_xticklabels(wang_tpt.net_labels(True), rotation=90)
    ax.legend(handles=legend_handles_pchange)
    txt0 = ax.text(-0.09, 0.5, "LZC change from REST (%) \u00B1 95% CI", rotation=90,
                   transform=ax.transAxes, ha='center', va='center')
    txt1 = ax.text(-0.06, 0.75, "Rest > Task", rotation=90, transform=ax.transAxes, ha='center', va='center')
    txt2 = ax.text(-0.06, 0.25, "Rest < Task", rotation=90, transform=ax.transAxes, ha='center', va='center')
    fig.subplots_adjust(hspace=0.05)
    print(savefig(fig, "rev1.lzc.visual.net", low=False, extra_artists=(txt0, txt1, txt2)))


def mf_task_net():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, True), rest_task_colors)]
    df = mf.gen_long_data(tpt).groupby(["task", "region", "network"]).mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    pt.RainCloud(data=df, x="network", y="metric", hue="task", order=tpt.net_order, alpha=.65, dodge=True,
                 move=.1, bw=.2, width_viol=.7, hue_order=task_order(tpt.space, True), ax=ax, offset=0.1,
                 pointplot=True, palette=rest_task_colors, scale="width")
    ax.set(xlabel="", ylabel=f"MF (Hz)")
    ax.set_xticklabels(tpt.net_labels(True), rotation=90)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.118, -0.035, 0.79, 1))
    print(savefig(fig, "rev1.mf.task.network", low=False, extra_artists=(lgn,)))


def mf_task_lh():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=task)
                      for task, color in zip(task_order(tpt.space, True), rest_task_colors)]
    df = mf.gen_long_data(tpt) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="metric", hue="task", order=["L", "H"], palette=rest_task_colors, move=.2,
                 hue_order=task_order(tpt.space, True), ax=ax, offset=0.1, alpha=.65, dodge=True, bw=.2, width_viol=.7,
                 pointplot=True, orient="h")
    ax.set(ylabel="", xlabel=f"MF (Hz)", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=3, mode="expand",
                     bbox_to_anchor=(0.108, 0.05, 0.81, 1))
    print(savefig(fig, "rev1.mf.task.lh", low=False, extra_artists=(lgn,)))


def mf_definition():
    from neuro_helper.hcp.fmri.storage import FMRILocalStorage, load_raw_file
    from neuro_helper.statistics import fir_filter, welch_psd
    from neuro_helper.measurement import calc_mf
    dt, fs = load_raw_file(FMRILocalStorage(FMRI_RAW_DATA_ROOT_DIR, "REST", "3", tpt.space).get_all()[0], tpt.space)
    dt, freq_l, freq_h = fir_filter(dt, fs, max_freq_low=0.01, pass_type="hp")
    data = dt[1].reshape(1, -1)
    freq, psd = welch_psd(data, freq_l, fs)
    filt = freq >= freq_l
    freq = freq[filt]
    psd = psd[:, filt]
    mf_val = calc_mf(freq, psd).item()
    psd = (psd / psd.max()).squeeze()

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    sns.lineplot(freq, psd, ax=ax)
    plt.fill_between(freq[freq <= mf_val], psd[freq <= mf_val], color="red")
    plt.fill_between(freq[freq > mf_val], psd[freq > mf_val], color="blue")
    ticks = [0.01, 0.1, mf_val, 0.3, 0.4, 0.5]
    ax.set_xticklabels(["0.01", "0.1", f"{mf_val:.3f}", "0.3", "0.4", "0.5"])
    ax.set_xticks(ticks)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"])
    ax.set_ylabel("Normalized Power")
    ax.set_xlabel("Frequency (Hz)")
    ax.grid(False)
    ax.axvline(mf_val, color='#000000', clip_on=False, linewidth=1, ymin=0, ymax=0.8)
    ax.arrow(mf_val, 0, 0, 0.8, head_width=0.02, head_length=0.04, fc='k', ec='k', color='#000000', linewidth=1)
    sns.set(font_scale=0.6, style="whitegrid")
    ax.text(0.5, 0.85, "Median Frequency: where AUC of red area = AUC of blue area",
            transform=ax.transAxes, ha='center', va='bottom')
    sns.set(font_scale=font_scale, style="whitegrid")
    print(savefig(fig, "rev1.mf.definition", low=False))


def mf_rest_net_3():
    df = mf.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    pt.RainCloud(data=df, x="network", y="metric", order=tpt.net_order, ax=ax, offset=0.1,
                 pointplot=True, palette=tpt.net_colors, scale="width")
    ax.set(xlabel="", ylabel=f"MF (Hz)")
    ax.set_xticklabels(tpt.net_labels(True), rotation=90)
    print(savefig(fig, "rev1.mf.rest.network.3", low=False))


def mf_rest_lh_3():
    df = mf.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(cole_tpt_32.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    fig, ax = plt.subplots(1, 1, figsize=(5, 2.5))
    pt.RainCloud(data=df, x="net_meta", y="metric", order=["L", "H"], ax=ax, offset=0.1,
                 pointplot=True, palette=LH_colors, orient="h")
    ax.set(ylabel="", xlabel=f"MF (Hz)", yticklabels=HierarchyName.LOWER_HIGHER_ORDER.labels)
    print(savefig(fig, "rev1.mf.rest.lh.3", low=False))


def relation():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=net)
                      for color, net in zip(tpt.net_colors, tpt.net_labels(False))]

    df = pd.merge(
        lzc.gen_long_data(tpt).rename(columns={"metric": "lzc"})
            .groupby(["task", "region", "network", "subject"]).mean().reset_index(),
        mf.gen_long_data(tpt).rename(columns={"metric": "mf"})
            .groupby(["task", "region", "network", "subject"]).mean().reset_index(),
        on=["task", "region", "network", "subject"]
    )

    # correlation
    for func, func_lbl in zip([stats.pearsonr, stats.spearmanr], ["pearson", "spearman"]):
        maps = []
        for task in task_order(tpt.space):
            corr = df.and_filter(task=task).sort_values("subject").reset_index(drop=True).groupby("region").apply(
                lambda x: pd.DataFrame(
                    np.asarray(func(x.lzc, x.mf)).reshape(1, -1), columns=["a", "p"])) \
                .reset_index().drop("level_1", 1)
            rejected, _, _, _ = multipletests(corr.p, method="fdr_bh")
            corr["a_sig"] = corr.a.copy()
            corr.loc[~rejected, "a_sig"] = 0
            maps.append(corr[["region", "a"]].build_single_topo_map(tpt))
            maps.append(corr[["region", "a_sig"]].build_single_topo_map(tpt))

        topo, brain, series = combine_topo_map(maps)
        savemap(f"relation.{func_lbl}", topo, brain, series)

    # scatter
    fig, axs = plt.subplots(2, 3, figsize=(14, 10), sharex='row')
    for ti, task in enumerate(task_order(tpt.space)):
        ax = axs[0, ti]
        dft = df.and_filter(task=task).groupby(["region", "network"]).mean().reset_index()
        sns.scatterplot(data=dft, x="mf", y="lzc", ax=ax, hue="network", hue_order=tpt.net_order,
                        palette=tpt.net_colors)
        ax.set(title=task, xlabel=f"MF {task} (Hz)", ylabel=f"LZC {task}", ylim=[0.6, 1.2])
        ax.get_legend().remove()

        ax = axs[1, ti]
        dft = df.and_filter(task=task).groupby(["subject"]).mean().reset_index()
        sns.scatterplot(data=dft, x="mf", y="lzc", ax=ax, color='black')
        ax.set(title=task, xlabel=f"MF {task} (Hz)", ylabel=f"LZC {task}", ylim=[0.88, 1.12])

    fig.subplots_adjust(wspace=0.3, hspace=0.4)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=6, mode="expand", bbox_to_anchor=(0.118, 0.01, 0.79, 1))
    sns.set(font_scale=1.2, style="whitegrid")
    txt1 = fig.text(0.5, 0.92, "LZC-MF Relationship of Regions", ha='center', va='center', clip_on=False)
    txt2 = fig.text(0.5, 0.47, "LZC-MF Relationship of Subjects", ha='center', va='center', clip_on=False)
    sns.set(font_scale=font_scale, style="whitegrid")
    print(savefig(fig, "rev1.relation.scatter", low=False, extra_artists=(lgn, txt1, txt2)))


def relation_simulation():
    from neuro_helper.statistics import welch_psd
    mat = loadmat("simulation/signals.mat")
    signals = mat["signals"]
    mf_values = mat["MF"]
    lzc_values = mat["LZC"]
    types = mat["types"].astype(bool)
    names = ["Pink", "White", "Sine"]

    psds = []
    for t in signals:
        freq, psd = welch_psd(t, 0.01, 1)
        psds.append(psd)

    fig, axs = plt.subplots(7, 1, figsize=(2.5, 14), sharex='all', sharey='all')
    for index in range(7):
        ax = axs[index]
        title = " & ".join(names[x] for x in range(len(names)) if types[index, x])
        sns.scatterplot(x=mf_values[index], y=lzc_values[index], ax=ax, color="darkorange")
        ax.set(title=title, xlabel="MF", ylabel="LZC")
    txt = fig.text(0.5, 0.92, "LZC-MF Relationship\nof Simulated Signals", ha='center', va='center', clip_on=False)
    print(savefig(fig, "rev1.relation.simulation", low=False, extra_artists=(txt,)))


def dynamic_range():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=net)
                      for color, net in zip(tpt.net_colors, tpt.net_labels(False))]

    df_lzc = lzc.gen_long_data(tpt) \
        .groupby(["task", "region", "network", "subject"]).mean().reset_index() \
        .groupby(["region", "network", "subject"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_3", 1).groupby(["task", "region", "network"]).mean().reset_index()
    df_mf_rest = mf.gen_long_data(tpt).and_filter(task="REST").groupby(["region", "network"]).mean().reset_index()

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for ti, task in enumerate(task_order(tpt.space, False)):
        dft = pd.merge(df_lzc.and_filter(task=task), df_mf_rest, on=["region", "network"])
        ax = axs[ti]
        sns.scatterplot(data=dft, x="metric", y="pchange", hue="network",
                        hue_order=tpt.net_order, ax=ax, palette=tpt.net_colors)
        ax.set(title=f"MF During REST ~ Change in LZC from REST to {task}", ylabel=f"LZC Change from REST (%)",
               xlabel=f"MF REST (Hz)")
        ax.get_legend().remove()

    fig.subplots_adjust(wspace=0.3)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=4, mode="expand", bbox_to_anchor=(0.118, 0.12, 0.79, 1))
    print(savefig(fig, "rev1.drlzc.scatter", low=False, extra_artists=(lgn,)))


def dynamic_range_split():
    df_lzc = lzc.gen_long_data(tpt) \
        .groupby(["task", "region", "network", "subject"]).mean().reset_index() \
        .groupby(["region", "network", "subject"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_3", 1).groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))
    df_mf_rest = mf.gen_long_data(tpt).and_filter(task="REST").groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER)) \
        .add_median_lh("metric")

    fig, axs = plt.subplots(1, 2, figsize=(12, 2.5))
    for ti, task in enumerate(task_order(tpt.space, False)):
        dft = pd.merge(df_lzc.and_filter(task=task), df_mf_rest, on=["region", "network"])
        ax = axs[ti]
        sns.boxplot(data=dft, y="metric_split", x="pchange", ax=ax, order=["L", "H"], showfliers=False, orient="h")
        ax.set(xlabel="", ylabel="", yticklabels=["Low MF REST", "High MF REST"],
               xticks=[], title=f"Change in LZC {task} from REST")
    fig.subplots_adjust(wspace=0.35, hspace=0.3)
    print(savefig(fig, "rev1.drlzc.split", low=False))


def t1t2():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=net)
                      for color, net in zip(tpt.net_colors, tpt.net_labels(False))]

    df = pd.merge(
        lzc.gen_long_data(tpt).rename(columns={"metric": "lzc"})
            .groupby(["task", "region", "network"]).mean().reset_index(),
        mf.gen_long_data(tpt).rename(columns={"metric": "mf"})
            .groupby(["task", "region", "network"]).mean().reset_index(),
        on=["task", "region", "network"]
    ).add_topo(t1t2_topo)

    t1t2_text = r'$\frac{T_1w}{T_2w}$'
    fig, axs = plt.subplots(2, 3, figsize=(14, 10))
    for ti, task in enumerate(task_order(tpt.space)):
        dft = df.and_filter(task=task)
        ax = sns.scatterplot(data=dft, x="t1t2", y="lzc", ax=axs[0, ti], hue="network", hue_order=tpt.net_order,
                             palette=tpt.net_colors)
        ax.set(title=task, xlabel=t1t2_text, ylabel=f"LZC {task}", ylim=[0.65, 1.15])
        ax.get_legend().remove()

        ax = sns.scatterplot(data=dft, x="t1t2", y="mf", ax=axs[1, ti], hue="network", hue_order=tpt.net_order,
                             palette=tpt.net_colors)
        ax.set(title=task, xlabel=t1t2_text, ylabel=f"MF {task} (Hz)", ylim=[0.035, 0.235])
        ax.set_yticklabels([f'{x:.2f}' for x in ax.get_yticks()])
        ax.get_legend().remove()

    fig.subplots_adjust(wspace=0.3, hspace=0.5)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=6, mode="expand", bbox_to_anchor=(0.118, 0.01, 0.79, 1))
    sns.set(font_scale=1.2, style="whitegrid")
    txt1 = fig.text(0.5, 0.92, f"LZC-{t1t2_text} Relationship of Regions", ha='center', va='center', clip_on=False)
    txt2 = fig.text(0.5, 0.455, f"MF-{t1t2_text} Relationship of Regions", ha='center', va='center', clip_on=False)
    sns.set(font_scale=font_scale, style="whitegrid")
    print(savefig(fig, "rev1.t1t2.scatter", low=False, extra_artists=(lgn, txt1, txt2)))


def relation_3():
    legend_handles = [Patch(facecolor=color, edgecolor=color, label=net)
                      for color, net in zip(cole_tpt_32.net_colors, cole_tpt_32.net_labels(False))]

    df = pd.merge(
        lzc.gen_long_data(cole_tpt_32).rename(columns={"metric": "lzc"})
            .groupby(["region", "network", "subject"]).mean().reset_index(),
        mf.gen_long_data(cole_tpt_32).rename(columns={"metric": "mf"})
            .groupby(["region", "network", "subject"]).mean().reset_index(),
        on=["region", "network", "subject"]
    )

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    dft = df.groupby(["region", "network"]).mean().reset_index()
    ax = sns.scatterplot(data=dft, x="mf", y="lzc", ax=axs[0, 0], hue="network", hue_order=cole_tpt_32.net_order,
                         palette=cole_tpt_32.net_colors)
    ax.set(title="LZC-MF Relationship of Regions", xlabel=f"MF REST", ylabel=f"LZC REST")
    ax.get_legend().remove()

    dft = df.groupby(["subject"]).mean().reset_index()
    ax = sns.scatterplot(data=dft, x="mf", y="lzc", ax=axs[0, 1], color="#000000")
    ax.set(title="LZC-MF Relationship of Subjects", xlabel="MF REST", ylabel="LZC REST")
    ax.set_yticklabels([f'{x:.2f}' for x in ax.get_yticks()])

    t1t2_text = r'$\frac{T_1w}{T_2w}$'
    dft = df.groupby(["region", "network"]).mean().reset_index().add_topo(t1t2_topo)
    ax = sns.scatterplot(data=dft, x="t1t2", y="lzc", ax=axs[1, 0], hue="network", hue_order=cole_tpt_32.net_order,
                         palette=cole_tpt_32.net_colors)
    ax.set(title=f"LZC-{t1t2_text} Relationship of Regions", xlabel=t1t2_text, ylabel=f"LZC REST")
    ax.get_legend().remove()

    ax = sns.scatterplot(data=dft, x="t1t2", y="mf", ax=axs[1, 1], hue="network", hue_order=cole_tpt_32.net_order,
                         palette=cole_tpt_32.net_colors)
    ax.set(title=f"MF-{t1t2_text} Relationship of Regions", xlabel=t1t2_text, ylabel="MF REST (Hz)")
    ax.set_yticklabels([f'{x:.2f}' for x in ax.get_yticks()])
    ax.get_legend().remove()

    fig.subplots_adjust(wspace=0.3, hspace=0.4)
    lgn = fig.legend(handles=legend_handles, loc=2, ncol=4, mode="expand", bbox_to_anchor=(0.118, 0.02, 0.79, 1))
    print(savefig(fig, "rev1.relation.scatter.3", low=False, extra_artists=(lgn, )))


if __name__ == "__main__":
    pass

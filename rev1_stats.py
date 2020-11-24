from neuro_helper.abstract.map import HierarchyName
from neuro_helper.hcp.fmri.generic import task_order
import hcp_mf as mf
from statsmodels.formula.api import ols
import statsmodels.api as sm
import hcp_lzc as lzc
from neuro_helper.dataframe import calc_percentage_change
from scipy import stats
from neuro_helper.statistics import cohend, anova_table
import numpy as np
import pandas as pd
import pingouin as pg
from config import *

tpt = cole_tpt


def print_ttest(label, d1, d2):
    t, p = stats.ttest_ind(d1, d2)
    d = cohend(d1, d2)
    print(f"{label}: T = {t:.2f}, Cohen's D = {d:.2f}, P = {p:.3f}")
    return t, d, p


def lzc_rest_net():
    df = lzc.gen_long_data(tpt) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    model = ols('metric ~ network', data=df)
    result = model.fit()
    result.summary()
    robust = None if result.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(result, typ=2, robust=robust))
    print(aov_table.to_string())
    # mc = MultiComparison(df.metric, df.network)
    # mc_results = mc.tukeyhsd()
    # print(mc_results)


def lzc_rest_lh():
    df = lzc.gen_long_data(tpt) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .and_filter(task="REST") \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    print_ttest(f"LZC REST",
                df.and_filter(net_meta="L").metric,
                df.and_filter(net_meta="H").metric)


def lzc_rest_net_3():
    df = lzc.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    model = ols('metric ~ network', data=df)
    result = model.fit()
    result.summary()
    robust = None if result.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(result, typ=2, robust=robust))
    print(aov_table.to_string())


def lzc_rest_lh_3():
    df = lzc.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST", ) \
        .groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(cole_tpt_32.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    print_ttest(f"LZC REST",
                df.and_filter(net_meta="L").metric,
                df.and_filter(net_meta="H").metric)


def lzc_task_net():
    df = lzc.gen_long_data(tpt) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .and_filter(NOTtask="Rest")

    model = ols('metric ~ C(task) + C(network) + C(task):C(network)', data=df).fit()
    model.summary()
    robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
    print(aov_table.to_string())

    for task in task_order(tpt.space, False):
        dft = df.and_filter(task=task)
        print(f"\n\n####### {task} #######")
        model = ols('metric ~ C(network)', data=dft).fit()
        model.summary()
        robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
        aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
        print(aov_table.to_string())


def lzc_task_lh():
    df = lzc.gen_long_data(tpt) \
        .and_filter(NOTtask="REST") \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    model = ols('metric ~ C(task) + C(net_meta) + C(task):C(net_meta)', data=df).fit()
    model.summary()
    robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
    print(aov_table.to_string())

    df["comb"] = pd.Series(df.task + "/" + df.net_meta, df.index, str)
    result = df.pairwise_tukey(dv="metric", between="comb", effsize="cohen")
    left = result.A.str.split("/", expand=True)
    right = result.B.str.split("/", expand=True)
    for task in task_order(tpt.space, False):
        print(result[(left[0] == task) & (right[0] == task)].to_string())


def lzc_pchange():
    shared_subjs = lzc.find_shared_subjects(tpt, task_order(tpt.space, True))
    df = lzc.gen_long_data(tpt).and_filter(subject=shared_subjs) \
        .groupby(["task", "subject", "region", "network"]).mean().reset_index() \
        .groupby(["subject", "network", "region"]).apply(calc_percentage_change, from_="REST") \
        .reset_index().drop("level_3", 1) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    model = ols('pchange ~ C(task) + C(network) + C(task):C(network)', data=df).fit()
    model.summary()
    robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
    print(aov_table.to_string())

    for task in task_order(tpt.space, False):
        dft = df.and_filter(task=task)
        print(f"\n\n####### {task} #######")
        model = ols('pchange ~ + C(network)', data=dft).fit()
        model.summary()
        robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
        aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
        print(aov_table.to_string())

    print("######################## LH ########################")
    model = ols('pchange ~ C(task) + C(net_meta) + C(task):C(net_meta)', data=df).fit()
    model.summary()
    robust = None if model.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(model, typ=2, robust=robust))
    print(aov_table.to_string())

    df["comb"] = pd.Series(df.task + "/" + df.net_meta, df.index, str)
    result = df.pairwise_tukey(dv="pchange", between="comb", effsize="cohen")
    left = result.A.str.split("/", expand=True)
    right = result.B.str.split("/", expand=True)
    for task in task_order(tpt.space, False):
        print(result[(left[0] == task) & (right[0] == task)].to_string())


def mf_rest_net_3():
    df = mf.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST") \
        .groupby(["region", "network"]).mean().reset_index()

    model = ols('metric ~ network', data=df)
    result = model.fit()
    result.summary()
    robust = None if result.diagn["omnipv"] > 0.05 else "hc3"
    aov_table = anova_table(sm.stats.anova_lm(result, typ=2, robust=robust))
    print(aov_table.to_string())


def mf_rest_lh_3():
    df = mf.gen_long_data(cole_tpt_32) \
        .and_filter(task="REST", ) \
        .groupby(["region", "network"]).mean().reset_index() \
        .add_net_meta(cole_tpt_32.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    print_ttest(f"MF REST",
                df.and_filter(net_meta="L").metric,
                df.and_filter(net_meta="H").metric)


def mf_lh():
    df = mf.gen_long_data(tpt) \
        .groupby(["task", "region", "network"]).mean().reset_index() \
        .add_net_meta(tpt.net_hierarchy(HierarchyName.LOWER_HIGHER_ORDER))

    for task in task_order(tpt.space):
        print_ttest(f"LZC {task}",
                    df.and_filter(task=task, net_meta="L").metric,
                    df.and_filter(task=task, net_meta="H").metric)


def mf_net():
    df = mf.gen_long_data(tpt) \
        .groupby(["task", "region", "network"]).mean().reset_index()

    for task in task_order(tpt.space):
        print(task)
        model = ols('metric ~ network', data=df.and_filter(task=task))
        result = model.fit()
        result.summary()
        robust = None if result.diagn["omnipv"] > 0.05 else "hc3"
        aov_table = anova_table(sm.stats.anova_lm(result, typ=2, robust=robust))
        print(aov_table.to_string())


def relation_range():
    def calc_range(x):
        med = x.metric.median()
        return pd.Series(index=["L", "H"], data=[
            np.ptp(x[x.metric < med].pchange.values),
            np.ptp(x[x.metric >= med].pchange.values)])

    df_lzc = lzc.gen_long_data(tpt).sort_values("subject") \
        .groupby(["task", "region", "subject"]).mean().reset_index() \
        .groupby(["region", "subject"]).apply(calc_percentage_change, from_="REST").reset_index().drop("level_2", 1)

    df_mf_rest = mf.gen_long_data(tpt) \
        .and_filter(task="REST").sort_values("subject") \
        .groupby(["region", "subject"]).mean().reset_index()

    for task in task_order(tpt.space, False):
        dft = pd.merge(
            df_mf_rest, df_lzc.and_filter(task=task),
            on=["region", "subject"]
        ).groupby("subject").apply(calc_range).reset_index()

        t, p = stats.ttest_rel(dft.L.values, dft.H.values)
        d = cohend(dft.L.values, dft.H.values)
        print(f"{task}: T = {t:.2f}, Cohen's D = {d:.2f}, P = {p:.3f}")

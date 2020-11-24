from neuro_helper.dataframe import add_median_lh, long_column_to_wide
from neuro_helper.hcp.fmri.generic import task_order

from config import cole_tpt
import pandas as pd
import hcp_lzc as lzc
import hcp_mf as mf
import pyarrow.feather as feather


tpt = cole_tpt


def med_model_1():
    # X = LZC REST, Y = LZC Task, M = MF REST
    tasks = task_order(tpt.space)
    df_lzc = lzc.gen_long_data(tpt).sort_values("subject") \
        .groupby(["task", "region", "subject"]).mean().reset_index() \
        .rename(columns={"metric": "lzc"}) \
        .long_column_to_wide("lzc", "task") \
        .rename(columns={x: x.lower() for x in tasks})

    df_mf_rest = mf.gen_long_data(tpt).sort_values("subject").and_filter(task="REST") \
        .groupby(["region", "subject"]).mean().reset_index().rename(columns={"metric": "mf_rest"})
    df = pd.merge(df_lzc, df_mf_rest, on=["region", "subject"]).groupby("subject").apply(add_median_lh, "mf_rest", (0, 1))

    # WARNING: THIS IS MED-SPLIT FOR EACH SUBJECT. IT IS DIFFERENT FROM A GLOBAL MED-SPLIT
    feather.write_feather(df, 'r/med_model_1.feather')


def med_model_2():
    # X = LZC REST, Y = LZC Task, M = MF REST
    tasks = task_order(tpt.space)

    df_lzc = lzc.gen_long_data(tpt).sort_values("subject") \
        .groupby(["task", "region", "subject"]).mean().reset_index() \
        .rename(columns={"metric": "lzc"}) \
        .long_column_to_wide("lzc", "task") \
        .rename(columns={x: x.lower() for x in tasks})

    df_mf_rest = mf.gen_long_data(tpt).sort_values("subject").and_filter(task="REST") \
        .groupby(["region", "subject"]).mean().reset_index().rename(columns={"metric": "mf_rest"})
    df_mf_rest = pd.merge(
        df_mf_rest,
        add_median_lh(df_mf_rest.groupby("region").mean().reset_index(), "mf_rest", (0, 1)).drop("mf_rest", 1),
        on=("region", )
    )
    df = pd.merge(df_lzc, df_mf_rest, on=["region", "subject"])
    feather.write_feather(df, 'r/med_model_2.feather')


def med_model_3():
    # X = LZC REST, Y = LZC Task, M = MF REST
    tasks = task_order(tpt.space)

    df_lzc = lzc.gen_long_data(tpt).sort_values("subject") \
        .groupby(["task", "region", "subject"]).mean().reset_index() \
        .rename(columns={"metric": "lzc"}) \
        .long_column_to_wide("lzc", "task") \
        .rename(columns={x: x.lower() for x in tasks}) \
        .groupby("region").mean().reset_index()

    df_mf_rest = mf.gen_long_data(tpt).sort_values("subject").and_filter(task="REST") \
        .groupby(["region"]).mean().reset_index() \
        .rename(columns={"metric": "mf_rest"})
    df = pd.merge(df_lzc, df_mf_rest, on=["region"]).add_median_lh("mf_rest", (0, 1))
    feather.write_feather(df, 'r/med_model_3.feather')

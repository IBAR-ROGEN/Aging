"""Global sidebar filters applied across all dashboard tabs."""

from __future__ import annotations

import math
from dataclasses import dataclass

import pandas as pd
import streamlit as st

from rogen_aging.eda_dashboard.schema import resolve_column


@dataclass(frozen=True, slots=True)
class GlobalFilters:
    age_min: float
    age_max: float
    sex_values: list[str]
    disease_values: list[str]
    col_age: str
    col_sex: str
    col_disease: str


def render_global_sidebar(df: pd.DataFrame) -> GlobalFilters:
    """Render persistent cohort filters; returns resolved column names and selections."""
    st.sidebar.header("Cohort filters")
    st.sidebar.caption("Selections apply to every tab.")

    col_age = resolve_column(df, ("Chronological_Age", "chronological_age", "Age", "age"), label="chronological age")
    col_sex = resolve_column(df, ("Sex", "sex", "SEX", "gender", "Gender"), label="sex")
    col_disease = resolve_column(
        df,
        ("Disease_Status", "disease_status", "Disease", "disease", "Case_Control", "case_control"),
        label="disease status",
    )

    a_min = float(pd.to_numeric(df[col_age], errors="coerce").min())
    a_max = float(pd.to_numeric(df[col_age], errors="coerce").max())
    if not np_finite_pair(a_min, a_max):
        a_min, a_max = 0.0, 100.0

    lo, hi = st.sidebar.slider(
        "Age range (chronological)",
        min_value=a_min,
        max_value=a_max,
        value=(a_min, a_max),
        help="Restrict all plots to participants within this age window.",
    )

    sex_opts = sorted({str(x) for x in df[col_sex].dropna().unique().tolist()})
    dis_opts = sorted({str(x) for x in df[col_disease].dropna().unique().tolist()})

    sex_pick = st.sidebar.multiselect("Sex", options=sex_opts, default=sex_opts)
    dis_pick = st.sidebar.multiselect("Disease status", options=dis_opts, default=dis_opts)

    return GlobalFilters(
        age_min=lo,
        age_max=hi,
        sex_values=sex_pick,
        disease_values=dis_pick,
        col_age=col_age,
        col_sex=col_sex,
        col_disease=col_disease,
    )


def np_finite_pair(a: float, b: float) -> bool:
    return math.isfinite(a) and math.isfinite(b)


def apply_global_filters(df: pd.DataFrame, filt: GlobalFilters) -> pd.DataFrame:
    """Return a filtered view of the cohort (copy)."""
    out = df.copy()
    age_num = pd.to_numeric(out[filt.col_age], errors="coerce")
    mask = (age_num >= filt.age_min) & (age_num <= filt.age_max)
    if filt.sex_values:
        mask &= out[filt.col_sex].astype(str).isin(filt.sex_values)
    if filt.disease_values:
        mask &= out[filt.col_disease].astype(str).isin(filt.disease_values)
    return out.loc[mask].reset_index(drop=True)

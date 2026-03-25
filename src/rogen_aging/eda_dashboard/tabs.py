"""Modular Streamlit tab renderers for the multi-omics aging EDA dashboard."""

from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import streamlit as st
from scipy.stats import kruskal, pearsonr

from rogen_aging.eda_dashboard.schema import clinical_numeric_columns_for_heatmap, list_snp_columns
from rogen_aging.eda_dashboard.sidebar import GlobalFilters


def render_tab_clinical_overview(df: pd.DataFrame, filt: GlobalFilters) -> None:
    """Tab 1 — clinical distributions and inter-correlation structure."""
    st.subheader("Clinical & phenotypic overview")
    st.markdown(
        """
This tab summarizes **demographic and phenotypic structure** of the filtered cohort. The stacked
age histogram highlights how disease labels partition chronological age, while the correlation
heatmap surfaces redundant or orthogonal clinical signals prior to multi-omic modeling.
        """.strip()
    )

    plot_df = df.dropna(subset=[filt.col_age, filt.col_disease]).copy()
    plot_df[filt.col_age] = pd.to_numeric(plot_df[filt.col_age], errors="coerce")
    plot_df = plot_df.dropna(subset=[filt.col_age])

    if plot_df.empty:
        st.warning("No rows left after filtering for age and disease status.")
        return

    hist = px.histogram(
        plot_df,
        x=filt.col_age,
        color=filt.col_disease,
        barmode="stack",
        nbins=30,
        color_discrete_sequence=px.colors.qualitative.Set2,
        labels={filt.col_age: "Chronological age (years)", filt.col_disease: "Disease status"},
        title="Age distribution (stacked by disease status)",
    )
    hist.update_layout(legend_title_text="Disease status", bargap=0.05)
    st.plotly_chart(hist, use_container_width=True)

    st.markdown("#### Clinical variable correlation structure")
    heat_cols = clinical_numeric_columns_for_heatmap(df)
    if len(heat_cols) < 2:
        st.info("Need at least two numeric clinical columns with variance for a heatmap.")
        return

    corr = df[heat_cols].corr(method="pearson", numeric_only=True)
    fig_corr, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(
        corr,
        ax=ax,
        cmap="vlag",
        center=0.0,
        square=True,
        linewidths=0.5,
        cbar_kws={"label": "Pearson r"},
    )
    ax.set_title("Pairwise Pearson correlation (continuous clinical features)")
    st.pyplot(fig_corr, clear_figure=True)
    plt.close(fig_corr)


def render_tab_epigenetic_clock_validation(df: pd.DataFrame, filt: GlobalFilters) -> None:
    """Tab 2 — chronological vs epigenetic age with OLS fit and accuracy metrics."""
    st.subheader("Epigenetic clock validation")
    st.markdown(
        """
We compare **chronological age** to **epigenetic age** estimated from the merged methylation layer.
A linear fit quantifies systematic bias; **MAE** and **Pearson r** summarize agreement for the
currently filtered cohort.
        """.strip()
    )

    if "Epigenetic_Age" not in df.columns:
        st.error("Merged cohort lacks `Epigenetic_Age` (or an accepted alias).")
        return

    sub = df[[filt.col_age, "Epigenetic_Age"]].copy()
    sub[filt.col_age] = pd.to_numeric(sub[filt.col_age], errors="coerce")
    sub["Epigenetic_Age"] = pd.to_numeric(sub["Epigenetic_Age"], errors="coerce")
    sub = sub.dropna()
    if len(sub) < 3:
        st.warning("Insufficient paired age observations under current filters.")
        return

    x = sub[filt.col_age].to_numpy(dtype=np.float64)
    y = sub["Epigenetic_Age"].to_numpy(dtype=np.float64)
    slope, intercept = np.polyfit(x, y, 1)
    y_hat = slope * x + intercept
    mae = float(np.mean(np.abs(y - y_hat)))
    r_value, _p_r = pearsonr(x, y)

    c1, c2, c3 = st.columns(3)
    c1.metric("Mean absolute error (years)", f"{mae:.2f}")
    c2.metric("Pearson r", f"{r_value:.3f}")
    c3.metric("N (paired)", f"{len(sub):,}")

    line_x = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    line_y = slope * line_x + intercept

    scatter = px.scatter(
        sub,
        x=filt.col_age,
        y="Epigenetic_Age",
        template="plotly_white",
        labels={
            filt.col_age: "Chronological age (years)",
            "Epigenetic_Age": "Epigenetic age (years)",
        },
        title="Chronological age vs epigenetic age",
    )
    scatter.update_traces(marker=dict(size=8, opacity=0.65), name="Samples")
    scatter.add_trace(
        go.Scatter(
            x=line_x,
            y=line_y,
            mode="lines",
            name=f"OLS fit (y = {slope:.3f}x + {intercept:.2f})",
            line=dict(color="#DC2626", width=3),
        )
    )
    scatter.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    st.plotly_chart(scatter, use_container_width=True)


def _genotype_labels(series: pd.Series) -> pd.Series:
    """Map dosage 0/1/2 or strings to display labels for stratified boxplots."""
    out: list[str] = []
    for v in series.tolist():
        if pd.isna(v):
            out.append("Missing")
            continue
        if isinstance(v, (int, np.integer)):
            if int(v) == 0:
                out.append("0/0")
            elif int(v) == 1:
                out.append("0/1")
            elif int(v) == 2:
                out.append("1/1")
            else:
                out.append(str(v))
            continue
        s = str(v).strip()
        if re.fullmatch(r"[01]/[01]", s):
            out.append(s)
        else:
            out.append(s)
    return pd.Series(out, index=series.index, dtype=object)


def render_tab_la_snp_genotype_phenotype(df: pd.DataFrame, _filt: GlobalFilters) -> None:
    """Tab 3 — local ancestry–informed SNP genotypes vs clinical traits."""
    st.subheader("Genotype–phenotype impact (LA-SNPs)")
    st.markdown(
        """
Select a **SNP column** from the merged integration table and a **continuous outcome** to inspect
dosage-stratified distributions. Genotypes are shown as **0/0**, **0/1**, and **1/1** when encoded
as dosage; a **Kruskal–Wallis** test provides a non-parametric omnibus p-value across groups.
        """.strip()
    )

    snp_cols = list_snp_columns(df.columns)
    if not snp_cols:
        st.warning("No SNP-like columns (matching `rs####`) were found in the merged cohort.")
        return

    trait_candidates = [
        c
        for c in df.select_dtypes(include=["number"]).columns.tolist()
        if c not in (filt.col_age,) and not re.search(r"rs\d+", str(c), flags=re.IGNORECASE)
    ]
    trait_candidates = sorted(set(trait_candidates))
    if not trait_candidates:
        st.warning("No numeric trait columns available for stratified plots.")
        return

    c1, c2 = st.columns(2)
    with c1:
        snp = st.selectbox("SNP / variant column", options=snp_cols, index=0)
    with c2:
        default_trait = "HDL_Cholesterol" if "HDL_Cholesterol" in trait_candidates else trait_candidates[0]
        trait = st.selectbox("Clinical / molecular trait", options=trait_candidates, index=trait_candidates.index(default_trait))

    plot_df = df[[snp, trait]].copy()
    plot_df = plot_df.rename(columns={snp: "_geno", trait: "_trait"})
    plot_df["_trait"] = pd.to_numeric(plot_df["_trait"], errors="coerce")
    plot_df = plot_df.dropna(subset=["_trait"])
    plot_df["Genotype"] = _genotype_labels(plot_df["_geno"])
    plot_df = plot_df[plot_df["Genotype"] != "Missing"].reset_index(drop=True)

    groups = [g["_trait"].to_numpy(dtype=np.float64) for _, g in plot_df.groupby("Genotype", sort=True) if len(g) >= 2]
    p_text: str
    if len(groups) < 2:
        p_text = "Not enough genotype groups with ≥2 samples for Kruskal–Wallis."
    else:
        h_stat, p_kw = kruskal(*groups)
        p_text = f"Kruskal–Wallis H = {h_stat:.3f}, p = {p_kw:.3g} (omnibus across genotype groups)."

    box = px.box(
        plot_df,
        x="Genotype",
        y="_trait",
        color="Genotype",
        points="outliers",
        title=f"{trait} by genotype at {snp}",
        labels={"_trait": trait, "Genotype": "Genotype"},
    )
    box.update_layout(showlegend=False, template="plotly_white")
    st.plotly_chart(box, use_container_width=True)
    st.caption(p_text)

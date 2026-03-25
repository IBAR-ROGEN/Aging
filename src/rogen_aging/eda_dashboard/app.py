"""Official Streamlit EDA dashboard for the ROGEN multi-omics aging merged cohort."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from rogen_aging.eda_dashboard.data import default_merged_parquet_path, load_merged_parquet, load_synthetic_cohort
from rogen_aging.eda_dashboard.sidebar import apply_global_filters, render_global_sidebar
from rogen_aging.eda_dashboard.tabs import (
    render_tab_clinical_overview,
    render_tab_epigenetic_clock_validation,
    render_tab_la_snp_genotype_phenotype,
)

_USE_SYNTHETIC = "eda_use_synthetic_cohort"
_CUSTOM_PATH = "eda_custom_parquet_path"


def _init_session_state(default_path: Path) -> None:
    if _USE_SYNTHETIC not in st.session_state:
        st.session_state[_USE_SYNTHETIC] = False
    if _CUSTOM_PATH not in st.session_state:
        st.session_state[_CUSTOM_PATH] = str(default_path)


def main() -> None:
    st.set_page_config(
        page_title="ROGEN Multi-Omics Aging — EDA",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    default_path = default_merged_parquet_path()
    _init_session_state(default_path)

    st.title("ROGEN multi-omics aging cohort — exploratory data analysis")
    st.markdown(
        """
This dashboard is the **primary interactive QA layer** for the merged integration Parquet produced
by the ROGEN multi-omics pipeline. Use the **sidebar** to define a global analytic subset; every
visualization below respects those filters.
        """.strip()
    )

    st.sidebar.divider()
    path_input = st.sidebar.text_input(
        "Merged cohort Parquet path",
        value=st.session_state[_CUSTOM_PATH],
        help="Override default location or set env `ROGEN_MERGED_COHORT_PARQUET`.",
    )
    st.session_state[_CUSTOM_PATH] = path_input.strip() or str(default_path)

    cohort_path = Path(st.session_state[_CUSTOM_PATH]).expanduser()
    parquet_exists = cohort_path.is_file()

    df: pd.DataFrame | None = None

    if st.session_state[_USE_SYNTHETIC]:
        df = load_synthetic_cohort()
    elif parquet_exists:
        df = load_merged_parquet(str(cohort_path.resolve()))
    else:
        st.error(
            f"Merged Parquet not found at `{cohort_path}`. "
            "Confirm the integration pipeline output or paste a valid path in the sidebar."
        )
        if st.button("Generate in-memory mock cohort", type="primary"):
            st.session_state[_USE_SYNTHETIC] = True
            st.rerun()
        st.stop()

    assert df is not None

    st.sidebar.success(
        "Using **synthetic in-memory cohort**."
        if st.session_state[_USE_SYNTHETIC]
        else f"Loaded **`{cohort_path.name}`** ({len(df):,} rows)."
    )
    if st.session_state[_USE_SYNTHETIC] and st.sidebar.button("Switch back to Parquet (if available)"):
        st.session_state[_USE_SYNTHETIC] = False
        st.rerun()

    filt = render_global_sidebar(df)
    filtered = apply_global_filters(df, filt)

    st.caption(f"Filtered cohort: **{len(filtered):,}** / {len(df):,} participants.")

    tab_clin, tab_clock, tab_snp = st.tabs(
        [
            "Clinical & phenotypic overview",
            "Epigenetic clock validation",
            "Genotype–phenotype (LA-SNPs)",
        ]
    )

    with tab_clin:
        render_tab_clinical_overview(filtered, filt)
    with tab_clock:
        render_tab_epigenetic_clock_validation(filtered, filt)
    with tab_snp:
        render_tab_la_snp_genotype_phenotype(filtered, filt)


if __name__ == "__main__":
    main()

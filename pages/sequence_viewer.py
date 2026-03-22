"""
Sequence Viewer — ProteoClinView
Uses OpenMS-Insight SequenceView + LinePlot + StateManager for cross-component linking.
"""
import sys
import hashlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import polars as pl
import numpy as np

from openms_insight import SequenceView, LinePlot, StateManager

from utils.fragments import (
    calculate_fragments, generate_demo_spectrum,
    DEMO_PROTEINS, sequence_coverage_pct, AA_MASSES,
)
from utils.plots import make_sequence_map

CACHE_DIR = Path("/tmp/proteoclinview_cache")
CACHE_DIR.mkdir(exist_ok=True)

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## ⚙️ Settings")

    demo_mode     = st.session_state.get("demo_mode", True)
    demo_protein  = st.session_state.get("demo_protein", "Human Insulin B chain")

    # Step 1
    st.markdown("### Step 1 · Protein sequence")
    input_mode = st.radio("Input mode", ["Demo protein", "Paste sequence"],
                          index=0 if demo_mode else 1)

    if input_mode == "Demo protein":
        protein_name = st.selectbox(
            "Select protein", list(DEMO_PROTEINS.keys()),
            index=list(DEMO_PROTEINS.keys()).index(demo_protein)
                  if demo_protein in DEMO_PROTEINS else 0,
        )
        info     = DEMO_PROTEINS[protein_name]
        sequence = info["sequence"]
        st.caption(f"**{info['description']}**  \nUniProt: `{info['uniprot']}`")
        st.info(f"🏥 **Clinical context:** {info['clinical_note']}")
    else:
        sequence = st.text_area(
            "Amino acid sequence (single-letter codes)",
            value="PEPTIDE", height=90,
        ).upper().replace(" ", "").replace("\n", "")

    # Step 2
    st.markdown("### Step 2 · Fragment ion types")
    st.caption(
        "**c / z** → ETD/ECD (top-down proteomics)  \n"
        "**b / y** → CID/HCD (bottom-up proteomics)"
    )
    col_l, col_r = st.columns(2)
    with col_l:
        use_c = st.checkbox("c-ions 🔵", value=True)
        use_b = st.checkbox("b-ions 🟣", value=False)
    with col_r:
        use_z = st.checkbox("z-ions 🟠", value=True)
        use_y = st.checkbox("y-ions 🟢", value=False)

    ion_types = (["c"] if use_c else []) + (["z"] if use_z else []) + \
                (["b"] if use_b else []) + (["y"] if use_y else [])
    if not ion_types:
        ion_types = ["c", "z"]

    # Step 3
    st.markdown("### Step 3 · Charge states")
    max_charge = st.slider("Maximum charge (z)", 1, 5, 2)
    charges    = list(range(1, max_charge + 1))

    # Step 4
    st.markdown("### Step 4 · Mass tolerance")
    tol_ppm = st.slider("Tolerance (ppm)", 1, 50, 10)

    st.markdown("---")
    st.caption("ProteoClinView · OpenMS GSoC 2026")

# ── Validate ──────────────────────────────────────────────────────────────────
bad = [aa for aa in sequence if aa not in AA_MASSES]
if bad:
    st.error(f"Unknown amino acid(s): {sorted(set(bad))}")
    st.stop()
if len(sequence) < 2:
    st.warning("Please enter a sequence with at least 2 residues.")
    st.stop()

# ── Compute theoretical + synthetic experimental ───────────────────────────────
try:
    theoretical = calculate_fragments(sequence, ion_types=ion_types, charges=charges)
except Exception as e:
    st.error(f"Fragment calculation error: {e}")
    st.stop()

experimental = generate_demo_spectrum(theoretical, seed=42)

# ── Build polars DataFrames for OpenMS-Insight ────────────────────────────────
# peaks_data: peak_id (int), mass (float = m/z), intensity (float)
exp_signal = experimental[experimental["ion_type"] != "noise"].copy()
peaks_df = pl.DataFrame({
    "peak_id":   list(range(len(exp_signal))),
    "mass":      exp_signal["mz"].tolist(),
    "intensity": exp_signal["intensity"].tolist(),
})

# ── Cache-key: regenerate when sequence/params change ─────────────────────────
cache_key = hashlib.md5(
    f"{sequence}|{''.join(ion_types)}|{max_charge}|{tol_ppm}".encode()
).hexdigest()[:8]

cache_path = str(CACHE_DIR)

# ── OpenMS-Insight components ─────────────────────────────────────────────────
state_manager = StateManager("pcv_state")

sequence_view = SequenceView(
    cache_id=f"pcv_sv_{cache_key}",
    sequence_data=sequence,          # plain string — simplest format
    peaks_data=peaks_df.lazy(),
    interactivity={"peak": "peak_id"},
    annotation_config={
        "ion_types": ion_types,
        "tolerance": float(tol_ppm),
        "tolerance_ppm": True,
        "neutral_losses": False,
    },
    cache_path=cache_path,
    title=f"Sequence coverage — {sequence[:15]}{'...' if len(sequence) > 15 else ''}",
    height=max(250, (len(sequence) // 10 + 1) * 70 + 80),
)

annotated_plot = LinePlot.from_sequence_view(
    sequence_view,
    cache_id=f"pcv_lp_{cache_key}",
    cache_path=cache_path,
    title="Annotated spectrum (click a peak or residue to cross-highlight)",
    x_label="m/z",
    y_label="Intensity",
    styling={
        "unhighlightedColor": "#94A3B8",
        "highlightColor":     "#3B82F6",
        "selectedColor":      "#F59E0B",
    },
)

# ── Page header ───────────────────────────────────────────────────────────────
st.markdown("# 🧬 Sequence Viewer")
st.caption(
    "Built on **OpenMS-Insight** `SequenceView` + `LinePlot` with `StateManager` cross-linking.  \n"
    "Click a residue in the sequence map → the corresponding peak is highlighted in the spectrum.  \n"
    "Click a peak in the spectrum → the covered residues are highlighted in the sequence."
)

# ── Metrics ───────────────────────────────────────────────────────────────────
from utils.fragments import match_peaks
matched = match_peaks(theoretical, experimental, tolerance_ppm=tol_ppm)
n_matched = matched["matched"].sum()
cov       = sequence_coverage_pct(sequence, matched)

m1, m2, m3, m4 = st.columns(4)
m1.metric("Sequence length",       f"{len(sequence)} AA")
m2.metric("Theoretical ions",      str(len(theoretical)))
m3.metric("Matched in demo spectrum", f"{n_matched} / {len(theoretical)}")
m4.metric("Sequence coverage",     f"{cov:.0%}")

st.markdown("---")

# ── OpenMS-Insight rendering ──────────────────────────────────────────────────
st.markdown("### Sequence coverage map  ·  OpenMS-Insight SequenceView")
st.caption(
    "Coloured residues = covered by matched fragment ions.  "
    "Ion type colours configured via `annotation_config`.  "
    "Cross-linked to the annotated spectrum below via **StateManager**."
)

sv_result = sequence_view(key="pcv_sv", state_manager=state_manager)

st.markdown("### Annotated spectrum  ·  OpenMS-Insight LinePlot")
st.caption(
    "Peaks coloured by ion type annotation from SequenceView.  "
    "Gold peak = currently selected.  "
    "Powered by `LinePlot.from_sequence_view()`."
)

annotated_plot(key="pcv_lp", state_manager=state_manager, sequence_view_key="pcv_sv")

# ── Selected peak info ────────────────────────────────────────────────────────
selected_peak = state_manager.get_selection("peak")
if selected_peak is not None:
    row = peaks_df.filter(pl.col("peak_id") == int(selected_peak))
    if not row.is_empty():
        mz  = row["mass"][0]
        inten = row["intensity"][0]
        # Find matching theoretical ion
        match_row = matched[
            abs(matched["mz"] - mz) / mz * 1e6 < tol_ppm * 2
        ]
        st.info(
            f"**Selected peak:** m/z = {mz:.4f}  ·  "
            f"Intensity = {inten:.2%}  ·  "
            f"Matched ion: **{match_row['ion_label'].iloc[0] if len(match_row) > 0 else 'unassigned'}**"
        )

st.markdown("---")

# ── Supplementary: Plotly coverage map ───────────────────────────────────────
with st.expander("📊 Supplementary: Plotly coverage map (for reference)", expanded=False):
    st.caption(
        "This coverage map uses our custom Plotly renderer. "
        "The OpenMS-Insight SequenceView above is the primary view."
    )
    fig_map = make_sequence_map(sequence, matched[matched["matched"]], residues_per_row=15)
    st.plotly_chart(fig_map, use_container_width=True, config={"displayModeBar": False})

# ── Fragment table ────────────────────────────────────────────────────────────
st.markdown("---")
with st.expander("🔍 Fragment ion table", expanded=False):
    display = matched[["ion_label", "ion_type", "charge", "mz",
                        "matched", "position_from", "position_to", "sequence_covered"]].copy()
    display.columns = ["Ion", "Type", "z+", "Theoretical m/z",
                       "Matched", "From (AA)", "To (AA)", "Covered sequence"]
    st.dataframe(display.reset_index(drop=True), use_container_width=True, height=280)

    csv = display.to_csv(index=False).encode("utf-8")
    st.download_button("⬇️ Download as CSV", data=csv,
                       file_name="proteoclinview_fragments.csv", mime="text/csv")

# ── Navigation ────────────────────────────────────────────────────────────────
st.markdown("---")
col1, col2 = st.columns([2, 1])
with col1:
    st.markdown("**Next:** compare experimental vs theoretical in the Mirror Plot.")
with col2:
    if st.button("→ Mirror Plot", type="primary", use_container_width=True):
        st.session_state["mirror_sequence"]  = sequence
        st.session_state["mirror_ion_types"] = ion_types
        st.session_state["mirror_charges"]   = charges
        st.switch_page("pages/mirror_plot.py")

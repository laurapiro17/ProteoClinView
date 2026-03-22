"""
Sequence Viewer — ProteoClinView
Uses OpenMS-Insight SequenceView + LinePlot + StateManager for cross-component linking.
Supports demo mode (synthetic spectrum) and real mzML file upload.
"""
import sys
import hashlib
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import polars as pl
import numpy as np
try:
    import pyopenms as oms
    OMS_AVAILABLE = True
except Exception:
    OMS_AVAILABLE = False

from openms_insight import SequenceView, LinePlot, StateManager

from utils.fragments import (
    calculate_fragments, generate_demo_spectrum,
    DEMO_PROTEINS, sequence_coverage_pct, AA_MASSES,
)
from utils.plots import make_sequence_map

CACHE_DIR = Path("/tmp/proteoclinview_cache")
CACHE_DIR.mkdir(exist_ok=True)

EXAMPLE_MZML = Path(__file__).parent.parent / "example_data" / "insulin_b_chain_topdown.mzML"


def load_mzml(path: str) -> list:
    """Load mzML and return list of MS2 spectra as dicts."""
    if not OMS_AVAILABLE:
        return []
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    spectra = []
    for i, s in enumerate(exp.getSpectra()):
        if s.getMSLevel() == 2 and s.size() > 0:
            mz_arr, int_arr = s.get_peaks()
            precs = s.getPrecursors()
            precursor_mz = precs[0].getMZ() if precs else 0.0
            spectra.append({
                "index": i,
                "scan_id": s.getNativeID() or f"scan={i}",
                "rt": round(s.getRT(), 2),
                "precursor_mz": round(precursor_mz, 4),
                "n_peaks": s.size(),
                "mz": mz_arr,
                "intensity": int_arr,
            })
    return spectra


# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## Settings")

    demo_mode    = st.session_state.get("demo_mode", True)
    demo_protein = st.session_state.get("demo_protein", "Human Insulin B chain")

    # Step 1 — Protein sequence
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
        st.info(f"**Clinical context:** {info['clinical_note']}")
    else:
        sequence = st.text_area(
            "Amino acid sequence (single-letter codes)",
            value="PEPTIDE", height=90,
        ).upper().replace(" ", "").replace("\n", "")

    # Step 2 — Ion types
    st.markdown("### Step 2 · Fragment ion types")
    st.caption(
        "**c / z** — ETD/ECD (top-down proteomics)  \n"
        "**b / y** — CID/HCD (bottom-up proteomics)"
    )
    col_l, col_r = st.columns(2)
    with col_l:
        use_c = st.checkbox("c-ions", value=True)
        use_b = st.checkbox("b-ions", value=False)
    with col_r:
        use_z = st.checkbox("z-ions", value=True)
        use_y = st.checkbox("y-ions", value=False)

    ion_types = (["c"] if use_c else []) + (["z"] if use_z else []) + \
                (["b"] if use_b else []) + (["y"] if use_y else [])
    if not ion_types:
        ion_types = ["c", "z"]

    # Step 3 — Charge states
    st.markdown("### Step 3 · Charge states")
    max_charge = st.slider("Maximum charge (z)", 1, 5, 2)
    charges    = list(range(1, max_charge + 1))

    # Step 4 — Spectrum source
    st.markdown("### Step 4 · Spectrum")
    spectrum_source = st.radio(
        "Spectrum source",
        ["Synthetic demo", "Load example mzML", "Upload mzML"],
    )

    selected_scan = None
    mzml_spectra  = []

    if spectrum_source in ("Load example mzML", "Upload mzML") and not OMS_AVAILABLE:
        st.info("mzML loading requires pyopenms (not available in this environment). Use Synthetic demo mode.")
        spectrum_source = "Synthetic demo"

    if spectrum_source == "Load example mzML":
        if EXAMPLE_MZML.exists():
            mzml_spectra = load_mzml(str(EXAMPLE_MZML))
            st.caption(f"Loaded: `insulin_b_chain_topdown.mzML`  \n{len(mzml_spectra)} MS2 spectrum/a")
        else:
            st.warning("Example file not found.")

    elif spectrum_source == "Upload mzML":
        uploaded = st.file_uploader("Upload mzML file", type=["mzML"])
        if uploaded:
            with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False) as tmp:
                tmp.write(uploaded.read())
                tmp_path = tmp.name
            mzml_spectra = load_mzml(tmp_path)
            st.caption(f"{len(mzml_spectra)} MS2 spectra found")

    if mzml_spectra:
        scan_labels = [f"{s['scan_id']}  RT={s['rt']}  prec={s['precursor_mz']}  peaks={s['n_peaks']}"
                       for s in mzml_spectra]
        sel_idx = st.selectbox("Select MS2 spectrum", range(len(scan_labels)),
                               format_func=lambda i: scan_labels[i])
        selected_scan = mzml_spectra[sel_idx]

    # Step 5 — Tolerance
    st.markdown("### Step 5 · Mass tolerance")
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

# ── Compute theoretical fragments ─────────────────────────────────────────────
try:
    theoretical = calculate_fragments(sequence, ion_types=ion_types, charges=charges)
except Exception as e:
    st.error(f"Fragment calculation error: {e}")
    st.stop()

# ── Build experimental peaks DataFrame ────────────────────────────────────────
if selected_scan is not None:
    # Real mzML data
    mz_arr  = np.array(selected_scan["mz"],       dtype=float)
    int_arr = np.array(selected_scan["intensity"], dtype=float)
    int_norm = int_arr / int_arr.max() if int_arr.max() > 0 else int_arr
    peaks_df = pl.DataFrame({
        "peak_id":   list(range(len(mz_arr))),
        "mass":      mz_arr.tolist(),
        "intensity": int_norm.tolist(),
    })
    spectrum_label = f"{selected_scan['scan_id']} (mzML)"
else:
    # Synthetic demo
    experimental = generate_demo_spectrum(theoretical, seed=42)
    exp_signal   = experimental[experimental["ion_type"] != "noise"].copy()
    peaks_df = pl.DataFrame({
        "peak_id":   list(range(len(exp_signal))),
        "mass":      exp_signal["mz"].tolist(),
        "intensity": exp_signal["intensity"].tolist(),
    })
    spectrum_label = "synthetic demo"

# ── Cache key ─────────────────────────────────────────────────────────────────
cache_key = hashlib.md5(
    f"{sequence}|{''.join(ion_types)}|{max_charge}|{tol_ppm}|{spectrum_label}".encode()
).hexdigest()[:8]

# ── OpenMS-Insight components ─────────────────────────────────────────────────
state_manager = StateManager("pcv_state")

sequence_view = SequenceView(
    cache_id=f"pcv_sv_{cache_key}",
    sequence_data=sequence,
    peaks_data=peaks_df.lazy(),
    interactivity={"peak": "peak_id"},
    annotation_config={
        "ion_types": ion_types,
        "tolerance": float(tol_ppm),
        "tolerance_ppm": True,
        "neutral_losses": False,
    },
    cache_path=str(CACHE_DIR),
    title=f"Sequence coverage — {sequence[:15]}{'...' if len(sequence) > 15 else ''}",
    height=max(250, (len(sequence) // 10 + 1) * 70 + 80),
)

annotated_plot = LinePlot.from_sequence_view(
    sequence_view,
    cache_id=f"pcv_lp_{cache_key}",
    cache_path=str(CACHE_DIR),
    title=f"Annotated spectrum — {spectrum_label}",
    x_label="m/z",
    y_label="Intensity",
    styling={
        "unhighlightedColor": "#94A3B8",
        "highlightColor":     "#3B82F6",
        "selectedColor":      "#F59E0B",
    },
)

# ── Page ──────────────────────────────────────────────────────────────────────
st.markdown("# Sequence Viewer")
st.caption(
    "Built on **OpenMS-Insight** `SequenceView` + `LinePlot` with `StateManager` cross-linking.  \n"
    "Click a residue to highlight the corresponding peak — click a peak to highlight the covered residues."
)

# Metrics
from utils.fragments import match_peaks
if selected_scan is not None:
    exp_for_match = pl.DataFrame({
        "ion_type": ["signal"] * len(peaks_df),
        "mz": peaks_df["mass"].to_list(),
        "intensity": peaks_df["intensity"].to_list(),
        "is_noise": [False] * len(peaks_df),
    })
    import pandas as pd
    exp_pd = pd.DataFrame({"ion_type": "signal", "mz": peaks_df["mass"].to_list(),
                            "intensity": peaks_df["intensity"].to_list(), "is_noise": False})
    matched = match_peaks(theoretical, exp_pd, tolerance_ppm=tol_ppm)
else:
    experimental = generate_demo_spectrum(theoretical, seed=42)
    matched = match_peaks(theoretical, experimental, tolerance_ppm=tol_ppm)

n_matched = matched["matched"].sum()
cov       = sequence_coverage_pct(sequence, matched)

m1, m2, m3, m4 = st.columns(4)
m1.metric("Sequence length",    f"{len(sequence)} AA")
m2.metric("Theoretical ions",   str(len(theoretical)))
m3.metric("Matched ions",       f"{n_matched} / {len(theoretical)}")
m4.metric("Sequence coverage",  f"{cov:.0%}")

if selected_scan is None:
    st.caption("Spectrum: synthetic demo — switch to 'Load example mzML' or upload your own in the sidebar.")

st.markdown("---")

# ── OpenMS-Insight rendering ──────────────────────────────────────────────────
st.markdown("### Sequence coverage map")
sv_result = sequence_view(key="pcv_sv", state_manager=state_manager)

st.markdown("### Annotated spectrum")
annotated_plot(key="pcv_lp", state_manager=state_manager, sequence_view_key="pcv_sv")

# Selected peak info
selected_peak = state_manager.get_selection("peak")
if selected_peak is not None:
    row = peaks_df.filter(pl.col("peak_id") == int(selected_peak))
    if not row.is_empty():
        mz_sel = row["mass"][0]
        match_row = matched[abs(matched["mz"] - mz_sel) / mz_sel * 1e6 < tol_ppm * 2]
        ion_label = match_row["ion_label"].iloc[0] if len(match_row) > 0 else "unassigned"
        st.info(
            f"**Selected peak:** m/z = {mz_sel:.4f}  ·  "
            f"Intensity = {row['intensity'][0]:.2%}  ·  "
            f"Matched ion: **{ion_label}**"
        )

st.markdown("---")

# Supplementary Plotly map
with st.expander("Supplementary: Plotly coverage map", expanded=False):
    fig_map = make_sequence_map(sequence, matched[matched["matched"]], residues_per_row=15)
    st.plotly_chart(fig_map, use_container_width=True, config={"displayModeBar": False})

# Fragment table
with st.expander("Fragment ion table", expanded=False):
    display = matched[["ion_label", "ion_type", "charge", "mz",
                        "matched", "position_from", "position_to", "sequence_covered"]].copy()
    display.columns = ["Ion", "Type", "z+", "Theoretical m/z",
                       "Matched", "From (AA)", "To (AA)", "Covered sequence"]
    st.dataframe(display.reset_index(drop=True), use_container_width=True, height=280)
    csv = display.to_csv(index=False).encode("utf-8")
    st.download_button("Download as CSV", data=csv,
                       file_name="proteoclinview_fragments.csv", mime="text/csv")

st.markdown("---")
col1, col2 = st.columns([2, 1])
with col1:
    st.markdown("**Next:** compare experimental vs theoretical in the Mirror Plot.")
with col2:
    if st.button("Mirror Plot ->", type="primary", use_container_width=True):
        st.session_state["mirror_sequence"]  = sequence
        st.session_state["mirror_ion_types"] = ion_types
        st.session_state["mirror_charges"]   = charges
        st.switch_page("pages/mirror_plot.py")

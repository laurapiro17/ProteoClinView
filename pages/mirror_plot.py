"""
Mirror Plot — ProteoClinView
Compare experimental spectrum vs theoretical fragment ions.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import streamlit as st
import pandas as pd
import numpy as np

from utils.fragments import (
    calculate_fragments, generate_demo_spectrum, match_peaks,
    DEMO_PROTEINS, sequence_coverage_pct,
)
from utils.plots import make_mirror_plot

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## Settings")

    # Carry over from sequence viewer if navigated from there
    default_seq   = st.session_state.get("mirror_sequence", "FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
    default_ions  = st.session_state.get("mirror_ion_types", ["c", "z"])
    default_chrgs = st.session_state.get("mirror_charges", [1, 2])

    st.markdown("### Step 1 · Protein sequence")
    input_mode = st.radio(
        "Input mode",
        ["Demo protein", "Paste sequence"],
        index=0,
    )

    if input_mode == "Demo protein":
        protein_name = st.selectbox("Select protein", list(DEMO_PROTEINS.keys()))
        info = DEMO_PROTEINS[protein_name]
        sequence = info["sequence"]
        st.caption(f"**{info['description']}**")
        st.info(info['clinical_note'])
    else:
        sequence = st.text_area(
            "Amino acid sequence",
            value=default_seq, height=90,
        ).upper().replace(" ", "").replace("\n", "")

    st.markdown("### Step 2 · Ion types & charges")
    col_c, col_z = st.columns(2)
    with col_c:
        use_c = st.checkbox("c-ions", value="c" in default_ions)
        use_b = st.checkbox("b-ions", value="b" in default_ions)
    with col_z:
        use_z = st.checkbox("z-ions", value="z" in default_ions)
        use_y = st.checkbox("y-ions", value="y" in default_ions)

    ion_types = (["c"] if use_c else []) + (["z"] if use_z else []) + \
                (["b"] if use_b else []) + (["y"] if use_y else [])
    if not ion_types:
        ion_types = ["c", "z"]

    max_charge = st.slider("Max charge (z)", 1, 5,
                           value=max(default_chrgs) if default_chrgs else 2)
    charges = list(range(1, max_charge + 1))

    st.markdown("### Step 3 · Experimental spectrum")
    exp_source = st.radio(
        "Spectrum source",
        ["Synthetic demo (auto-generated)", "Paste peak list (m/z, intensity)"],
        help="Demo mode generates a realistic synthetic spectrum from the theoretical ions.",
    )

    if exp_source == "Paste peak list (m/z, intensity)":
        raw_peaks = st.text_area(
            "Peak list (one peak per line: m/z intensity)",
            placeholder="200.123 0.85\n201.456 0.32\n...",
            height=140,
        )
    else:
        raw_peaks = None

    st.markdown("### Step 4 · Matching parameters")
    tolerance_ppm = st.slider(
        "Mass tolerance (ppm)",
        min_value=1, max_value=50, value=10,
        help="Maximum allowed mass error between experimental and theoretical m/z.",
    )
    st.caption(
        f"At 10 ppm, a peak at m/z 500 matches theoretical ions within ±{500*10*1e-6*1000:.2f} mDa."
    )

    st.markdown("---")
    st.caption("ProteoClinView · OpenMS GSoC 2026")

# ── Main ──────────────────────────────────────────────────────────────────────
st.markdown("# Mirror Plot")
st.caption(
    "**Top half:** experimental spectrum — coloured peaks are matched fragment ions.  \n"
    "**Bottom half:** theoretical ions (mirrored) — solid = matched, dotted = unmatched.  \n"
    "Hover any peak to see its label, m/z, and mass error in ppm."
)

# Validate
from utils.fragments import AA_MASSES
bad = [aa for aa in sequence if aa not in AA_MASSES]
if bad:
    st.error(f"Unknown amino acid(s): {sorted(set(bad))}")
    st.stop()
if len(sequence) < 2:
    st.warning("Please enter a sequence with at least 2 residues.")
    st.stop()

# Compute theoretical
try:
    theoretical = calculate_fragments(sequence, ion_types=ion_types, charges=charges)
except Exception as e:
    st.error(f"Fragment calculation error: {e}")
    st.stop()

# Experimental spectrum
if raw_peaks:
    try:
        lines = [l.strip() for l in raw_peaks.strip().splitlines() if l.strip()]
        mz_list, int_list = [], []
        for line in lines:
            parts = line.split()
            mz_list.append(float(parts[0]))
            int_list.append(float(parts[1]) if len(parts) > 1 else 1.0)
        max_int = max(int_list) if int_list else 1.0
        experimental = pd.DataFrame({
            "ion_type": "signal",
            "ion_label": [f"peak_{i}" for i in range(len(mz_list))],
            "position_from": 0, "position_to": 0, "sequence_covered": "",
            "mz": mz_list,
            "charge": 1, "neutral_mass": 0.0, "color": "#9CA3AF",
            "intensity": [v / max_int for v in int_list],
            "is_noise": False,
        })
    except Exception as e:
        st.error(f"Could not parse peak list: {e}")
        st.stop()
else:
    experimental = generate_demo_spectrum(theoretical, seed=42)

# Match peaks
matched = match_peaks(theoretical, experimental, tolerance_ppm=tolerance_ppm)

# ── Metrics ───────────────────────────────────────────────────────────────────
n_matched   = matched["matched"].sum()
n_total     = len(theoretical)
match_pct   = n_matched / n_total * 100 if n_total > 0 else 0
cov         = sequence_coverage_pct(sequence, matched)
n_exp_peaks = len(experimental[~experimental.get("is_noise", pd.Series(True, index=experimental.index))])

m1, m2, m3, m4 = st.columns(4)
m1.metric("Matched ions", f"{n_matched} / {n_total}", f"{match_pct:.0f}%")
m2.metric("Sequence coverage", f"{cov:.0%}")
m3.metric("Experimental peaks", str(n_exp_peaks))
m4.metric("Mass tolerance", f"{tolerance_ppm} ppm")

if exp_source == "Synthetic demo (auto-generated)":
    st.info(
        "**Demo mode** — the experimental spectrum is synthetically generated from the theoretical ions "
        "(with ~25% random missing peaks, <5 ppm m/z noise, and background noise peaks). "
        "Upload your own peak list in the sidebar for real data.",
        icon=None,
    )

st.markdown("---")

# ── Mirror plot ───────────────────────────────────────────────────────────────
fig = make_mirror_plot(
    theoretical_df=theoretical,
    experimental_df=experimental,
    matched_df=matched,
    title=f"Mirror Plot — {sequence[:12]}{'...' if len(sequence)>12 else ''}",
    tolerance_ppm=tolerance_ppm,
)
st.plotly_chart(fig, use_container_width=True)

# ── Matched ion table ─────────────────────────────────────────────────────────
st.markdown("---")
st.markdown("### Matched fragment ions")
st.caption(
    "Each row is a theoretical fragment ion that was matched to an experimental peak within the "
    "specified mass tolerance. Δm/z (ppm) shows the measurement accuracy."
)

matched_display = matched[matched["matched"]].copy()
matched_display["delta_ppm"] = (
    (matched_display["exp_mz"] - matched_display["mz"]) / matched_display["mz"] * 1e6
).round(2)
matched_display["exp_intensity_pct"] = (matched_display["exp_intensity"] * 100).round(1)

cols_show = ["ion_label", "ion_type", "charge", "mz", "exp_mz",
             "delta_ppm", "exp_intensity_pct", "sequence_covered"]
col_names  = ["Ion", "Type", "z+", "Theoretical m/z", "Experimental m/z",
              "Δ (ppm)", "Intensity (%)", "Covered sequence"]

st.dataframe(
    matched_display[cols_show].rename(columns=dict(zip(cols_show, col_names))).reset_index(drop=True),
    use_container_width=True, height=280,
)

# Export
csv = matched_display[cols_show].rename(columns=dict(zip(cols_show, col_names))).to_csv(index=False).encode("utf-8")
st.download_button("Download matched ions as CSV", data=csv,
                   file_name="proteoclinview_matched.csv", mime="text/csv")

# ── Navigation ────────────────────────────────────────────────────────────────
st.markdown("---")
col_nav1, col_nav2 = st.columns([2, 1])
with col_nav1:
    st.markdown("**Go back** to the Sequence Viewer to adjust ion types or explore a different protein.")
with col_nav2:
    if st.button("<- Sequence Viewer", use_container_width=True):
        st.switch_page("pages/sequence_viewer.py")

import streamlit as st

# ── Hero ──────────────────────────────────────────────────────────────────────
col_hero, col_info = st.columns([2, 1], gap="large")

with col_hero:
    st.markdown("# 🔬 ProteoClinView")
    st.markdown("### Interactive Proteoform Visualization for Top-Down Proteomics")
    st.markdown(
        """
        **ProteoClinView** bridges the gap between computational proteomics and clinical research.
        Visualise intact proteoforms, annotate fragment ions, and validate top-down MS/MS results —
        designed to be readable by wet-lab scientists and clinical researchers, not just bioinformaticians.
        """
    )
    st.markdown("---")

    col_a, col_b = st.columns(2)
    with col_a:
        if st.button(
            "🚀 Launch Demo  ·  no upload needed",
            use_container_width=True,
            type="primary",
            help="Loads Human Insulin B chain with synthetic MS/MS data. No files required.",
        ):
            st.session_state["demo_mode"] = True
            st.session_state["demo_protein"] = "Human Insulin B chain"
            st.switch_page("pages/sequence_viewer.py")

    with col_b:
        if st.button(
            "📂 Upload my mzML / sequence",
            use_container_width=True,
            help="Paste your own amino acid sequence and explore fragment ions.",
        ):
            st.session_state["demo_mode"] = False
            st.switch_page("pages/sequence_viewer.py")

with col_info:
    st.markdown("### Modules")
    st.markdown(
        """
| Module | Status |
|--------|--------|
| 🧬 Sequence Viewer | ✅ Live |
| 📊 Mirror Plot | ✅ Live |
| 🔍 FLASHTnT Search | 🔜 Coming |
| 🗺️ Feature Map | 🔜 Coming |
        """
    )
    st.caption("Built on OpenMS-WebApps · pyOpenMS · Streamlit")

# ── Differentiators ───────────────────────────────────────────────────────────
st.markdown("---")
st.markdown("### Designed for every researcher in the lab")

c1, c2, c3 = st.columns(3)
with c1:
    st.markdown("**👩‍🔬 Wet-lab scientists**")
    st.markdown(
        "Plain-language labels, guided sidebar workflow, one-click demo with real proteins. "
        "No command line, no configuration files."
    )
with c2:
    st.markdown("**🖥️ Computational researchers**")
    st.markdown(
        "Full control over ion types (c, z, b, y), charge states, and mass tolerance. "
        "Export-ready Plotly figures."
    )
with c3:
    st.markdown("**🏥 Clinical researchers**")
    st.markdown(
        "Clinically annotated demo proteins: hemoglobin (sickle cell variants), insulin (diabetes). "
        "Designed for method development and result interpretation."
    )

# ── Quick explanation ─────────────────────────────────────────────────────────
with st.expander("ℹ️  What is top-down proteomics? (plain language)"):
    st.markdown(
        """
        In **top-down proteomics**, intact proteins (proteoforms) are fragmented inside the mass spectrometer
        and the resulting pieces — called **fragment ions** — are measured.

        By matching those pieces back to a known protein sequence, we can:
        - Pinpoint exactly **which form** of the protein is present (e.g. phosphorylated? truncated?)
        - Detect **disease-linked variants** (e.g. the E6V mutation in sickle-cell haemoglobin)
        - Quantify multiple proteoforms **in a single experiment**

        ProteoClinView shows you a colour-coded map of which parts of the protein were detected,
        and compares your experimental spectrum directly against theoretical predictions.
        """
    )

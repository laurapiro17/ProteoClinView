# 🔬 ProteoClinView

**Interactive proteoform-centric visualization for top-down proteomics**

Live demo → **https://laurapiro17-proteoclinview-app-tv7oqq.streamlit.app/**

GSoC 2026 proposal prototype · Organization: OpenMS · Project C1

---

## What it does

ProteoClinView bridges the gap between computational proteomics and clinical research.
Visualise intact proteoforms, annotate fragment ions, and validate top-down MS/MS results —
designed to be readable by wet-lab scientists and clinical researchers, not just bioinformaticians.

## Modules

| Module | Status | Description |
|--------|--------|-------------|
| 🧬 Sequence Viewer | ✅ Live | OpenMS-Insight `SequenceView` + `LinePlot` + `StateManager` cross-linking |
| 📊 Mirror Plot | ✅ Live | Experimental vs theoretical spectrum comparison with Δppm per peak |
| 🔍 FLASHTnT Search | 🔜 Planned | On-demand targeted search via pyOpenMS |
| 🗺️ Feature Map | 🔜 Planned | m/z vs RT proteoform feature visualization |

## Architecture

Built on the **OpenMS-WebApps** framework (same stack as TOPPView-Lite and FLASHApp):

- **Frontend:** Streamlit + OpenMS-Insight components (`SequenceView`, `LinePlot`, `StateManager`)
- **Fragment matching:** OpenMS-Insight Vue-side annotation engine
- **Data layer:** Polars LazyFrames; parquet caching for fast reload
- **Computation:** Monoisotopic fragment ion calculation (c, z, b, y ions; multi-charge)
- **Export:** CSV download of matched ions and fragment tables

## Demo proteins

| Protein | Length | Clinical relevance |
|---------|--------|-------------------|
| Human Insulin B chain | 30 AA | Diabetes research; neonatal diabetes mutations |
| Human Ubiquitin | 76 AA | Neurodegeneration; cancer (ubiquitin pathway) |
| Hemoglobin α N-terminal | 30 AA | Sickle cell disease (E6V variant); thalassemia |

## Run locally

```bash
git clone https://github.com/laurapiro17/ProteoClinView.git
cd ProteoClinView
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Author

Laura Piñero Roig · 3rd year Medicine, Universitat de Barcelona
GitHub: [@laurapiro17](https://github.com/laurapiro17)
GSoC 2026 proposal: OpenMS Project C1 — Interactive Proteoform-Centric Visualization

import streamlit as st
from pathlib import Path

st.set_page_config(
    page_title="ProteoClinView",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

pages = {
    "ProteoClinView": [
        st.Page(Path("pages", "welcome.py"),          title="Home",            icon="🏠", default=True),
        st.Page(Path("pages", "sequence_viewer.py"),  title="Sequence Viewer", icon="🧬"),
        st.Page(Path("pages", "mirror_plot.py"),      title="Mirror Plot",     icon="📊"),
    ],
}

pg = st.navigation(pages)
pg.run()

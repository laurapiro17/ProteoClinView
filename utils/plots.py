"""
Plotly figure builders for ProteoClinView.
"""
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import Optional

from utils.fragments import ION_COLORS

ION_NAMES = {"c": "c-ions (N-terminal)", "z": "z-ions (C-terminal)",
             "b": "b-ions", "y": "y-ions"}

# Coverage colors
COV_NONE  = "#F3F4F6"   # gray-100
COV_C     = "#BFDBFE"   # blue-200
COV_Z     = "#FED7AA"   # orange-200
COV_BOTH  = "#DDD6FE"   # violet-200
COV_TEXT  = {"none": "#6B7280", "c": "#1D4ED8", "z": "#C2410C", "both": "#5B21B6"}


def _stems(mzs, intensities, color, width=1.5, dash="solid"):
    """Build x/y arrays for a stem plot (list of [mz, mz, None] segments)."""
    x, y = [], []
    for m, i in zip(mzs, intensities):
        x += [m, m, None]
        y += [0, i, None]
    return go.Scatter(x=x, y=y, mode="lines",
                      line=dict(color=color, width=width, dash=dash),
                      showlegend=False, hoverinfo="skip")


# ─────────────────────────────────────────────────────────────────────────────
# Sequence coverage map
# ─────────────────────────────────────────────────────────────────────────────

def make_sequence_map(
    sequence: str,
    fragments_df: pd.DataFrame,
    residues_per_row: int = 10,
    highlight_label: Optional[str] = None,
) -> go.Figure:
    """
    ProSight Lite–style sequence coverage figure.
    Each residue is a coloured box. Cleavage ticks shown between residues.
    highlight_label: an ion_label (e.g. 'c12') to highlight with a border.
    """
    seq = sequence.upper()
    n = len(seq)

    # Build per-residue coverage
    c_pos, z_pos = set(), set()
    c_cleave, z_cleave = set(), set()   # cleavage site index (between residue i and i+1)

    for _, row in fragments_df.iterrows():
        pf, pt = int(row["position_from"]), int(row["position_to"])
        if row["ion_type"] == "c":
            for p in range(pf - 1, pt):
                c_pos.add(p)
            c_cleave.add(pt)          # cleavage after residue pt (1-based)
        elif row["ion_type"] == "z":
            for p in range(pf - 1, pt):
                z_pos.add(p)
            z_cleave.add(pf - 1)      # cleavage before residue pf (0-based)

    # Highlight positions for a selected ion
    highlight_pos = set()
    if highlight_label and not fragments_df.empty:
        sel = fragments_df[fragments_df["ion_label"] == highlight_label]
        if not sel.empty:
            pf, pt = int(sel.iloc[0]["position_from"]), int(sel.iloc[0]["position_to"])
            highlight_pos = set(range(pf - 1, pt))

    n_rows = (n + residues_per_row - 1) // residues_per_row
    BOX_W, BOX_H = 0.88, 0.80
    ROW_GAP = 2.2

    shapes, annotations = [], []

    for i, aa in enumerate(seq):
        col = i % residues_per_row
        row = i // residues_per_row
        x = col
        y = -row * ROW_GAP

        in_c, in_z = i in c_pos, i in z_pos
        if in_c and in_z:
            fill, tc = COV_BOTH, COV_TEXT["both"]
        elif in_c:
            fill, tc = COV_C, COV_TEXT["c"]
        elif in_z:
            fill, tc = COV_Z, COV_TEXT["z"]
        else:
            fill, tc = COV_NONE, COV_TEXT["none"]

        border_color = "#F59E0B" if i in highlight_pos else "#D1D5DB"
        border_width = 3 if i in highlight_pos else 1

        shapes.append(dict(
            type="rect",
            x0=x - BOX_W / 2, y0=y - BOX_H / 2,
            x1=x + BOX_W / 2, y1=y + BOX_H / 2,
            fillcolor=fill,
            line=dict(color=border_color, width=border_width),
            layer="below",
        ))

        annotations.append(dict(
            x=x, y=y, text=f"<b>{aa}</b>", showarrow=False,
            font=dict(size=13, color=tc, family="monospace"),
            xanchor="center", yanchor="middle",
        ))

        # Position labels every 5 residues
        if (i + 1) % 5 == 0 or i == 0:
            annotations.append(dict(
                x=x, y=y - BOX_H / 2 - 0.28,
                text=str(i + 1), showarrow=False,
                font=dict(size=8, color="#9CA3AF"),
                xanchor="center", yanchor="top",
            ))

        # Cleavage ticks
        next_col = (i + 1) % residues_per_row
        if (i + 1) in c_cleave and i < n - 1 and next_col != 0:
            shapes.append(dict(type="line",
                x0=x + BOX_W / 2, y0=y + BOX_H / 2 + 0.05,
                x1=x + BOX_W / 2, y1=y + BOX_H / 2 + 0.40,
                line=dict(color=ION_COLORS["c"], width=2.5)))

        if i in z_cleave and col > 0:
            shapes.append(dict(type="line",
                x0=x - BOX_W / 2, y0=y - BOX_H / 2 - 0.05,
                x1=x - BOX_W / 2, y1=y - BOX_H / 2 - 0.40,
                line=dict(color=ION_COLORS["z"], width=2.5)))

    fig = go.Figure()
    # Invisible hover scatter
    fig.add_trace(go.Scatter(
        x=[i % residues_per_row for i in range(n)],
        y=[-(i // residues_per_row) * ROW_GAP for i in range(n)],
        mode="markers",
        marker=dict(size=32, opacity=0),
        text=[f"Position {i+1} · {seq[i]}" + (
            f"<br>Covered by c-ions" if i in c_pos else "") + (
            f"<br>Covered by z-ions" if i in z_pos else "") for i in range(n)],
        hovertemplate="%{text}<extra></extra>",
        showlegend=False,
    ))

    fig.update_layout(
        shapes=shapes, annotations=annotations,
        xaxis=dict(range=[-0.65, residues_per_row - 0.35], visible=False, fixedrange=True),
        yaxis=dict(range=[-(n_rows - 1) * ROW_GAP - 1.1, 1.1], visible=False, fixedrange=True),
        height=max(120, n_rows * 95 + 50),
        margin=dict(l=8, r=8, t=8, b=8),
        plot_bgcolor="white", paper_bgcolor="white",
        showlegend=False, hovermode="closest",
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Theoretical spectrum (stem plot)
# ─────────────────────────────────────────────────────────────────────────────

def make_theoretical_spectrum(
    fragments_df: pd.DataFrame,
    title: str = "Theoretical Fragment Ions",
) -> go.Figure:
    fig = go.Figure()

    for ion_type, group in fragments_df.groupby("ion_type"):
        color = ION_COLORS.get(ion_type, "#9CA3AF")
        fig.add_trace(_stems(group["mz"], [1.0] * len(group), color))
        fig.add_trace(go.Scatter(
            x=group["mz"], y=[1.0] * len(group),
            mode="markers", marker=dict(size=6, color=color),
            name=ION_NAMES.get(ion_type, ion_type),
            text=group.apply(
                lambda r: f"<b>{r['ion_label']}</b> (z={r['charge']}+)<br>"
                          f"m/z = {r['mz']:.4f}<br>"
                          f"covers: {r['sequence_covered']}", axis=1),
            hovertemplate="%{text}<extra></extra>",
        ))

    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title="m/z",
        yaxis=dict(range=[0, 1.35], tickformat=".0%", title="Relative Intensity"),
        height=300, margin=dict(l=55, r=15, t=45, b=50),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="top", y=1.18, xanchor="right", x=1, font=dict(size=11)),
    )
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Mirror plot
# ─────────────────────────────────────────────────────────────────────────────

def make_mirror_plot(
    theoretical_df: pd.DataFrame,
    experimental_df: pd.DataFrame,
    matched_df: pd.DataFrame,
    title: str = "Mirror Plot — Experimental vs Theoretical",
    tolerance_ppm: float = 10.0,
) -> go.Figure:
    """
    Top half  (+y): experimental spectrum, coloured where matched.
    Bottom half (−y): theoretical fragment ions, mirrored downward.
    A horizontal zero-line separates the two halves.
    """
    fig = go.Figure()

    matched_labels = set(matched_df.loc[matched_df["matched"], "ion_label"])

    # ── Experimental ────────────────────────────────────────────────────────
    noise = experimental_df[experimental_df.get("is_noise", experimental_df["ion_type"] == "noise")]
    signal = experimental_df[~experimental_df.index.isin(noise.index)]

    # Noise (light gray, thin)
    if not noise.empty:
        fig.add_trace(_stems(noise["mz"], noise.get("intensity", pd.Series(0.2, index=noise.index)) * 0.5,
                             "#E5E7EB", width=0.8))

    # Unmatched signal (mid gray)
    matched_exp_mzs = set(matched_df.loc[matched_df["matched"], "exp_mz"].dropna())
    unmatched = signal[~signal["mz"].round(4).isin({round(m, 4) for m in matched_exp_mzs})]
    if not unmatched.empty:
        fig.add_trace(_stems(unmatched["mz"], unmatched.get("intensity", pd.Series(0.6, index=unmatched.index)),
                             "#9CA3AF", width=1.2))

    # Matched signal (coloured by ion type)
    for ion_type, group in matched_df[matched_df["matched"]].groupby("ion_type"):
        color = ION_COLORS.get(ion_type, "#9CA3AF")
        exp_mzs  = group["exp_mz"].values
        exp_ints = group["exp_intensity"].fillna(0.8).values
        fig.add_trace(_stems(exp_mzs, exp_ints, color, width=2.2))
        fig.add_trace(go.Scatter(
            x=exp_mzs, y=exp_ints,
            mode="markers", marker=dict(size=8, color=color),
            name=f"Matched {ION_NAMES.get(ion_type, ion_type)}",
            text=[
                f"<b>{row['ion_label']}</b> matched<br>"
                f"exp m/z = {row['exp_mz']:.4f}<br>"
                f"Δ = {(row['exp_mz'] - row['mz']) / row['mz'] * 1e6:.1f} ppm<br>"
                f"covers: {row['sequence_covered']}"
                for _, row in group.iterrows()
            ],
            hovertemplate="%{text}<extra></extra>",
        ))

    # ── Theoretical (mirrored, negative y) ──────────────────────────────────
    for ion_type, group in theoretical_df.groupby("ion_type"):
        color = ION_COLORS.get(ion_type, "#9CA3AF")
        is_matched = group["ion_label"].isin(matched_labels)

        # Matched theoretical: full intensity, solid
        m_group = group[is_matched]
        if not m_group.empty:
            fig.add_trace(_stems(m_group["mz"], [-1.0] * len(m_group), color, width=2.0))
            fig.add_trace(go.Scatter(
                x=m_group["mz"], y=[-1.0] * len(m_group),
                mode="markers", marker=dict(size=7, color=color, symbol="triangle-down"),
                showlegend=False,
                text=[f"<b>{r['ion_label']}</b> (theoretical, matched)<br>m/z = {r['mz']:.4f}"
                      for _, r in m_group.iterrows()],
                hovertemplate="%{text}<extra></extra>",
            ))

        # Unmatched theoretical: dimmer, dashed
        u_group = group[~is_matched]
        if not u_group.empty:
            fig.add_trace(_stems(u_group["mz"], [-0.5] * len(u_group),
                                 color, width=1.0, dash="dot"))
            fig.add_trace(go.Scatter(
                x=u_group["mz"], y=[-0.5] * len(u_group),
                mode="markers", marker=dict(size=5, color=color, symbol="triangle-down", opacity=0.4),
                showlegend=False,
                text=[f"<b>{r['ion_label']}</b> (theoretical, unmatched)<br>m/z = {r['mz']:.4f}"
                      for _, r in u_group.iterrows()],
                hovertemplate="%{text}<extra></extra>",
            ))

    # ── Layout ──────────────────────────────────────────────────────────────
    fig.add_hline(y=0, line=dict(color="#374151", width=1.5))

    fig.add_annotation(x=0.01, y=0.97, xref="paper", yref="paper",
                       text="<b>Experimental</b>", showarrow=False,
                       font=dict(size=12, color="#374151"), xanchor="left")
    fig.add_annotation(x=0.01, y=0.03, xref="paper", yref="paper",
                       text="<b>Theoretical</b>", showarrow=False,
                       font=dict(size=12, color="#374151"), xanchor="left")

    fig.update_layout(
        title=dict(text=title, font=dict(size=14)),
        xaxis_title="m/z",
        yaxis=dict(
            range=[-1.35, 1.35],
            tickvals=[-1.0, -0.5, 0, 0.5, 1.0],
            ticktext=["100%", "50%", "0%", "50%", "100%"],
            title="Relative Intensity",
        ),
        height=520,
        margin=dict(l=65, r=20, t=55, b=70),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=-0.22,
                    xanchor="center", x=0.5, font=dict(size=10)),
        hovermode="closest",
    )
    return fig

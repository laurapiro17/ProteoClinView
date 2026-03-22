"""
Fragment ion calculation for top-down proteomics.
Uses monoisotopic residue masses — no pyOpenMS dependency for core math.
"""
import pandas as pd
import numpy as np

# Monoisotopic residue masses (Da)
AA_MASSES = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
    'C': 103.00919, 'E': 129.04259, 'Q': 128.05858, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
    'M': 131.04049, 'F': 147.06841, 'P': 97.05276,  'S': 87.03203,
    'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841,
}

PROTON   = 1.007276
H2O      = 18.010565
NH3      = 17.026549
NH2_RAD  = 16.018724   # NH2 radical (z-ion correction)

ION_COLORS = {
    "c": "#3B82F6",   # blue
    "z": "#F97316",   # orange
    "b": "#8B5CF6",   # purple
    "y": "#10B981",   # green
    "noise": "#D1D5DB",
}

DEMO_PROTEINS = {
    "Human Insulin B chain": {
        "sequence": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
        "description": "Insulin B chain (30 AA) — regulates blood glucose. Central to diabetes research.",
        "uniprot": "P01308",
        "clinical_note": "Mutations in insulin are causative in neonatal diabetes and MODY.",
    },
    "Human Ubiquitin": {
        "sequence": "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        "description": "Ubiquitin (76 AA) — protein degradation tag. Highly conserved.",
        "uniprot": "P0CG48",
        "clinical_note": "Ubiquitin pathway dysregulation is implicated in neurodegeneration and cancer.",
    },
    "Hemoglobin α N-terminal": {
        "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
        "description": "Hemoglobin alpha subunit N-terminal (30 AA) — oxygen transport.",
        "uniprot": "P69905",
        "clinical_note": "Top-down proteomics detects sickle cell (E6V) and thalassemia variants directly.",
    },
}


def calculate_fragments(
    sequence: str,
    ion_types: list = ("c", "z"),
    charges: list = (1,),
) -> pd.DataFrame:
    """
    Calculate theoretical fragment ions (c, z, b, y) for a protein/peptide sequence.

    Parameters
    ----------
    sequence : str
        Amino acid sequence (single-letter codes, uppercase).
    ion_types : list of str
        Which ion series to compute. Any subset of ['c', 'z', 'b', 'y'].
    charges : list of int
        Charge states to compute.

    Returns
    -------
    pd.DataFrame with columns:
        ion_type, ion_label, position_from, position_to,
        sequence_covered, mz, charge, neutral_mass, color
    """
    seq = sequence.upper().replace(" ", "")
    n = len(seq)

    bad = [aa for aa in seq if aa not in AA_MASSES]
    if bad:
        raise ValueError(f"Unknown amino acid(s): {sorted(set(bad))}")

    masses = [AA_MASSES[aa] for aa in seq]
    cs = np.cumsum(masses)  # cs[i] = sum of residues 0..i

    rows = []

    for z in charges:
        if "c" in ion_types:
            for i in range(1, n):
                # c[i]: N-terminal fragment of length i, ends with -NH2
                # neutral mass = sum(residues[0:i]) + H (N-term) + NH3
                neutral = cs[i - 1] + 1.007825 + NH3
                mz = (neutral + z * PROTON) / z
                rows.append({
                    "ion_type": "c",
                    "ion_label": f"c{i}",
                    "position_from": 1,
                    "position_to": i,
                    "sequence_covered": seq[:i],
                    "mz": round(mz, 5),
                    "charge": z,
                    "neutral_mass": round(neutral, 5),
                    "color": ION_COLORS["c"],
                })

        if "z" in ion_types:
            for i in range(1, n):
                # z•[i]: C-terminal radical fragment of length i
                # neutral mass = sum(residues[n-i:n]) + OH (C-term) - NH2•
                suffix_mass = cs[n - 1] - (cs[n - i - 1] if i < n else 0)
                neutral = suffix_mass + 17.002740 - NH2_RAD  # OH - NH2•
                mz = (neutral + z * PROTON) / z
                rows.append({
                    "ion_type": "z",
                    "ion_label": f"z{i}",
                    "position_from": n - i + 1,
                    "position_to": n,
                    "sequence_covered": seq[n - i:],
                    "mz": round(mz, 5),
                    "charge": z,
                    "neutral_mass": round(neutral, 5),
                    "color": ION_COLORS["z"],
                })

        if "b" in ion_types:
            for i in range(1, n):
                # b[i] = sum(residues[0:i]) + H (N-term)
                neutral = cs[i - 1] + 1.007825
                mz = (neutral + z * PROTON) / z
                rows.append({
                    "ion_type": "b",
                    "ion_label": f"b{i}",
                    "position_from": 1,
                    "position_to": i,
                    "sequence_covered": seq[:i],
                    "mz": round(mz, 5),
                    "charge": z,
                    "neutral_mass": round(neutral, 5),
                    "color": ION_COLORS["b"],
                })

        if "y" in ion_types:
            for i in range(1, n):
                # y[i] = sum(residues[n-i:n]) + H2O + H (N-term)
                suffix_mass = cs[n - 1] - (cs[n - i - 1] if i < n else 0)
                neutral = suffix_mass + H2O + 1.007825
                mz = (neutral + z * PROTON) / z
                rows.append({
                    "ion_type": "y",
                    "ion_label": f"y{i}",
                    "position_from": n - i + 1,
                    "position_to": n,
                    "sequence_covered": seq[n - i:],
                    "mz": round(mz, 5),
                    "charge": z,
                    "neutral_mass": round(neutral, 5),
                    "color": ION_COLORS["y"],
                })

    return pd.DataFrame(rows)


def generate_demo_spectrum(
    fragments_df: pd.DataFrame,
    missing_fraction: float = 0.25,
    noise_fraction: float = 0.30,
    mz_error_ppm: float = 4.0,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Simulate a realistic experimental spectrum from theoretical fragments.
    - Randomly drops ~25% of fragments (incomplete coverage is expected)
    - Adds sub-5 ppm m/z measurement error
    - Adds low-intensity noise peaks
    """
    rng = np.random.default_rng(seed)

    keep = rng.random(len(fragments_df)) > missing_fraction
    exp = fragments_df[keep].copy()

    # m/z measurement error
    errors = rng.normal(0, mz_error_ppm / 3, len(exp))
    exp["mz"] = exp["mz"] * (1 + errors * 1e-6)

    # Intensities: exponential distribution (most peaks moderate, few dominant)
    base = rng.exponential(0.5, len(exp)) + 0.1
    exp["intensity"] = base / base.max()
    exp["is_noise"] = False

    # Noise peaks
    n_noise = max(1, int(len(fragments_df) * noise_fraction))
    mz_min = fragments_df["mz"].min() * 0.9 if len(fragments_df) > 0 else 100
    mz_max = fragments_df["mz"].max() * 1.1 if len(fragments_df) > 0 else 2000

    noise = pd.DataFrame({
        "ion_type": ["noise"] * n_noise,
        "ion_label": [f"noise_{i}" for i in range(n_noise)],
        "position_from": 0,
        "position_to": 0,
        "sequence_covered": "",
        "mz": rng.uniform(mz_min, mz_max, n_noise),
        "charge": 1,
        "neutral_mass": 0.0,
        "color": ION_COLORS["noise"],
        "intensity": rng.exponential(0.08, n_noise),
        "is_noise": True,
    })

    result = pd.concat([exp, noise], ignore_index=True)
    return result.sort_values("mz").reset_index(drop=True)


def match_peaks(
    theoretical: pd.DataFrame,
    experimental: pd.DataFrame,
    tolerance_ppm: float = 10.0,
) -> pd.DataFrame:
    """
    Match theoretical fragment ions to experimental peaks within tolerance_ppm.
    Returns theoretical df enriched with 'matched', 'exp_mz', 'exp_intensity'.
    """
    result = theoretical.copy()
    result["matched"] = False
    result["exp_mz"] = np.nan
    result["exp_intensity"] = np.nan

    exp_signal = experimental[~experimental.get("is_noise", pd.Series(False, index=experimental.index))]

    for idx, row in result.iterrows():
        tol = row["mz"] * tolerance_ppm * 1e-6
        candidates = exp_signal[
            (exp_signal["mz"] >= row["mz"] - tol) &
            (exp_signal["mz"] <= row["mz"] + tol)
        ]
        if not candidates.empty:
            best = candidates.loc[(candidates["mz"] - row["mz"]).abs().idxmin()]
            result.at[idx, "matched"] = True
            result.at[idx, "exp_mz"] = best["mz"]
            result.at[idx, "exp_intensity"] = best.get("intensity", 1.0)

    return result


def sequence_coverage_pct(sequence: str, matched_df: pd.DataFrame) -> float:
    """Fraction of residues covered by at least one matched fragment ion."""
    n = len(sequence)
    covered = set()
    for _, row in matched_df[matched_df["matched"]].iterrows():
        for pos in range(row["position_from"] - 1, row["position_to"]):
            covered.add(pos)
    return len(covered) / n if n > 0 else 0.0

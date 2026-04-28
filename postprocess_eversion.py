#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
postprocess_eversion.py

Post-processing for Drosophila wing disc eversion simulation.
Reads VTK output frames + simulation log, produces three figures:

  1. Cross-section plots (along-DV X-axis cut, across-DV Y-axis cut)
     showing apical and basal surface profiles.
  2. Curvature profiles along each cross-section, separately for
     apical and basal surfaces.
  3. Volume trajectory: simulated V_curr and V_target from the log,
     overlaid with biological Liu/Carthew 2024 measurements.

Usage:
    python3 postprocess_eversion.py \
        --vtk-dir /path/to/vtk/frames \
        --log /path/to/run.log \
        --num-layers 5 \
        --out-dir /path/to/figures

Notes
-----
Mesh ordering assumption: XMLParser adds apical nodes first, then body
layers descending, then basal last. So:
    apical = node indices [0,         N_per_layer)
    basal  = node indices [N - N_per_layer, N)
where N_per_layer = total_nodes / num_layers.

Cross-section logic: along-DV uses |Y| < band_y to pick nodes near
the DV midline; across-DV uses |X| < band_x to pick nodes near X=0.
The bands are tunable via CLI flags.

Curvature is computed by fitting a cubic spline to (s, Z) where s is
arc length along the cross-section, then evaluating
    kappa(s) = |Z''(s)| / (1 + Z'(s)**2)**1.5

Volume plots parse the simulation log for the "V_target = ... | V_curr = ..."
pattern emitted in solveSystem(). The log MUST be the stdout from a run
that uses the 3-stage volume-ramp loop.

Sources for biological volumes:
    Liu, O'Connell, Wall & Carthew (2024) eLife 12:RP91572
    DOI 10.7554/eLife.91572.3
    Animal_Properties.xlsx, Pouch Measurements sheet, wildtype rows.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
# Force a non-interactive backend BEFORE importing pyplot. On headless
# clusters (no DISPLAY) matplotlib otherwise tries to load Tk and may hit
# a version mismatch in the Tk backend; Agg is pure-software, writes PNGs,
# and works on every system.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# VTK parsing
# ---------------------------------------------------------------------------

def parse_vtk(path):
    """Return (points Nx3 ndarray, edges Mx2 ndarray, strain M-array or None).

    Handles the legacy ASCII VTK format emitted by Storage::print_VTK_File:
        DATASET UNSTRUCTURED_GRID
        POINTS N float
        ...
        CELLS M total
        2 i j   (VTK_LINE for edges)
        ...
        CELL_TYPES M
        3
        ...
        CELL_DATA M
        SCALARS Strain double 1
        LOOKUP_TABLE default
        ...
    """
    with open(path) as f:
        lines = f.readlines()

    # Find POINTS section
    n_points = None
    points_start = None
    for i, line in enumerate(lines):
        if line.startswith('POINTS'):
            n_points = int(line.split()[1])
            points_start = i + 1
            break
    if n_points is None:
        raise ValueError(f"No POINTS section found in {path}")

    pts = np.zeros((n_points, 3))
    for j in range(n_points):
        toks = lines[points_start + j].split()
        pts[j] = [float(toks[0]), float(toks[1]), float(toks[2])]

    # Find CELLS section, parse only edges (size==2)
    edges = []
    cells_start = None
    n_cells = None
    for i, line in enumerate(lines):
        if line.startswith('CELLS'):
            n_cells = int(line.split()[1])
            cells_start = i + 1
            break
    if cells_start is not None:
        for j in range(n_cells):
            toks = lines[cells_start + j].split()
            if len(toks) >= 3 and toks[0] == '2':
                edges.append([int(toks[1]), int(toks[2])])
    edges = np.array(edges, dtype=int) if edges else np.zeros((0, 2), dtype=int)

    # Find SCALARS Strain (per-cell)
    strain = None
    for i, line in enumerate(lines):
        if line.startswith('SCALARS') and 'Strain' in line:
            # next line is LOOKUP_TABLE; data starts after
            data_start = i + 2
            strain = np.array([float(lines[data_start + k].split()[0])
                               for k in range(n_cells)])
            break

    return pts, edges, strain


def find_frames(vtk_dir, pattern=r'_(\d+)\.vtk$', name_prefix=None):
    """Return list of (frame_index, path) sorted by frame index.

    If name_prefix is provided, only files whose name starts with that
    prefix are included. Useful when the directory contains output from
    multiple runs or simulations with different mesh sizes.
    """
    rx = re.compile(pattern)
    out = []
    for p in sorted(Path(vtk_dir).iterdir()):
        if p.suffix != '.vtk':
            continue
        if name_prefix and not p.name.startswith(name_prefix):
            continue
        m = rx.search(p.name)
        if m:
            out.append((int(m.group(1)), p))
    out.sort(key=lambda t: t[0])
    return out


def filter_frames_by_node_count(frames, target_n=None, num_layers=None):
    """Filter `frames` to only those whose VTK has the right node count.

    If target_n is given, keeps only frames matching that count exactly
    (with optional 1-point dummy trim).

    If target_n is None, infers it from the *most common* node count
    among the frames (after trimming) -- robust to a few broken frames.

    Returns (filtered_frames, target_n, removed_node_counts).
    """
    counts = {}  # node_count -> [(fi, path), ...]
    for fi, path in frames:
        try:
            with open(path) as f:
                # Just read until POINTS line to avoid full parse
                for line in f:
                    if line.startswith('POINTS'):
                        n = int(line.split()[1])
                        # Pick best candidate after dummy trim
                        if num_layers:
                            for trim in (0, 1):
                                if (n - trim) % num_layers == 0 and (n - trim) > 0:
                                    counts.setdefault(n - trim, []).append((fi, path))
                                    break
                            else:
                                counts.setdefault(n, []).append((fi, path))
                        else:
                            counts.setdefault(n, []).append((fi, path))
                        break
        except Exception:
            continue

    if not counts:
        return frames, None, set()

    if target_n is None:
        # Infer the most common count
        target_n = max(counts.keys(), key=lambda k: len(counts[k]))

    kept = sorted(counts.get(target_n, []), key=lambda t: t[0])
    removed = set(k for k in counts if k != target_n)
    return kept, target_n, removed


# ---------------------------------------------------------------------------
# Layer extraction
# ---------------------------------------------------------------------------

def split_by_layer(pts, num_layers):
    """Return list of point-index arrays, [layer_max, ..., layer_0].

    Assumes XMLParser ordering: apical first, layers descend.
    apical = layer (num_layers - 1)
    basal  = layer 0

    The simulation's nodes_in_upperhem flag stores the layer number,
    but we don't have that in the VTK so we rely on the index-based
    correspondence guaranteed by BuildPrismsFromLayerFlags.
    """
    n = len(pts)
    if n % num_layers != 0:
        raise ValueError(f"Total nodes {n} not divisible by num_layers {num_layers}. "
                         f"VTK may include a trailing dummy point; pass --vtk-trims-dummy 1 "
                         f"or check the mesh.")
    per = n // num_layers
    layers = []
    for L in range(num_layers):
        # L=0 is apical (top of XML), L=num_layers-1 is basal (bottom of XML)
        idx = np.arange(L * per, (L + 1) * per)
        layers.append(idx)
    return layers  # layers[0] = apical, layers[-1] = basal


# ---------------------------------------------------------------------------
# Cross-section extraction & plotting
# ---------------------------------------------------------------------------

def cross_section_along_axis(pts, indices, axis, band, perp_axis):
    """Return (s, Z) profile of points in `indices` whose `perp_axis` value
    lies within +/- band, sorted by `axis` value.

    `axis`: 0 (X) or 1 (Y) -- the axis we're cutting along.
    `perp_axis`: the axis we filter on (the orthogonal one).
    Returns axis_values, Z_values, both sorted by axis.
    """
    sub = pts[indices]
    mask = np.abs(sub[:, perp_axis]) < band
    if mask.sum() < 3:
        return None, None
    sel = sub[mask]
    order = np.argsort(sel[:, axis])
    return sel[order, axis], sel[order, 2]


def plot_cross_sections(pts, layers, out_path, band_along=3.0, band_across=3.0,
                        start_pts=None, start_layers=None):
    """4-panel figure: along-DV apical/basal, across-DV apical/basal.

    If start_pts/start_layers are provided, the starting (frame 0) cross-
    sections are plotted as faint dashed lines behind the final frame so
    the deformation is immediately visible.
    """
    apical_idx = layers[0]
    basal_idx = layers[-1]

    fig, axes = plt.subplots(2, 1, figsize=(9, 8))

    has_start = start_pts is not None and start_layers is not None
    if has_start:
        s_apical_idx = start_layers[0]
        s_basal_idx = start_layers[-1]

    # ---- Along DV (X-axis cut, slice at |Y| < band_along) ----
    ax = axes[0]
    if has_start:
        s_xa, s_za = cross_section_along_axis(start_pts, s_apical_idx, axis=0,
                                              band=band_along, perp_axis=1)
        s_xb, s_zb = cross_section_along_axis(start_pts, s_basal_idx, axis=0,
                                              band=band_along, perp_axis=1)
        if s_xa is not None:
            ax.plot(s_xa, s_za, ':', color='tab:red', alpha=0.5,
                    linewidth=1.5, label='apical (start)')
        if s_xb is not None:
            ax.plot(s_xb, s_zb, ':', color='tab:blue', alpha=0.5,
                    linewidth=1.5, label='basal (start)')

    xa, za = cross_section_along_axis(pts, apical_idx, axis=0, band=band_along, perp_axis=1)
    xb, zb = cross_section_along_axis(pts, basal_idx, axis=0, band=band_along, perp_axis=1)
    if xa is not None:
        ax.plot(xa, za, 'o:', color='tab:red',
                label='apical (final)' if has_start else 'apical',
                markersize=4, linewidth=1.5)
    if xb is not None:
        ax.plot(xb, zb, 's:', color='tab:blue',
                label='basal (final)' if has_start else 'basal',
                markersize=4, linewidth=1.5)
    ax.set_xlabel('X (along DV midline, $\\mu$m)')
    ax.set_ylabel('Z ($\\mu$m)')
    ax.set_title(f'Cross-section along DV midline  (|Y| < {band_along})')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='datalim')

    # ---- Across DV (Y-axis cut, slice at |X| < band_across) ----
    ax = axes[1]
    if has_start:
        s_ya, s_za = cross_section_along_axis(start_pts, s_apical_idx, axis=1,
                                              band=band_across, perp_axis=0)
        s_yb, s_zb = cross_section_along_axis(start_pts, s_basal_idx, axis=1,
                                              band=band_across, perp_axis=0)
        if s_ya is not None:
            ax.plot(s_ya, s_za, ':', color='tab:red', alpha=0.5,
                    linewidth=1.5, label='apical (start)')
        if s_yb is not None:
            ax.plot(s_yb, s_zb, ':', color='tab:blue', alpha=0.5,
                    linewidth=1.5, label='basal (start)')

    ya, za = cross_section_along_axis(pts, apical_idx, axis=1, band=band_across, perp_axis=0)
    yb, zb = cross_section_along_axis(pts, basal_idx, axis=1, band=band_across, perp_axis=0)
    if ya is not None:
        ax.plot(ya, za, 'o-', color='tab:red',
                label='apical (final)' if has_start else 'apical',
                markersize=4, linewidth=1.5)
    if yb is not None:
        ax.plot(yb, zb, 's-', color='tab:blue',
                label='basal (final)' if has_start else 'basal',
                markersize=4, linewidth=1.5)
    ax.set_xlabel('Y (across DV midline, $\\mu$m)')
    ax.set_ylabel('Z ($\\mu$m)')
    ax.set_title(f'Cross-section across DV midline  (|X| < {band_across})')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='datalim')

    suptitle = 'Apical and basal surface cross-sections'
    if has_start:
        suptitle += '  (dashed = wL3 start, solid = final frame)'
    fig.suptitle(suptitle, fontsize=12)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Curvature
# ---------------------------------------------------------------------------

def _is_curve_degenerate(s, z, full_span=None, min_span_frac=0.3,
                         max_bin_share=0.30, n_bins=10):
    """Decide whether the (s, z) cross-section is well-behaved enough for
    a meaningful curvature plot.

    Returns (is_bad, reason_string). Reasons we consider the curve unsuitable:
      - Fewer than 5 points
      - All points clustered in a span shorter than min_span_frac of the
        full layer footprint (means this layer has collapsed and the
        cross-section isn't representative)
      - Points severely clustered: more than max_bin_share of points fall
        into a single bin out of n_bins evenly spaced bins along s
      - Multiple z-values at the same s (folded curve, not a function)

    full_span is the maximum span along this axis seen at any layer; we
    use it so a basal layer that shrank to 10um is flagged when the apical
    is at 100um.
    """
    s = np.asarray(s)
    z = np.asarray(z)
    n = len(s)
    if n < 5:
        return True, f"too few points (n={n})"
    span = s.max() - s.min()
    if span < 1e-6:
        return True, "zero span along cut axis"
    if full_span is not None and span < min_span_frac * full_span:
        return True, (f"span {span:.1f}um is less than {int(min_span_frac*100)}% "
                      f"of full layer span {full_span:.1f}um (layer collapsed)")
    # Bin uniformly and check the largest bin's share
    counts, _ = np.histogram(s, bins=n_bins, range=(s.min(), s.max()))
    if counts.max() / n > max_bin_share:
        return True, (f"clustered: {counts.max()}/{n} points in 1 of {n_bins} bins "
                      f"(>{int(max_bin_share*100)}% threshold)")
    # Check if values fold back on themselves (multi-valued in s)
    s_sorted = np.sort(s)
    duplicate_thresh = 0.005 * span
    if np.any(np.diff(s_sorted) < duplicate_thresh):
        # Many near-duplicates means the curve isn't a function of s
        n_dups = int(np.sum(np.diff(s_sorted) < duplicate_thresh))
        if n_dups > n // 5:
            return True, f"curve folds back on itself ({n_dups} near-duplicate s values)"
    return False, ""


def signed_curvature(s, z, smoothing=None):
    """Signed curvature of curve z(s) using local quadratic fits.

    For each interior point i, fit a quadratic z = a*s^2 + b*s + c to a
    moving window of W neighbors (default W = max(7, len(s)//4)). Then
        z'(s_i)  = 2*a*s_i + b
        z''(s_i) = 2*a
        kappa(s_i) = z''(s_i) / (1 + z'(s_i)^2)^(3/2)

    Sign convention: positive when curve is concave upward (smile-shape).
    A central upward bump in apical Z(X) gives NEGATIVE central curvature
    because the curve bends downward on each side of the peak.

    Local-fit is preferred over global spline here because:
      - Cross-sections often have <30 points clustered non-uniformly
      - Spline second derivatives blow up at endpoints when data is sparse
      - Local fits degrade gracefully when nodes are clustered
    """
    s = np.asarray(s, dtype=float)
    z = np.asarray(z, dtype=float)
    n = len(s)
    if n < 5:
        return s, np.zeros_like(s)

    # Window: at least 7, scale with N, capped to N
    W = min(n, max(7, n // 4))
    half = W // 2

    kappa = np.zeros(n)
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        # Pull data into window-local coordinates centered on s[i]
        sw = s[lo:hi] - s[i]
        zw = z[lo:hi]
        if len(sw) < 3:
            kappa[i] = 0.0
            continue
        # Quadratic fit: zw = a*sw^2 + b*sw + c
        coeffs = np.polyfit(sw, zw, 2)
        a, b, _ = coeffs
        zp = b
        zpp = 2 * a
        kappa[i] = zpp / (1.0 + zp ** 2) ** 1.5

    return s, kappa


def _load_fuhrmann_curve(path, devstage=None, genotype=None):
    """Load a Fuhrmann reference curvature curve from CSV.

    Expected columns:
      - arclength column: 'arclength' or 'arclength_scaled'
      - curvature column: 'curvature' or 'curvature_scaled'
      - optional: curvature_sd (uncertainty band)
      - optional: crosssection ('DV' / 'AP' / 'Along_DV' / 'Across_DV' / etc.)
      - optional: devstage ('wL3' / '0hAPF' / '2hAPF' / '4hAPF') -- if
        present, the file is filtered to rows matching the devstage arg
      - optional: genotype -- same idea, if multiple genotypes are present

    Returns a dict {crosssection_name: (arclength, curvature, curvature_sd)}
    or {None: (a, k, sd)} if no crosssection column exists. Returns None on
    failure (with a warning printed to stderr).
    """
    try:
        # Lazy import: pandas is only needed if user passes --fuhrmann-curve
        import pandas as pd
    except ImportError:
        print("  WARN: pandas not installed; --fuhrmann-curve requires "
              "'pip install --user pandas'")
        return None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"  WARN: failed to read {path}: {e}")
        return None

    # Filter on devstage / genotype if those columns exist and an arg was given
    if devstage and 'devstage' in df.columns:
        avail = sorted(df['devstage'].unique())
        df = df[df['devstage'] == devstage]
        if df.empty:
            print(f"  WARN: no rows in {path} match devstage={devstage}; "
                  f"available: {avail}")
            return None
    if genotype and 'genotype' in df.columns:
        avail = sorted(df['genotype'].unique())
        df = df[df['genotype'] == genotype]
        if df.empty:
            print(f"  WARN: no rows in {path} match genotype={genotype}; "
                  f"available: {avail}")
            return None

    # Resolve column names
    arc_col = None
    for c in ('arclength', 'arclength_scaled'):
        if c in df.columns:
            arc_col = c; break
    kappa_col = None
    for c in ('curvature', 'curvature_scaled'):
        if c in df.columns:
            kappa_col = c; break
    if arc_col is None or kappa_col is None:
        print(f"  WARN: {path} missing arclength/curvature columns; "
              f"got {list(df.columns)}")
        return None
    sd_col = 'curvature_sd' if 'curvature_sd' in df.columns else None
    cs_col = 'crosssection' if 'crosssection' in df.columns else None

    out = {}
    if cs_col:
        for cs_val, sub in df.groupby(cs_col):
            sub = sub.sort_values(arc_col)
            a = sub[arc_col].to_numpy()
            k = sub[kappa_col].to_numpy()
            sd = sub[sd_col].to_numpy() if sd_col else None
            out[str(cs_val)] = (a, k, sd)
    else:
        df = df.sort_values(arc_col)
        a = df[arc_col].to_numpy()
        k = df[kappa_col].to_numpy()
        sd = df[sd_col].to_numpy() if sd_col else None
        out[None] = (a, k, sd)
    return out


def _resolve_fuhrmann_curve_for_panel(fuhrmann_data, panel):
    """Pick the right Fuhrmann curve for a given panel ('along' or 'across').

    Matches case-insensitively. Common labels:
      - along DV (X axis cut)  -> 'Along_DV', 'DV', 'along_DV', 'along'
      - across DV (Y axis cut) -> 'Across_DV', 'AP', 'PD', 'across_DV', 'across'
    If no crosssection column exists (single-curve file), use it for both.
    Returns (arclength, curvature, sd) or None.
    """
    if fuhrmann_data is None:
        return None
    if None in fuhrmann_data:
        return fuhrmann_data[None]
    # case-insensitive lookup
    keys_lower = {str(k).lower(): k for k in fuhrmann_data.keys()}
    if panel == 'along':
        for cand in ('along_dv', 'dv', 'along'):
            if cand in keys_lower:
                return fuhrmann_data[keys_lower[cand]]
    elif panel == 'across':
        for cand in ('across_dv', 'ap', 'pd', 'across'):
            if cand in keys_lower:
                return fuhrmann_data[keys_lower[cand]]
    # Fall back: single key, use it
    if len(fuhrmann_data) == 1:
        return next(iter(fuhrmann_data.values()))
    return None


def plot_curvature(pts, layers, out_path, band_along=3.0, band_across=3.0,
                   fuhrmann_data=None, fuhrmann_label='Fuhrmann et al. 2024',
                   flip_sim_curvature=False):
    """4-panel figure: curvature for apical/basal x along/across DV.

    Each panel sets its own y-range based on a robust 2-98 percentile so
    apical and basal can show meaningful detail at their own scales.
    Apical and basal x-ranges are tied per row so they line up vertically.

    A panel is replaced with an explanatory message if the data is too
    clustered, sparse, or folded for a curvature plot to be meaningful.
    This is the typical situation for the basal layer when the apical
    has fanned out and the basal has stayed compact: trying to compute
    curvature on a tiny clustered point set produces visual noise that
    is worse than no plot at all.

    If flip_sim_curvature=True, the simulation curvature is negated before
    plotting. This is provided as a comparison aid against published data
    (e.g. Fuhrmann 2024) where coordinate-system conventions differ. When
    set, the figure is annotated so the flip is visible in the rendered
    output, not just in the code.
    """
    apical_idx = layers[0]
    basal_idx = layers[-1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Pre-compute the four cross-sections.
    along_apical = cross_section_along_axis(pts, apical_idx, axis=0,
                                            band=band_along, perp_axis=1)
    along_basal = cross_section_along_axis(pts, basal_idx, axis=0,
                                           band=band_along, perp_axis=1)
    across_apical = cross_section_along_axis(pts, apical_idx, axis=1,
                                             band=band_across, perp_axis=0)
    across_basal = cross_section_along_axis(pts, basal_idx, axis=1,
                                            band=band_across, perp_axis=0)

    # Layer-wide spans, used as the reference scale for "is this layer
    # too small to compare to its sibling?". Use ALL apical points (not
    # just the slice) since we want the layer footprint.
    apical_pts = pts[apical_idx]
    full_span_x = apical_pts[:, 0].max() - apical_pts[:, 0].min()
    full_span_y = apical_pts[:, 1].max() - apical_pts[:, 1].min()

    sim_sign = -1.0 if flip_sim_curvature else 1.0

    def safe_curve(xs_zs, full_span):
        xs, zs = xs_zs
        if xs is None or len(xs) < 5:
            return None, None, "no data"
        is_bad, reason = _is_curve_degenerate(xs, zs, full_span=full_span)
        if is_bad:
            return None, None, reason
        x_d, kappa = signed_curvature(xs, zs)
        return x_d, sim_sign * kappa, ""

    def shared_xrange(*pairs):
        mins, maxs = [], []
        for p in pairs:
            if p[0] is not None:
                mins.append(np.min(p[0])); maxs.append(np.max(p[0]))
        if not mins:
            return None
        return min(mins), max(maxs)

    def robust_ylim(k):
        if k is None or len(k) == 0:
            return None
        lo, hi = np.percentile(k, [2, 98])
        m = max(abs(lo), abs(hi), 1e-3)
        return -m * 1.2, m * 1.2

    def render_panel(ax, x_d, kappa, reason, color, title):
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        if reason:
            ax.text(0.5, 0.5,
                    'curvature plot skipped:\n' + reason,
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=10, color='gray', wrap=True,
                    bbox=dict(facecolor='whitesmoke', edgecolor='lightgray',
                              boxstyle='round,pad=0.5'))
            ax.set_xlabel('')
            ax.set_ylabel('')
            return
        ax.plot(x_d, kappa, color=color, linewidth=1.5)
        ax.axhline(0, color='gray', linewidth=0.5)

    # ---- Top row: along DV ----
    sa_x, sa_k, sa_reason = safe_curve(along_apical, full_span_x)
    sb_x, sb_k, sb_reason = safe_curve(along_basal, full_span_x)
    n_along_apical = 0 if along_apical[0] is None else len(along_apical[0])
    n_along_basal = 0 if along_basal[0] is None else len(along_basal[0])
    xr = shared_xrange(along_apical, along_basal)

    def overlay_fuhrmann(ax, panel):
        if fuhrmann_data is None:
            return
        curve = _resolve_fuhrmann_curve_for_panel(fuhrmann_data, panel)
        if curve is None:
            return
        a, k, sd = curve
        ax.plot(a, k, '-', color='black', linewidth=1.5,
                label=fuhrmann_label, zorder=3, alpha=0.8)
        if sd is not None:
            ax.fill_between(a, k - sd, k + sd, color='black', alpha=0.15, zorder=2)
        # Expand the y-range so the Fuhrmann curve (+SD band) is fully
        # visible alongside the simulation curve, with 10% padding.
        cur_lo, cur_hi = ax.get_ylim()
        f_lo = float(np.min(k - (sd if sd is not None else 0)))
        f_hi = float(np.max(k + (sd if sd is not None else 0)))
        new_lo = min(cur_lo, f_lo)
        new_hi = max(cur_hi, f_hi)
        pad = 0.1 * (new_hi - new_lo) if new_hi > new_lo else 1e-3
        ax.set_ylim(new_lo - pad, new_hi + pad)
        ax.legend(loc='best', fontsize=8)

    ax = axes[0, 0]
    render_panel(ax, sa_x, sa_k, sa_reason, 'tab:red',
                 f'Along DV  -  apical  (n={n_along_apical})')
    if not sa_reason:
        ax.set_xlabel('X ($\\mu$m)')
        ax.set_ylabel('Curvature $\\kappa$ (1/$\\mu$m)')
        if xr: ax.set_xlim(xr)
        yr = robust_ylim(sa_k)
        if yr: ax.set_ylim(yr)
        overlay_fuhrmann(ax, 'along')

    ax = axes[0, 1]
    render_panel(ax, sb_x, sb_k, sb_reason, 'tab:blue',
                 f'Along DV  -  basal  (n={n_along_basal})')
    if not sb_reason:
        ax.set_xlabel('X ($\\mu$m)')
        ax.set_ylabel('Curvature $\\kappa$ (1/$\\mu$m)')
        if xr: ax.set_xlim(xr)
        yr = robust_ylim(sb_k)
        if yr: ax.set_ylim(yr)

    # ---- Bottom row: across DV ----
    ca_x, ca_k, ca_reason = safe_curve(across_apical, full_span_y)
    cb_x, cb_k, cb_reason = safe_curve(across_basal, full_span_y)
    n_across_apical = 0 if across_apical[0] is None else len(across_apical[0])
    n_across_basal = 0 if across_basal[0] is None else len(across_basal[0])
    xr = shared_xrange(across_apical, across_basal)

    ax = axes[1, 0]
    render_panel(ax, ca_x, ca_k, ca_reason, 'tab:red',
                 f'Across DV  -  apical  (n={n_across_apical})')
    if not ca_reason:
        ax.set_xlabel('Y ($\\mu$m)')
        ax.set_ylabel('Curvature $\\kappa$ (1/$\\mu$m)')
        if xr: ax.set_xlim(xr)
        yr = robust_ylim(ca_k)
        if yr: ax.set_ylim(yr)
        overlay_fuhrmann(ax, 'across')

    ax = axes[1, 1]
    render_panel(ax, cb_x, cb_k, cb_reason, 'tab:blue',
                 f'Across DV  -  basal  (n={n_across_basal})')
    if not cb_reason:
        ax.set_xlabel('Y ($\\mu$m)')
        ax.set_ylabel('Curvature $\\kappa$ (1/$\\mu$m)')
        if xr: ax.set_xlim(xr)
        yr = robust_ylim(cb_k)
        if yr: ax.set_ylim(yr)

    suptitle = ('Surface curvature along and across DV  '
                '(positive = concave up; each panel y-axis at 2-98 percentile)')
    fig.suptitle(suptitle, fontsize=11)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Volume trajectory parsing
# ---------------------------------------------------------------------------

# Matches lines like:
#   Stage 0/3 step 0/20 | E = 3.28 | Linear E = 2.49 | V_target = 64211.4 | V_curr = 63323.8 | Mov = ... | Relax steps = 446
STAGE_LINE_RX = re.compile(
    r'Stage\s+(\d+)/(\d+)\s+step\s+(\d+)/(\d+).*?'
    r'V_target\s*=\s*([0-9.eE+\-]+).*?'
    r'V_curr\s*=\s*([0-9.eE+\-]+)',
    re.IGNORECASE)

# Initial volume line: "Initial volume: 62340.6"
INIT_VOL_RX = re.compile(r'Initial volume:\s*([0-9.eE+\-]+)')

# Sim-to-bio scaling: "sim-to-bio scaling factor = 0.087208"
SCALE_RX = re.compile(r'sim-to-bio scaling factor\s*=\s*([0-9.eE+\-]+)')


def parse_log_volumes(log_path):
    """Return dict with keys: step (1-D arr), v_target, v_curr, stage, scale, v_init."""
    rows = []
    scale = None
    v_init = None
    with open(log_path) as f:
        for line in f:
            m = STAGE_LINE_RX.search(line)
            if m:
                stage = int(m.group(1))
                stages_total = int(m.group(2))
                step = int(m.group(3))
                steps_per_stage = int(m.group(4))
                v_t = float(m.group(5))
                v_c = float(m.group(6))
                rows.append((stage, step, steps_per_stage, v_t, v_c))
                continue
            m = SCALE_RX.search(line)
            if m and scale is None:
                scale = float(m.group(1))
                continue
            m = INIT_VOL_RX.search(line)
            if m and v_init is None:
                v_init = float(m.group(1))
                continue
    if not rows:
        return None

    # ----------------------------------------------------------------
    # Auto-detect 1-indexed stage labels.
    #
    # The simulation's stage loop variable is 0-indexed (0,1,2 for a
    # 3-stage run), but a known older System.cu prints "Stage 1/3"
    # for the first stage instead of "Stage 0/3". If the smallest
    # observed stage label in the log equals 1 (and stages_total = 3),
    # we shift everything down by 1 so global_step starts at 0.
    # ----------------------------------------------------------------
    stages_seen = set(r[0] for r in rows)
    min_stage = min(stages_seen)
    if min_stage == 1:
        # 1-indexed log -- shift down
        rows = [(s - 1, st, sps, vt, vc) for (s, st, sps, vt, vc) in rows]
        # tag this so we can mention it on the plot
        offset_applied = True
    else:
        offset_applied = False

    # Now compute global step
    rows_with_global = [
        (st_idx * sps + step, st_idx, step, vt, vc)
        for (st_idx, step, sps, vt, vc) in rows
    ]
    arr = np.array(rows_with_global, dtype=float)
    return {
        'global_step': arr[:, 0].astype(int),
        'stage': arr[:, 1].astype(int),
        'step_in_stage': arr[:, 2].astype(int),
        'v_target': arr[:, 3],
        'v_curr': arr[:, 4],
        'scale': scale,
        'v_init': v_init,
        'offset_applied': offset_applied,
    }


def plot_volume(log_data, out_path):
    """2-panel: simulation units (top), biological um^3 (bottom).

    Bottom panel overlays Liu/Carthew 2024 measurements.
    """
    if log_data is None:
        print("  WARN: no volume data parsed from log; skipping volume plot")
        return None

    gs = log_data['global_step']
    vt = log_data['v_target']
    vc = log_data['v_curr']
    scale = log_data['scale']  # sim_per_bio (so bio = sim / scale)
    stage_changes = np.where(np.diff(log_data['stage']) > 0)[0] + 1

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Top: simulation units
    ax = axes[0]
    ax.plot(gs, vt, '-', color='gray', linewidth=2, label='target (Liu/Carthew interpolated)')
    ax.plot(gs, vc, 'o-', color='tab:red', markersize=4, label='simulated')
    for sc in stage_changes:
        ax.axvline(gs[sc] - 0.5, color='black', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.set_xlabel('Simulation substep')
    ax.set_ylabel('Volume (sim units)')
    ax.set_title('Simulated vs target volume (simulation units)')
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

    # Bottom: biological micrometers^3
    ax = axes[1]
    if scale is not None and scale > 0:
        # bio = sim / scale  (since scale = sim_per_bio = V_sim / V_bio)
        ax.plot(gs, vt / scale, '-', color='gray', linewidth=2, label='target')
        ax.plot(gs, vc / scale, 'o-', color='tab:red', markersize=4, label='simulated')

        # Overlay Liu/Carthew 2024 wildtype pouch measurements, mapped to
        # estimated substep based on Fuhrmann's 3-stage convention:
        #   stage 0 ends at step 20 -> 0hAPF
        #   stage 1 ends at step 40 -> 2hAPF
        #   stage 2 ends at step 60 -> 4hAPF
        # Liu/Carthew measured these timepoints (um^3, all wildtype):
        #   wL3        :  714849   (interpolated, -4hAPF)
        #   WPP (0hAPF):  1143897
        #   WPP+1      :  1117134  (~+1hAPF)
        # Linear interp gives ~+2hAPF at 1135478 and ~+4hAPF at 1117134.
        bio_points = [
            (0,  714849, 'wL3 (-4 hAPF)'),
            (20, 1143897, 'WPP (0 hAPF)'),
            (40, 1135478, '~2 hAPF'),
            (60, 1117134, '4 hAPF'),
        ]
        bx = [p[0] for p in bio_points]
        by = [p[1] for p in bio_points]
        ax.plot(bx, by, 'D', color='tab:green', markersize=10,
                markeredgecolor='black', label='Liu/Carthew 2024 measured')
        for x, y, lab in bio_points:
            ax.annotate(lab, (x, y), textcoords='offset points', xytext=(8, -8),
                        fontsize=8)

        for sc in stage_changes:
            ax.axvline(gs[sc] - 0.5, color='black', linewidth=0.5,
                       linestyle='--', alpha=0.5)
        ax.set_xlabel('Simulation substep')
        ax.set_ylabel('Volume ($\\mu$m$^3$)')
        ax.set_title('Volume converted to biological units, with Liu/Carthew 2024 overlay')
        ax.legend(loc='lower right')
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, 'No sim-to-bio scaling factor in log',
                ha='center', va='center', transform=ax.transAxes)

    title = 'Volume trajectory: predicted vs simulated'
    if log_data.get('offset_applied'):
        title += '  (note: log used 1-indexed stage labels, shifted to 0-indexed)'
    fig.suptitle(title, fontsize=12)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Cell height at the DV midline center
# ---------------------------------------------------------------------------

def cell_height_at_center(pts, layers, n_nearest=5):
    """Return apical-basal thickness at the disc center.

    Computes:
      - the n_nearest apical nodes to (X=0, Y=0), averaging their Z
      - the n_nearest basal nodes to (X=0, Y=0), averaging their Z
      - cell_height = mean(Z_apical_center) - mean(Z_basal_center)

    "Apical" is layer index 0 in our split (closest to top of XML, highest Z
    for a normal dome) and "basal" is layer index -1 (closest to bottom of
    XML, lowest Z).

    Returns (cell_height, z_apical_mean, z_basal_mean, n_apical_used,
    n_basal_used) so the caller can sanity-check.
    """
    apical_idx = layers[0]
    basal_idx = layers[-1]
    apical = pts[apical_idx]
    basal = pts[basal_idx]
    # Distance from origin in the XY plane only -- we want the center of
    # the disc, regardless of how high or low this layer sits in Z.
    d_apical = np.hypot(apical[:, 0], apical[:, 1])
    d_basal = np.hypot(basal[:, 0], basal[:, 1])
    # Pick nearest n_nearest from each layer
    n_a = min(n_nearest, len(apical))
    n_b = min(n_nearest, len(basal))
    sel_a = np.argsort(d_apical)[:n_a]
    sel_b = np.argsort(d_basal)[:n_b]
    z_a = apical[sel_a, 2].mean()
    z_b = basal[sel_b, 2].mean()
    return (z_a - z_b), z_a, z_b, n_a, n_b


def plot_cell_height(frames, num_layers, out_path, n_nearest=5,
                     log_data=None):
    """Read every VTK frame, compute cell height at center, plot vs frame.

    `frames` is a list of (frame_index, Path) pairs as returned by
    find_frames(). num_layers is needed to split each frame's points into
    layers. log_data, if provided, lets us mark stage boundaries on the
    same axis as the volume plot (substep numbers).

    Stages assumed at substeps 0, 20, 40, 60 (3-stage Fuhrmann schedule).
    """
    print(f"  Computing cell height across {len(frames)} frames...")
    rows = []
    skipped_count = 0
    skipped_node_counts = set()
    target_n = None  # node count to expect; locked in on first parseable frame
    for fi, path in frames:
        try:
            pts, _, _ = parse_vtk(path)
            n = len(pts)
            # First time: lock target_n to the first frame that's divisible
            # by num_layers (with optional 1-point dummy trim).
            trimmed = None
            for trim in (0, 1):
                if (n - trim) % num_layers == 0 and (n - trim) > 0:
                    trimmed = pts[:n - trim]
                    candidate_n = n - trim
                    break
            if trimmed is None:
                # No valid trim; skip with a count, not a per-frame WARN
                skipped_count += 1
                skipped_node_counts.add(n)
                continue
            # If we've already locked target_n, the candidate must match.
            if target_n is None:
                target_n = candidate_n
            elif candidate_n != target_n:
                skipped_count += 1
                skipped_node_counts.add(n)
                continue
            pts = trimmed
            layers = split_by_layer(pts, num_layers)
            h, z_a, z_b, n_a, n_b = cell_height_at_center(pts, layers, n_nearest)
            rows.append((fi, h, z_a, z_b))
        except Exception as e:
            skipped_count += 1
            continue

    if skipped_count > 0:
        print(f"    Skipped {skipped_count} frames whose node count did not "
              f"match the rest of the run (expected {target_n}; "
              f"saw {sorted(skipped_node_counts)}).")

    if not rows:
        print("  No frames successfully processed; cell-height plot skipped.")
        return None

    arr = np.array(rows, dtype=float)
    frame_idx = arr[:, 0].astype(int)
    height = arr[:, 1]
    z_apical = arr[:, 2]
    z_basal = arr[:, 3]

    # Print summary stats
    print(f"  Cell height summary (n_nearest={n_nearest} per layer):")
    print(f"    First frame {int(frame_idx[0])}: height = {height[0]:.3f} um")
    print(f"    Last frame  {int(frame_idx[-1])}: height = {height[-1]:.3f} um")
    if len(height) >= 2 and height[0] > 0:
        pct = 100.0 * (height[-1] - height[0]) / height[0]
        print(f"    Change: {pct:+.1f}%")

    # ---- Plot ----
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Panel 1: cell height vs frame
    ax = axes[0]
    ax.plot(frame_idx, height, 'o-', color='tab:purple', markersize=4,
            linewidth=1.5, label='cell height (Z_apical_center - Z_basal_center)')
    ax.axhline(height[0], color='gray', linestyle='--', linewidth=0.8,
               label=f'wL3 reference: {height[0]:.2f} um')

    # Mark Fuhrmann 3-stage boundaries if frame range covers them
    stage_boundaries = [(0, 'wL3'), (20, '0hAPF'), (40, '2hAPF'), (60, '4hAPF')]
    for sb_frame, sb_label in stage_boundaries:
        if frame_idx.min() <= sb_frame <= frame_idx.max():
            ax.axvline(sb_frame, color='black', linewidth=0.5,
                       linestyle=':', alpha=0.5)
            ax.annotate(sb_label, (sb_frame, ax.get_ylim()[1]),
                        textcoords='offset points', xytext=(2, -10),
                        fontsize=8, color='gray')

    ax.set_xlabel('Simulation frame')
    ax.set_ylabel('Cell height ($\\mu$m)')
    ax.set_title(f'Apical-basal cell height at DV midline center  '
                 f'(averaged over {n_nearest} nearest nodes per layer)')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Panel 2: apical and basal Z separately, so the user can see whether
    # height change comes from apical going up, basal going down, or both
    ax = axes[1]
    ax.plot(frame_idx, z_apical, 'o-', color='tab:red', markersize=4,
            linewidth=1.5, label='apical center (mean Z)')
    ax.plot(frame_idx, z_basal, 's-', color='tab:blue', markersize=4,
            linewidth=1.5, label='basal center (mean Z)')
    for sb_frame, _ in stage_boundaries:
        if frame_idx.min() <= sb_frame <= frame_idx.max():
            ax.axvline(sb_frame, color='black', linewidth=0.5,
                       linestyle=':', alpha=0.5)
    ax.set_xlabel('Simulation frame')
    ax.set_ylabel('Z position ($\\mu$m)')
    ax.set_title('Apical and basal center Z positions over time')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    fig.suptitle('Cell height (apical-basal thickness) at DV midline', fontsize=13)
    try:
        fig.tight_layout()
    except Exception:
        pass
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vtk-dir', type=str, help='Directory containing VTK frames')
    ap.add_argument('--vtk-prefix', type=str, default=None,
                    help='Only consider VTK files whose filename starts with '
                         'this prefix. Use this when --vtk-dir contains files '
                         'from multiple runs and you want to disambiguate.')
    ap.add_argument('--vtk-file', type=str,
                    help='Single VTK file to use as the final frame (alternative to --vtk-dir)')
    ap.add_argument('--log', type=str,
                    help='Simulation log file with V_target / V_curr lines')
    ap.add_argument('--num-layers', type=int, default=5,
                    help='Number of mesh layers (default 5; 12-layer res4 mesh uses 12)')
    ap.add_argument('--out-dir', type=str, default='./postprocess_out',
                    help='Output directory for figures')
    ap.add_argument('--band-along', type=float, default=3.0,
                    help='|Y| band for along-DV cross-section in micrometers (default 3.0)')
    ap.add_argument('--band-across', type=float, default=3.0,
                    help='|X| band for across-DV cross-section in micrometers (default 3.0)')
    ap.add_argument('--frame-index', type=int, default=-1,
                    help='Which frame index to plot from --vtk-dir (-1 = last; default -1)')
    ap.add_argument('--n-center-nodes', type=int, default=5,
                    help='How many nodes near (X=0,Y=0) to average per layer for cell-height (default 5)')
    ap.add_argument('--fuhrmann-curve', type=str, default=None,
                    help='Path to a Fuhrmann reference curvature CSV with columns '
                         "'arclength' (or 'arclength_scaled') and 'curvature' (or "
                         "'curvature_scaled'); optional 'curvature_sd', "
                         "'crosssection' (DV/AP/Along_DV/Across_DV), 'devstage', "
                         "'genotype'. When provided, overlays the reference curve "
                         "on the apical curvature panels.")
    ap.add_argument('--fuhrmann-stage', type=str, default='4hAPF',
                    help="Filter Fuhrmann CSV to one developmental stage "
                         "('wL3', '0hAPF', '2hAPF', '4hAPF'). Default '4hAPF' "
                         "to compare against the simulation endpoint.")
    ap.add_argument('--fuhrmann-genotype', type=str, default=None,
                    help='Filter Fuhrmann CSV to one genotype (default: '
                         'use whatever is in the file).')
    ap.add_argument('--flip-sim-curvature', action='store_true',
                    help='Negate the simulation curvature in the curvature '
                         'plot, to align coordinate-system conventions for '
                         'comparison with reference data (e.g. Fuhrmann 2024). '
                         'The sign flip is annotated on the rendered figure.')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_paths = []

    # ---- Load Fuhrmann reference curve, if provided ----
    fuhrmann_data = None
    if args.fuhrmann_curve:
        print(f"  Loading Fuhrmann reference curve from {args.fuhrmann_curve} "
              f"(stage={args.fuhrmann_stage}, genotype={args.fuhrmann_genotype}) ...")
        fuhrmann_data = _load_fuhrmann_curve(args.fuhrmann_curve,
                                             devstage=args.fuhrmann_stage,
                                             genotype=args.fuhrmann_genotype)
        if fuhrmann_data is not None:
            keys = list(fuhrmann_data.keys())
            print(f"    loaded {len(fuhrmann_data)} curve(s); cross-section keys: {keys}")

    # ---- Resolve final-frame VTK ----
    final_vtk = None
    all_frames = []  # full list, used for cell-height time series
    if args.vtk_file:
        final_vtk = Path(args.vtk_file)
    elif args.vtk_dir:
        all_frames = find_frames(args.vtk_dir, name_prefix=args.vtk_prefix)
        if not all_frames:
            extra = f" with prefix '{args.vtk_prefix}'" if args.vtk_prefix else ""
            print(f"  ERROR: no VTK frames matching _NNNNN.vtk in {args.vtk_dir}{extra}")
        else:
            n_total = len(all_frames)
            # Auto-detect the most common node count and filter to those frames.
            # This protects against directories that contain multiple runs with
            # different mesh sizes (e.g. 1505-node and 14537-node frames mixed).
            all_frames, target_n, removed_counts = filter_frames_by_node_count(
                all_frames, target_n=None, num_layers=args.num_layers)
            if removed_counts:
                print(f"  Filtered {n_total - len(all_frames)} frames with "
                      f"different node counts {sorted(removed_counts)}; "
                      f"kept {len(all_frames)} frames at node count {target_n}")
            if not all_frames:
                print(f"  ERROR: no frames remained after node-count filtering")
            else:
                idx = args.frame_index if args.frame_index >= 0 else len(all_frames) + args.frame_index
                idx = max(0, min(idx, len(all_frames) - 1))
                final_vtk = all_frames[idx][1]
                print(f"  Using frame {all_frames[idx][0]} ({final_vtk.name}) "
                      f"out of {len(all_frames)} frames")

    # ---- Cross-sections + curvature ----
    if final_vtk is not None and final_vtk.exists():
        print(f"  Parsing {final_vtk} ...")
        pts, edges, strain = parse_vtk(final_vtk)
        # The VTK output sometimes has a trailing dummy point; trim if so.
        n = len(pts)
        for trim in (0, 1):
            if (n - trim) % args.num_layers == 0:
                if trim:
                    print(f"  Trimming {trim} trailing dummy point(s) so "
                          f"{n - trim} % {args.num_layers} == 0")
                    pts = pts[:n - trim]
                break
        else:
            print(f"  WARN: {n} points not divisible by --num-layers={args.num_layers}; "
                  f"layer-split may be wrong")

        try:
            layers = split_by_layer(pts, args.num_layers)

            # Optional starting cross-section overlay: scan forward from
            # frame 0 to find the first frame whose node count matches the
            # final frame. Some early frames may have a different layout
            # (e.g. cumulative dumps from initialization phases) -- we skip
            # those and use the earliest *parseable* frame as the start.
            start_pts = None
            start_layers = None
            target_n = len(pts)  # already trimmed final-frame node count
            if all_frames and all_frames[0][0] != all_frames[idx][0]:
                attempted = 0
                skipped_n = []
                for cand_fi, cand_path in all_frames[:idx]:
                    attempted += 1
                    try:
                        c_pts, _, _ = parse_vtk(cand_path)
                        cn = len(c_pts)
                        # Try with and without trailing-dummy trim
                        for trim in (0, 1):
                            if cn - trim == target_n:
                                c_pts = c_pts[:cn - trim]
                                start_pts = c_pts
                                start_layers = split_by_layer(c_pts, args.num_layers)
                                break
                        if start_pts is not None:
                            print(f"  Overlaying starting frame {cand_fi} "
                                  f"({cand_path.name}) on cross-sections")
                            if attempted > 1:
                                print(f"    (skipped {attempted-1} earlier frames "
                                      f"with mismatched node counts: {skipped_n})")
                            break
                        skipped_n.append(cn)
                    except Exception as e:
                        skipped_n.append(f"err({e})")
                        continue
                if start_pts is None:
                    print(f"  WARN: could not find a parseable starting frame "
                          f"for overlay (target node count = {target_n}; "
                          f"early frames had: {skipped_n[:5]}{'...' if len(skipped_n)>5 else ''})")

            cs_path = out_dir / 'cross_sections.png'
            plot_cross_sections(pts, layers, cs_path,
                                band_along=args.band_along,
                                band_across=args.band_across,
                                start_pts=start_pts,
                                start_layers=start_layers)
            out_paths.append(cs_path)
            print(f"  Cross-sections -> {cs_path}")

            cv_path = out_dir / 'curvature.png'
            f_label = 'Fuhrmann et al. 2024'
            if args.fuhrmann_curve:
                f_label = f'Fuhrmann 2024 ({args.fuhrmann_stage}'
                if args.fuhrmann_genotype:
                    f_label += f', {args.fuhrmann_genotype}'
                f_label += ')'
            plot_curvature(pts, layers, cv_path,
                           band_along=args.band_along,
                           band_across=args.band_across,
                           fuhrmann_data=fuhrmann_data,
                           fuhrmann_label=f_label,
                           flip_sim_curvature=args.flip_sim_curvature)
            out_paths.append(cv_path)
            print(f"  Curvature -> {cv_path}")
        except Exception as e:
            print(f"  ERROR generating cross-section/curvature plots: {e}")
    else:
        print("  No VTK file specified or found; skipping cross-section/curvature plots.")

    # ---- Volume trajectory ----
    if args.log:
        print(f"  Parsing {args.log} ...")
        log_data = parse_log_volumes(args.log)
        if log_data is None:
            print(f"  WARN: no V_target/V_curr lines in {args.log}; skipping")
        else:
            print(f"  Found {len(log_data['global_step'])} substeps with volume data, "
                  f"scale={log_data['scale']}, v_init={log_data['v_init']}")
            v_path = out_dir / 'volume_trajectory.png'
            plot_volume(log_data, v_path)
            out_paths.append(v_path)
            print(f"  Volume -> {v_path}")

    # ---- Cell height across all frames ----
    if all_frames:
        ch_path = out_dir / 'cell_height.png'
        try:
            result = plot_cell_height(all_frames, args.num_layers, ch_path,
                                      n_nearest=args.n_center_nodes)
            if result is not None:
                out_paths.append(ch_path)
                print(f"  Cell height -> {ch_path}")
        except Exception as e:
            print(f"  ERROR generating cell-height plot: {e}")
    elif args.vtk_file:
        print("  Cell-height time series needs all frames; pass --vtk-dir instead "
              "of --vtk-file to enable it.")

    if not out_paths:
        print("  No outputs generated. Pass --vtk-dir or --vtk-file and/or --log.")
        return 1

    print('\nDone. Outputs:')
    for p in out_paths:
        print(f"  {p}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
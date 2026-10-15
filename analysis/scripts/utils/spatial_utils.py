"""
Shared utilities for spatial visualization of cell segmentation data.

Provides common functions for:
- Scale bars
- Cell contour plotting
- Legends
- Color palettes
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon, FancyBboxPatch
from pathlib import Path

from config import SAMPLE_DATA_DIR

# =============================================================================
# Constants
# =============================================================================

# Pixel scale for H&E slides (typical for 20x scanning)
UM_PER_PIXEL = 0.325

# ConSeP color palette - ColorBrewer Set1 (publication-quality)
CONSEP_COLORS = {
    "Malignant epithelial": "#e41a1c",  # Red
    "Healthy epithelial": "#377eb8",    # Blue
    "Stroma": "#4daf4a",                # Green
    "Muscle": "#984ea3",                # Purple
    "Inflammatory": "#ff7f00",          # Orange
    "Other": "#999999",                 # Gray
}

# Ordering for consistent display
CONSEP_ORDER = [
    "Malignant epithelial",
    "Healthy epithelial",
    "Inflammatory",
    "Stroma",
    "Muscle",
    "Other",
]


def get_project_root() -> Path:
    """Get the project root directory."""
    return Path(__file__).parent.parent.parent


def get_sample_path_info(sample_id: str) -> dict:
    """
    Get file paths for a sample's data files.

    For detailed visualizations, the full SpatialData zarrs and H&E images
    live under SAMPLE_DATA_DIR (see config.py).

    Parameters
    ----------
    sample_id : str
        Sample identifier (e.g., 'dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A')

    Returns
    -------
    dict
        Dictionary with paths to zarr, H&E tif, etc.
    """
    sample_dir = SAMPLE_DATA_DIR / sample_id

    return {
        "consep_zarr": sample_dir / f"{sample_id}_consep_spatialdata.zarr",
        "he_tif": sample_dir / f"{sample_id}.tif",
        "sample_id": sample_id,
    }


def add_scale_bar(
    ax,
    length_um: float,
    pixel_scale: float = UM_PER_PIXEL,
    position: str = "lower right",
    color: str = "black",
    fontsize: int = 8,
    bar_height: float = 0.02,
):
    """
    Add a scale bar to a matplotlib axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to add the scale bar to
    length_um : float
        Length of scale bar in micrometers
    pixel_scale : float
        Micrometers per pixel
    position : str
        Position of scale bar ('lower right', 'lower left', 'upper right', 'upper left')
    color : str
        Color of scale bar and text
    fontsize : int
        Font size for the label
    bar_height : float
        Height of scale bar as fraction of axis height
    """
    length_px = length_um / pixel_scale

    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Calculate position (account for potentially inverted y-axis)
    y_range = abs(ylim[1] - ylim[0])
    x_range = abs(xlim[1] - xlim[0])

    margin_x = x_range * 0.05
    margin_y = y_range * 0.08

    if "lower" in position:
        y_pos = max(ylim) - margin_y if ylim[0] < ylim[1] else min(ylim) + margin_y
    else:
        y_pos = min(ylim) + margin_y if ylim[0] < ylim[1] else max(ylim) - margin_y

    if "right" in position:
        x_pos = max(xlim) - margin_x - length_px
    else:
        x_pos = min(xlim) + margin_x

    # Draw scale bar
    ax.plot(
        [x_pos, x_pos + length_px],
        [y_pos, y_pos],
        color=color,
        linewidth=2,
        solid_capstyle="butt",
    )

    # Add label
    text_offset = y_range * 0.03
    if "lower" in position:
        text_y = y_pos - text_offset if ylim[0] < ylim[1] else y_pos + text_offset
    else:
        text_y = y_pos + text_offset if ylim[0] < ylim[1] else y_pos - text_offset

    ax.text(
        x_pos + length_px / 2,
        text_y,
        f"{length_um:.0f} µm",
        ha="center",
        va="top" if "lower" in position else "bottom",
        fontsize=fontsize,
        color=color,
        fontweight="bold",
    )


def plot_contours(
    ax,
    contours_gdf,
    cell_type_colors: dict,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    linewidth: float = 0.6,
    alpha: float = 0.9,
):
    """
    Plot cell contours for a region.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot on
    contours_gdf : geopandas.GeoDataFrame
        GeoDataFrame with cell contours and 'classification_name' column
    cell_type_colors : dict
        Mapping of cell type names to colors
    x_min, x_max, y_min, y_max : float
        Bounding box of region to plot
    linewidth : float
        Width of contour lines
    alpha : float
        Transparency of contours

    Returns
    -------
    tuple
        (cell_counts dict, total_cells int)
    """
    # Filter contours to region using centroid coordinates
    centroids = contours_gdf.geometry.centroid
    mask = (
        (centroids.x >= x_min) & (centroids.x <= x_max) &
        (centroids.y >= y_min) & (centroids.y <= y_max)
    )
    contours_in_region = contours_gdf[mask]

    cell_counts = {}
    for cell_type, color in cell_type_colors.items():
        ct_contours = contours_in_region[
            contours_in_region["classification_name"] == cell_type
        ]
        n_cells = len(ct_contours)
        cell_counts[cell_type] = n_cells

        if n_cells > 0:
            patches = []
            for geom in ct_contours.geometry:
                if geom.geom_type == "Polygon":
                    coords = np.array(geom.exterior.coords)
                    patches.append(Polygon(coords, closed=True))
                elif geom.geom_type == "MultiPolygon":
                    for poly in geom.geoms:
                        coords = np.array(poly.exterior.coords)
                        patches.append(Polygon(coords, closed=True))

            if patches:
                collection = PatchCollection(
                    patches,
                    facecolor="none",
                    edgecolor=color,
                    alpha=alpha,
                    linewidth=linewidth,
                )
                ax.add_collection(collection)

    return cell_counts, len(contours_in_region)


def create_legend(
    ax,
    colors_dict: dict,
    counts_dict: dict,
    title: str = None,
    ncol: int = 1,
    fontsize: int = 7,
    loc: str = "upper right",
    show_pct: bool = True,
):
    """
    Create a clean legend with cell counts and percentages.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to add legend to
    colors_dict : dict
        Mapping of cell type names to colors
    counts_dict : dict
        Mapping of cell type names to counts
    title : str, optional
        Legend title
    ncol : int
        Number of columns
    fontsize : int
        Font size
    loc : str
        Legend location
    show_pct : bool
        Whether to show percentages

    Returns
    -------
    matplotlib.legend.Legend
        The created legend
    """
    handles = []
    labels = []

    # Sort by count descending
    sorted_items = sorted(counts_dict.items(), key=lambda x: x[1], reverse=True)
    total = sum(counts_dict.values())

    for cell_type, count in sorted_items:
        if count > 0 and cell_type in colors_dict:
            handles.append(
                plt.Line2D(
                    [0], [0],
                    marker="o",
                    color="w",
                    markerfacecolor=colors_dict[cell_type],
                    markersize=6,
                )
            )
            if show_pct and total > 0:
                pct = 100 * count / total
                labels.append(f"{cell_type} ({pct:.0f}%)")
            else:
                labels.append(cell_type)

    legend = ax.legend(
        handles,
        labels,
        loc=loc,
        fontsize=fontsize,
        framealpha=0.9,
        title=title,
        title_fontsize=fontsize + 1 if title else None,
        ncol=ncol,
        handletextpad=0.3,
        columnspacing=0.5,
    )
    legend.get_frame().set_linewidth(0.5)
    return legend


def create_color_legend_panel(
    ax,
    colors_dict: dict,
    title: str,
    order: list = None,
    y_start: float = 0.95,
    y_step: float = 0.08,
    fontsize: int = 8,
):
    """
    Create a standalone color legend panel.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw legend on (should have axis='off')
    colors_dict : dict
        Mapping of cell type names to colors
    title : str
        Legend title
    order : list, optional
        Order of cell types to display
    y_start : float
        Starting y position (0-1)
    y_step : float
        Step between items
    fontsize : int
        Font size
    """
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    y_pos = y_start
    ax.text(
        0.05, y_pos, title,
        fontsize=fontsize + 1,
        fontweight="bold",
        transform=ax.transAxes,
    )
    y_pos -= y_step * 0.6

    items = order if order else list(colors_dict.keys())

    for ct in items:
        if ct in colors_dict:
            color = colors_dict[ct]
            ax.add_patch(
                FancyBboxPatch(
                    (0.05, y_pos - 0.02),
                    0.08,
                    0.04,
                    boxstyle="round,pad=0.01",
                    facecolor=color,
                    edgecolor="none",
                    transform=ax.transAxes,
                )
            )
            ax.text(
                0.16, y_pos,
                ct,
                fontsize=fontsize,
                va="center",
                transform=ax.transAxes,
            )
            y_pos -= y_step

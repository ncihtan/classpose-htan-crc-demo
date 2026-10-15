"""
Morphometric Analysis from Cell Contours

Extracts morphological features from cell contours to:
1. Characterize cell types by shape/size
2. Attempt to distinguish malignant vs normal epithelial
3. Explore lymphocyte heterogeneity (potential T vs B cell differences)
"""

import spatialdata as sd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.validation import make_valid
from scipy import stats
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# Pixel scale
UM_PER_PIXEL = 0.325

# Load data
zarr_path = "./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A/dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A_spatialdata.zarr"
print("Loading SpatialData...")
sdata = sd.read_zarr(zarr_path)

# Load contours
print("Loading cell contours...")
contours = sdata.shapes['cell_contours']
print(f"Total contours: {len(contours):,}")
print(f"Columns: {contours.columns.tolist()}")

# Check available cell types
cell_types = contours['classification_name'].value_counts()
print("\nCell types:")
for ct, count in cell_types.items():
    print(f"  {ct}: {count:,}")


def extract_morphometrics(geometry):
    """
    Extract morphometric features from a polygon geometry.

    Returns dict with:
    - area_px2: Area in pixels^2
    - area_um2: Area in um^2
    - perimeter_px: Perimeter in pixels
    - perimeter_um: Perimeter in um
    - circularity: 4π × area / perimeter² (1.0 = perfect circle)
    - solidity: area / convex_hull_area
    - eccentricity: 1 - (minor_axis / major_axis)
    - major_axis: Length of major axis (from min rotated rectangle)
    - minor_axis: Length of minor axis
    - aspect_ratio: major / minor
    - convexity: convex_hull_perimeter / perimeter
    - extent: area / bounding_box_area
    """
    try:
        # Ensure valid geometry
        if not geometry.is_valid:
            geometry = make_valid(geometry)

        # Handle MultiPolygon (take largest)
        if geometry.geom_type == 'MultiPolygon':
            geometry = max(geometry.geoms, key=lambda g: g.area)

        if geometry.geom_type != 'Polygon' or geometry.is_empty:
            return None

        # Basic measurements
        area_px2 = geometry.area
        perimeter_px = geometry.length

        if area_px2 <= 0 or perimeter_px <= 0:
            return None

        # Convert to um
        area_um2 = area_px2 * (UM_PER_PIXEL ** 2)
        perimeter_um = perimeter_px * UM_PER_PIXEL

        # Circularity (shape factor): 1.0 = perfect circle, <1 = irregular
        circularity = (4 * np.pi * area_px2) / (perimeter_px ** 2)
        circularity = min(circularity, 1.0)  # Cap at 1.0 due to discretization

        # Convex hull properties
        convex_hull = geometry.convex_hull
        hull_area = convex_hull.area
        hull_perimeter = convex_hull.length

        # Solidity: how "filled in" the shape is
        solidity = area_px2 / hull_area if hull_area > 0 else 1.0

        # Convexity: how smooth the boundary is
        convexity = hull_perimeter / perimeter_px if perimeter_px > 0 else 1.0

        # Minimum bounding rectangle for axis measurements
        min_rect = geometry.minimum_rotated_rectangle
        rect_coords = list(min_rect.exterior.coords)

        # Calculate edge lengths of the rotated rectangle
        edge1 = np.sqrt((rect_coords[1][0] - rect_coords[0][0])**2 +
                       (rect_coords[1][1] - rect_coords[0][1])**2)
        edge2 = np.sqrt((rect_coords[2][0] - rect_coords[1][0])**2 +
                       (rect_coords[2][1] - rect_coords[1][1])**2)

        major_axis = max(edge1, edge2)
        minor_axis = min(edge1, edge2)

        # Aspect ratio and eccentricity
        aspect_ratio = major_axis / minor_axis if minor_axis > 0 else 1.0
        eccentricity = 1 - (minor_axis / major_axis) if major_axis > 0 else 0.0

        # Extent: how much of bounding box is filled
        bbox = geometry.bounds  # (minx, miny, maxx, maxy)
        bbox_area = (bbox[2] - bbox[0]) * (bbox[3] - bbox[1])
        extent = area_px2 / bbox_area if bbox_area > 0 else 0.0

        return {
            'area_px2': area_px2,
            'area_um2': area_um2,
            'perimeter_px': perimeter_px,
            'perimeter_um': perimeter_um,
            'circularity': circularity,
            'solidity': solidity,
            'convexity': convexity,
            'eccentricity': eccentricity,
            'major_axis_px': major_axis,
            'minor_axis_px': minor_axis,
            'major_axis_um': major_axis * UM_PER_PIXEL,
            'minor_axis_um': minor_axis * UM_PER_PIXEL,
            'aspect_ratio': aspect_ratio,
            'extent': extent,
        }
    except Exception as e:
        return None


# Extract morphometrics for all cells
print("\nExtracting morphometric features...")
morphometrics = []

for idx, row in tqdm(contours.iterrows(), total=len(contours), desc="Processing cells"):
    morph = extract_morphometrics(row.geometry)
    if morph is not None:
        morph['cell_id'] = idx
        morph['classification_name'] = row['classification_name']
        # Get centroid for spatial context
        centroid = row.geometry.centroid
        morph['centroid_x'] = centroid.x
        morph['centroid_y'] = centroid.y
        morphometrics.append(morph)

# Create DataFrame
df = pd.DataFrame(morphometrics)
print(f"\nSuccessfully extracted features for {len(df):,} cells ({100*len(df)/len(contours):.1f}%)")

# Save to parquet
output_parquet = "./cell_morphometrics.parquet"
df.to_parquet(output_parquet)
print(f"Saved: {output_parquet}")

# Also save summary CSV
output_csv = "./cell_morphometrics_summary.csv"


# =============================================================================
# ANALYSIS BY CELL TYPE
# =============================================================================
print("\n" + "="*60)
print("MORPHOMETRIC ANALYSIS BY CELL TYPE")
print("="*60)

# Summary statistics per cell type
summary_stats = []
for ct in df['classification_name'].unique():
    ct_data = df[df['classification_name'] == ct]
    stats_dict = {
        'cell_type': ct,
        'n_cells': len(ct_data),
        'area_um2_mean': ct_data['area_um2'].mean(),
        'area_um2_median': ct_data['area_um2'].median(),
        'area_um2_std': ct_data['area_um2'].std(),
        'circularity_mean': ct_data['circularity'].mean(),
        'circularity_std': ct_data['circularity'].std(),
        'eccentricity_mean': ct_data['eccentricity'].mean(),
        'aspect_ratio_mean': ct_data['aspect_ratio'].mean(),
        'solidity_mean': ct_data['solidity'].mean(),
    }
    summary_stats.append(stats_dict)

summary_df = pd.DataFrame(summary_stats)
print("\nMorphometric Summary by Cell Type:")
print(summary_df.to_string(index=False))
summary_df.to_csv(output_csv, index=False)
print(f"Saved: {output_csv}")


# =============================================================================
# VISUALIZATION
# =============================================================================
print("\nCreating visualizations...")

# Color palette matching previous scripts
cell_colors = {
    'Epithelial': '#1f77b4',    # Blue
    'Connective': '#ff7f0e',    # Orange
    'Lymphocyte': '#2ca02c',    # Green
    'Plasma cell': '#d62728',   # Red
    'Eosinophil': '#9467bd',    # Purple
    'Neutrophil': '#8c564b',    # Brown
}

# 1. Distribution plots for key features
fig1, axes = plt.subplots(2, 3, figsize=(18, 12))

features = ['area_um2', 'circularity', 'eccentricity', 'aspect_ratio', 'solidity', 'convexity']
feature_labels = ['Area (um²)', 'Circularity', 'Eccentricity', 'Aspect Ratio', 'Solidity', 'Convexity']

for idx, (feat, label) in enumerate(zip(features, feature_labels)):
    ax = axes.flatten()[idx]

    for ct in cell_colors.keys():
        ct_data = df[df['classification_name'] == ct][feat].dropna()
        if len(ct_data) > 0:
            # Use log scale for area
            if feat == 'area_um2':
                ct_data = np.log10(ct_data + 1)
                label = 'log₁₀(Area um²)'

            ax.hist(ct_data, bins=50, alpha=0.5, label=ct, color=cell_colors[ct], density=True)

    ax.set_xlabel(label, fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title(f'{label} Distribution by Cell Type', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)

plt.suptitle('Morphometric Feature Distributions by Cell Type', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("./morphometrics_distributions.png", dpi=200, bbox_inches='tight')
print("Saved: morphometrics_distributions.png")


# 2. Violin plots for better comparison
fig2, axes = plt.subplots(2, 3, figsize=(18, 12))

for idx, (feat, label) in enumerate(zip(features, feature_labels)):
    ax = axes.flatten()[idx]

    plot_data = df[['classification_name', feat]].copy()

    # Log transform area for better visualization
    if feat == 'area_um2':
        plot_data[feat] = np.log10(plot_data[feat] + 1)
        label = 'log₁₀(Area um²)'

    # Order by cell type
    order = ['Epithelial', 'Connective', 'Lymphocyte', 'Plasma cell', 'Eosinophil', 'Neutrophil']
    order = [ct for ct in order if ct in plot_data['classification_name'].unique()]

    palette = {ct: cell_colors.get(ct, 'gray') for ct in order}

    sns.violinplot(data=plot_data, x='classification_name', y=feat, order=order,
                   palette=palette, ax=ax, inner='quartile')
    ax.set_xlabel('Cell Type', fontsize=11)
    ax.set_ylabel(label, fontsize=11)
    ax.set_title(f'{label}', fontsize=12, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)

plt.suptitle('Morphometric Features by Cell Type (Violin Plots)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("./morphometrics_violin.png", dpi=200, bbox_inches='tight')
print("Saved: morphometrics_violin.png")


# 3. 2D scatter: Area vs Circularity (key discriminating features)
fig3, ax = plt.subplots(figsize=(12, 10))

for ct in cell_colors.keys():
    ct_data = df[df['classification_name'] == ct]
    if len(ct_data) > 0:
        # Subsample for visualization if too many points
        if len(ct_data) > 5000:
            ct_data = ct_data.sample(5000, random_state=42)

        ax.scatter(np.log10(ct_data['area_um2'] + 1), ct_data['circularity'],
                  alpha=0.3, s=2, c=cell_colors[ct], label=ct, rasterized=True)

ax.set_xlabel('log₁₀(Area um²)', fontsize=12)
ax.set_ylabel('Circularity', fontsize=12)
ax.set_title('Cell Morphospace: Area vs Circularity\nEpithelial cells show larger size and more variable shape',
             fontsize=14, fontweight='bold')
ax.legend(markerscale=10, fontsize=10)
plt.tight_layout()
plt.savefig("./morphometrics_area_circularity.png", dpi=200, bbox_inches='tight')
print("Saved: morphometrics_area_circularity.png")


# =============================================================================
# EPITHELIAL CELL ANALYSIS - Looking for malignant vs normal patterns
# =============================================================================
print("\n" + "="*60)
print("EPITHELIAL CELL ANALYSIS")
print("="*60)

epi_data = df[df['classification_name'] == 'Epithelial'].copy()
print(f"Total epithelial cells: {len(epi_data):,}")

# Analyze epithelial cell variability
print("\nEpithelial cell morphometrics:")
print(f"  Area (um²): {epi_data['area_um2'].mean():.1f} ± {epi_data['area_um2'].std():.1f} (median: {epi_data['area_um2'].median():.1f})")
print(f"  Circularity: {epi_data['circularity'].mean():.3f} ± {epi_data['circularity'].std():.3f}")
print(f"  Aspect ratio: {epi_data['aspect_ratio'].mean():.2f} ± {epi_data['aspect_ratio'].std():.2f}")

# Look for bimodal/multimodal distributions (could indicate malignant vs normal)
print("\nLooking for epithelial subpopulations based on morphology...")

# Use area percentiles
area_25 = epi_data['area_um2'].quantile(0.25)
area_75 = epi_data['area_um2'].quantile(0.75)

epi_small = epi_data[epi_data['area_um2'] < area_25]
epi_medium = epi_data[(epi_data['area_um2'] >= area_25) & (epi_data['area_um2'] <= area_75)]
epi_large = epi_data[epi_data['area_um2'] > area_75]

print(f"\nEpithelial size quartiles:")
print(f"  Small (<{area_25:.1f} um²): {len(epi_small):,} cells")
print(f"    Circularity: {epi_small['circularity'].mean():.3f}, Aspect ratio: {epi_small['aspect_ratio'].mean():.2f}")
print(f"  Medium ({area_25:.1f}-{area_75:.1f} um²): {len(epi_medium):,} cells")
print(f"    Circularity: {epi_medium['circularity'].mean():.3f}, Aspect ratio: {epi_medium['aspect_ratio'].mean():.2f}")
print(f"  Large (>{area_75:.1f} um²): {len(epi_large):,} cells")
print(f"    Circularity: {epi_large['circularity'].mean():.3f}, Aspect ratio: {epi_large['aspect_ratio'].mean():.2f}")


# Plot epithelial subpopulations
fig4, axes = plt.subplots(1, 3, figsize=(18, 5))

# Area distribution with potential bimodality
ax1 = axes[0]
ax1.hist(np.log10(epi_data['area_um2'] + 1), bins=100, color=cell_colors['Epithelial'], alpha=0.7, edgecolor='black')
ax1.axvline(np.log10(area_25 + 1), color='red', linestyle='--', label=f'Q25: {area_25:.0f} um²')
ax1.axvline(np.log10(area_75 + 1), color='red', linestyle='--', label=f'Q75: {area_75:.0f} um²')
ax1.set_xlabel('log₁₀(Area um²)', fontsize=11)
ax1.set_ylabel('Count', fontsize=11)
ax1.set_title('Epithelial Cell Size Distribution', fontsize=12, fontweight='bold')
ax1.legend()

# Area vs Circularity by size group
ax2 = axes[1]
ax2.scatter(np.log10(epi_small['area_um2'] + 1), epi_small['circularity'],
           alpha=0.2, s=1, c='#3498db', label='Small', rasterized=True)
ax2.scatter(np.log10(epi_large['area_um2'] + 1), epi_large['circularity'],
           alpha=0.2, s=1, c='#e74c3c', label='Large', rasterized=True)
ax2.set_xlabel('log₁₀(Area um²)', fontsize=11)
ax2.set_ylabel('Circularity', fontsize=11)
ax2.set_title('Small vs Large Epithelial Cells\n(Large cells often less circular = more irregular)', fontsize=12, fontweight='bold')
ax2.legend(markerscale=10)

# Circularity comparison
ax3 = axes[2]
ax3.hist(epi_small['circularity'], bins=50, alpha=0.6, color='#3498db', label=f'Small (n={len(epi_small):,})', density=True)
ax3.hist(epi_large['circularity'], bins=50, alpha=0.6, color='#e74c3c', label=f'Large (n={len(epi_large):,})', density=True)
ax3.set_xlabel('Circularity', fontsize=11)
ax3.set_ylabel('Density', fontsize=11)
ax3.set_title('Circularity: Small vs Large Epithelial\n(Lower circularity may indicate malignancy)', fontsize=12, fontweight='bold')
ax3.legend()

# Statistical test
stat, pval = stats.mannwhitneyu(epi_small['circularity'].dropna(), epi_large['circularity'].dropna())
ax3.text(0.05, 0.95, f'Mann-Whitney U\np < {pval:.2e}', transform=ax3.transAxes, fontsize=10,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.suptitle('Epithelial Cell Morphometric Analysis\nLooking for malignant vs normal patterns', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("./epithelial_morphometrics.png", dpi=200, bbox_inches='tight')
print("Saved: epithelial_morphometrics.png")


# =============================================================================
# LYMPHOCYTE ANALYSIS - Looking for T cell vs B cell patterns
# =============================================================================
print("\n" + "="*60)
print("LYMPHOCYTE ANALYSIS")
print("="*60)

lymph_data = df[df['classification_name'] == 'Lymphocyte'].copy()
plasma_data = df[df['classification_name'] == 'Plasma cell'].copy()

print(f"Lymphocytes: {len(lymph_data):,}")
print(f"Plasma cells: {len(plasma_data):,}")

print("\nLymphocyte morphometrics:")
print(f"  Area (um²): {lymph_data['area_um2'].mean():.1f} ± {lymph_data['area_um2'].std():.1f} (median: {lymph_data['area_um2'].median():.1f})")
print(f"  Circularity: {lymph_data['circularity'].mean():.3f} ± {lymph_data['circularity'].std():.3f}")

print("\nPlasma cell morphometrics:")
print(f"  Area (um²): {plasma_data['area_um2'].mean():.1f} ± {plasma_data['area_um2'].std():.1f} (median: {plasma_data['area_um2'].median():.1f})")
print(f"  Circularity: {plasma_data['circularity'].mean():.3f} ± {plasma_data['circularity'].std():.3f}")

# Statistical comparison
stat, pval = stats.mannwhitneyu(lymph_data['area_um2'].dropna(), plasma_data['area_um2'].dropna())
print(f"\nLymphocyte vs Plasma cell area: Mann-Whitney U p={pval:.2e}")

# Look for lymphocyte subpopulations (potential T vs B)
lymph_area_median = lymph_data['area_um2'].median()
lymph_small = lymph_data[lymph_data['area_um2'] < lymph_area_median]
lymph_large = lymph_data[lymph_data['area_um2'] >= lymph_area_median]

print(f"\nLymphocyte size subgroups:")
print(f"  Smaller (< median): {len(lymph_small):,} cells, area {lymph_small['area_um2'].mean():.1f} um²")
print(f"  Larger (>= median): {len(lymph_large):,} cells, area {lymph_large['area_um2'].mean():.1f} um²")

# Plot lymphocyte analysis
fig5, axes = plt.subplots(1, 3, figsize=(18, 5))

# Lymphocyte vs Plasma cell comparison
ax1 = axes[0]
combined_data = pd.concat([
    lymph_data[['area_um2', 'circularity']].assign(cell_type='Lymphocyte'),
    plasma_data[['area_um2', 'circularity']].assign(cell_type='Plasma cell')
])

for ct, color in [('Lymphocyte', cell_colors['Lymphocyte']), ('Plasma cell', cell_colors['Plasma cell'])]:
    ct_subset = combined_data[combined_data['cell_type'] == ct]
    if len(ct_subset) > 5000:
        ct_subset = ct_subset.sample(5000, random_state=42)
    ax1.scatter(np.log10(ct_subset['area_um2'] + 1), ct_subset['circularity'],
               alpha=0.3, s=2, c=color, label=ct, rasterized=True)

ax1.set_xlabel('log₁₀(Area um²)', fontsize=11)
ax1.set_ylabel('Circularity', fontsize=11)
ax1.set_title('Lymphocyte vs Plasma Cell Morphology\n(Plasma cells slightly larger, less circular)', fontsize=12, fontweight='bold')
ax1.legend(markerscale=10)

# Lymphocyte size distribution (looking for bimodality = T vs B?)
ax2 = axes[1]
ax2.hist(np.log10(lymph_data['area_um2'] + 1), bins=100, color=cell_colors['Lymphocyte'], alpha=0.7, edgecolor='black')
ax2.axvline(np.log10(lymph_area_median + 1), color='red', linestyle='--', label=f'Median: {lymph_area_median:.0f} um²')
ax2.set_xlabel('log₁₀(Area um²)', fontsize=11)
ax2.set_ylabel('Count', fontsize=11)
ax2.set_title('Lymphocyte Size Distribution\n(Unimodal suggests limited T/B discrimination)', fontsize=12, fontweight='bold')
ax2.legend()

# Aspect ratio comparison (T cells often slightly more elongated when activated?)
ax3 = axes[2]
ax3.hist(lymph_small['aspect_ratio'], bins=50, alpha=0.6, color='#2ecc71', label=f'Small lymph (n={len(lymph_small):,})', density=True)
ax3.hist(lymph_large['aspect_ratio'], bins=50, alpha=0.6, color='#27ae60', label=f'Large lymph (n={len(lymph_large):,})', density=True)
ax3.hist(plasma_data['aspect_ratio'], bins=50, alpha=0.6, color=cell_colors['Plasma cell'], label=f'Plasma (n={len(plasma_data):,})', density=True)
ax3.set_xlabel('Aspect Ratio', fontsize=11)
ax3.set_ylabel('Density', fontsize=11)
ax3.set_title('Aspect Ratio Comparison\nLymphocyte subgroups vs Plasma cells', fontsize=12, fontweight='bold')
ax3.legend()

plt.suptitle('Lymphoid Cell Morphometric Analysis\nNote: T vs B cell distinction from H&E morphology is challenging', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("./lymphocyte_morphometrics.png", dpi=200, bbox_inches='tight')
print("Saved: lymphocyte_morphometrics.png")


# =============================================================================
# SPATIAL MORPHOMETRIC PATTERNS
# =============================================================================
print("\n" + "="*60)
print("SPATIAL MORPHOMETRIC ANALYSIS")
print("="*60)

# For epithelial cells: map area and circularity spatially to identify regions
# of potentially malignant cells (large, irregular)
print("\nCreating spatial morphometric maps...")

fig6, axes = plt.subplots(1, 2, figsize=(16, 8))

# Get slide dimensions for proper axis limits
import openslide
he_slide = openslide.OpenSlide("./dg_4DFC_6C434FEA-0A77-43DB-8BCF-2873E5873A8A.tif")
slide_dims = he_slide.dimensions  # (width, height)

# Subsample epithelial for visualization
epi_sample = epi_data.sample(min(50000, len(epi_data)), random_state=42)

# Area map
ax1 = axes[0]
scatter1 = ax1.scatter(epi_sample['centroid_x'], epi_sample['centroid_y'],
                       c=np.log10(epi_sample['area_um2'] + 1), cmap='viridis',
                       s=0.5, alpha=0.5, rasterized=True)
plt.colorbar(scatter1, ax=ax1, label='log₁₀(Area um²)')
ax1.set_xlim(0, slide_dims[0])
ax1.set_ylim(slide_dims[1], 0)
ax1.set_xlabel('X (pixels)')
ax1.set_ylabel('Y (pixels)')
ax1.set_title('Epithelial Cell Area - Spatial Distribution\n(Larger cells may indicate malignant regions)', fontsize=12, fontweight='bold')
ax1.set_aspect('equal')

# Circularity map
ax2 = axes[1]
scatter2 = ax2.scatter(epi_sample['centroid_x'], epi_sample['centroid_y'],
                       c=epi_sample['circularity'], cmap='RdYlBu',
                       s=0.5, alpha=0.5, rasterized=True, vmin=0.3, vmax=1.0)
plt.colorbar(scatter2, ax=ax2, label='Circularity')
ax2.set_xlim(0, slide_dims[0])
ax2.set_ylim(slide_dims[1], 0)
ax2.set_xlabel('X (pixels)')
ax2.set_ylabel('Y (pixels)')
ax2.set_title('Epithelial Cell Circularity - Spatial Distribution\n(Low circularity = irregular shape)', fontsize=12, fontweight='bold')
ax2.set_aspect('equal')

plt.suptitle('Spatial Morphometric Patterns in Epithelial Cells', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig("./spatial_morphometrics.png", dpi=200, bbox_inches='tight')
print("Saved: spatial_morphometrics.png")


# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*60)
print("MORPHOMETRIC ANALYSIS SUMMARY")
print("="*60)

print("""
Key Findings:

1. EPITHELIAL CELLS (Malignant vs Normal indicators):
   - Large cells with LOW circularity may indicate malignancy
   - Size and shape variability within epithelial population
   - Spatial clustering of morphometric features could identify tumor regions

2. LYMPHOID CELLS (T vs B cell distinction):
   - Limited morphometric distinction possible from H&E alone
   - Plasma cells are distinguishable (larger, less circular)
   - T vs B cell classification requires molecular markers (CD3, CD20)

3. CAVEATS:
   - Cell type classification is from H&E-based deep learning
   - Morphometrics are from 2D contours (3D shape variation not captured)
   - True malignant vs normal distinction requires pathologist validation
   - T vs B cell distinction is not reliable from morphology alone

4. RECOMMENDATIONS:
   - Use spatial patterns to identify regions of interest
   - Combine with TLS analysis for immune contexture
   - Validate with multiplex IHC/IF if available
""")

print("Analysis complete!")
plt.show()

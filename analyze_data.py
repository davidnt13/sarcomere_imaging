import os
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from skimage import io
import seaborn as sns
from src.nuclei_count import find_and_count_nuclei
from scipy.stats import linregress
from src.tissue_volume import find_myofibril_density_vol

# ---------- HELPER METHODS FOR ANGLE PLOTTING AND CALCULATIONS ----------

def find_sarcgraph_angle_mean(file_path):
    sg_df = pd.read_csv(file_path)
    angles = sg_df['angle'].dropna().to_numpy() * (180/np.pi)
    return np.mean(angles), np.std(angles)

# ---------- HELPER METHODS FOR LENGTH PLOTTING AND CALCULATIONS ----------

def find_sarcgraph_length_mean(file_path):
    sg_df = pd.read_csv(file_path)
    sg_df['Length_microns'] = sg_df['length'] * pixel_size
    sg_df.to_csv(file_path)
    return np.mean(sg_df['Length_microns']), np.max(sg_df['sarc_id'])#, len(sg_df['Length_microns'].dropna())

# ---------- PLOTTING LENGTHS AND ANGLES ----------

pixel_size = 0.312 # microns
original_date = '6-11'
# stain = 'vimentin'
stain = 'cx43'

title = 'both_cx43_vimentin_0cf'
# title = 'all_cx43'
save_stats_fig = False

data_dir = f'../../Data/{stain}_CZI_OG'
data_dir_vim = f'../../Data/vimentin_CZI_OG'
data_dir_path = Path(data_dir)
data_dir_vim_path = Path(data_dir_vim)
# files = [f.name for f in data_dir_path.glob('IPSC-CF_fibers*MAX.tif')]
# files = [f.name for f in data_dir_path.glob('IPSC-CF_fibers*.tif') if 'MAX' not in f.name]
files = [f.name for f in data_dir_path.glob('0CF*.tif') if 'MAX' not in f.name]
files.extend([f.name for f in data_dir_vim_path.glob('0CF*.tif') if 'MAX' not in f.name])
# files = [f.name for f in data_dir_path.glob('*.tif') if 'MAX' not in f.name]
# files.extend([f.name for f in data_dir_vim_path.glob('*.tif') if 'MAX' not in f.name])

info_cols = ['Filename', 'Cell Type', 'Fibers', 'Stain', 'Sarcomere Count', 'Mean Length (µm)'
             , 'Mean Angle (Degrees)', 'Stdev Angle (Degrees)', 'Nuclei Count', 'Nuclei Count 3D',
             'Myofibril Density', 'Myofibril Density (Vol)']
image_info = pd.DataFrame(columns=info_cols)

for img in files:
    if 'vimentin' in img: data_dir = data_dir_vim

    img_name, _ = os.path.splitext(img)
    img_dir = f'./Image_Tests/{original_date}-Work/{img_name}'
    sarcgraph_csv_path = f'{img_dir}/{img_name}_sg_sarcgraph_output/sarcomeres.csv'

    parts = img_name.split('_')
    if len(parts) >= 3:
        cell_type = parts[0]
        structure = parts[1]
        marker = parts[2]
    else:
        cell_type = structure = marker = None  # or skip this img if preferred

    # Sarcgraph angles
    sg_mean_angle, sg_angle_std = find_sarcgraph_angle_mean(sarcgraph_csv_path)

    # Sarcgraph mean length
    sg_mean_length, sg_count = find_sarcgraph_length_mean(sarcgraph_csv_path)

    # Nuclei Count
    thresh_mult = 0.5                    
    min_distance_2d = 15 // 4
    min_distance_3d = 10  

    max_file_path = f'{data_dir}/{img_name}_MAX.tif'
    if os.path.exists(max_file_path):
        max_img = io.imread(max_file_path)
        nuclei_count = find_and_count_nuclei(max_img, thresh_mult=thresh_mult, 
                                             min_distance=min_distance_2d, show_image=False)
    else:
        print(f"File not found: {max_file_path}")
        nuclei_count = np.nan

    img = io.imread(f'{data_dir}/{img}')
    nuclei_count_3d = find_and_count_nuclei(img, thresh_mult=thresh_mult, min_distance=min_distance_3d,
                                            three_dim=True)
    print(f'{img_name} 3D Count: {nuclei_count_3d}')

    # Z Bands Density
    if os.path.exists(max_file_path):
        max_img = io.imread(max_file_path)
        sarc_mask = io.imread(f'./Image_Tests/6-11-Work/{img_name}/{img_name}_sg_filtered/sarcomere_mask.tif')
        myofibril_density = find_myofibril_density(max_img, sarc_mask)
    else:
        print(f"File not found: {max_file_path}")
        myofibril_density = np.nan

    sarc_mask = io.imread(f'./Image_Tests/6-11-Work/{img_name}/{img_name}_sg_filtered/sarcomere_mask.tif')
    myofibril_density_vol = find_myofibril_density_vol(img, sarc_mask)

    # Store Info
    img_details = [img_name, cell_type, structure, marker, sg_count, 
                   sg_mean_length, sg_mean_angle, sg_angle_std, nuclei_count, 
                   nuclei_count_3d, myofibril_density, myofibril_density_vol]

    # Append using pd.concat
    new_row = pd.DataFrame([img_details], columns=image_info.columns)
    image_info = pd.concat([image_info, new_row], ignore_index=True)

image_info.to_csv(f"./Image_Analysis_Results/image_analysis_details_{title}.csv")

#################### PLOTTING ###################################
image_info['Group'] = image_info['Cell Type'] + ' - ' + image_info['Fibers']
agg_df = image_info.groupby('Group').agg({
    'Mean Length (µm)': ['mean', 'std'],
    'Mean Angle (Degrees)': ['mean', 'std'],
    'Stdev Angle (Degrees)': ['mean', 'std'],
    'Nuclei Count 3D': ['mean', 'std'],
    'Myofibril Density': ['mean', 'std']
}).reset_index()

# Flatten multi-index columns
agg_df.columns = ['Group',
                  'Mean Length Mean', 'Mean Length Std',
                  'Mean Angle Mean', 'Mean Angle Std',
                  'Stdev Angle Mean', 'Stdev Angle Std',
                  'Nuclei Count 3D Mean', 'Nuclei Count 3D Std',
                  'Myofibril Density Mean', 'Myofibril Density Std']

palette = sns.color_palette('Dark2', n_colors=len(agg_df))
fig, axes = plt.subplots(1, 5, figsize=(24, 6))

def barplot_with_error(ax, x, y_mean, y_err, title, ylabel, colors):
    ax.bar(x, y_mean, yerr=y_err, capsize=5, color=colors)
    ax.set_title(title, fontsize=14)
    ax.set_xticklabels(x, rotation=45, ha='right', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis='x', rotation=45)

# 1. Mean Length
barplot_with_error(
    axes[0],
    agg_df['Group'],
    agg_df['Mean Length Mean'],
    agg_df['Mean Length Std'],
    'Mean Length (µm)',
    'Length (µm)',
    palette
)

# 2. Mean Angle
barplot_with_error(
    axes[1],
    agg_df['Group'],
    agg_df['Mean Angle Mean'],
    agg_df['Mean Angle Std'],
    'Mean Angle (Degrees)',
    'Angle (Degrees)',
    palette
)

# 3. Stdev Angle
barplot_with_error(
    axes[2],
    agg_df['Group'],
    agg_df['Stdev Angle Mean'],
    agg_df['Stdev Angle Std'],
    'Stdev Angle (Degrees)',
    'Angle Std Dev (Degrees)',
    palette
)

# 4. Stdev Angle
barplot_with_error(
    axes[3],
    agg_df['Group'],
    agg_df['Nuclei Count 3D Mean'],
    agg_df['Nuclei Count 3D Std'],
    'Nuclei Count 3D Mean',
    'Nuclei Count 3D',
    palette
)

# 5. Stdev Angle
barplot_with_error(
    axes[4],
    agg_df['Group'],
    agg_df['Myofibril Density Mean'],
    agg_df['Myofibril Density Std'],
    'Myofibril Density Mean',
    'Myofibril Density',
    palette
)

plt.tight_layout()
if save_stats_fig: plt.savefig(f"./Image_Analysis_Results/image_analysis_details_{title}.png", dpi=300)
plt.show()

# Comparing Nuclei Counting with MAX vs 3D
image_info['Nuclei Count'] = pd.to_numeric(image_info['Nuclei Count'], errors='coerce')
image_info['Nuclei Count 3D'] = pd.to_numeric(image_info['Nuclei Count 3D'], errors='coerce')
valid = image_info.dropna(subset=['Nuclei Count', 'Nuclei Count 3D'])

r = np.corrcoef(valid['Nuclei Count'], valid['Nuclei Count 3D'])[0,1]
slope, intercept, r_value, p_value, std_err = linregress(valid['Nuclei Count'], valid['Nuclei Count 3D'])

plt.figure(figsize=(6,6), dpi=120)
plt.scatter(valid['Nuclei Count'], valid['Nuclei Count 3D'], alpha=0.7, edgecolor='k', label='Data')

x_vals = np.array([0, 100])
y_vals = slope * x_vals + intercept
plt.plot(x_vals, y_vals, 'b-', label=f'Best fit: y = {slope:.2f}x + {intercept:.2f}')
plt.plot(x_vals, x_vals, 'r--', label='y = x')

plt.xlabel('Nuclei Count (2D)')
plt.ylabel('Nuclei Count (3D)')
plt.title(f'Comparison of 2D vs 3D Counting\nPearson r = {r:.3f}')
plt.xlim([0, 100])
plt.ylim([0, 100])
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Comparing Sarcomere Count to Myofibril Density
image_info['Sarcomere Count'] = pd.to_numeric(image_info['Sarcomere Count'], errors='coerce')
image_info['Myofibril Density'] = pd.to_numeric(image_info['Myofibril Density'], errors='coerce')
valid = image_info.dropna(subset=['Sarcomere Count', 'Myofibril Density'])

r = np.corrcoef(valid['Sarcomere Count'], valid['Myofibril Density'])[0,1]
slope, intercept, r_value, p_value, std_err = linregress(valid['Sarcomere Count'], valid['Myofibril Density'])

plt.figure(figsize=(6,6), dpi=120)
plt.scatter(valid['Sarcomere Count'], valid['Myofibril Density'], alpha=0.7, edgecolor='k', label='Data')

x_vals = np.array([np.min(valid['Sarcomere Count']), np.max(valid['Sarcomere Count'])])
y_vals = slope * x_vals + intercept
plt.plot(x_vals, y_vals, 'b-', label=f'Best fit: y = {slope:.4f}x + {intercept:.4f}')

plt.xlabel('Sarcomere Count')
plt.ylabel('Myofibril Density')
plt.title(f'Correlation of Sarcomere Count and Myofibril Density\nPearson r = {r:.3f}')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Comparing Sarcomere Count to Myofibril Density (Vol)
image_info['Sarcomere Count'] = pd.to_numeric(image_info['Sarcomere Count'], errors='coerce')
image_info['Myofibril Density (Vol)'] = pd.to_numeric(image_info['Myofibril Density (Vol)'], errors='coerce')
valid = image_info.dropna(subset=['Sarcomere Count', 'Myofibril Density (Vol)'])

r = np.corrcoef(valid['Sarcomere Count'], valid['Myofibril Density (Vol)'])[0,1]
slope, intercept, r_value, p_value, std_err = linregress(valid['Sarcomere Count'], valid['Myofibril Density (Vol)'])

plt.figure(figsize=(6,6), dpi=120)
plt.scatter(valid['Sarcomere Count'], valid['Myofibril Density (Vol)'], alpha=0.7, edgecolor='k', label='Data')

x_vals = np.array([np.min(valid['Sarcomere Count']), np.max(valid['Sarcomere Count'])])
y_vals = slope * x_vals + intercept
plt.plot(x_vals, y_vals, 'b-', label=f'Best fit: y = {slope:.4f}x + {intercept:.4f}')

plt.xlabel('Sarcomere Count')
plt.ylabel('Myofibril Density')
plt.title(f'0CF - Correlation of Sarcomere Count and Myofibril Density (Vol)\nPearson r = {r:.3f}')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
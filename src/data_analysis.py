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


def create_statistics_df(original_image_dir, folders, pixel_size=0.312):

    info_cols = ['Filename', 'Cell Type', 'Fibers', 'Stain', 'Sarcomere Count', 'Mean Length (µm)'
             , 'Mean Angle (Degrees)', 'Stdev Angle (Degrees)', 'Nuclei Count 3D', 'Myofibril Density (Vol)']
    
    image_info = pd.DataFrame(columns=info_cols)

    def find_sarcgraph_angle_mean(file_path):
        sg_df = pd.read_csv(file_path)
        angles = sg_df['angle'].dropna().to_numpy() * (180/np.pi)
        return np.mean(angles), np.std(angles)

    def find_sarcgraph_length_mean(file_path):
        sg_df = pd.read_csv(file_path)
        sg_df['Length_microns'] = sg_df['length'] * pixel_size
        sg_df.to_csv(file_path)
        return np.mean(sg_df['Length_microns']), np.max(sg_df['sarc_id'])

    for folder in folders:

        img_name = os.path.basename(os.path.normpath(folder))
        sarcgraph_csv_path = f'{folder}/sarcomeres.csv'

        parts = img_name.split('_')
        if len(parts) >= 3:
            cell_type = parts[0]
            structure = parts[1]
            marker = parts[2]
        else:
            cell_type = structure = marker = None

        # Sarcgraph angles
        sg_mean_angle, sg_angle_std = find_sarcgraph_angle_mean(sarcgraph_csv_path)

        # Sarcgraph mean length
        sg_mean_length, sg_count = find_sarcgraph_length_mean(sarcgraph_csv_path)

        # Nuclei Count
        thresh_mult = 0.5
        min_distance_3d = 10  

        img = io.imread(f'{original_image_dir}/{img_name}.tif')
        nuclei_count_3d = find_and_count_nuclei(img, thresh_mult=thresh_mult, min_distance=min_distance_3d,
                                                three_dim=True)
        # print(f'{img_name} 3D Count: {nuclei_count_3d}')

        # Z Bands Density
        sarc_mask = io.imread(f'{folder}/sarcomere_mask.tif')
        myofibril_density_vol = find_myofibril_density_vol(img, sarc_mask)

        # Store Info
        img_details = [img_name, cell_type, structure, marker, sg_count, 
                    sg_mean_length, sg_mean_angle, sg_angle_std, 
                    nuclei_count_3d, myofibril_density_vol]

        # Append using pd.concat
        new_row = pd.DataFrame([img_details], columns=image_info.columns)
        image_info = pd.concat([image_info, new_row], ignore_index=True)

    return image_info

def plot_values_from_df(image_info, cols_to_plot, save_dir, title):
    image_info['Group'] = image_info['Cell Type'] + ' - ' + image_info['Fibers']

    # Dynamic aggregation
    agg_dict = {col: ['mean', 'std'] for col in cols_to_plot}
    agg_df = image_info.groupby('Group').agg(agg_dict)
    agg_df.columns = [f'{col} {stat.title()}' for col, stat in agg_df.columns]
    agg_df = agg_df.reset_index()

    palette = sns.color_palette('Dark2', n_colors=len(agg_df))
    fig, axes = plt.subplots(1, len(cols_to_plot), figsize=(6 * len(cols_to_plot), 6))

    if len(cols_to_plot) == 1:
        axes = [axes]

    def barplot_with_error(ax, x, y_mean, y_err, title, ylabel, colors):
        ax.bar(x, y_mean, yerr=y_err, capsize=5, color=colors)
        ax.set_xticks(range(len(x)))
        ax.set_xticklabels(x, rotation=45, ha='right', fontsize=10)
        ax.set_title(title)
        ax.set_ylabel(ylabel)

    for i, col in enumerate(cols_to_plot):
        barplot_with_error(
            axes[i],
            agg_df['Group'],
            agg_df[f'{col} Mean'],
            agg_df[f'{col} Std'],
            title=col,
            ylabel=col,
            colors=palette
        )

    plt.tight_layout()
    plt.savefig(f'{save_dir}/{title}.png')
    plt.show()

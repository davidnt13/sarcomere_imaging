import os
import matplotlib.pyplot as plt
from sarcasm import *
from skimage import io

def compute_sarcasm(filepath, save_path, filters_used, visualize=True, pixel_size=0.312):
    sarc = Structure(filepath, pixelsize=pixel_size, auto_save=True)
    sarc.detect_sarcomeres(max_patch_size=(512, 512), clip_thres=(0,100))
    sarc.analyze_z_bands()
    sarc.analyze_sarcomere_vectors()#linewidth=0.9)

    frame = 0 if 'MAX' in filepath else 8
    if visualize: visualize_sarcasm(sarc, frame)

    save_sarcasm_output(sarc, save_path, filters_used)

    zbands_dir, _ = os.path.splitext(filepath)
    zbands_img = io.imread(f'{zbands_dir}/zbands.tif')

    return sarc, zbands_img

def visualize_sarcasm(sarc, frame):
    
    fig, axs = plt.subplots(ncols=3, figsize=(12, 4), dpi=300)
    Plots.plot_image(axs[0], sarc, frame=frame, title='Image')
    Plots.plot_z_bands(axs[1], sarc, frame=frame, title='Z-bands', cmap='Greys_r')
    Plots.plot_cell_mask(axs[2], sarc, frame=frame, title='Cell Mask')
    plt.tight_layout()
    plt.show()
    
    fig, axs = plt.subplots(ncols=2, figsize=(8, 4), dpi=300)
    Plots.plot_image(axs[0], sarc, frame=frame, title='Image')
    Plots.plot_z_segmentation(axs[1], sarc, frame=frame, title='Z-band Segmentation')
    plt.tight_layout()
    plt.show()
    
    fig, axs = plt.subplots(ncols=3, figsize=(12, 4), dpi=300)
    Plots.plot_image(axs[0], sarc, frame=frame, title='Image')
    Plots.plot_sarcomere_vectors(axs[1], sarc, frame=frame, s_points=10, title='Sarcomere Vectors')
    Plots.plot_sarcomere_mask(axs[2], sarc, frame=frame, title='Sarcomere Mask')
    plt.tight_layout()
    plt.show()

def save_sarcasm_output(sarc, save_path, filters):
    save_path = f'{save_path}_{filters}_sarcasm_output.txt'
    structure_dict = export.Export.get_structure_dict(sarc)
    with open(save_path, 'w') as f:
        f.write(str(structure_dict))
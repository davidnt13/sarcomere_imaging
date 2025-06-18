#!/usr/bin/env python

# Imports
import os
import tkinter as tk
from tkinter import filedialog

from skimage import io
from pathlib import Path
import pandas as pd
import shutil

from src.sarcgraph_analysis import compute_sarcgraph, visualize_sarcgraph
from src.sarcasm_analysis import compute_sarcasm
from src.thresholding_pipeline import separate_sarc_channel_and_save, final_filtering

def get_folder_and_date():
    result = {"folder": None, "date": None}

    def select_folder():
        folder = filedialog.askdirectory(title="Select a Folder")
        folder_var.set(folder)

    def submit():
        result["folder"] = folder_var.get()
        result["date"] = date_var.get()
        window.destroy()

    # Set up the main window
    window = tk.Tk()
    window.title("Select Folder and Enter Date")
    window.geometry("400x200")

    # Folder selection
    folder_var = tk.StringVar()
    tk.Label(window, text="Select folder:").pack(pady=(10, 0))
    tk.Button(window, text="Browse", command=select_folder).pack()
    tk.Label(window, textvariable=folder_var).pack()

    # Date input
    date_var = tk.StringVar()
    tk.Label(window, text="Enter date (MM-DD):").pack(pady=(10, 0))
    tk.Entry(window, textvariable=date_var).pack()

    # Submit button
    tk.Button(window, text="Submit", command=submit).pack(pady=10)
    window.mainloop()

    return result["folder"], result["date"]


def save_specific_files(file_names, original_directory, new_directory):
    original_path = Path(original_directory)
    new_path = Path(new_directory)
    new_path.mkdir(parents=True, exist_ok=True)

    # Convert list to a set for faster lookup
    file_names_set = set(file_names)

    # Move files that match by stem (name without extension)
    for file in original_path.iterdir():
        if file.is_file() and file.stem in file_names_set:
            shutil.move(str(file), new_path / file.name)

    # Delete the original directory and all remaining contents
    shutil.rmtree(original_path)


if __name__ == "__main__":

    data_dir, date = get_folder_and_date()

    data_dir_path = Path(data_dir)
    files = [f.name for f in data_dir_path.glob('*.tif') if 'MAX' not in f.name]

    for img in files:  
        try:
            img_name, _ = os.path.splitext(img)
            data = io.imread(f'{data_dir}/{img_name}.tif')

            save_dir = f'./{date}-Work/{img_name}'
            os.makedirs(save_dir, exist_ok=True)
            image_save_path = f'{save_dir}/{img_name}'

            both_filters = False
            filter_string = 'clust_sg' if both_filters else 'sg'

            # Generating Unfiltered Channel 1 (Just Sarcomeres)
            unfiltered_input = separate_sarc_channel_and_save(data, image_save_path)

            # Filtering Images with Custom & SarcGraph
            _, filter_image_path = final_filtering(unfiltered_input, image_save_path, both_filters=both_filters)

            # Getting the Sarcasm Z-Bands Mask
            _, zbands_mask = compute_sarcasm(filepath=filter_image_path, save_path=image_save_path, \
                                            filters_used=filter_string, visualize=False)
            
            # Computing the SarcGraph with the Z-Bands Mask from SarcAsM
            sg, sarcs, zdiscs = compute_sarcgraph(input_image=zbands_mask, save_path=image_save_path, 
                                                filters_used=filter_string, visualize=False)
            save_specific_files(['sarcomeres'], f'{image_save_path}_{filter_string}_sarcgraph_output/', save_dir)

        except:
            print(f"Error processing {img}. Skipping this file.")
            continue
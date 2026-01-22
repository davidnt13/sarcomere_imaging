#!/usr/bin/env python

# Imports
import os
import tkinter as tk
from tkinter import filedialog

from skimage import io
from pathlib import Path
import pandas as pd
import shutil
import traceback

from src.sarcgraph_analysis import compute_sarcgraph, visualize_sarcgraph
from src.sarcasm_analysis import compute_sarcasm
from src.thresholding_pipeline import separate_sarc_channel_and_save, final_filtering
from src.data_analysis import *

def get_folder_and_date():
    result = {
        "folder": None,
        "date": None,
        "save_folder": None,
        "stat_title": None,
        "plot": False,
        "cols": []
    }

    def select_folder(var):
        folder = filedialog.askdirectory(title="Select a Folder")
        var.set(folder)

    def submit():
        result["folder"] = folder_var.get()
        result["date"] = date_var.get()
        result["save_folder"] = save_folder_var.get()
        result["num_img_channels"] = num_channels_var.get()
        result["stat_title"] = title_var.get()
        result["plot"] = plot_choice.get() == "Yes"
        result["cols"] = [col for col, var in checkbox_vars.items() if var.get()]
        window.destroy()

    def toggle_checkboxes():
        if plot_choice.get() == "Yes":
            plot_label.pack(pady=(10, 0))
            for cb in checkboxes:
                cb.pack(anchor="w", padx=20)
            submit_btn.pack(pady=20)
        else:
            plot_label.pack_forget()
            for cb in checkboxes:
                cb.pack_forget()
            submit_btn.pack(pady=20)

    window = tk.Tk()
    window.title("Select Folder and Enter Date")
    window.geometry("500x600")

    # Folder selection
    folder_var = tk.StringVar()
    tk.Label(window, text="Select Folder with Data:").pack(pady=(10, 0))
    tk.Button(window, text="Browse", command=lambda: select_folder(folder_var)).pack()
    tk.Label(window, textvariable=folder_var).pack()

    # Date input
    date_var = tk.StringVar()
    tk.Label(window, text="Enter date (MM-DD) for organization purposes:").pack(pady=(10, 0))
    tk.Entry(window, textvariable=date_var).pack()

    # Save folder selection
    save_folder_var = tk.StringVar()
    tk.Label(window, text="Select Save Destination:").pack(pady=(10, 0))
    tk.Button(window, text="Browse", command=lambda: select_folder(save_folder_var)).pack()
    tk.Label(window, textvariable=save_folder_var).pack()
    
    # Number Channels input
    num_channels_var = tk.IntVar(value=4)
    tk.Label(window, text="Enter Number of Image Channels:").pack(pady=(10, 0))
    tk.Entry(window, textvariable=num_channels_var).pack()

    # Stats sheet title
    title_var = tk.StringVar()
    tk.Label(window, text="Summary Stats Spreadsheet Title:").pack(pady=(10, 0))
    tk.Entry(window, textvariable=title_var).pack()

    # Plotting choice
    plot_choice = tk.StringVar(value="No")
    tk.Label(window, text="Would you like to plot the data?").pack(pady=(15, 0))
    tk.Radiobutton(window, text="Yes", variable=plot_choice, value="Yes", command=toggle_checkboxes).pack()
    tk.Radiobutton(window, text="No", variable=plot_choice, value="No", command=toggle_checkboxes).pack()

    # Label and checkboxes (initially hidden)
    plot_label = tk.Label(window, text="Which columns would you like to plot?")
    checkbox_options = [
        'Sarcomere Count',
        'Mean Length (µm)',
        'Mean Angle (Degrees)',
        'Stdev Angle (Degrees)',
        'Nuclei Count 3D',
        'Myofibril Density (Vol)'
    ]
    checkbox_vars = {col: tk.BooleanVar() for col in checkbox_options}
    checkboxes = [tk.Checkbutton(window, text=col, variable=checkbox_vars[col]) for col in checkbox_options]

    # Submit button (always at the bottom)
    submit_btn = tk.Button(window, text="Submit", command=submit)
    submit_btn.pack(pady=20)

    window.mainloop()

    return (result["folder"], result["date"], result["save_folder"],
            result["num_img_channels"], result["stat_title"], result["plot"], result["cols"])


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

    data_dir, date, save_folder, num_img_channels, stat_title, plot_data, cols_to_plot = get_folder_and_date()

    if date is None: date = 'Output'

    data_dir_path = Path(data_dir)
    files = [f.name for f in data_dir_path.glob('*.tif') if 'MAX' not in f.name]

    for img in files:  
        try:
            img_name, _ = os.path.splitext(img)
            data = io.imread(f'{data_dir}/{img_name}.tif')

            save_dir = f'{save_folder}/{date}-Work/{img_name}'
            os.makedirs(save_dir, exist_ok=True)
            image_save_path = f'{save_dir}/{img_name}'

            both_filters = False
            filter_string = 'clust_sg' if both_filters else 'sg'

            # Generating Unfiltered Channel 1 (Just Sarcomeres)
            unfiltered_input = separate_sarc_channel_and_save(data, image_save_path, num_channels=num_img_channels)

            # Filtering Images with Custom & SarcGraph
            _, filter_image_path = final_filtering(unfiltered_input, image_save_path, both_filters=both_filters)

            # Getting the Sarcasm Z-Bands Mask
            _, zbands_mask = compute_sarcasm(filepath=filter_image_path, save_path=image_save_path, \
                                            filters_used=filter_string, visualize=False)
            
            # Computing the SarcGraph with the Z-Bands Mask from SarcAsM
            sg, sarcs, zdiscs = compute_sarcgraph(input_image=zbands_mask, save_path=image_save_path, 
                                                filters_used=filter_string, visualize=False)
            save_specific_files(['sarcomeres'], f'{image_save_path}_{filter_string}_sarcgraph_output/', save_dir)
            save_specific_files(['sarcomere_mask', 'zbands'], f'{image_save_path}_{filter_string}_filtered/', save_dir)

        except Exception as e:
            print(f"Error processing {img}")
            traceback.print_exc()
            continue
    
    overall_save_directory = f'{save_folder}/{date}-Work'
    data_folders = [os.path.join(overall_save_directory, name) for name in os.listdir(overall_save_directory)
              if os.path.isdir(os.path.join(overall_save_directory, name))]
    print(data_folders)
    data_statistics = create_statistics_df(original_image_dir=data_dir, folders=data_folders, num_channels=num_img_channels)
    data_statistics.to_csv(f'{overall_save_directory}/image_analysis_{stat_title}.csv')

    if plot_data:
        plot_values_from_df(image_info=data_statistics, cols_to_plot=cols_to_plot,
                            save_dir=overall_save_directory, title=f'image_analysis_plots_{stat_title}')

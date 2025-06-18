#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from skimage import io, feature, filters, morphology, draw
from scipy import ndimage
import nibabel as nib

def find_nuclei(image, cleaning_size, thresh_mult, cell_radius, min_distance, 
                show_mask=False, show_template_match=False):
    # Threshold image
    mask = image > filters.threshold_otsu(image) * thresh_mult
    # print(f'Otsu: {filters.threshold_otsu(image)}')
    # print(f'Yen: {filters.threshold_yen(image)}')
    mask = morphology.binary_opening(mask, footprint=morphology.disk(cleaning_size))
    dist = ndimage.distance_transform_edt(mask)

    # Template image
    template = morphology.disk(cell_radius)
    temp_dist = ndimage.distance_transform_edt(template)

    # Match template and find peaks
    result = feature.match_template(dist, temp_dist)
    result_pad = np.zeros_like(dist)
    result_pad[cell_radius:-cell_radius, cell_radius:-cell_radius] = result
    centers = result_pad > 0.1
    result_pad = result_pad*centers
    coords = feature.peak_local_max(result_pad, min_distance=min_distance, threshold_abs=0.1)

    # Figures if asked
    if show_mask:
        plt.figure()
        plt.imshow(image, cmap='binary_r')
        plt.imshow(mask, alpha=0.5, cmap='viridis', vmin=0, vmax=1)
        plt.plot(coords[:,1], coords[:,0], 'r.', ms=2)
        plt.axis('off')
        plt.show()

    if show_template_match:
        plt.figure()
        plt.imshow(result_pad, cmap='binary_r')
        plt.plot(coords[:,1], coords[:,0], 'r.', ms=2)
        plt.axis('off')
        plt.show()

    return coords

def max_project(image):
    image = np.sum(image, axis=0)  
    image = np.clip(image, 0, 1) 
    image = image > 0.5
    return image

def find_nuclei_3d(image, cleaning_size, thresh_mult, cell_radius, min_distance, 
                show_mask=True, show_template_match=True):
    # Threshold image
    mask = image > filters.threshold_otsu(image) * thresh_mult
    # print(f'Otsu: {filters.threshold_otsu(image)}')
    # print(f'Yen: {filters.threshold_yen(image)}')
    mask = morphology.binary_opening(mask, footprint=draw.ellipsoid(cleaning_size//3, cleaning_size, cleaning_size))
    dist = ndimage.distance_transform_edt(mask)

    # Template image
    template = draw.ellipsoid(cell_radius // 4, cell_radius, cell_radius) 
    temp_dist = ndimage.distance_transform_edt(template)

    # Match template and find peaks
    result = feature.match_template(dist, temp_dist, pad_input=True)
    result_pad = np.pad(result, pad_width=cell_radius, mode='constant', constant_values=0)
    centers = result_pad > 0.1
    result_pad = result_pad*centers
    coords = feature.peak_local_max(result_pad, min_distance=min_distance, threshold_abs=0.25)

    # Figures if asked
    if show_mask:
        plt.figure()
        plt.imshow(np.sum(image, axis=0), cmap='binary_r')
        plt.imshow(np.sum(mask, axis=0), alpha=0.5, cmap='viridis', vmin=0, vmax=1)
        plt.plot(coords[:,2], coords[:,1], 'r.', ms=2)
        plt.axis('off')
        plt.show()

    if show_template_match:
        # pixdim = 0.312, 0.312, 2.56
        # rp = np.swapaxes(result_pad[cell_radius:-cell_radius, cell_radius:-cell_radius, cell_radius:-cell_radius], 0, 2)
        # affine = np.eye(4)
        # affine[:3,:3] = np.diag(pixdim)
        # nii_image = nib.Nifti1Image(rp, affine)
        # nib.save(nii_image, 'data/res_pad_0cf_10.nii')
        plt.figure()
        plt.imshow(np.sum(result_pad, axis = 0), cmap='binary_r')
        plt.plot(coords[:,2], coords[:,1], 'r.', ms=2)
        plt.axis('off')
        plt.show()

    return coords

def separate_nuclei_channel_and_save(data, save_path=''):
    if data.ndim != 4:
        data = data[np.newaxis, ...]
        stack_size = 1
    else:
        stack_size = data.shape[0]
        
    channels_stack = []
    for channel_index in range(4):
        channel_data = data[:, :, :, channel_index]
        channel_stack = np.zeros((512, 512, stack_size), dtype=np.uint16)
        for i in range(stack_size):
            channel_stack[:, :, i] = channel_data[i]
        channels_stack.append(channel_stack)

    channel_0 = np.transpose(channels_stack[0], (2, 0, 1))
    
    # io.imsave(f'data/{save_path}_channel_1.tif', channel_1)

    return channel_0

def find_and_count_nuclei(image, min_distance=10, save_file='', three_dim = False, thresh_mult=1.5, cell_radius=15, show_image=True, save_coords = False):
    image = separate_nuclei_channel_and_save(image, '')

    cleaning_size = cell_radius // 3

    if three_dim:
        count = 0
        max_attempts = 4
        attempts = 0
        while count < 10:
            coords = find_nuclei_3d(image, cleaning_size, thresh_mult, cell_radius, min_distance, 
                    show_mask=False, show_template_match=False)
            count = len(coords)
            thresh_mult -= 0.1
        while count < 10 and attempts < max_attempts:
            coords = find_nuclei_3d(image, cleaning_size, thresh_mult, cell_radius, min_distance, 
                    show_mask=False, show_template_match=False)
            thresh_mult -= 0.1
            count = len(coords)

            attempts += 1
        if attempts == max_attempts:
            print("Max attempts reached, skipping this image.")
            return np.nan  # or np.nan
    else:
        image = image[0]
        coords = find_nuclei(image, cleaning_size, thresh_mult, cell_radius, min_distance, 
                show_mask=False, show_template_match=False)

        if show_image:
            plt.figure()
            plt.imshow(image, cmap='binary_r')
            plt.plot(coords[:,1], coords[:,0], 'r.', ms=2)
            plt.axis('off')
            plt.show()

    if save_coords: np.savetxt(save_file, coords)

    return len(coords)

# data_dir = f'../../../Data/cx43_CZI_OG'
# img_name = 'iPSC-CF_fibers_cx43_40x-10'
# image = io.imread(f'{data_dir}/{img_name}.tif')
# image = separate_nuclei_channel_and_save(image, '')
# # image = image[0]

# out_file = f'data/nuclei_centers_{img_name}.txt'
# thresh_mult = 0.5                   # If the thresholding didnt do a good job, you can modify the threshold value with this
# cell_radius = 15                    # Aprox cell radius in pixels
# cleaning_size = cell_radius // 3    # This defines the disk size for cleaning the mask (i.e get rid of isolated pixels and others)
# min_distance_2d = cell_radius // 4
# min_distance_3d = 10                # Minimum distance between two cell nuclei

# coords = find_nuclei_3d(image, cleaning_size, thresh_mult, cell_radius, min_distance_3d, 
#                 show_mask=True, show_template_match=True)         # This two options are to visualize intermediate steps
# print(len(coords))

# np.savetxt(out_file, coords)

# Show result
# img_name = 'iPSC-CF_fibers_cx43_40x-10_MAX'
# max_image = io.imread(f'{data_dir}/{img_name}.tif')
# max_image = separate_nuclei_channel_and_save(max_image, '')

# plt.figure()
# plt.imshow(image, cmap='binary_r')
# plt.plot(coords[:,1], coords[:,0], 'r.', ms=2)
# plt.axis('off')
# plt.show()

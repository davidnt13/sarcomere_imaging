from skimage import io, filters, morphology
from scipy.ndimage import binary_fill_holes, gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib

def separate_channels_and_save(data, save_path=''):

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

    # channel_0 = np.transpose(channels_stack[0], (2, 0, 1))
    
    # io.imsave(f'data/{save_path}_channel_1.tif', channel_1)

    return channels_stack

def find_tissue_vol(image, thresh_mult=0.1):
    channels_stack = separate_channels_and_save(image, '')

    masks = []
    thresh_mults = [0.15, 0.15, 0.1, 0.15]

    for i in range(len(channels_stack)):
        mask = channels_stack[i] > filters.threshold_isodata(channels_stack[i]) * thresh_mults[i]
        print(filters.threshold_isodata(channels_stack[i]))
        masks.append(mask)
        
    combined_mask = np.logical_or.reduce(masks)
    cleaned = morphology.remove_small_objects(combined_mask, min_size=15000)
    combined_mask = gaussian_filter(cleaned.astype(float), sigma=2)

    # cm = np.swapaxes(combined_mask, 0, 2)
    # cm = combined_mask
    # print(cm.shape)
    # pixdim = 0.312, 0.312, 2.56
    # affine = np.eye(4)
    # affine[:3,:3] = np.diag(pixdim)
    # nii_image = nib.Nifti1Image(cm.astype(float), affine)
    # nib.save(nii_image, './10_combined_mask2.nii')

    pixel_vol_count = np.count_nonzero(combined_mask)
    pixel_size = 0.312
    pixel_size_z = 2.56
    vol_microns = pixel_vol_count * (pixel_size ** 2) * pixel_size_z
    
    print(f"Combined Volume: {vol_microns:.2f} µm³")
    return vol_microns

def find_myofibril_volume(myofibril_mask):

    myofibril_mask = np.clip(myofibril_mask, 0, 1)
    myofibril_mask = myofibril_mask > 0.1

    # mm = np.swapaxes(myofibril_mask, 0, 2)
    # mm = np.swapaxes(mm, 0, 1)
    # print(mm.shape)
    # pixdim = 0.312, 0.312, 2.56
    # affine = np.eye(4)
    # affine[:3,:3] = np.diag(pixdim)
    # nii_image = nib.Nifti1Image(mm.astype(float), affine)
    # nib.save(nii_image, './05_myofibril_mask.nii')
    
    myofibril_pixel_area_count = np.count_nonzero(myofibril_mask)

    pixel_size = 0.312
    pixel_size_z = 2.56
    myofibril_vol_microns = myofibril_pixel_area_count * (pixel_size ** 2) * (pixel_size_z)
    
    print(f"Myofibril volume: {myofibril_vol_microns:.2f} µm³")
    return myofibril_vol_microns

def find_myofibril_density_vol(image, myofibril_mask):
    total_vol_in_microns = find_tissue_vol(image)
    myofibril_vol_in_microns = find_myofibril_volume(myofibril_mask)
    density = myofibril_vol_in_microns / total_vol_in_microns

    print(f"Density: {density:.2f} µm³/µm³")
    return density

# data_dir = f'../../../Data/cx43_CZI_OG'
# img_name = 'iPSC-CF_fibers_cx43_40x-10'
# # img_name = 'iPSC-CF_nofibers_cx43_40x-01'
# image = io.imread(f'{data_dir}/{img_name}.tif')
# myofibril_mask = io.imread(f'../Image_Tests/6-11-Work/{img_name}/{img_name}_sg_filtered/sarcomere_mask.tif')

# area_in_microns = find_tissue_vol(image)
# myofibril_area_in_microns = find_myofibril_volume(myofibril_mask)
# density = find_myofibril_density_vol(image, myofibril_mask)

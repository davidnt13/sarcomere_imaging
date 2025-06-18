import numpy as np
from skimage import io
from sklearn.cluster import KMeans
from sarcgraph import SarcGraph

from skimage.morphology import binary_dilation, disk
from skimage.filters import threshold_local, threshold_otsu,\
                            gaussian
from skimage.measure import label, regionprops
from scipy.ndimage import binary_dilation

# ---------- ISOLATING SARCOMERE CHANNEL ----------

def separate_sarc_channel_and_save(data, save_path):
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

    channel_1 = np.transpose(channels_stack[1], (2, 0, 1))
    
    io.imsave(f'{save_path}_channel_1.tif', channel_1)

    return channel_1

# ---------- MAIN PIPELINE INVOLVING SARCGRAPH ----------

# Using My Custom Algorithm, Then SarcGraph
def final_filtering(input_image, save_path, both_filters=False):
    if both_filters: 
        input_image = clustering_filtering(input_image)
        filter_string = 'clust_sg'
    else:
        filter_string = 'sg'
    post_sg = sg_filtering(input_image)
    ret_path = f'{save_path}_{filter_string}_filtered.tif'
    io.imsave(ret_path, post_sg)
    return post_sg, ret_path

# Custom Clustering Filtering Algorithm
def clustering_filtering(input_image):
    _, local_mask = two_stage_threshold(input_image)
    _, expanded_mask = dilate_mask_and_apply(local_mask, input_image, radius=5)
    clustered, _ = cluster_filter_regions(expanded_mask, input_image)
    return clustered

# SarcGraph Filtering
def sg_filtering(input_image):
    sg = SarcGraph(save_output=False, input_type='image')
    ret_frames = sg.filter_frames(input_image)
    return ret_frames

# ---------- HELPER METHODS FOR FILTERING ----------

def two_stage_threshold(image, local_radius=5, local_sigma=1.0, otsu_sigma=1.0, otsu_offset=0.07):
    # Local adaptive threshold
    blurred_local = gaussian(image, sigma=local_sigma)
    local_thresh = threshold_local(blurred_local, block_size=local_radius)
    local_mask = blurred_local > local_thresh
    local_selected = np.zeros_like(image)
    local_selected[local_mask] = image[local_mask]

    # Global Otsu threshold
    blurred_global = gaussian(local_selected, sigma=otsu_sigma)
    otsu_thresh = threshold_otsu(blurred_global) + otsu_offset
    final_mask = blurred_global > otsu_thresh
    final_selected = np.zeros_like(image)
    final_selected[final_mask] = image[final_mask]

    return final_selected, final_mask

def dilate_mask_and_apply(mask, original, radius=1):
    expanded_mask = np.zeros_like(mask)
    selem = disk(radius)
    for i in range(mask.shape[2]):
        expanded_mask[:, :, i] = binary_dilation(mask[:, :, i], structure=selem)
    result = np.zeros_like(original)
    result[expanded_mask] = original[expanded_mask]
    return result, expanded_mask

def cluster_filter_regions(mask, image, n_clusters=2):
    labeled, _ = label(mask, return_num=True)
    regions = regionprops(labeled, intensity_image=image)
    features = np.array([[r.area, r.mean_intensity] for r in regions])
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(features)
    labels = kmeans.labels_
    cluster_idx = np.argmax([features[labels == i][:, 1].mean() for i in range(n_clusters)])
    final_mask = np.zeros_like(mask, dtype=bool)
    for region, label_id in zip(regions, kmeans.labels_):
        if label_id == cluster_idx:
            final_mask[labeled == region.label] = True
    final_selected = np.zeros_like(image)
    final_selected[final_mask] = image[final_mask]
    return final_selected, final_mask
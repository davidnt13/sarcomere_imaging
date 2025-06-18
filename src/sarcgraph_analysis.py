import os
import numpy as np
import matplotlib.pyplot as plt
from sarcgraph import SarcGraph
from skimage import io

PIXEL_SIZE = 0.312
MAX_LEN = 2.2/PIXEL_SIZE

def compute_sarcgraph(input_image, save_path, filters_used, visualize=True):
    if input_image.ndim == 2:
        input_image = input_image[np.newaxis, ...]

    save_dir = f'{save_path}_{filters_used}_sarcgraph_output'
    sg = SarcGraph(save_output=True, input_type='image', output_dir=save_dir, \
                   avg_sarc_length = MAX_LEN/2, max_sarc_length = MAX_LEN)

    sarcomeres, _ = sg.sarcomere_detection(raw_frames=input_image, sigma=1.0)
    zdiscs = sg.zdisc_segmentation(raw_frames=input_image)
    frame_num = 0 if input_image.shape[0] == 1 else 5

    if visualize: visualize_sarcgraph(input_image, sarcomeres, zdiscs, frame_num)

    return sg, sarcomeres, zdiscs

def visualize_sarcgraph(input_image, sarcomeres, zdiscs, frame_num):
    sarcs_frame = sarcomeres[(sarcomeres.frame == frame_num)] #& (sarcomeres.length > (3/0.312))]
    zdiscs_frame = zdiscs[zdiscs.frame == frame_num]
    plt.figure(figsize=(7,7))
    plt.axis('off')
    plt.imshow(input_image[frame_num, :, :], cmap='gray')
    plt.plot(sarcs_frame.y, sarcs_frame.x, 'bo', ms=1, label='Sarcomeres')
    plt.plot(zdiscs_frame.y, zdiscs_frame.x, 'ro', ms=1, label='Z Discs')
    plt.show()
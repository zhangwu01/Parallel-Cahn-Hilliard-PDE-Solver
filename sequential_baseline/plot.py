#!/usr/bin/env python3

import pandas as pd
import glob
import numpy as np
from matplotlib import pyplot as plt
import imageio

file_list = glob.glob('./test_data/raw/*.csv')

for file in file_list:
    
    snap = pd.read_csv(file).to_numpy()

    fig, ax = plt.subplots()
    pcm = ax.imshow(snap , cmap='RdBu' ,vmin=0, vmax=1)
    fig.colorbar(pcm, ax=ax)
    timepoint = float(file.split('_')[2].replace('.csv',''))
    ax.set_title(f'elapsed time = {timepoint}')
    plt.imshow(snap,cmap='RdBu')
    plt.savefig(f'test_data/frames/timepoint_{timepoint}.png')
    plt.close()
    print('done')

total_duration = 5
images = []
timepoints = [float(file.split('_')[2].replace('.csv','')) for file in file_list]
for timepoint in sorted(timepoints):
    images.append(imageio.imread(f'test_data/frames/timepoint_{timepoint}.png'))

imageio.mimsave('./test_data/animation.gif', images, fps=len(images)//total_duration)

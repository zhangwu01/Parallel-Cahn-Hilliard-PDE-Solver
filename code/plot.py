
#!/usr/bin/env python3

import os, glob
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import imageio

frames = []
files = glob.glob('steps/*.bin')
files = sorted(files, key = lambda x: int(x.split("_")[1].replace(".bin", "")))

plt.rcParams["figure.figsize"] = (10,10)
plt.rcParams["figure.dpi"] = 100

for filename in files:    
    with open(filename, 'r') as f:
        dim = np.fromfile(f, dtype=np.int32, count=1)[0]
        snap = np.fromfile(
            f, dtype=np.float64, count=dim * dim
        ).reshape(dim, dim)

    fig, ax = plt.subplots()
    pcm = ax.imshow(snap , cmap='RdBu' ,vmin=0, vmax=1)
    ax.set_title(f'elapsed time = {float(filename.replace("steps/step_", "").replace(".bin", ""))/1000}')
    fig.colorbar(pcm, ax=ax)
    plt.imshow(snap, cmap='RdBu')

    temp_image = 'temp.png'
    plt.savefig(temp_image, dpi=100)
    plt.close()
    
    # Read the temporary image file and append it to the frames list
    frames.append(imageio.imread(temp_image))

    # Delete the temporary image file
    os.remove(temp_image)

gif_file = 'output.gif'
imageio.mimsave(gif_file, frames, fps=len(frames)//5)

print('done')
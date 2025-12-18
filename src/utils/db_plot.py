from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import dbdicom as db


def db_mosaic(series_to_display, pngfile, title="Center slices"):

    # Build list of center slices
    center_slices = []
    for series in tqdm(series_to_display, desc='Reading center images'):
        vol = db.volume(series, verbose=0)
        center_slice = vol.values[:,:,round(vol.shape[-1]/2)]
        center_slices.append(center_slice)

    # Display center slices as mosaic
    n_imgs = len(center_slices)
    aspect_ratio = 16/9
    nrows = int(np.ceil(np.sqrt(n_imgs / aspect_ratio)))
    if nrows==0:
        nrows = 1
    ncols = int(np.ceil(aspect_ratio * nrows))
    if ncols==0:
        ncols=1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(ncols, nrows), dpi=300)

    i=0
    for row in tqdm(ax, desc='Building png'):
        if nrows==1:
            _db_moisaic_tile(row, center_slices, series_to_display, i)
            i += 1 
        else:
            for col in row:
                _db_moisaic_tile(col, center_slices, series_to_display, i)
                i += 1 

    fig.suptitle(title, fontsize=14)
    fig.savefig(pngfile)
    plt.close()


def _db_moisaic_tile(col, center_slices, series_to_display, i):

    col.set_xticklabels([])
    col.set_yticklabels([])
    col.set_aspect('equal')
    col.axis("off")

    if i < len(center_slices):

        # Show center image
        col.imshow(
            center_slices[i].T, 
            cmap='gray', 
            interpolation='none', 
            vmin=0, 
            vmax=np.mean(center_slices[i]) + 2 * np.std(center_slices[i])
        )
        # Add white text with black background in upper-left corner
        patient_id = series_to_display[i][1]
        series_desc = series_to_display[i][-1][0]
        study = series_to_display[i][2][0]
        col.text(
            0.01, 0.99,                   
            f'{patient_id}_{study}\n{series_desc}',   
            color='white',
            fontsize=2,
            ha='left',
            va='top',
            transform=col.transAxes,     # Use axes coordinates
            bbox=dict(facecolor='black', alpha=0.7, boxstyle='round,pad=0.3')
        )
import os
import logging
import argparse

import numpy as np
from tqdm import tqdm
import dbdicom as db
from miblab_plot import mosaic_overlay

from utils.total_segmentator_class_maps import class_map


def run(build):
    run_site(build, 'Patients', site='Leeds')
    return
    for site in ['Leeds', 'Exeter', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site=site)
    run_site(build, 'Controls')


def run_site(build, group, site=None):

    source = 'stage_1_fat_fraction_maps'
    masks = 'stage_3_aligned_masks'
    mosaics = 'stage_4_check_alignment'

    build_pdff = os.path.join(build, 'dixon-pdff')
    if group == 'Controls':
        source_db = os.path.join(build_pdff, source, group) 
        masks_db = os.path.join(build_pdff, masks, group)
        mosaics_db = os.path.join(build_pdff, mosaics, group)
    else:
        source_db = os.path.join(build_pdff, source, group, site)
        masks_db = os.path.join(build_pdff, masks, group, site)
        mosaics_db = os.path.join(build_pdff, mosaics, group, site)

    os.makedirs(mosaics_db, exist_ok=True)

    kidney_class_map = {1: "kidney_left", 2: "kidney_right"}
    total_mr_class_map = class_map['total_mr']
    tissue_types_mr_class_map = class_map['tissue_types_mr']

    # Get all water series in the source database
    all_source_series = db.series(source_db)
    all_pdff_series = [s for s in all_source_series if s[3][0][-len('pdff'):]=='pdff']

    # Loop over the masks
    for pdff_series in tqdm(all_pdff_series, 'Displaying masks..'):

        try:

            # Patient and study IDs
            patient_id = pdff_series[1]
            study_desc = pdff_series[2][0]
            series_desc = pdff_series[3][0]
            sequence = series_desc[:-len('_pdff')]
            
            # No need to align to post-contrast
            if 'post_contrast' in series_desc:
                continue

            # Skip if file exists
            prefix = f'{patient_id}_{study_desc}_{sequence}'
            kidney_png = os.path.join(mosaics_db, f'{prefix}_kidney_masks.png')
            total_mr_png = os.path.join(mosaics_db, f'{prefix}_total_mr.png')
            tissue_types_mr_png = os.path.join(mosaics_db, f'{prefix}_tissue_types_mr.png')
            if os.path.exists(kidney_png):
                continue
            
            # Get mask and fat seriesx
            mask_study = [masks_db, patient_id, (study_desc, 0)]
            series_kidney_mask = mask_study + [(f'{sequence}_kidney_masks', 0)]
            series_total_mr = mask_study + [(f'{sequence}_total_mr', 0)]
            series_tissue_types_mr = mask_study + [(f'{sequence}_tissue_types_mr', 0)]
            
            # Read mask volume
            vol_pdff = db.volume(pdff_series)
            vol_kidney_mask = db.volume(series_kidney_mask)#.slice_like(vol_pdff)
            vol_total_mr = db.volume(series_total_mr)#.slice_like(vol_pdff)
            vol_tissue_types_mr = db.volume(series_tissue_types_mr)#.slice_like(vol_pdff)

            # Build mosaics
            rois = {roi: (vol_kidney_mask.values==idx).astype(np.int16) for idx, roi in kidney_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, kidney_png, vmin=0, vmax=0.3, margin=[16,16,2], opacity=0.4)

            rois = {roi: (vol_total_mr.values==idx).astype(np.int16) for idx, roi in total_mr_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, total_mr_png, vmin=0, vmax=0.3, opacity=0.4)

            rois = {roi: (vol_tissue_types_mr.values==idx).astype(np.int16) for idx, roi in tissue_types_mr_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, tissue_types_mr_png, vmin=0, vmax=0.3, opacity=0.4)

        except:
            logging.exception(f"{patient_id}_{study_desc}_{sequence}: error building mosaic")


if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    dixon_build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(dixon_build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(dixon_build, 'stage_4_check_alignment.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build)

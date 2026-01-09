import os
import logging
import argparse

import numpy as np
from tqdm import tqdm
import dbdicom as db
from miblab_plot import mosaic_overlay, mosaic_checkerboard

from utils.total_segmentator_class_maps import class_map


def run(build):
    for site in ['Leeds', 'Exeter', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site=site)
    run_site(build, 'Controls')
    

def run_site(build, group, site=None):

    build_dixon = os.path.join(build, 'dixon', 'stage_5_clean_dixon_data')
    build_pdff = os.path.join(build, 'dixon-pdff')

    pdff = 'stage_1_fat_fraction_maps'
    masks = 'stage_3_aligned_masks'
    mosaics = 'stage_4_check_alignment'

    if group == 'Controls':
        dixon_db = os.path.join(build_dixon, group) 
        pdff_db = os.path.join(build_pdff, pdff, group) 
        masks_db = os.path.join(build_pdff, masks, group)
        mosaics_db = os.path.join(build_pdff, mosaics, group)
    else:
        dixon_db = os.path.join(build_dixon, group, site)
        pdff_db = os.path.join(build_pdff, pdff, group, site)
        masks_db = os.path.join(build_pdff, masks, group, site)
        mosaics_db = os.path.join(build_pdff, mosaics, group, site)

    os.makedirs(mosaics_db, exist_ok=True)

    kidney_class_map = {1: "kidney_left", 2: "kidney_right"}
    total_mr_class_map = class_map['total_mr']
    tissue_types_mr_class_map = class_map['tissue_types_mr']

    # Get all water series in the source database
    all_mask_series = db.series(masks_db)
    all_kidney_mask_series = [s for s in all_mask_series if s[3][0][-len('kidney_masks'):]=='kidney_masks']

    # Loop over the masks
    for series_kidney_mask in tqdm(all_kidney_mask_series, 'Displaying masks..'):

        try:

            # Patient and study IDs
            patient_id = series_kidney_mask[1]
            study_desc = series_kidney_mask[2][0]
            series_desc = series_kidney_mask[3][0]
            sequence = series_desc[:-len('_kidney_masks')]

            # Skip if file exists
            prefix = f'{patient_id}_{study_desc}_{sequence}'
            kidney_png = os.path.join(mosaics_db, f'{prefix}_kidney_masks.png')
            total_mr_png = os.path.join(mosaics_db, f'{prefix}_total_mr.png')
            tissue_types_mr_png = os.path.join(mosaics_db, f'{prefix}_tissue_types_mr.png')
            water_png = os.path.join(mosaics_db, f'{prefix}_water.png')
            if os.path.exists(water_png):
                continue
            
            # Get source studies
            dixon_study = [dixon_db, patient_id, (study_desc, 0)]
            pdff_study = [pdff_db, patient_id, (study_desc, 0)]
            mask_study = [masks_db, patient_id, (study_desc, 0)]

            # Get mask and fat series
            series_fixed_water = dixon_study + [(f'{sequence}_water', 0)]
            series_coreg_water = mask_study + [(f'{sequence}_water_ref_coreg', 0)]
            series_total_mr = mask_study + [(f'{sequence}_total_mr', 0)]
            series_tissue_types_mr = mask_study + [(f'{sequence}_tissue_types_mr', 0)]
            series_pdff = pdff_study + [(f'{sequence}_pdff', 0)]
            
            # Read mask volume
            vol_pdff = db.volume(series_pdff, verbose=0)
            vol_fixed_water = db.volume(series_fixed_water, verbose=0)
            vol_coreg_water = db.volume(series_coreg_water, verbose=0)
            vol_kidney_mask = db.volume(series_kidney_mask, verbose=0)
            vol_total_mr = db.volume(series_total_mr, verbose=0)
            vol_tissue_types_mr = db.volume(series_tissue_types_mr, verbose=0)

            # Build mosaics
            mosaic_checkerboard(vol_fixed_water.values, vol_coreg_water.values, water_png, normalize=True)

            rois = {roi: (vol_kidney_mask.values==idx).astype(np.int16) for idx, roi in kidney_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, kidney_png, vmin=0, vmax=0.3, margin=[16,16,2], opacity=0.4)

            rois = {roi: (vol_total_mr.values==idx).astype(np.int16) for idx, roi in total_mr_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, total_mr_png, vmin=0, vmax=0.3, opacity=0.7)

            rois = {roi: (vol_tissue_types_mr.values==idx).astype(np.int16) for idx, roi in tissue_types_mr_class_map.items()}
            mosaic_overlay(vol_pdff.values, rois, tissue_types_mr_png, vmin=0, vmax=0.3, opacity=0.7)

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


import os
import logging
import argparse

import numpy as np
# from scipy.ndimage import label
from tqdm import tqdm
import vreg
import dbdicom as db
# from mdreg.ants import coreg, transform
from mdreg.elastix import coreg, transform
# from mdreg.skimage import coreg, transform

import utils.data


def run(build):
    for site in ['Leeds', 'Exeter', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site=site)
    run_site(build, 'Controls')
    

def run_site(build, group, site=None):

    # Define database paths
    source = 'stage_5_clean_dixon_data'
    results = 'stage_3_aligned_masks'
    pdff_results = 'stage_1_fat_fraction_maps'

    if group == 'Controls':
        dixon_data = os.path.join(build, 'dixon', source, group) 
        pdff_db = os.path.join(build, 'dixon-pdff', pdff_results, group)
        coreg_db = os.path.join(build, 'dixon-pdff', results, group)
        kidney_masks_db = os.path.join(build, 'kidneyvol', 'stage_3_edit', group)
        total_masks_db = os.path.join(build, 'totseg', 'stage_1_segment', group)
    else:
        dixon_data = os.path.join(build, 'dixon', source, group, site)
        pdff_db = os.path.join(build, 'dixon-pdff', pdff_results, group, site)
        coreg_db = os.path.join(build, 'dixon-pdff', results, group, site)
        kidney_masks_db = os.path.join(build, 'kidneyvol', 'stage_3_edit', group, site)
        total_masks_db = os.path.join(build, 'totseg', 'stage_1_segment', group, site)

    os.makedirs(coreg_db, exist_ok=True)

    # Get all precontrast water series in the source database
    all_source_series = db.series(dixon_data)
    all_water_series = [s for s in all_source_series if s[3][0][-5:]=='water']
    all_water_series = [s for s in all_water_series if 'post_contrast' not in s[3][0]]

    # List existing series so they can be skipped in the loop
    existing_result_series = db.series(coreg_db)

    # List of reference dixon series
    src = os.path.dirname(os.path.abspath(__file__))
    record = utils.data.dixon_record(src)

    # Loop over the fat series in the source database
    for fixed_series_water in tqdm(all_water_series, desc='Aligning masks'):

        # Patient and study IDs
        patient_id = fixed_series_water[1]
        study_desc = fixed_series_water[2][0]
        series_desc = fixed_series_water[3][0]
        sequence = series_desc[:-len('_water')]

        # Define the results series
        coreg_study = [coreg_db, patient_id, (study_desc, 0)]
        coreg_series_water = coreg_study + [(f'{sequence}_water_ref_coreg', 0)]
        coreg_series_kidney_masks = coreg_study + [(f'{sequence}_kidney_masks', 0)]
        coreg_series_total_mr = coreg_study + [(f'{sequence}_total_mr', 0)]
        coreg_series_tissue_types_mr = coreg_study + [(f'{sequence}_tissue_types_mr', 0)]
        
        # Skip computation if the results are already there
        if coreg_series_water in existing_result_series:
            continue

        # Get the original masks on the reference series
        kidney_mask_study = [kidney_masks_db, patient_id, (study_desc, 0)]
        totseg_mask_study = [total_masks_db, patient_id, (study_desc, 0)]
        moving_series_kidney_masks = kidney_mask_study + [(f'kidney_masks', 0)]
        moving_series_total_mr = totseg_mask_study + [(f'total_mr', 0)]
        moving_series_tissue_types_mr = totseg_mask_study + [(f'tissue_types_mr', 0)]

        # Get the refence sequence and moving water series
        source_study = [dixon_data, patient_id, (study_desc, 0)]
        reference_sequence = utils.data.dixon_series_desc(record, patient_id, study_desc)
        moving_series_water = source_study + [(f'{reference_sequence}_water', 0)]

        # Get the pdff series
        fixed_series_pdff = [pdff_db, patient_id, (study_desc, 0), (f'{sequence}_pdff', 0)]

        # Any error happens in the computations, log and move on to the next
        try:

            # The fixed sequence does not need coregistering
            if sequence == reference_sequence:
                db.copy(moving_series_water, coreg_series_water)
                db.copy(moving_series_kidney_masks, coreg_series_kidney_masks)
                db.copy(moving_series_total_mr, coreg_series_total_mr)
                db.copy(moving_series_tissue_types_mr, coreg_series_tissue_types_mr)
                continue
            
            # Coregister fixed and moving water volumes
            fixed_vol_water = db.volume(fixed_series_water, verbose=0)
            moving_vol_water = db.volume(moving_series_water, verbose=0)
            coreg_values_water, transfo = coreg(
                moving_vol_water.slice_like(fixed_vol_water).values.astype(float), 
                fixed_vol_water.values.astype(float), 
                # Elastix options
                spacing=fixed_vol_water.spacing.tolist(),
                FinalGridSpacingInPhysicalUnits="16",
                # skimage options
                # attachment=30,
            )

            # Read the other volumes and reslice to target geometry
            fixed_vol_pdff = db.volume(fixed_series_pdff, verbose=0)
            moving_vol_kidney_masks = db.volume(moving_series_kidney_masks, verbose=0)
            moving_vol_total_mr = db.volume(moving_series_total_mr, verbose=0)
            moving_vol_tissue_types_mr = db.volume(moving_series_tissue_types_mr, verbose=0)

            # Apply the same transformation to the other series
            coreg_values_kidney_masks = transform_label(moving_vol_kidney_masks, transfo, fixed_vol_pdff, 'total_mr')
            coreg_values_total_mr = transform_label(moving_vol_total_mr, transfo, fixed_vol_pdff, 'total_mr')
            coreg_values_tissue_types_mr = transform_label(moving_vol_tissue_types_mr, transfo, fixed_vol_pdff, 'tissue_types_mr')
            
            # Save results 
            db.write_volume((coreg_values_water, fixed_vol_water.affine), coreg_series_water, ref=moving_series_water, verbose=0)
            db.write_volume((coreg_values_kidney_masks, fixed_vol_water.affine), coreg_series_kidney_masks, ref=moving_series_kidney_masks, verbose=0)
            db.write_volume((coreg_values_total_mr, fixed_vol_water.affine), coreg_series_total_mr, ref=moving_series_total_mr, verbose=0)
            db.write_volume((coreg_values_tissue_types_mr, fixed_vol_water.affine), coreg_series_tissue_types_mr, ref=moving_series_tissue_types_mr, verbose=0)

        except:

            logging.exception(f"Error aligning {fixed_series_water}")



def transform_label(vol_label, transfo, vol_pdff, task):

    # The labels have been interpolated
    scl = 1000
    values = np.unique(vol_label.values)
    transformed = np.zeros(vol_pdff.shape, dtype=vol_label.values.dtype)
    for v in values:

        # First create mask and reslice if needed
        mask_values = (vol_label.values == v).astype(np.float32)
        mask_vol = vreg.volume(mask_values, vol_label.affine).slice_like(vol_pdff)
        mask = mask_vol.values.astype(np.float32)

        # Transform and turn into binary
        mask_transformed = transform(scl * mask, transfo, spacing=vol_pdff.spacing.tolist())
        
        # Remove pixels with unphysiological pdff values (usually at the boundaries)
        min_pdff = 0.0
        max_pdff = 0.4
        if task=='total_mr':
            if v in [5,19]: # liver, vertebrae
                max_pdff = 1.0
        elif task=='tissue_types_mr':
            if v in [1,2]:  # visceral and subcutaneous fat
                min_pdff = 0.0
                max_pdff = np.inf
        mask_transformed[vol_pdff.values < min_pdff] = 0
        mask_transformed[vol_pdff.values > max_pdff] = 0

        transformed[mask_transformed > scl / 2] = v

    return transformed


# def snap_to_values(arr, values):
#     return values[np.abs(arr[..., None] - values).argmin(axis=-1)]


# def keep_largest_cc_per_label(lbl):
#     """
#     Keep only the largest connected component for each non-zero label.

#     Parameters
#     ----------
#     lbl : ndarray
#         Label image (int), any dimensionality.

#     Returns
#     -------
#     out : ndarray
#         Label image with only largest CC per label retained.
#     """
#     out = np.zeros_like(lbl)

#     for lab in np.unique(lbl):
#         if lab == 0:
#             continue  # skip background

#         mask = lbl == lab
#         cc, ncc = label(mask)

#         if ncc == 0:
#             continue

#         sizes = np.bincount(cc.ravel())
#         sizes[0] = 0  # ignore background
#         largest = sizes.argmax()

#         out[cc == largest] = lab

#     return out





if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    dixon_build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(dixon_build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(dixon_build, 'stage_3_aligned_masks.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build)



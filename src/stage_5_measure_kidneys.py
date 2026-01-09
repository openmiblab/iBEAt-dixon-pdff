import os
import logging
import shutil
from pathlib import Path
from datetime import date
import argparse

from tqdm import tqdm
import numpy as np
import dbdicom as db
import vreg
import pydmr

from utils import radiomics



def run(build):
    for site in ['Exeter', 'Leeds', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site)
    run_site(build, 'Controls')
    combine(build)


def run_site(build, group, site=None):

    build_pdff = os.path.join(build, 'dixon-pdff')
    results = os.path.join(build_pdff, 'stage_5_measure_kidneys')

    pdff = 'stage_1_fat_fraction_maps'
    masks = 'stage_3_aligned_masks'

    if group == 'Controls':
        pdff_db = os.path.join(build_pdff, pdff, group) 
        masks_db = os.path.join(build_pdff, masks, group)
        results_db = os.path.join(results, group)         
    else:
        pdff_db = os.path.join(build_pdff, pdff, group, site)
        masks_db = os.path.join(build_pdff, masks, group, site)
        results_db = os.path.join(results, group, site)
    
    all_pdff_series = [s for s in db.series(pdff_db) if s[3][0][-len('pdff'):]=='pdff']
    all_mask_series = db.series(masks_db)

    for pdff_series in tqdm(all_pdff_series, desc='Extracting metrics'):

        patient, study, sequence = pdff_series[1], pdff_series[2][0], pdff_series[3][0][:-len('_pdff')]

        # Get input series and check completeness
        mask_series = [masks_db, patient, (study, 0), (f'{sequence}_kidney_masks', 0)]
        wat_series = [pdff_db, patient, (study, 0), (f"{sequence}_water", 0)]
        fat_series = [pdff_db, patient, (study, 0), (f"{sequence}_fat", 0)]
        if mask_series not in all_mask_series:
            logging.info(f"Patient {patient}, Study {study}, Sequence {sequence} - no kidney masks available.")
            continue

        # Skip if results exist
        dmr_file = os.path.join(results_db, f"{patient}-{study}.dmr.zip")
        if os.path.exists(dmr_file):
            continue
  
        # Read data
        pdff_vol = db.volume(pdff_series, verbose=0)
        label_arr = db.volume(mask_series, verbose=0).values
        fat_arr = db.volume(fat_series, verbose=0).values   
        wat_arr = db.volume(wat_series, verbose=0).values   
        
        # Compute results
        dmr = {'data':{}, 'pars':{}}

        # Single kidneys
        class_map = {1: "kidney_left", 2: "kidney_right"}
        for idx, roi in class_map.items():
            _compute_roi_vals(label_arr == idx, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr)

        # Both kidneys
        roi = 'kidneys_both'
        _compute_roi_vals(label_arr != 0, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr)

        # Write results to file
        pydmr.write(dmr_file, dmr)


def _compute_roi_vals(mask_bool, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr):

    mask = mask_bool.astype(np.float32)
    if np.sum(mask) == 0:
        return
    
    try:
        roi_vol = vreg.volume(mask, pdff_vol.affine)
        results = radiomics.texture_features(roi_vol, pdff_vol, roi, 'dixon_pdff', binWidth=0.001)
    except:
        logging.exception(f"Patient {patient}, Study {study}, Sequence {sequence} - error computing {roi}")
        return

    dmr['data'] = dmr['data'] | {p:v[1:] for p, v in results.items()}
    dmr['pars'] = dmr['pars'] | {(patient, study, p): v[0] for p, v in results.items()}

    # Compute ROI-based
    fat = fat_arr[mask_bool].mean()
    tot = (wat_arr + fat_arr)[mask_bool].mean()
    pdff = fat / tot if tot > 0 else 0
    p = f"{roi}-dixon_pdff-roi"
    dmr['data'][p] = [p, '', 'float']
    dmr['pars'][(patient, study, p)] = pdff


def combine(build):
    """
    Concatenate all dmri files in a folder into a single dmr file. 
    Create long and wide format csv's for export.
    """
    results = os.path.join(build, 'dixon-pdff', 'stage_5_measure_kidneys')
    for group in ['Controls', 'Patients']:

        # Combine all dmr files into one
        folder = os.path.join(results, group) 
        folder = Path(folder)
        dmr_files = list(folder.rglob("*.dmr.zip"))  
        if dmr_files == []:
            continue
        dmr_files = [str(f) for f in dmr_files]
        dmr_file = os.path.join(results, f'{group}_kidneys_pdff')
        pydmr.concat(dmr_files, dmr_file)



if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    dixon_build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(dixon_build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(dixon_build, 'stage_5_measure_kidneys.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build)

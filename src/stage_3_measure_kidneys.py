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
import pandas as pd

from utils import data, radiomics



def run(build):
    run_site(build, 'Controls')
    for site in ['Exeter', 'Leeds', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site)


def run_site(build_path, group, site=None):

    maskpath = os.path.join(build_path, 'kidneyvol', 'stage_3_edit')
    datapath = os.path.join(build_path, 'dixon-pdff', 'stage_1_fat_fraction_maps')
    resultspath = os.path.join(build_path, 'dixon-pdff', 'stage_3_measure')

    if group == 'Controls':
        datapath = os.path.join(datapath, "Controls") 
        maskpath = os.path.join(maskpath, "Controls")
        resultspath = os.path.join(resultspath, "Controls")         
    else:
        datapath = os.path.join(datapath, "Patients", site) 
        maskpath = os.path.join(maskpath, "Patients", site)
        resultspath = os.path.join(resultspath, "Patients", site)

    src = os.path.dirname(os.path.abspath(__file__))
    record = data.dixon_record(src)
    
    all_pdff_series = [s for s in db.series(datapath) if s[3][0][-4:]=='pdff']
    all_mask_series = db.series(maskpath)

    for pdff_series in tqdm(all_pdff_series, desc='Extracting metrics'):

        patient, study, series = pdff_series[1], pdff_series[2][0], pdff_series[3][0]

        # Skip if it is not the right sequence
        sequence = series[:-5]  # remove '_pdff' suffix
        selected_sequence = data.dixon_series_desc(record, patient, study)
        if sequence != selected_sequence:
            continue

        # Get mask series
        mask_series = [maskpath, patient, (study, 0), ('kidney_masks', 0)]
        if mask_series not in all_mask_series:
            logging.error(f"Patient {patient}, Study {study}, Sequence {sequence} - no kidney masks available.")
            continue

        # Skip if results exist
        dmr_file = os.path.join(resultspath, f"{patient}-{study}.dmr.zip")
        if os.path.exists(dmr_file):
            continue

        # Get in an opposed-phase series
        wat_series = [datapath, patient, (study, 0), (f"{sequence}_water", 0)]
        fat_series = [datapath, patient, (study, 0), (f"{sequence}_fat", 0)]
        
        # Read data
        fat_arr = db.volume(fat_series).values   
        wat_arr = db.volume(wat_series).values   
        pdff_vol = db.volume(pdff_series)           
        mask_arr = db.volume(mask_series).values
        
        # Compute results
        dmr = {'data':{}, 'pars':{}}

        # Single kidneys
        class_map = {1: "kidney_left", 2: "kidney_right"}
        for idx, roi in class_map.items():
            _compute_roi_vals(mask_arr == idx, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr)

        # Both kidneys
        roi = 'kidneys_both'
        _compute_roi_vals(mask_arr != 0, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr)

        # Write results to file
        pydmr.write(dmr_file, dmr)

    #concat_patient(measurepath, site)


def _compute_roi_vals(mask_bool, fat_arr, wat_arr, pdff_vol, patient, study, sequence, roi, dmr):

    mask = mask_bool.astype(np.float32)
    if np.sum(mask) == 0:
        return
    
    try:
        roi_vol = vreg.volume(mask, pdff_vol.affine)
        results = radiomics.texture_features(roi_vol, pdff_vol, roi, 'dixon_pdff', binWidth=0.01)
    except Exception as e:
        logging.error(f"Patient {patient}, Study {study}, Sequence {sequence} - error computing {roi}: {e}")
        return

    dmr['data'] = dmr['data'] | {p:v[1:] for p, v in results.items()}
    dmr['pars'] = dmr['pars'] | {(patient, study, p): v[0] for p, v in results.items()}

    # Compute ROI-based
    fat = fat_arr[mask_bool].mean()
    tot = (wat_arr + fat_arr)[mask_bool].mean()
    pdff = fat / tot if tot > 0 else 0
    p = f"{roi}-dixon_pdff"
    dmr['data'][p] = ['Dixon proton density fat fraction', '', 'float']
    dmr['pars'][(patient, study, p)] = pdff


def combine(build_path):
    """
    Concatenate all dmri files in a folder into a single dmr file. 
    Create long and wide format csv's for export.
    """
    measurepath = os.path.join(build_path, 'kidneyvol', 'stage_5_measure')
    for group in ['Controls', 'Patients']:

        # Combine all dmr files into one
        folder = os.path.join(measurepath, group) 
        folder = Path(folder)
        dmr_files = list(folder.rglob("*.dmr.zip"))  
        if dmr_files == []:
            continue
        dmr_files = [str(f) for f in dmr_files]
        dmr_file = os.path.join(measurepath, f'{group}_kidneyvol')
        pydmr.concat(dmr_files, dmr_file)



if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    os.makedirs(args.build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(args.build, 'dixon-pdff', 'stage_3_measure.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build)

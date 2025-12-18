
import os
import logging
import argparse
import subprocess
import time

from tqdm import tqdm
import numpy as np
import dbdicom as db
from miblab_dl import fatwater
#import utils.data


def run(build, cache=None):
    compute(build, 'Controls', cache=cache)
    for site in ['Exeter', 'Leeds', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        compute(build, 'Patients', site=site, cache=cache)

def compute(build, group, site=None, cache=None):

    # Define site paths
    if group == 'Controls':
        datapath = os.path.join(build, 'dixon', 'stage_5_clean_dixon_data', group) 
        resultspath = os.path.join(build, 'dixon-pdff', 'stage_2_fat_fraction_maps', group)
    else:
        datapath = os.path.join(build, 'dixon', 'stage_5_clean_dixon_data', group, site)
        resultspath = os.path.join(build, 'dixon-pdff', 'stage_2_fat_fraction_maps', group, site)
    os.makedirs(resultspath, exist_ok=True)

    # Get all out_phase series
    series = db.series(datapath)
    series_out_phase = [s for s in series if s[3][0][-9:]=='out_phase']

    # List existing series so they can be skipped in the loop
    existing_series = db.series(resultspath)

    # TODO: This requires a selected PREcontrast series
    # # List of selected dixon series
    # src = os.path.dirname(os.path.abspath(__file__))
    # record = utils.data.dixon_record(src)

    # Loop over the out_phase series
    for series_op in tqdm(series_out_phase, desc='Computing fat-water maps'):

        # Patient and study IDs
        patient = series_op[1]
        study = series_op[2][0]

        # Remove '_out_phase' suffix from series description
        sequence = series_op[3][0][:-10] 

        # # Skip if it is not the right sequence
        # selected_sequence = utils.data.dixon_series_desc(record, patient, study)
        # if sequence != selected_sequence:
        #     continue

        # Define output series
        fat_series = [resultspath, patient, (study, 0), (f'{sequence}_fat', 0)]
        water_series = [resultspath, patient, (study, 0), (f'{sequence}_water', 0)]
        pdff_series = [resultspath, patient, (study, 0), (f'{sequence}_pdff', 0)]

        # Skip if the output already exists
        if pdff_series in existing_series:
            continue

        # Get in_phase series
        series_ip = series_op[:3] + [(f'{sequence}_in_phase', 0)]

        # Read the in-phase and opposed-phase data
        try:
            op = db.volume(series_op)
            ip = db.volume(series_ip)
        except:
            logging.exception(f"Patient {patient} - error reading images for {sequence}")
            continue

        # Read the echo times for in-phase and opposed-phase data
        try:
            pars_op = db.unique(['EchoTime', 'FlipAngle', 'RepetitionTime'], series_op)
            pars_ip = db.unique(['EchoTime', 'FlipAngle', 'RepetitionTime'], series_ip)
        except:
            logging.exception(f"Patient {patient} - error reading echo times for {sequence}")
            continue     

        # Compute fat and water images with T2* and T1 correction
        try:
            fat, water = fatwater(
                op.values, 
                ip.values, 
                te_o=pars_op['EchoTime'][0], # msec
                te_i=pars_ip['EchoTime'][0], # msec
                t2s_w=30, 
                t2s_f=15, 
                tr=pars_op['RepetitionTime'][0], # msec
                fa=pars_op['FlipAngle'][0],  # deg
                t1_w=1400, # Muscle, kidney cortex at 3T
                t1_f=350,
                cache=cache,
            )
        except:
            logging.exception(f"Error predicting fat-water maps for {patient} {sequence}")
            continue

        # Compute pdff
        pdff = np.divide(fat, fat + water, out=np.zeros_like(fat, dtype=float), where=(fat + water) != 0)
        
        # Save results in DICOM
        db.write_volume((fat, op.affine), fat_series, ref=series_op, Manufacturer='miblab')
        db.write_volume((water, op.affine), water_series, ref=series_op, Manufacturer='miblab')
        db.write_volume((pdff, op.affine), pdff_series, ref=series_op, Manufacturer='miblab')




if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    parser.add_argument("--cache", type=str, default=None, help="Cache folder")
    args = parser.parse_args()

    pdff_build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(pdff_build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(pdff_build, 'stage_2_fat_fraction_maps.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build, args.cache)



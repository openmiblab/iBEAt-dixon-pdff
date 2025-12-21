import os
import argparse
import logging

import dbdicom as db
from utils.db_plot import db_mosaic


def run(build):
    for site in ['Leeds', 'Exeter', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        run_site(build, 'Patients', site)
    run_site(build, 'Controls')


def run_site(build, group, site=None):

    input = 'stage_1_fat_fraction_maps'
    output = 'stage_2_check_pdff'

    build = os.path.join(build, 'dixon-pdff')
    resultspath = os.path.join(build, output)
    os.makedirs(resultspath, exist_ok=True)
    summary = os.path.join(resultspath, 'summary.txt')

    if group == 'Controls':
        datapath = os.path.join(build, input, group) 
        csv = os.path.join(resultspath, f'{group}.csv')
        png = os.path.join(resultspath, f'{group}.png')
    else:
        datapath = os.path.join(build, input, group, site)
        csv = os.path.join(resultspath, f'{group}_{site}.csv')
        png = os.path.join(resultspath, f'{group}_{site}.png')
    
    # Build csv

    db.to_csv(datapath, csv)
    _write_summary(datapath, summary, group, site)
    _build_mosaic(datapath, png)



def _write_summary(datapath, summary, group, site):

    patients = db.patients(datapath)
    studies = db.studies(datapath)
    series = db.series(datapath)

    nr_studies = {}
    for patient in patients:
        n = len([s for s in studies if s[:2]==patient[:2]])
        if n in nr_studies:
            nr_studies[n] += 1
        else:
            nr_studies[n] = 1

    nr_series = {}
    for study in studies:
        n = len([s for s in series if s[:3]==study[:3]])
        if n in nr_series:
            nr_series[n] += 1
        else:
            nr_series[n] = 1

    with open(summary, 'a') as f:

        f.write('\n')
        f.write(f"{group} ({'all sites' if site is None else site})\n")

        f.write('\n')
        f.write(f"  {len(patients)} patients with {len(studies)} studies\n")
        f.write('\n')
        for n in sorted(nr_studies.keys()):
            f.write(f"    {nr_studies[n]} with {n} studies\n")

        f.write('\n')
        f.write(f"  {len(studies)} studies with {len(series)} series\n")
        f.write('\n')
        for n in sorted(nr_series.keys()):
            f.write(f"    {nr_series[n]} with {n} series\n")




def _build_mosaic(datapath, png):

    series = db.series(datapath)
    series_desc = [s[-1][0] for s in series]
    series_pdff = [s for i, s in enumerate(series) if series_desc[i][-4:]=='pdff']

    if os.path.exists(png):
        return
    
    db_mosaic(series_pdff, png, title="PDFF maps", vmin=0, vmax=0.3)



if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(build, 'stage_2_check_pdff.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.build)
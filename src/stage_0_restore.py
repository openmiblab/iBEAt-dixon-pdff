import os
import logging
import argparse

import dbdicom as db


def run(archive, build):
    # dixons(archive, build, 'Controls', site=None)
    # kidney_segmentations(archive, build, 'Controls', site=None)
    totseg(archive, build, 'Controls', site=None)
    for site in ['Exeter', 'Leeds', 'Bari', 'Bordeaux', 'Sheffield', 'Turku']:
        # dixons(archive, build, 'Patients', site)
        kidney_segmentations(archive, build, 'Patients', site)
        totseg(archive, build, 'Patients', site)


def dixons(archivepath, buildpath, group, site=None):
    local_buildpath = os.path.join(buildpath, 'dixon', 'stage_5_clean_dixon_data')
    local_archivepath = os.path.join(archivepath, "dixon", "stage_5_clean_dixon_data")
    if group == 'Controls':
        sitebuildpath = os.path.join(local_buildpath, 'Controls')
        sitearchivepath = os.path.join(local_archivepath, 'Controls')
    else:
        sitebuildpath = os.path.join(local_buildpath, 'Patients', site)
        sitearchivepath = os.path.join(local_archivepath, 'Patients', site)
    db.restore(sitearchivepath, sitebuildpath)


def kidney_segmentations(archivepath, buildpath, group, site=None):
    local_buildpath = os.path.join(buildpath, 'kidneyvol', 'stage_3_edit')
    local_archivepath = os.path.join(archivepath, "kidneyvol", "stage_3_edit")
    if group == 'Controls':
        sitebuildpath = os.path.join(local_buildpath, 'Controls')
        sitearchivepath = os.path.join(local_archivepath, 'Controls')
    else:
        sitebuildpath = os.path.join(local_buildpath, 'Patients', site)
        sitearchivepath = os.path.join(local_archivepath, 'Patients', site)
    db.restore(sitearchivepath, sitebuildpath)


def totseg(archivepath, buildpath, group, site=None):
    local_buildpath = os.path.join(buildpath, 'totseg', 'stage_1_segment')
    local_archivepath = os.path.join(archivepath, 'totseg', 'stage_1_segment')
    if group == 'Controls':
        sitebuildpath = os.path.join(local_buildpath, 'Controls')
        sitearchivepath = os.path.join(local_archivepath, 'Controls')
    else:
        sitebuildpath = os.path.join(local_buildpath, 'Patients', site)
        sitearchivepath = os.path.join(local_archivepath, 'Patients', site)
    db.restore(sitearchivepath, sitebuildpath)



if __name__=='__main__':

    ARCHIVE = r"G:\Shared drives\iBEAt_Build"
    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=str, default=ARCHIVE, help="Archive folder")
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    os.makedirs(args.build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(args.build, 'stage_0_restore.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    run(args.archive, args.build)




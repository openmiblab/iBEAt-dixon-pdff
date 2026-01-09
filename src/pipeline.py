import argparse
import os
import logging

# import stage_1_fat_fraction_maps
# import stage_2_check_pdff
import stage_3_align_masks
import stage_4_check_alignment
import stage_5_measure_kidneys


if __name__=='__main__':

    BUILD = r'C:\Users\md1spsx\Documents\Data\iBEAt_Build'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=str, default=BUILD, help="Build folder")
    args = parser.parse_args()

    dixon_build = os.path.join(args.build, 'dixon-pdff')
    os.makedirs(dixon_build, exist_ok=True)

    logging.basicConfig(
        filename=os.path.join(dixon_build, 'pipeline.log'),
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # stage_1_fat_fraction_maps.run(args.build)
    # stage_2_check_pdff.run(args.build)
    # stage_3_align_masks.run(args.build)
    # stage_4_check_alignment.run(args.build)
    stage_5_measure_kidneys.run(args.build)
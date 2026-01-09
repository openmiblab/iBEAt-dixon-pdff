# import os
# import logging
# import tempfile

import numpy as np
# import scipy.ndimage as ndi
# from scipy.interpolate import interpn
# import skimage
# import vreg
# from radiomics import featureextractor # Replace with numpyradiomics
import numpyradiomics as npr


def texture_features(roi_vol, img_vol, roi, img, binWidth=25):

    result = npr.firstorder(
        img_vol.values, 
        roi_vol.values, 
        voxelVolume=np.prod(roi_vol.spacing), 
        binWidth=binWidth,
    )
    unit = npr.firstorder_units('', 'mm^3')

    # Format return value
    rval = {}
    for p, v in result.items():
        name = f"{roi}-{img}-histogram-{p}"
        vals = [float(v), name, unit[p], 'float']
        rval[name] = vals
    return rval


# def texture_features(roi_vol, img_vol, roi, img, binWidth=25):

#     with tempfile.TemporaryDirectory() as tmp:
#         roi_file = os.path.join(tmp, 'roi.nii.gz')
#         img_file = os.path.join(tmp, 'img.nii.gz')
#         print('radiomics texture ', roi)
#         # Downsample large ROIs to avoid memory over
#         # TODO: Not enough for large regions - still RAM issue
#         roi_vol_box = roi_vol
#         img_vol_box = img_vol
#         # roi_vol_box = roi_vol.crop(mask=roi_vol) # some edits in vreg. Check
#         # img_vol_box = img_vol.crop(mask=roi_vol)
#         # roi_vol_box = roi_vol_box.resample(5.0)
#         # img_vol_box = img_vol_box.resample(5.0)
#         vreg.write_nifti(roi_vol_box, roi_file)
#         vreg.write_nifti(img_vol_box, img_file)
#         # All features for water
#         extractor = featureextractor.RadiomicsFeatureExtractor(binWdith=binWidth)
#         extractor.disableAllFeatures()
#         # TODO: try without first order
#         classes = ['firstorder', 'glcm', 'glrlm', 'glszm', 'gldm', 'ngtdm'] # glcm seems to fail a lot
#         # classes = ['firstorder', 'glrlm', 'glszm', 'gldm', 'ngtdm']
#         for cl in classes:
#             extractor.enableFeatureClassByName(cl)
#         # extractor.enableImageTypeByName('Wavelet')
#         # extractor.enableImageTypeByName('LoG', {'sigma': [1.0, 1.0, 1.0]}) 
#         # extractor.enableImageTypeByName('Gradient')
#         result = extractor.execute(img_file, roi_file)

#     # Format return value
#     rval = {}
#     for p, v in result.items():
#         if p[:8] == 'original':
#             name = roi + '-' + img + '-' 
#             for cl in classes:
#                 if cl in p:
#                     if cl=='firstorder':
#                         name += p.replace(f'original_{cl}_', f'histogram-')
#                     else:
#                         name += p.replace(f'original_{cl}_', f'texture-')
#                     break
#             vals = [float(v), name, 'unit', 'float']
#             rval[name] = vals
#     return rval
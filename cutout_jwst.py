
"""
Python functions for making cutout images of COSMOS-Web at given sky locations.
These functions assume that the COSMOS-Web images are stored in your local machine.
If not, please download the images from the COSMOS-Web website (https://exchg.calet.org/cosmosweb-public/).

Before using this function, please make sure that the path to the COSMOS-Web images is correctly set.

Created by    : Sang Hyeok Im (tkdgur0117@gmail.com)
Last modified : 19th October 2024

"""



import glob

import numpy as np
import pandas as pd

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord


def make_COSMOWeb_cutouts(coords, jwst_filter, size=3*u.arcsec, dr=0.5, pixScale=30):

    """
    This function makes cutouts of both SCIENCE and ERROR images COSMOS-Web survey.
    Check also the pd.DataFrame object that contains the information of the objects and their cutouts.

    Parameters
    ----------
        coords: astropy.coordinates.SkyCoord
            Sky positions of the objects
        jwst_filter: str
            The filter name of the JWST image
            Should be one of the ['f115w', 'f150w', 'f277w', 'f444w'] for COSMOS-Web
        size: astropy.units.Quantity
            The size of the cutout images (default: 3 arcsec)
        dr: float
            The DR number of the public COSMOS-Web images (default: 0.5)
            Should be one of the [0.2, 0.5] for COSMOS-Web (as of Oct. 2024)
        pixScale: int
            The pixel scale of the public COSMOS-Web images in [mas] 
            Should be one of the [30, 60] for COSMOS-Web
    
    Returns
    -------
        sciCutouts: list of :class:'~astropy.nddata.Cutout2D'
            List of cutouts of COSMOS-Web SCIENCE images
        errCutouts: list of :class:'~astropy.nddata.Cutout2D'
            List of cutouts of COSMOS-Web ERROR images
        objData: pandas.DataFrame
            DataFrame containing the information of the objects
                * RA: Right Ascension of the objects [deg]
                * DEC: Declination of the objects [deg]
                * inJWST: Whether the object is in the public COSMOS-Web images or not
                * mosaic: The name of the COSMOS-Web mosaic that contains the object
                            ('none' if not in the images)
                * warnNaN: True if the SCIENCE cutout contains only NaN values
    """

    jwstPwd    = f"/md/imsang/ODIN_obs_LSS/JWST_images/COSMOS-Web/dr{dr}/NIRCam/Apr23/"
    jwstFnames = f"mosaic_nircam_{jwst_filter}_COSMOS-Web_{pixScale}mas_A*_v0_5_i2d.fits"
    jwstFlist  = sorted(glob.glob(jwstPwd + jwstFnames))
    
    mosaicStr = find_COSMOSWeb_mosaic(coords, jwstFlist)

    inJWST  = np.invert(mosaicStr == 'none')
    warnNaN = np.full(len(coords), False)

    sciCutouts = []
    errCutouts = []

    for i in range(len(coords)):
        if (inJWST[i] == False):
            sciCutout = None
            errCutout = None
        else:
            jwstFname = f"mosaic_nircam_{jwst_filter}_COSMOS-Web_{pixScale}mas_{mosaicStr[i]}_v0_5_i2d.fits"
            with fits.open(jwstPwd + jwstFname, mode='denywrite') as hdul:
                sciData, sciHdr = hdul[1].data, hdul[1].header
                errData, errHdr = hdul[2].data, hdul[2].header
                wcs             = WCS(sciHdr, naxis=2, fix=False)
            
            sciCutout = Cutout2D(sciData, position=coords[i], size=size, wcs=wcs)
            errCutout = Cutout2D(errData, position=coords[i], size=size, wcs=wcs)

            # << Check for NaN >> #
            # Even if the object is in the JWST image, 
            # the SCIENCE cutout may contain only NaN values for some cases (e.g., outskirts of the mosaic).
            nNaN = np.sum(np.isnan(sciCutout.data))
            nPix = sciCutout.data.size
            if (nNaN == nPix):
                warnNaN[i] = True

        sciCutouts.append(sciCutout)
        errCutouts.append(errCutout)


    objData = pd.DataFrame({'RA (deg)': coords.ra.deg, 
                            'DEC (deg)': coords.dec.deg, 
                            'inJWST': inJWST, 
                            'mosaic': mosaicStr, 
                            'warnNaN': warnNaN})

    return sciCutouts, errCutouts, objData



def find_COSMOSWeb_mosaic(coords, jwstFlist):
    
    """
    The whole public COSMOS-Web image is divide into several mosaics (due to its size).
    This function finds the mosaic in which the objects are located.

    Parameters
    ----------
        coords: astropy.coordinates.SkyCoord
            Sky positions of the objects
        jwstFlist: list
            List of the filenames of the COSMOS-Web images 
    
    Returns
    -------
        mosaicStr: numpy.ndarray
            Array containing the name of the mosaic that contains each object
            ('none' if the object is not in any COSMOS-Web mosaics)
    """

    mosaicStr = np.full(len(coords), 'none')

    for i in range(len(jwstFlist)):

        fileName = np.char.split(jwstFlist[i], '/').tolist()[-1]
        mosaic   = np.char.split(fileName, sep='_').tolist()[5]

        with fits.open(jwstFlist[i], mode='denywrite') as hdul:
            sciData = hdul[1].data
            sciHdr  = hdul[1].header
            wcs = WCS(sciHdr, naxis=2, fix=False)
        
        obj_Xpix, obj_Ypix = wcs.world_to_pixel(coords)

        Xpix_mask  = (0 <= obj_Xpix) & (obj_Xpix <= sciData.shape[1])
        Ypix_mask  = (0 <= obj_Ypix) & (obj_Ypix <= sciData.shape[0])
        pixLocMask = Xpix_mask & Ypix_mask

        mosaicStr[pixLocMask] = mosaic

    return mosaicStr



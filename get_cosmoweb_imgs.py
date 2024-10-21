
"""
A Python function for downloading COSMOS-Web images from the website (https://exchg.calet.org/cosmosweb-public/).
The whole image of COSMOS-Web in each filter is divided into 10 different mosaic images (A1 to A10).

Created by    : Sang Hyeok Im (tkdgur0117@gmail.com)
Last modified : 21th October 2024

"""

import os
import numpy as np

import requests
from bs4 import BeautifulSoup


def get_cosmoweb_imgs(dr, sav_pwd, jwst_filter=None, pixelScale=None, mosaic=None, make_shell=False, print_progress=True):
    
    """
    Get the list of the file names for the COSMOS-Web images with the given JWST filter, pixel scale, and mosaic.
    By default, all the images with the desired DR number will be downloaded (i.e., for all available filters, pixel sacles, mosaics).
    The user can also give the JWST filter, pixel scale, and mosaic to be downloaded.

    Parameters
    ----------
        dr: float
            The DR number of the public COSMOS-Web images
            Should be one of the [0.2, 0.5] for COSMOS-Web (as of Oct. 2024)
        sav_pwd: str
            The path where the COSMOS-Web images are saved
        jwst_filter: str
            The filter name of the JWST image
            Should be one of the ['f115w', 'f150w', 'f277w', 'f444w'] for COSMOS-Web
            If None, the files for all the filters will be downloaded.
        pixelScale: int
            The pixel scale of the public COSMOS-Web images in [mas] 
            Should be one of the [30, 60] for COSMOS-Web
            If None, the files for all the pixel scales will be downloaded.
        mosaic: str
            The mosaic name of the COSMOS-Web images
            Should be one of the ['A1', 'A2', ... , 'A9', 'A10'] for COSMOS-Web
            If None, the files for all the mosaics will be downloaded.
        make_shell: bool (default: False)
            If True, a shell script file will be created to download the images.
            If False, the images will be downloaded directly.
        print_progress: bool (default: True)
            If True, the progress of will be printed.
    
    Returns
    -------
    
    """

    webb_url = f"https://exchg.calet.org/cosmosweb-public/DR{dr}/NIRCam/Apr23/"

    response  = requests.get(webb_url)
    soup      = BeautifulSoup(response.text, "html5lib")
    fileNames = np.array(soup.select("a"))[4:]


    ## Find the files with the given JWST filter, pixel scale, and mosaic ##
    if jwst_filter is not None:
        filterMask = np.char.find(fileNames, jwst_filter) != -1
        fileNames  = fileNames[filterMask]
    
    if pixelScale is not None:
        pixScaleMask = np.char.find(fileNames, f"{pixelScale}mas") != -1
        fileNames    = fileNames[pixScaleMask]
    
    if mosaic is not None:
        mosaicMask = np.char.find(fileNames, mosaic) != -1
        fileNames  = fileNames[mosaicMask]
    

    ## Download the Files or Make a Shell File ##
    if not os.path.exists(sav_pwd):
        os.makedirs(sav_pwd)

    if not make_shell:
        for i in range(len(fileNames)):
            command = f"wget {webb_url}{fileNames[i]} -P {sav_pwd}"
            if print_progress: print(f"Start downloading.. {fileNames[i]}")        
            os.system(command)
            if print_progress: print(f"-> Completed \n")
    else:
        shellFile = f"{sav_pwd}/download_COSMOSWeb_imgs.sh"
        with open(shellFile, 'w') as shFile:
            shFile.write("#!/bin/bash\n")
            for i in range(len(fileNames)):
                shFile.write(f"wget {webb_url}{fileNames[i]} -P {sav_pwd}\n")
        print(f" The shell script file is created: {shellFile}")
    
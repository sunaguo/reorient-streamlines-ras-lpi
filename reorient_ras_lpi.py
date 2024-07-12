#!/usr/bin/env python3

"""
sunaguo 2023.09.03
to fix ras/lpi inconsistency between subj in MDLFang, MDLFslp, and Uncinate

sunaguo 2024.07.12
Note: 
xyz follow dwi volume orientation.
(orientation: small - big)
x: R - L
y: S - I (vertical)
z: P - A (horizontal)
"""
## TODO: define major orientation for all tracts for consistency

def relabel(fdir):

    import glob
    import numpy as np

    import nibabel as nib

    def get_center(img):
        d = img.get_fdata()
        xs, ys, zs = np.where(d)
        return {"x": xs.mean(), "y": ys.mean(), "z": zs.mean()}
    
    def swap_fnames(rasfn, lpifn):
        print("swapping ras/lpi")

        tempfn = f"{fdir}/temp"
        print("reveresing fnames")
        os.rename(rasfn, tempfn)
        os.rename(lpifn, rasfn)
        os.rename(tempfn, lpifn)

    
    ## define the major orientation of the tracts
    tract_orientations = {
        "MDLFang":"y", 
        "MDLFspl":"y", 
        "Uncinate":"z", 
        "Aslant":"y"
    }

    for tname, orientation in tract_orientations.items():
        print(f"====={tname}=====")
        fnames = {}
        imgs = {}

        for fn in glob.glob(f"{fdir}/*"):
            if ((f"left{tname}" in fn) and ("RAS" in fn)):
                # print(fn)
                imgs["lras"] = nib.load(fn)
                fnames["lras"] = fn
            elif ((f"left{tname}" in fn) and ("LPI" in fn)):
                # print(fn)
                imgs["llpi"] = nib.load(fn)
                fnames["llpi"] = fn
            elif ((f"right{tname}" in fn) and ("RAS" in fn)):
                # print(fn)
                imgs["rras"] = nib.load(fn)
                fnames["rras"] = fn
            elif ((f"right{tname}" in fn) and ("LPI" in fn)):
                # print(fn)
                imgs["rlpi"] = nib.load(fn)
                fnames["rlpi"] = fn
        if len(imgs) < 4: 
            raise Exception(f"missing data for {tname} (only loaded {list(imgs.keys())}). Aborted.")

        centers = {tlab: get_center(img) for tlab, img in imgs.items()}
        for lab, cs in centers.items():
            print(lab, cs)

        ## sanity check: left x > right x
        print("=== Sanity check: left x > right x for all left/right pairs ===")
        for llab in ["lras", "llpi"]:
            for rlab in ["rras", "rlpi"]:
                if centers[llab]["x"] < centers[rlab]["x"]:
                    raise Exception(f"{llab} {rlab} failed: data in different coordinates. Aborted.")
        print("passed")

        ## check each major orientation
        print(f"=== Checking ras {orientation} > lpi {orientation} ===")
        if orientation == "y": 
            ## y: S - I (vertical) --> ras < lpi
            
            for raslab, lpilab in [["lras", "llpi"], ["rras", "rlpi"]]:
                print(raslab, lpilab)
                if centers[raslab][orientation] < centers[lpilab][orientation]:
                    print("True")
                else: 
                    print("swapping ras/lpi ")
                    rasfn = fnames[raslab]
                    lpifn = fnames[lpilab]
                    swap_fnames(rasfn, lpifn)
            
        if orientation == "z": 
            ## z: P - A (horizontal) --> ras > lpi
            
            for raslab, lpilab in [["lras", "llpi"], ["rras", "rlpi"]]:
                print(raslab, lpilab)
                if centers[raslab][orientation] > centers[lpilab][orientation]:
                    print("True")
                else: 
                    print("swapping ras/lpi ")
                    rasfn = fnames[raslab]
                    lpifn = fnames[lpilab]
                    swap_fnames(rasfn, lpifn)
    return True
    
	
if __name__ == '__main__':

    import os, shutil
    import json

    # load config
    with open('config.json','r') as config_f:
        config = json.load(config_f)

    # define paths and variables
    roisdir = config['rois']

    # make output directory
    out_path = './rois_relabeled/'
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    out_path += "rois/"
    shutil.copytree(roisdir, out_path)

    res = relabel(out_path)

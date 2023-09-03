#!/usr/bin/env python3

"""
sunaguo 2023.09.03
to fix ras/lpi inconsistency between subj in MDLFang, MDLFslp, and Uncinate
"""


def relabel(fdir):

    import glob
    import numpy as np

    import nibabel as nib

    def get_center(img):
        d = img.get_fdata()
        xs, ys, zs = np.where(d)
        return xs.mean(), ys.mean(), zs.mean()

    for tname in ["MDLFang", "MDLFspl", "Uncinate"]:
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
        print("=== Sanity check: left x > right x ===")
        for llab in ["lras", "llpi"]:
            for rlab in ["rras", "rlpi"]:
                if centers[llab][0] < centers[rlab][0]:
                    raise Exception(f"{llab} {rlab} failed: data in different coordinates. Aborted.")
        print("passed")

        ## check y: ras y > lpi y
        print(f"=== Checking ras y > lpi y & swapping ras/lpi ===")
        for raslab, lpilab in [["lras", "llpi"], ["rras", "rlpi"]]:
            print(raslab, lpilab)
            if centers[raslab][1] > centers[lpilab][1]:
                print("True")
            else: 
                rasfn = fnames[raslab]
                lpifn = fnames[lpilab]
                tempfn = f"{fdir}/temp"
                print("reveresing fnames")
                os.rename(rasfn, tempfn)
                os.rename(lpifn, rasfn)
                os.rename(tempfn, lpifn)

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
    out_path = './rois_relabeled'
    if not os.path.exists(out_path):
        shutil.copytree(roisdir, out_path)

    res = relabel(out_path)

#!/usr/bin/env python3

"""
sunaguo 2023.08.21
adapted from https://github.com/brainlife/app-tractanalysisprofiles/blob/dipy-1.0/tractAnalysisProfilesDipy.py
"""

def generateRAScentroid(centroid_cluster):

	import numpy as np

	[x_diff, y_diff, z_diff] = np.diff([centroid_cluster[:,0],centroid_cluster[:,1],centroid_cluster[:,2]])
	[x_diff, y_diff, z_diff] = [np.append(x_diff,[centroid_cluster[0,0] - centroid_cluster[-1,0]]),np.append(y_diff,[centroid_cluster[0,1] - centroid_cluster[-1,1]]),np.append(z_diff,[centroid_cluster[0,2] - centroid_cluster[-1,2]])]
	[x_diff_max, y_diff_max, z_diff_max] = [np.max(np.absolute(x_diff)),np.max(np.absolute(y_diff)),np.max(np.absolute(z_diff))]

	max_dim = np.where([x_diff_max,y_diff_max,z_diff_max] == np.max([x_diff_max,y_diff_max,z_diff_max]))[0][0]

	if centroid_cluster[0,max_dim] < centroid_cluster[-1,max_dim]:
		centroid_cluster = np.flip(centroid_cluster,0)

	return centroid_cluster

def reorient(reference_anat_path,streamlines_path,classification_path,n_points,new_streamlines_outpath):

    import numpy as np
    from copy import deepcopy as dc

    import nibabel as nib
    import scipy.io as sio
    from dipy.io.streamline import load_tractogram, save_tractogram
    from dipy.io.stateful_tractogram import StatefulTractogram as sft
    import dipy.tracking.streamline as dts
    from dipy.segment.clustering import QuickBundles
    from dipy.segment.metric import AveragePointwiseEuclideanMetric
    from dipy.segment.metric import ResampleFeature
    # ## udpated for dipy>=1.5.0: https://dipy.org/documentation/1.7.0/api_changes/
    # from dipy.segment.featurespeed import ResampleFeature

    ## ===== load stuff for init
    # load reference anatomy (dwi)
    print('loading reference anatomy')
    ref_anat = nib.load(reference_anat_path)

    # load tractogram
    print('loading tractogram')
    streamlines = load_tractogram(streamlines_path,ref_anat)

    # load classification
    print('loading classification')
    classification = sio.loadmat(classification_path)

    # extract names and indices from classification
    names = list(np.ravel(list(classification['classification'][0]['names'][0][0])))
    indices = classification['classification'][0]['index'][0][0]

    # define metrics to use for reorienting streamlines using quickbundles 
    feature = ResampleFeature(nb_points=n_points)
    metric = AveragePointwiseEuclideanMetric(feature)

    ## to save new/potentially flipped streamlines
    new_streamlines = dc(streamlines.streamlines)

    for tii, tname in enumerate(names): 
        tract_indices = np.where(indices==(tii+1))[0]
        fg = streamlines.streamlines[tract_indices]

        ## ===== reorient for RAS/LPI
        # reorient streamlines to match orientation of first streamline. then compute centroid
        fg_oriented = dts.orient_by_streamline(fg,fg[0])

        # run quickbundles, find centroid, and reorient streamlines
        qb = QuickBundles(np.inf,metric=metric)
        tract_cluster = qb.cluster(fg_oriented)
        centroid_cluster = tract_cluster.centroids[0]
        centroid_cluster_ras = generateRAScentroid(centroid_cluster)
        oriented_tract = dts.orient_by_streamline(fg,centroid_cluster_ras)
	
        ## store
        new_streamlines[tract_indices] = oriented_tract

    ## ===== save things back to orignal format
    new_tck = sft(new_streamlines, space=streamlines.space, reference=ref_anat)
    save_tractogram(new_tck, new_streamlines_outpath)

def main():
	
	import os
	import json

	# load config
	with open('config.json','r') as config_f:
		config = json.load(config_f)

	# make output directory
	out_path = './tract_oriented'
	if not os.path.exists(out_path):
		os.mkdir(out_path)
	out_path += '/tract.tck'

	# define paths and variables
	# subjectID = config['_inputs'][0]['meta']['subject']
	streamlines_path = config['track']
	classification_path = config['classification']
	reference_anat_path = config['dwi']
	n_points = 100

	reorient(reference_anat_path,streamlines_path,classification_path,n_points,out_path)


if __name__ == '__main__':
    main()

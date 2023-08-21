"""
sunaguo 2023.08.21
adapted from https://github.com/brainlife/app-tractanalysisprofiles/blob/dipy-1.0/tractAnalysisProfilesDipy.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import sklearn
import dipy
import nibabel as nib
import scipy.io as sio
from dipy.io.streamline import load_tractogram
import dipy.tracking.streamline as dts
from dipy.segment.clustering import QuickBundles
from dipy.segment.metric import AveragePointwiseEuclideanMetric
## udpated in dipy 1.5.0: https://dipy.org/documentation/1.7.0/api_changes/
from dipy.segment.featurespeed import ResampleFeature

def generateRAScentroid(centroid_cluster):

	import numpy as np

	[x_diff, y_diff, z_diff] = np.diff([centroid_cluster[:,0],centroid_cluster[:,1],centroid_cluster[:,2]])
	[x_diff, y_diff, z_diff] = [np.append(x_diff,[centroid_cluster[0,0] - centroid_cluster[-1,0]]),np.append(y_diff,[centroid_cluster[0,1] - centroid_cluster[-1,1]]),np.append(z_diff,[centroid_cluster[0,2] - centroid_cluster[-1,2]])]
	[x_diff_max, y_diff_max, z_diff_max] = [np.max(np.absolute(x_diff)),np.max(np.absolute(y_diff)),np.max(np.absolute(z_diff))]

	max_dim = np.where([x_diff_max,y_diff_max,z_diff_max] == np.max([x_diff_max,y_diff_max,z_diff_max]))[0][0]

	if centroid_cluster[0,max_dim] < centroid_cluster[-1,max_dim]:
		centroid_cluster = np.flip(centroid_cluster,0)

	return centroid_cluster

def reorient():

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

    # tracts
    tracts = {}
    images_json = {}
    images_json['images'] = []

    # error messages
    failed_tracts = np.array([])
    failed_tracts_lows = np.array([])

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

        ## ===== save thing back to orignal format

    pass



if __name__ == '__main__':

	import os,sys
	import json

	# load config
	with open('config.json','r') as config_f:
		config = json.load(config_f)

	# make output directories
	if not os.path.exists('./images'):
		os.mkdir('./images')
	if not os.path.exists('./profiles'):
		os.mkdir('./profiles')
	if not os.path.exists('./tractmeasures'):
		os.mkdir('./tractmeasures')

	# define paths and variables
	subjectID = config['_inputs'][0]['meta']['subject']
	streamlines_path = config['track']
	classification_path = config['classification']
	reference_anat_path = config['dwi']
	n_points = config['num_nodes']

	out_path = './'

	# loop through measures to create measures_path
	df_measures = []
	
	# dti
	fa = config['fa']
	if os.path.isfile(fa):
		df_measures = df_measures+['ad','fa','md','rd']

	# dki
	ga = config['ga']
	if os.path.isfile(ga):
		df_measures = df_measures+['ga','ak','mk','rk']

	# noddi
	if 'odi' in config:
		odi = config['odi']
		if os.path.isfile(odi):
			df_measures = df_measures+['ndi','odi','isovf']
		
	# myelin-map
	if 'myelin' in config:
		myelin = config['myelin']
		if os.path.isfile(myelin):
			df_measures = df_measures+['myelin']
		
	# qmri
	if 'T1' in config:
		qmri = ["T1","R1","M0","PD","MTV","VIP","SIR","WF"]
		for i in qmri:
			test_met = config[i]
			if os.path.isfile(test_met):
				df_measures = df_measures+[i]
	
	measure_path = [ config[f] for f in df_measures ]

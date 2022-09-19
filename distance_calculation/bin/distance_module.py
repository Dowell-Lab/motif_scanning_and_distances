######################################### Imports #########################################
import os
from os import path
import sys
import glob
import time
import datetime
import multiprocessing
from functools import partial
import numpy as np
import pandas as pd
######################################### Distance Main #########################################
def run_distance_calculation(verbose, outdir, sample, annotation, window, cpus, seq_type, pre_scan):
    if verbose==True:
        print('--------------Pulling in Annotation and Getting List of Motifs---------------')
    md_dirs(verbose, outdir=outdir, seq_type=seq_type)
    tf_list=get_scanned_tfs(verbose=verbose, outdir=outdir, sample=sample, 
                            pre_scan=pre_scan, seq_type=seq_type)
    annotation_df, chr_list = read_annotation(verbose=verbose, sample=sample, 
                                                 outdir=outdir, window=window, 
                                              pre_scan=pre_scan,
                                              annotation=annotation, seq_type=seq_type)
    if verbose == True: 
        print('--------------Beginning Distance Calculation---------------')
        print('Initializing ' + str(cpus) + ' threads to calculate distances from mu.')
        start_time = int(time.time())
        print('Start time: %s' % str(datetime.datetime.now()))

    for tf in tf_list:
        motif_df = read_motif(verbose=verbose, outdir=outdir, pre_scan=pre_scan,
                              chr_list=chr_list, seq_type=seq_type, tf=tf)  ##handle prescan here
        pool = multiprocessing.Pool(cpus)
        distance_dfs = pool.map(partial(distance_calculation, 
                               inputs=[window, annotation_df, motif_df]), chr_list)
        pool.close()
        pool.join()
        distance_df = pd.concat(distance_dfs, axis=0)
        distance_df.to_csv(outdir + '/distances/' + seq_type + '/' + tf +'_distances.txt',
                          sep='\t', index=False)
    if verbose == True:
        print("---------Distance Calculation Complete----------")
        stop_time = int(time.time())
        print('Stop time: %s' % str(datetime.datetime.now()))
        print('Total Run time :', (stop_time-start_time)/60, ' minutes')
######################################### Motif Distance Calculation Functions #########################################
def md_dirs(verbose, outdir, seq_type):
    if (path.exists(outdir + '/distances') == False):
        os.system('mkdir -p ' + outdir + '/distances') 
    try:
        os.system('mkdir -p ' + outdir + '/distances/' + seq_type)
    except OSError:
        print ('Creation of the directory %s failed' % outdir + '/distances/' + seq_type)
        sys.exit(1)
    else:
        if verbose == True:
            print(outdir + '/distances/' + seq_type + ' exists.')

def get_scanned_tfs(verbose, outdir, sample, pre_scan, seq_type):
    if seq_type == 'experimental' and pre_scan is not None:
        tf_motif_path = pre_scan + '/*'
    else:
        tf_motif_path = outdir + '/motifs/' + seq_type + '/*'
        
    tf_list = []
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    if verbose == True:
        print("Processing %s motif files in %s" % (motif_count, tf_motif_path))
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        filename_no_path=filename_no_path.replace('.sorted.bed','')
        tf_list.append(filename_no_path)
    if verbose == True:
        print('There are ' + str(len(tf_list)) + ' motifs with hits in this dataset.')
    return tf_list

def read_annotation(verbose, sample, outdir, window, pre_scan, annotation, seq_type):
    if verbose == True:
        print('Reading annotation file and centering regions...')
    if seq_type == 'experimental' and pre_scan is not None:
        annotation_df = prescan_annotation_windower(bed=annotation, outdir=outdir, 
                                                    sample=sample, window=window, seq_type=seq_type)
    else:
        annotation_file = outdir + '/annotations/' + sample + '_' + seq_type + '_window.bed'
        annotation_df = pd.read_csv(annotation_file, sep='\t', header=None)
    
    annotation_df.columns = ['chr', 'start', 'stop', 'region_id']
    annotation_df['center'] = round((annotation_df['stop'] + annotation_df['start'])/2)
    annotation_df.center = annotation_df.center.astype(int)
    chr_list = list(annotation_df['chr'].unique())
    return annotation_df, chr_list

def prescan_annotation_windower(bed, outdir, sample, window, seq_type):
    print(bed)
    df = pd.read_csv(bed, sep='\t', header=None)
    df.columns = ['chr', 'start', 'stop', 'region_id']
    df['start_new'] = df.apply(lambda x: round((x['start'] + x['stop'])/2), axis=1)
    df['stop_new'] = df.apply(lambda x: x['start_new'] + 1, axis = 1)

    ##the -1500 position from 'origin'
    df['start'] = df.apply(lambda x: x['start_new'] - int(window), axis=1)
    df['stop'] = df.apply(lambda x: x['stop_new'] + int(window), axis=1)
    df=df[['chr', 'start', 'stop', 'region_id']]
    
    os.system('mkdir -p ' + outdir + '/annotations/')
    df.to_csv(outdir + '/annotations/' + sample + '_' + seq_type + '_prescan_windowed.bed', header=None, index=False, sep='\t')
    return df
    
def read_motif(verbose, outdir, pre_scan, chr_list, seq_type, tf):
    if seq_type == 'experimental' and pre_scan is not None:
        motif_file = pre_scan + '/' + tf + '.sorted.bed'
    else:      
        motif_file = outdir + '/motifs/' + seq_type + '/' + tf + '.sorted.bed'
       
    motif_df = pd.read_csv(motif_file, sep='\t', header=None)
    motif_df.columns = ['chr','start', 'stop', 'score', 'strand', 'motif_id']
    motif_df['center'] = round((motif_df['stop'] + motif_df['start'])/2)
    motif_df = motif_df.sort_values(by='score', ascending=False)
    motif_df = motif_df.drop_duplicates(subset=['center'], keep='first')
    motif_df=motif_df[motif_df['chr'].isin(chr_list)]
    motif_df = motif_df.sort_values(by='center', ascending=True)
    motif_df.center = motif_df.center.astype(int)
    return motif_df

def distance_calculation(chr, inputs):
    window, annotation_df, motif_df = inputs
    
    single_chr_annotation_df = annotation_df[annotation_df['chr'] == chr]
    annotation_array = single_chr_annotation_df['center'].to_numpy(dtype=int)
    
    single_chr_motif_df = motif_df[motif_df['chr'] == chr]
    motif_array = single_chr_motif_df['center'].to_numpy(dtype=int)

    all_distances = annotation_array[None, :] - motif_array[:, None]
    hit = np.where((all_distances >= -1*window) & (all_distances <= window), True, False)
    
    annotation_shape = annotation_array.shape[0]
    motif_shape = motif_array.shape[0]
    motif_positions=np.repeat(motif_array[None, :],annotation_shape , axis=1).reshape(motif_shape, annotation_shape)
    annotation_positions = np.repeat(annotation_array[None,:],motif_shape , axis=0).reshape(motif_shape, annotation_shape)
    
    distances = all_distances[hit]
    motif_hits = motif_positions[hit]
    annotation_hits= annotation_positions[hit]
    distance_array = np.column_stack((annotation_hits[:,None],motif_hits[:,None], distances[:,None]))
    
    distance_df = pd.DataFrame(distance_array, 
                               columns=['annotation_center', 'motif_center', 'distance'])
    distance_df=distance_df.merge(single_chr_annotation_df, right_on='center', left_on='annotation_center')
    distance_df=distance_df[['region_id', 'distance', 'motif_center']]
    distance_df=distance_df.merge(single_chr_motif_df, right_on='center', left_on='motif_center')
    distance_df['ABS_distance'] = abs(distance_df['distance'])
    distance_df['distance_rank']=distance_df.groupby(['region_id'])['ABS_distance'].rank(method ='dense')
    distance_df['quality_rank']=distance_df.groupby(['region_id'])['score'].rank(method ='dense')
    distance_df=distance_df[['region_id', 'motif_id', 'distance', 'distance_rank', 'quality_rank']]
    return distance_df
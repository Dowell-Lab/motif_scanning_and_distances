######################################### Imports #########################################
import argparse
from argparse import RawTextHelpFormatter
import os
from os import path
import sys
import time
import datetime
import glob
import pandas as pd
import numpy as np

######################################### Arguments #########################################
if __name__ == "__main__":
    
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')    

    parser = argparse.ArgumentParser(description = 'This is a script to submit a whole genome fimo scan.', formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    optional.add_argument('-v', '--verbose',dest="verbose", type=str2bool, nargs='?', const=True, default=False, help = 'will generate verbose output file', metavar="True/False", required=False)
    required.add_argument('-g', '--genome', dest="genome", help = 'reference genome in fasta format. Genome index must be in the same directory (genome.fa.fai).', metavar="genome.fa")
    required.add_argument('-o', '--outdir', dest="outdir", help = 'directory for output', metavar="/full/path/to/output")
    required.add_argument('-m', '--motifs', dest="motifs", help = 'meme file for TFs', metavar="motif_database.meme")
    optional.add_argument('-t', '--threshold_fimo',dest="threshold_fimo",type=float, default='1e-5', help = 'threshold for motifs called by fimo. Default: 1e-5', metavar="float", required=False)
    optional.add_argument('-b', '--background_file', dest="background_file", help = 'background base composition of a given genome. This flag is HIGHLY recommended otherwise a uniform background distribution is assumed ie A/T/G/C = 25/25/25/25', metavar="background.csv", default=None, required=False)    

    args = parser.parse_args()

######################################### FIMO Scanner Main #########################################
    def run_fimo_scanner(outdir, genome, motifs, background_file=None, threshold_fimo='1e-5',
                         verbose=False):
        if verbose == True: 
            print("---------Calculating GC Content of Motifs from Motif Database (.meme file)----------")
        motif_list = fimo_motif_names(verbose=verbose,motifs=motifs)
        get_gc(verbose, outdir=outdir, motifs=motifs, motif_list=motif_list, alphabet=['A','C','G','T'])

        if verbose == True: 
            print("---------FIMO Scan: Identifying Locations of Motif Hits in Simulated Genome----------")
            print('Start time: %s' % str(datetime.datetime.now()))
            start_time = int(time.time())
        fimo_dirs(verbose, outdir=outdir)
        for i,motif in enumerate(motif_list):
            if ((i+100)%200 == 0):
                scanner(motif=motif, verbose=verbose, outdir=outdir,
                        motifs=motifs, threshold_fimo=threshold_fimo, 
                        background_file=background_file, genome=genome)
                if verbose == True:
                    print('200 jobs submitted. Sleep for 2hr')
                time.sleep(7200)
            else:
                scanner(motif=motif, verbose=verbose, outdir=outdir,
                        motifs=motifs, threshold_fimo=threshold_fimo, 
                        background_file=background_file, genome=genome)
        if verbose == True: 
            stop_time = int(time.time())
            print('Stop time: %s' % str(datetime.datetime.now()))
            print("Total Run time :", (stop_time-start_time)/3600, " hrs")
            print("----------------Convert FIMO to Bed Format-----------------------")
        fimotobed(verbose=verbose, outdir=outdir, motif_list=motif_list)
        if verbose == True: 
            print("Bed file columns: ['chr','start', 'stop','score', 'strand', 'motif_region_name']")
            print('Scan complete.')

######################################### FIMO Functions #########################################
    def fimo_dirs(verbose, outdir):
        if verbose == True:
            print('Making output directories.')
        os.system('mkdir -p ' + outdir + '/temp') 
        os.system('mkdir -p ' + outdir + '/temp/fimo_out')
        os.system('mkdir -p ' + outdir + '/motifs')    
        if verbose == True:
            if ((path.exists(outdir + '/temp') == True) & (path.exists(outdir + '/temp/fimo_out') == True) & (path.exists(outdir + '/motifs') == True)):
                print('Output directories exisit.')

    def fimo_motif_names(verbose, motifs):
        '''Extracts motif names from a MEME formatted motif database'''
        motif_list = list()
        with open(motifs) as F:
            for line in F:
                if 'MOTIF' in line:
                    line = line.strip('\n').split()
                    motif_name = line[-1]
                    motif_list.append(motif_name)
        first = [motif_list[0:1]]
        last =  [motif_list[-1:]]
        if verbose == True:
            print('There are ' + str(len(motif_list)) + " motifs in this meme file. " + str(first) + '...' + str(last))
        return motif_list

    def get_gc(verbose, outdir, motifs, motif_list, alphabet=['A','C','G','T']):
        '''Obtain a pssm model from a meme formatted database file.'''
        gc_out={}
        for motif in motif_list:
            motif_hit = False
            PSSM = []
            with open(motifs,'r') as F:
                for line in F:
                    if 'MOTIF' in line:
                        if motif in line:
                            motif_hit = True
                        else:
                            motif_hit = False
                    elif motif_hit and 'URL' not in line and 'letter-probability' not in line and line != '\n':
                        acgt_probabilities = [float(x) for x in line.strip('\n').split()]
                        total_prob = sum(acgt_probabilities)
                        acgt_probabilities = [x/total_prob for x in acgt_probabilities] #Convert to probabilities
                        PSSM.append(acgt_probabilities)

                gc = 0
                for base in PSSM:
                    gc += base[alphabet.index('C')]
                    gc += base[alphabet.index('G')]

                gc = gc/float(len(PSSM))
                if verbose == True:
                    print('Percent GC content for ' + motif + ' is ' + str(gc))
            gc_out[motif] = motif,gc,(len(PSSM))

            df_gc_out = pd.DataFrame.from_dict(gc_out)
            df_gc_out = df_gc_out.transpose()
            df_gc_out.columns =['motif', 'percent_gc','motif_length']
            df_gc_out.to_csv(outdir + "/motif_gc_percentage.txt", header=None, index=False, sep='\t')

    def scanner(motif, verbose, outdir, motifs, threshold_fimo, background_file, genome):
        if verbose == True:
            fimo_verbosity = '2'
        else:
            fimo_verbosity = '1'
        full_outdir= str(outdir + '/temp/fimo_out/' + motif)
        os.system('sbatch /Users/tajo5912/fimo_scanner/fimo_scanner/fimo_scan.sbatch ' + fimo_verbosity + ' ' + str(threshold_fimo) + ' ' + background_file + ' ' + motif + ' ' + full_outdir + ' ' + motifs + ' ' + genome)


    def fimotobed(verbose, outdir, motif_list):
        for motif in motif_list:
            df = pd.read_csv(outdir+'/temp/fimo_out/'+motif+'/fimo.tsv', sep ='\t')
            if verbose == True:
                print("Post-processing FIMO output for " + motif + '.')
            df.drop(df.tail(3).index,inplace=True)
            if df.empty == True:
                if verbose == True:
                    print('Skipping ' + motif + ' -no motif hits.')
            else:
                df = df.sort_values(by=['sequence_name', 'start']).reset_index()
                df['start'] = df['start'].astype(int) 
                df['stop'] = df['stop'].astype(int)
                df['count'] = (np.arange(len(df)))
                df['count'] = (df['count']+1).apply(str)
                df['motif_region_name'] = df['sequence_name'] + ';motif_' + df['count']
                df.drop(['count'], axis=1, inplace=True)
                df = df[['sequence_name','start', 'stop','score', 'strand', 'motif_region_name']]
                df.to_csv(outdir + '/motifs/' + motif + '.sorted.bed', sep='\t', header=None, index=False)          
######################################### Call function #########################################
run_fimo_scanner(args.outdir, args.genome, args.motifs, args.background_file, args.threshold_fimo, args.verbose)
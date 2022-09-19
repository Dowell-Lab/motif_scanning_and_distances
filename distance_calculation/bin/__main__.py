import argparse
from argparse import RawTextHelpFormatter
from multiprocessing import Pool
import main
import sys

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
    
    p = argparse.ArgumentParser(description = 'This work flow functions to calculate the distances from the center of a given motif to the center of a single bidirectional.', formatter_class=RawTextHelpFormatter)    
    
    general = p.add_argument_group('Arguments')

    #### General options
    general.add_argument('-v', '--verbose',dest="verbose", type=str2bool, nargs='?', const=True, default=False, help = 'will generate verbose output file', metavar="True/False", required=False)
 
    general.add_argument('-a', '--annotation', dest="annotation", help = 'input bed file of bidirectionals or ATAC peaks, ends with .bed or .sorted.bed', metavar="annotation.bed", required=True)
    general.add_argument('-o', '--outdir', dest="outdir", help = 'directory for output', metavar="/full/path/to/output", required=True)
    general.add_argument('-s', '--sample', dest="sample", help = 'name of the sample to be run (str)', metavar="name_of_sample", required=True)
    general.add_argument('-c', '--cpus', type=int, dest='cpus', metavar='int', help='number of CPUs for multiprocessing. Default=1', default=1, required=False)
    general.add_argument('-w', '--window', dest="window",type=int, default=1500, help = 'window to extract sequences. Default=1500', metavar="int", required=False)
    general.add_argument('-q', '--pre_scan', dest="pre_scan", default=None, help = 'directory containing pre-scanned motif hits in bed format over the whole genome obtained using the same fimo parameters (background and threshold) as the simulated dataset. If path is set experimental genome scan will be skipped and pre-scanned motif hits will be used instead.', metavar="/full/path/to/pre-scanned/motifs", required=False)                   
    general.add_argument('-x', '--experimental_fimo', dest="experimental_fimo", type=str2bool, nargs='?', const=True, default=False, help = 'will run fimo over only the annotated regions from the experimental dataset. True will increase run time. Recommended if you are only looking at a single dataset. If False, provide destination of the pre-scanned genome using the "--pre_scan" flag. Default: False.', metavar="True/False", required=False)

    args = p.parse_args()

main.run(args.verbose, args.outdir, args.sample, args.annotation, 
        args.cpus, args.window, args.experimental_fimo, args.pre_scan)


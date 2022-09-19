# Motif scanning and distances
The code used to generate the motif information using in the DBNascent repository paper for TRN TF anchoring.

There are two main steps in generating the needed motif information for this project.
1. We needed to scan the whole genome for all motif instances in a given motif database.
2. We needed to associate these motif hits with our bidirectional information in a meaningful way. To do this we calculate the distance from the motif to the center of the bidirectional as it has been shown that the center of the bidirectional transcript colocalizes with TF motifs (doi: 10.1101/gr.225755.117).

## Motif Scanning
To perform motif scanning we used the scripts located in fimo_scan. We used the meme suite (meme v5.0.3) to perform the scan. (Charles E. Grant, Timothy L. Bailey and William Stafford Noble, "FIMO: Scanning for occurrences of a given motif", Bioinformatics 27(7):1017-1018, 2011.)
- /fimo_scan/bin contains the main python script used to generate the motif bed files. This script calls fimo_scan/bin/fimo_scan.sbatch
- /fimo_scan/files contains the motif databases and background file used.
  * The motifs used were pulled from HOCOMOCOv11
  * The background assumption is an equal distribution of A/T/C/G
- /fimo_scan/scripts contains the sbatch scripts used to initiate the scanning processes. 
  * Note that a 1e-5 cutoff was used for fimo.

## Distance Calculation
To perform the distance calculation we used the scripts located in distance_calculation.
- distance_calculation/split_master_eRNA_file_and_recombine_distances.ipynb was used to first split the master eRNA files for both hg38 and mm10 into 10 separate, smaller bed files that were more practical to work with. It was then used again after completing the distance calculation to re-combine the distance_calculation output.
- distance_calculation/bin contains the python scripts used to calculate the distances from the center of each bidirectional to the center of each motif provided. 
  * Each TF motif has its own distance table.
- distance_calculation/scripts contains the sbatch scripts used to initiate the scanning processes.
```
Example distance output:

| region_id      | motif_id     | distance | distance_rank | quality_rank |
| -------------- | ------------ | -------- | ------------- | ------------ |
| chr1;region_16 | chr1;motif_1 | 56       | 1             | 1            |
| chr1;region_20 | chr1;motif_4 | 1303     | 3             | 1            | 
| chr1;region_20 | chr1;motif_5 | -742     | 1             | 3            |
| chr1;region_20 | chr1;motif_6 | -751     | 2             | 2            |

```
- The region_id matches the ids in the provided annotation file.
- The motif_id matches the 6th column of the motif bedfiles for whichever TF you are looking at from fimo_scan.
- The distance is the distance from the center of the annotated region (ie the bidir) to the center of the motif for whatever TF you're looking at. All motif distances within +/-1500bp of the center of the bidirectional were calculated. A negative distance indicates upstream of the center of the bidirectional, whereas a positive distance indicates downstream of the bidirectional. Typically, we consider a motif hit within +/-150bp of the center of the bidirectional as "active."
- Distance rank is for the case where there are 2+ motif hits within ONE bidirectional. The example here is chr1;region_20... For the distance rank 1 means closest to the center.
- Quality rank is for the case where there are 2+ motif hits within ONE bidirectional. The example here is chr1;region_20... For the quality rank 1 means the highest confidence fimo call.

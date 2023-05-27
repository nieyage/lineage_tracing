
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk --version
# https://github.com/caleblareau/mgatk/wiki/Modes
mgatk tenx --help
  -bt, --barcode-tag TEXT         Read tag (generally two letters) to separate
                                  single cells; valid and required only in
                                  `bcall` mode.
  -b, --barcodes TEXT             File path to barcodes that will be
                                  extracted; useful only in `bcall` mode. If
                                  none supplied, mgatk will learn abundant
                                  barcodes from the bam file (threshold
                                  defined by the -mb tag).
  -mb, --min-barcode-reads INTEGER
                                  Minimum number of mitochondrial reads for a
                                  barcode to be genotyped; useful only in
                                  `bcall` mode; will not overwrite the
                                  `--barcodes` logic. Default = 1000.
  --NHmax INTEGER                 Maximum number of read alignments allowed as
                                  governed by the NH flag. Default = 1.
  --NMmax INTEGER                 Maximum number of paired mismatches allowed
                                  represented by the NM/nM tags. Default = 4.
  -kd, --keep-duplicates          Retained dupliate (presumably PCR) reads
  -ub, --umi-barcode TEXT         Read tag (generally two letters) to specify
                                  the UMI tag when removing duplicates for
                                  genotyping.
  -ho, --handle-overlap           Only count each base in the overlap region
                                  between a pair of reads once
  -lc, --low-coverage-threshold INTEGER
                                  Variant count for each cell will be ignored
                                  below this when calculating VMR
  -jm, --max-javamem TEXT         Maximum memory for java for running
                                  duplicate removal per core. Default = 8000m.
  -pp, --proper-pairs             Require reads to be properly paired.
  -q, --base-qual INTEGER         Minimum base quality for inclusion in the
                                  genotype count. Default = 0.
  -aq, --alignment-quality INTEGER
                                  Minimum alignment quality to include read in
                                  genotype. Default = 0.
  -eb, --emit-base-qualities      Output mean base quality per alt allele as
                                  part of the final output.
  -ns, --nsamples INTEGER         The number of samples / cells to be
                                  processed per iteration; Default = 7000.
                                  Supply 0 to try all.
  -k, --keep-samples TEXT         Comma separated list of sample names to
                                  keep; ALL (special string) by default.
                                  Sample refers to basename of .bam file
  -x, --ignore-samples TEXT       Comma separated list of sample names to
                                  ignore; NONE (special string) by default.
                                  Sample refers to basename of .bam file
  -z, --keep-temp-files           Add this flag to keep all intermediate
                                  files.
  -qc, --keep-qc-bams             Add this flag to keep the quality-controlled
                                  bams after processing.
  -sr, --skip-R                   Generate plain-text only output. Otherwise,
                                  this generates a .rds obejct that can be
                                  immediately read into R for downstream
                                  analysis.
  -so, --snake-stdout             Write snakemake log to sdout rather than a
                                  file. May be necessary for certain HPC
                                  environments.
  -nfg, --ncells_fg INTEGER       number of "foreground" cells to use for
                                  CellBender. Default = 1000.
  -nbg, --ncells_bg INTEGER       number of "background" cells to use for
                                  CellBender. Default = 20000.


samtools faidx /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa chrM > assembly_hardmask_chrM.fasta


## PBS configure 
#PBS -N mutaion
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=32
#PBS -l mem=32G

source /public/home/nieyg/.bash_profile
cd /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test
source /md01/nieyg/ori/biosoft/conda/etc/profile.d/conda.sh
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk tenx -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/ \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/barcodes.tsv 

#######################
######## test #########
#######################

nohup mgatk tenx -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv &

nohup mgatk check -i /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test4 \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa \
-c 32 -bt CB -b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv & 





/md01/nieyg/software/subset-bam_linux --bam /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/outs/atac_possorted_bam.bam --cell-barcodes test.tsv --cores 10 --out-bam test.bam
samtools sort -@ 12 -O bam -o test_sorted.bam test.bam
samtools index test_sorted.bam 

mgatk check -i test_sorted.bam \
-o /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test4 \
-n plogF1met5d_masked -g /md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask_chrM.fasta -c 32 -bt CB \
-b /md01/nieyg/project/lineage_tracing/jointpolgF1met5d/plogF1met5d_masked/mgatk/test.tsv


Given the .fa or .fasta file that was part of your reference genome for whichever alignment tool, 
one can extract only the mitochondrial genome using the following command for chrM 
(update your mitochondrial genome name accordingly-- other examples may be MT).
 For example:


mgatk tenx -i test_barcode.bam -n bc1 -o test -bt CB -b test_barcodes.txt -c 32

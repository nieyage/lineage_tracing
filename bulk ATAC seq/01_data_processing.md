1. Install mgatk
```
python3 -m venv venv3
source venv3/bin/activate
pip3 install mgatk
mgatk --version
```
2. run mgatk by default

```
conda activate r4-base
source /md01/nieyg/software/venv3/bin/activate
mgatk --version
# https://github.com/caleblareau/mgatk/wiki/Modes
(venv3) (base) mgatk support
Tue Apr 18 15:52:28 CST 2023: mgatk v0.6.7
Tue Apr 18 15:52:28 CST 2023: List of built-in genomes supported in mgatk:
Tue Apr 18 15:52:28 CST 2023: ['GRCh37', 'GRCh38', 'GRCm38', 'GRCz10', 'NC_012920', 'hg19', 'hg19_chrM', 'hg38', 'mm10', 'mm9', 'rCRS']
Tue Apr 18 15:52:28 CST 2023: Specify one of these genomes or provide your own .fasta file with the --mito-genome flag

mgatk call --help

mgatk call -i /md01/nieyg/project/lineage_tracing/bulk_ATAC/2_bam/sort_bam \
-o /data/R02/nieyg/project/lineage_tracing/bulk_ATAC/7_mgatk \
-n P7_bulk_ATAC \
-g mm10 \
-c 12
```
3. Increasing coverage (Increasing mito reads / uniform coverage)
### “Masked reference genome and NUMT comparison.” ([Lareau 等, 2021, p. 463]
#https://github.com/caleblareau/mitoblacklist/



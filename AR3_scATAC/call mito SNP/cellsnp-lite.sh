conda activate cellsnp-lite

nohup cellsnp-lite -s /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/outs/possorted_bam.bam \
 -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C4_scATAC_add500G/mgatk/barcode.tsv \
 -O /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/cellsnp-lite/AR3_C4 \
 -p 20 --minCOUNT 20 --gzip --chrom chrM -f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa &


nohup cellsnp-lite -s /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/outs/possorted_bam.bam \
 -b /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/AR3_C5_scATAC_add500G/mgatk/barcode.tsv \
 -O /md01/nieyg/project/lineage_tracing/heart_regeneration/00_data/AR3_data/scATAC/cellsnp-lite/AR3_C5 \
 -p 20 --minCOUNT 20 --gzip --chrom chrM -f /md01/nieyg/ref/hard-mask/mm10_hard_masked/fasta/genome.fa &

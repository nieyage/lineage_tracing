# Updating CellRanger to have a hard-masked nuclear genome
### 1) generate the hard-masked .fa file and then 
### 2) execute the cellranger-atac mkref command as suggested.

```
cp ~/ref/10X/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa .
mv genome.fa old_genome.fa
bedtools maskfasta -fi old_genome.fa -bed /md01/nieyg/ref/mitoblacklist/mitoblacklist-master/combinedBlacklist/mm10.full.blacklist.bed  -fo assembly_hardmask.fa

```

* cellranger-arc create a reference 
#the config file 
```
{
    organism: "mouse" # same as before
    genome: ["mm10_hard_masked"] # same as before
    input_fasta: ["/md01/nieyg/ref/hard-mask/genome_modify/assembly_hardmask.fa"] # updated from bedtools maskfasta with masked regions
    input_gtf: ["/md01/nieyg/ref/10X/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz"] # same as before
    non_nuclear_contigs: ["chrM"] # same as before
    input_motifs: "/md01/nieyg/ref/hard-mask/genome_modify/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt" # same as before
}

```

`nohup cellranger-arc mkref --config=mm10_hard_masked.config &`


* bowtie build index for bulk ATACseq 
`bowtie2-build ../genome_modify/assembly_hardmask.fa mm10_hard_masked`



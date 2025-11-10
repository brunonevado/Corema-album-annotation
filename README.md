# Corema-album-annotation

## Bruno Nevado, 2025


# Notes: 
- Following (with a few changes) the annotation protocol here: https://begibmcp.notion.site/Genome-Annotation-Protocol-df7a2c4b05444290b879556efcd3b1e2, kindly provided by Aureliano Bombarely.  

- Annotation covered only the 13 largest scaffolds of the assembly, which represent the 13 chromosomes of *C. album*.  

- RNAseq data used from 16 samples, representing 2 individuals (1 male, 8 females), 4 tissues and 2 sampling times (early morning and late afternoon).  


## Repeat annotation

Software used:  
- repeatmodeler v 2.0.6  
- repeatmasker v 4.1.8  
- TEsorter v 1.5.0  

```
module load repeatmodeler/2.0.6

## 1- Generate the RepeatModeler sequence database
BuildDatabase --name coremaV1 Calbum.V1.fas

## 2- Run RepeatModeler2
RepeatModeler -database coremaV1 -threads 48 -LTRStruct

# Output:
#Program Time: 29:19:15 (hh:mm:ss) Elapsed Time
#Working directory:  /home/bnevado/corema/RM_3997133.MonMar171313582025
#may be deleted unless there were problems with the run.
#
#The results have been saved to:
#  /home/bnevado/corema/
#    coremaV1-families.fa  - Consensus sequences for each family identified.
#    coremaV1-families.stk - Seed alignments for each family identified.
#    coremaV1-rmod.log     - Execution log.  Useful for reproducing results.
#
#


## 3- Run RepeatMasker

module load repeatmasker/4.1.8
# Loading requirement: python/3.10.11 hmmer/3.2.1 trf/4.10.0
RepeatMasker -lib coremaV1-families.fa -xsmall -gff -e ncbi -pa 4 Calbum.V1.fas

## 4- Run TESorted

module load TEsorter/1.5.0
#  Loading requirement: python/3.10.11 hmmer/3.3.2 blast+/2.13.0
  
TEsorter coremaV1-families.fa -p 48 -db rexdb-plant
```

<Details>
<summary><span style="color:darkred">Repeat classification summary</span></summary>

- Repeat classification summary done with R package repeatR  
- In total, ca. 51% of the genome assembly consists of repeat elements, with the majority (30% of the assembly) being unknown.  

- Repeat annotation results with all classes.  
![Repeat content per class, including unknown elements](./annotationResults/repearContentPerClass.1.png){#id .class width=60% height=60%}
- Repeat annotation results excluding unclassified repeats.  
![Repeat content per class, excluding unknown elements](./annotationResults/repearContentPerClass.2.png){#id .class width=60% height=60%}

</Details>

## Evidence-based annotation with RNAseq data

Software used:  
- trimgalore v 0.6.10  
- fastqc v 0.12.1  
- STAR v 2.7.11b  
- samtools/bcftools v1.17  
- stringtie v 3.0.1  
- transdecoder v 5.7.1
- blast+ v 2.13.0  
- diamond v 2.1.11  
- hmmer v 3.3.2
```
# 1. Read trimming and adaptor removal
module load trimgalore/0.6.10
module load fastqc/0.12.1


fastqc $FREADS -o /home/bnevado/Lupinus/resequencing/raw
fastqc $RREADS -o /home/bnevado/Lupinus/resequencing/raw

trim_galore -q 20 --length 50 --fastqc --output_dir $OUTPUTFOLDER \
	    --paired $FREADS $RREADS

# 2. QC of raw and trimmed reads
fastqc $FREADS $RREADS [...]

# 3. Genome preparation  

STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./genome --genomeFastaFiles ./genome/Calbum.V1.fas --genomeSAindexNbases 10

# 4. Mapping RNAseq data to genome

module load star/2.7.11b

STAR --runThreadN $NT --genomeDir $GENOMEFOLDER \
 --readFilesIn $FREADS $RREADS \
 --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate \
 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical \
 --outFileNamePrefix $OUTPUTFOLDER/$SAMPLE.bam --limitBAMsortRAM 10000000000

# 5. QC of mapped data  
 samtools flagstats $SAMPLE.bam
 
# 6. Run StringTie on each of the BAM files obtained from the mapping.

module load stringtie/3.0.1
stringtie -o bams/$_.stringtie.gtf -p 8 bams/$_.bamAligned.sortedByCoord.out.bam > out.stringtie.$_ 2> err.stringtie.$_\n";'

# 7. Merge all transcripts  
ls bams/*gtf > files.gtfs
stringtie --merge -o merged.STAR.my_genome.gtf files.gtfs

# 8. Identify CDS with TRansDecoder  

module load transdecoder/5.7.1

gtf_genome_to_cdna_fasta.pl\
  merged.STAR.my_genome.gtf genome/Calbum.V1.fas > merged.STAR.my_genome.transcripts.fasta

# 9. Extract the long ORFs
TransDecoder.LongOrfs -t merged.STAR.my_genome.transcripts.fasta

# 10. Search for protein evidences by protein domains and sequence homology
module load blast+/2.13.0
# downloaded uniprot database 
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz
# make blast DB
makeblastdb  -in uniprot_trembl.fasta -input_type fasta -dbtype prot -parse_seqids -out uniprot_trembl_250325

module load diamond/2.1.11

diamond blastp --query merged.STAR.my_genome.transcripts.fasta.transdecoder_dir/longest_orfs.pep \
  --db /home/bnevado/software/blast/db/uniprot_trembl_250325 --outfmt 6   --max-target-seqs 1 \
  --evalue 1e-25 --out longest_orfs.pep.dmd.diamond.o6.txt --threads 20 --sensitive


module load hmmer/3.3.2
hmmscan --cpu 48 --domtblout longest_orfs.pfam.domtblout \
  /home/bnevado/software/blast/db/Pfam-A.hmm \
  merged.STAR.my_genome.transcripts.fasta.transdecoder_dir/longest_orfs.pep

# 11- Prediction of the CDS using the evidences from protein and domain hits

TransDecoder.Predict -t merged.STAR.my_genome.transcripts.fasta --retain_pfam_hits longest_orfs.pfam.domtblout --retain_blastp_hits longest_orfs.pep.dmd.diamond.o6.txt 

# 12- Generate a genome-based coding region annotation file
gtf_to_alignment_gff3.pl merged.STAR.my_genome.gtf > all_reads.STAR.my_genome.gff3

cdna_alignment_orf_to_genome_orf.pl \
  merged.STAR.my_genome.transcripts.fasta.transdecoder.gff3 \
  all_reads.STAR.my_genome.gff3 \
  merged.STAR.my_genome.transcripts.fasta  > all_reads.STAR.my_genome.transcripts.genome.gff3 2> err.cdna_alignment_orf_to_genome_orf.pl
```




## Annotation using genomes from closely related species  

- For this step we will use the genome of *Rhododendron vialii*, accession GCA_030253575.1, which seems to be very good and has a NCBI RefSeq annotation. I blasted several random sequences from the genome and the closest meaningfull hit was always to *R. vialii*, thus should be the closest species available.  

Software used:  
- gemoma v 1.9

```
module load gemoma/1.9
gemoma CLI GeMoMaPipeline threads=12 outdir=gemoma_Rvialii GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=genome/Calbum.V1.fas a=closeSpecies/Rvialii/GCF_030253575.1_ASM3025357v1_genomic.gff g=closeSpecies/Rvialii/GCF_030253575.1_ASM3025357v1_genomic.fna
```

## Annotation ab initio with Augustus

Software used:  

- busco v 5.4.7  
- augustus v 3.5.0  
- braker v 3.0.8  
- agat v 1.4.2  

To run BRAKER and train Augustus, it is necessary the following input files:

1. **Softmasked reference genome in FASTA format**. This file can be obtained from RepeatMasker (see the previous sections).
```
./annotation/Calbum.V1.fas.masked
```

2. **RNA-Seq data mapped to the reference in BAM format**. This file can be obtained from the RNA-Seq mapping with START or Hisat2 (see the previous sections). The files can be merged with samtools merge.
```
samtools merge -o all_reads.bam --threads 26 \
Calb001.bamAligned.sortedByCoord.out.bam  Calb004.bamAligned.sortedByCoord.out.bam  Calb007.bamAligned.sortedByCoord.out.bam \
Calb010.bamAligned.sortedByCoord.out.bam  Calb013.bamAligned.sortedByCoord.out.bam  Calb016.bamAligned.sortedByCoord.out.bam \
Calb002.bamAligned.sortedByCoord.out.bam  Calb005.bamAligned.sortedByCoord.out.bam  Calb008.bamAligned.sortedByCoord.out.bam \
Calb011.bamAligned.sortedByCoord.out.bam  Calb014.bamAligned.sortedByCoord.out.bam Calb003.bamAligned.sortedByCoord.out.bam \
Calb006.bamAligned.sortedByCoord.out.bam  Calb009.bamAligned.sortedByCoord.out.bam  Calb012.bamAligned.sortedByCoord.out.bam \
Calb015.bamAligned.sortedByCoord.out.bam
```

3. **Protein sequences of a closed related species in FASTA format**.  

```
closeSpecies/Rvialii/GCF_030253575.1_ASM3025357v1_protein.faa
```

<Details>
<summary><span style="color:darkred">BUSCO assessment before continuing</span></summary>

**Let's run busco on each of these to see how much information we have to work with**

```
# C. album genome
busco -i Calbum.V1.fas.masked -o outBUSCO_Calbum.V1.fas.masked -m genome -l viridiplantae 
        --------------------------------------------------
        |Results from dataset viridiplantae_odb10         |
        --------------------------------------------------
        |C:98.6%[S:96.5%,D:2.1%],F:0.2%,M:1.2%,n:425      |
        |419    Complete BUSCOs (C)                       |
        |410    Complete and single-copy BUSCOs (S)       |
        |9      Complete and duplicated BUSCOs (D)        |
        |1      Fragmented BUSCOs (F)                     |
        |5      Missing BUSCOs (M)                        |
        |425    Total BUSCO groups searched               |
        --------------------------------------------------
        
```
```
# proteins from R. vialii
busco -i GCF_030253575.1_ASM3025357v1_protein.faa -o outBUSCO_GCF_030253575.1_ASM3025357v1_protein.faa -m proteins -l viridiplantae 
        --------------------------------------------------
        |Results from dataset viridiplantae_odb10         |
        --------------------------------------------------
        |C:99.5%[S:51.5%,D:48.0%],F:0.0%,M:0.5%,n:425     |
        |423    Complete BUSCOs (C)                       |
        |219    Complete and single-copy BUSCOs (S)       |
        |204    Complete and duplicated BUSCOs (D)        |
        |0      Fragmented BUSCOs (F)                     |
        |2      Missing BUSCOs (M)                        |
        |425    Total BUSCO groups searched               |
        --------------------------------------------------
```

```
# transcripts
busco -i merged.STAR.my_genome.transcripts.fasta.transdecoder.cds -o outBUSCO_merged.STAR.my_genome.transcripts.fasta.transdecoder.cds -m transcriptome -l viridiplantae 
        --------------------------------------------------
        |Results from dataset viridiplantae_odb10         |
        --------------------------------------------------
        |C:98.1%[S:55.3%,D:42.8%],F:0.9%,M:1.0%,n:425     |
        |417    Complete BUSCOs (C)                       |
        |235    Complete and single-copy BUSCOs (S)       |
        |182    Complete and duplicated BUSCOs (D)        |
        |4      Fragmented BUSCOs (F)                     |
        |4      Missing BUSCOs (M)                        |
        |425    Total BUSCO groups searched               |
        --------------------------------------------------
```

</Details>

- Pipeline used 

```
module load braker/3.0.8

braker.pl --genome Calbum.V1.fas.masked \
          --prot_seq ./closeSpecies/Rvialii/GCF_030253575.1_ASM3025357v1_protein.faa \
          --bam ./bams/all_reads.bam --softmasking --thread 8 --species=Rvialii

## Fix the GTF for braker (TSEBRA script)
fix_gtf_ids.py --gtf  braker/braker.gtf --out  braker/braker.fixed.gtf

## Add the prefix for the Augustus genes (TSEBRA script)
rename_gtf.py --gtf braker/braker.fixed.gtf --prefix Calbum_v1_Rvialii --translation_tab augustus_renamed_genes.txt --out braker/braker.fixed.gtf

## Add UTRs with ingeannot

module load agat/1.4.2
agat_convert_sp_gxf2gxf.pl -g merged.STAR.my_genome.gtf -o ./merged.STAR.my_genome.WGenes.gtf

## Generate a BAM file of files
ls | grep bam$ | awk '{ print $1"\tTrue\tTrue"}'  | grep -v all_reads | perl -pe 'print "bams/"' > bam.fof

## Run ingenannot clusterize
ingenannot clusterize -g all_reads.STAR.my_genome.transcripts.genome.gff3 --keep_atts merged.STAR.my_genome.WGenes.gtf clusterised_transcripts.gff

## Fix the format of the clusterised_transcripts.gff
agat_convert_sp_gxf2gxf.pl -g clusterised_transcripts.gff -o clusterised_transcripts.formatted.gff

## Rank the isoforms with IngenAnnot
ingenannot -v 2 -p 8 isoform_ranking -f bams/bam.fof clusterised_transcripts.formatted.gff

## Add the UTRs
ingenannot -v 2 -p 48 utr_refine all_reads.STAR.my_genome.transcripts.genome.gff3 isoforms.alternatives.gff --erase --utr_mode rank  braker.wutrs.gff3 

```

## Annotation ab initio with SNAP and maker

Software used:

- snap v 20060728  
- maker v 3.01.04  
  
### Round 1
```
## 01- Generate the configuration files
maker -CTL

## 02- Modify the following lines in the configuration file
##
# genome=./annotation/Calbum.V1.fas.masked #genome sequence (fasta file or fasta embeded in GFF3 file)
# est=braker/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.fasta
# protein=closeSpecies/Rvialii/GCF_030253575.1_ASM3025357v1_protein.faa
# est2genome=1
# protein2genome=1
# rmlib=../annotation/coremaV1-families.fa #provide an organism specific repeat library in fasta format for RepeatMasker

## 03 - Run MAKER
maker -cpus 80 > out 2> err

## 04- Collect the annotations in a GFF file with the FASTA sequences at the end
gff3_merge -o my_genome.maker.snap_t.r00.gff3 -d Calbum.V1.fas.maker.output/Calbum.V1.fas_master_datastore_index.log

## 05 - Check how many models have an AED < 0.25
grep -v "#" my_genome.maker.snap_t.r00.gff3 | awk '{ if ($3 ~ /mRNA/) print $0}' | sed -r 's/.+_AED=//' | sed -r 's/;.+//' | awk '{ if ($1 < 0.25) print $0}' | wc -l

## 06 - Transform the GFF maker file into ZFF filtering the models with a min. AED of 0.25 and a min. length of 50 bp
maker2zff -x 0.25 -l 50 -d Calbum.V1.fas.maker.output/Calbum.V1.fas_master_datastore_index.log
mkdir snap_round01
mv genome.* snap_round01
cd snap_round01

## 07- Get some gene stats and validate the models
fathom genome.ann genome.dna -gene-stats > gene_stats.txt
fathom genome.ann genome.dna -validate > validate.txt

## 08- Collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus

## 09- Generate the training parameters
forge export.ann export.dna

## 10- Generate the HMM training file
hmm-assembler.pl my_genome_SNAP_r01 ./ > ../my_genome_SNAP_r01.hmm
```

### ROUNDS 2-3

- Change the options below in maker ctl file and rerun twice, using the new models generated on each round.  
```
## snaphmm=my_genome_SNAP_r01.hmm
## est2genome=0
## protein2genome=0

maker -cpus 80 > out 2> err
```


## Final result: merging and tidying gff files
 
NOTE: Initial runs with repeat + gene annotation as the final step in the maker pipeline failed for several chromosomes despite multiple tries and different options used. I thus opted to report the repeat and gene coding annotation in separate gff files.  

Software used:  
- genometools v 1.6.2

```
gt gff3 -sort -sortlines yes -tidy Calbum.repeats.gff3 > Calbum.repeats.final.gff3
gt gff3validator Calbum.repeats.final.gff3

gt gff3 -sort -tidy Calbum.geneannot.allCHR.gff3 > Calbum.geneannot.allCHR.final.gff3
gt gff3validator Calbum.geneannot.allCHR.final.gff3

```

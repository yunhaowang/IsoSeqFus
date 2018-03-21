# IsoSeqFus (version 0.2)
Time-stamp: <2018-02-06 Yunhao Wang, Email: yunhaowang@126.com>


## Introduction

Full-length isoform sequencing (Iso-Seq) technology originally developed by Pacific Biosciences (PacBio) has been widely applied to transcriptome study. IsoSeqFus is a bioinformatics tool to identify fusion genes and gene isoforms using PacBio Iso-Seq data.


## Prerequisite

- Linux system

- python 2.7

- Numpy (tested with version 1.24.0)


## Install and Run

1. Download the package (e.g., `IsoSeqFus-0.2.tar.gz`) to a directory (e.g., `/home/`)

2. Unpack it using the command `tar -zxvf /home/IsoSeqFus-0.2.tar.gz`

3. Now, you can run IsoSeqFus by the executable file `/home/IsoSeqFus-0.2/bin/isoseqcon`. Optional, you can add IsoSeqFus into your PATH so that you can run IsoSeqFus without having to specify the entire path. For example, you can add one line `export PATH=/home/IsoSeqFus-0.2/bin:$PATH` to your `~/.bashrc`


## Input

### 1. Gene annotation file (GTF format)

Note: considering gene duplication phenomenon, wherein multiple genic loci (from same/different chromosomes) can be transcribed into the sequence-identical transcript, thus make sure the transcript ID is unique for each derived genic loci. Suggest to use the well-annotated gene annotation libraries (e.g., RefSeq, Ensembl and Gencode).

### 2. PacBio Iso-Seq long read alignment file(s) (SAM format)

Multiple sam files can be as inputs split by "space".

Please follow the processes below to generate the Iso-Seq long read alignment files. (Here, we take SIRV (Spike-In RNA Variant by Lexogen Inc.) Iso-Seq data as an example.)

#### Raw Iso-Seq dataset

After sequencing SIRV sample by PacBio RSII platform, we get a batch of `*.bax.h5` files. Now, the BAM format (typically, `*.subreads.bam`) is the standard output of PacBio sequencer (e.g., PacBio Sequel system).

#### Step1: extract ROI

Extract ROI (Reads of Insert, also historically called CCS) by PacBio SMRT® Analysis Software or SMRT® Link Software.

In our test, we use SMRT Analysis (v2.3.0): `ConsensusTools.sh CircularConsensus --minFullPasses 0 --minPredictedAccuracy 70` to extract ROI. Now, we get the `SIRV_ROI.fasta` file. 

Alternatively, if your raw Iso-Seq data is `*.subreads.bam` file, you can use SMRT Link (v5.0.0): `ccs --minPasses 0 --minPredictedAccuracy 0.7` to get `*.ROI.bam` file.

#### Step2: classifiy ROI

Classifiy ROI by PacBio SMRT® Analysis Software or SMRT® Link Software.

In our test, we use SMRT Analysis (v2.3.0): `pbtranscript.py classify --flnc SIRV_ROI.flnc.fasta --nfl SIRV_ROI.nfl.fasta -d ./output/ -p SIRV_ROI.primer_info.csv --detect_chimera_nfl SIRV_ROI.fasta SIRV_ROI.classify.fasta`. Now, we get the FLNC (full-length non-chimera) ROI (i.e., `SIRV_ROI.flnc.fasta`) and nFLNC (non-full-length non-chimera) ROI (i.e., `./output/nflnc.fasta`)

Alternatively, if your data is `*.ROI.bam` format, you can use SMRT Link (v5.0.0): `pbtranscript.py classify --flnc *.flnc.fasta --nfl *.nfl.fasta -d ./output/ -p *.primer_info.csv --detect_chimera_nfl *.ROI.bam *.ROI.classify.fasta`. Similarly, you can get the FLNC and nFLNC.

#### Step3: separate nFLNC ROI

Separate the nFLNC ROI by the script `./IsoSeqCon-0.2/utilities/py_isoseqcon_separate_nflnc_fasta.py -i ./output/nflnc.fasta -s ./output/SIRV_ROI_nflnc_sense.fasta -u ./output/SIRV_ROI_nflnc_unknown.fasta`. Now we have 3 fasta files: 1) `SIRV_ROI.flnc.fasta` (with strand information); 2) `./output/SIRV_ROI_nflnc_sense.fasta` (with strand information); and 3) `./output/SIRV_ROI_nflnc_unknown.fasta` (without strand information).

#### Step4: align ROI

Align 2 fasta files (`SIRV_ROI.flnc.fasta` and `./output/SIRV_ROI_nflnc_sense.fasta`) to reference genome by GMAP aligner (version 2016-06-09). We used the parameter `-z sense_force -f samse -n 0 --split-output ./SIRV_ROI_sense`. Now, we get 4 output files (SAM format) but only take 1 of them (`SIRV_ROI_sense.transloc`) as the input of isoseqcon.

Alternatively, you can use other aligners (e.g., BLAT) to align your Iso-Seq long read to reference genome. However, some points need to be noted: (1) Iso-Seq long read sequences in `FLNC.fasta` and `nFLNC_sense.fasta` files are sense strand; (2) For fusion long reads, must be split into two parts and at most one alignment path for each part is output; (3) The aligned strand should be marked using the Tag:Type:Value = 'XS:A:+/-/?' in the optional fields (>=12 colunm) of SAM output file.

We strongly suggest you to use GMAP which is best aligner for Iso-Seq long read transcriptome data (see the paper "Evaluation of tools for long read RNA-seq splice-aware alignment. bioRxiv (2017)").

### 4. PacBio Iso-Seq long read primer information file (CSV format)

In the Step2 of generating the PacBio Iso-Seq long read alignment file (i.e., classifiy ROI by PacBio SMRT® Analysis Software or SMRT® Link Software), you can output the primer information (CSV format) using the parameter `-p` (e.g., we use SMRT Analysis (v2.3.0): `pbtranscript.py classify --flnc SIRV_ROI.flnc.fasta --nfl SIRV_ROI.nfl.fasta -d ./output/ -p SIRV_ROI.primer_info.csv --detect_chimera_nfl SIRV_ROI.fasta SIRV_ROI.classify.fasta`). The output file (`SIRV_ROI.primer_info.csv`) contains the primer information for all ROIs.


## Output

### 1. Identified fusion gene and gene isoform (modified GPD format)

(1) Fusion transcript (gene) ID

The prefix of ID is "Fusion_transcript_"

(2) Number of supporting Iso-Seq full-length long reads

(3) Number of supporting Iso-Seq long reads (both full-length and non-full-length)

(4) If two parts of fusion transcript are from same chromosome

(5) If two parts of fusion transcript have overlap

(6) Derived annotated gene of first part of fusion transcript

If the prefix is "Novel_loci_", it means a novel genic locus

(7) Derived annotated isoform of first part of fusion transcript

If the prefix is "Novel_loci_", it means a novel genic locus

(8) Fusion site information of first part of fusion transcript. [min, max, mode, mode_frequency, median, mean, standard deviation]

(9) Derived annotated gene of second part of fusion transcript

If the prefix is "Novel_loci_", it means a novel genic locus

(10) Derived annotated isoform of second part of fusion transcript

If the prefix is "Novel_loci_", it means a novel genic locus

(11) Fusion site information of second part of fusion transcript. [min, max, mode, mode_frequency, median, mean, standard deviation]

(12-19) Genic structue of first part of fusion transcript:

(12) Chromosome

(13) Strand

(14) Transcription start site (for "+" strand)

(15) Transcription end site (for "+" strand)

(16) Sequence length

(17) Exon number

(18) Exon start site set

(19) Exon end site set

(20-27) Genic structue of second part of fusion transcript:

(20) Chromosome

(21) Strand

(22) Transcription start site (for "+" strand)

(23) Transcription end site (for "+" strand)

(24) Sequence length

(25) Exon number

(26) Exon start site set

(27) Exon end site set

### 2. Identified fusion gene and gene isoform (GTF format)

In the attribute field (9th column if split by 'tab'), there are some tag-value pairs corresponding to the modified GPD file above:

(1) gene_id

1st column of GPD file

(2) transcript_id

1st column of GPD file, but with the suffix "_First_Part" if from first part of fusion transcript; with the suffix "_Second_Part" if from second part of fusion transcript

(3) full_length_LR_count

2nd column of GPD file

(4) LR_count

3rd column of GPD file

(5) same_chr

4th column of GPD file

(6) overlap_two_parts

5th column of GPD file

(7) fst_derived_gene

6th column of GPD file

(8) fst_derived_isoform

7th column of GPD file

(9) fst_fusion_site_info

8th column of GPD file

(10) fst_sequence_length

16th column of GPD file

(11) sec_derived_gene

9th column of GPD file

(12) sec_derived_isoform

10th column of GPD file

(13) sec_fusion_site_info

11th column of GPD file

(14) sec_sequence_length

24th column of GPD file


## Usage and Example

`./bin/isoseqfus -a example/input/SIRV_gene_annotation.gtf -l example/input/SIRV_ROI_sense.transloc -c example/input/SIRV_ROI.primer_info.csv --tempdir example/temp --output_gpd example/SIRV_test.gpd --output_gtf example/SIRV_test.gtf`

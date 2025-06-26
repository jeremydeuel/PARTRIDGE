# PARTRIDGE

<img src="partridge_logo.png" alt="partridge logo" width="200"/>

Precise Analysis of retroTransposed replications of IAPs by deep genomic enrichment


by Jeremy Deuel <jeremy.deuel@usz.ch>

## Purpose
This program is used to analyse IAP driven insertions in hybridisation capture enriched paired-end sequencing data.

# Wetlab part

You need good quality bulk genomic DNA from the tissue of interest. 

## Materials

* High quality DNA extracted from the tissue of interest. We only tested DNA directly extracted from fresh RBC-lysed bone marrow
* Qiagen All-Prep Kit
* xGen DNA Library Prep EZ UNI Kit
* xGen UDI-UMI full length adapters
* xGen Blocking Primers for UDI-UMI full length adapters
* xGEN Hybridisation Kit
* Biotinylated hybridisation capture probes at 1µM final concentration (not concentration per probe but summed concentration of all probes), see below
* Human COT DNA
* Qbit HS DNA Assay
* Bioanalyser HS DNA Assay

## Equipment
* Nanodrop
* Two Thermocyclers
* Bioanalyser
* Qbit
* SpeedVac

## Biotinylated hybridisation capture probes
Order these Oligos with a 5' or 3' biotinylation (we used Sigma)
```
5'-TGTTGGGAGCCGCGCCCACATTCGCCGTTACAAGATGGCGCTGACAGCTGTGTTCTAAGTGGTAAACAAATAATCTGCGCATGTGCCAAGGGTATCTTATGACTACTTGTGCTCTGCCT-3'
5'-TGTGGGGAGCCGCCCCCACATTCGCCGTTGCAAGATGGCGCTGACATCCTGTGTTCTAAGTGGTAAACAAATAATCTGCGCATGTGCCAAGGGTAGTTCTCCACCCCATGTGCTCTGCC-3'
5'-TCTTGCTCTCTTGCTTCTTGCACTCTGGCTCCTGAAGATGTAAGTAATAAAGCTTTGCCGCAGAAGATTCTGGTCTGTGGTGTTCTTCCTGGCCGGTCGTGAGAACGCGTCTAATAACA-3'
5'-GCAATAGAGCTCTTGCTCTCTTGCTCTCTGGCTCCTGAAGATGTAAGCAATAAAGTTTTGCCGCAGAAGATTCCGGTTTGTTGCGTTCTTCCTGGCCGGTCGCGAGAACGCGTGTAAGA-3'
```

## Protocol
Generally, follow xGen Protocol. Here, only the deviations are noted.
Use a 1 min less than the fragmentation time recommended for 350bp fragment length to generate slightly longer fragments. After ligating the full length adapters, clean up by adding 48µl Ampure XP beads to each sample and directly generate pools with a total of 500ng-4µg of DNA per pool. Use a SpeedVac to dry the samples after elution of DNA. For hybridisation use 4µL of 1µM pooled hybridistion capture probes (4pmol) and follow the xGen protocol. After washing, amplify by a total of 11 PCR cycles. Don't forget to inform the sequencing facility about the special settings needed to sequence the UMI barcodes (see xGen protocol for details).

# Bioinformatic analysis

## Overview

The program is run in several steps

1. Generate bam files
2. Extract insertions isa_hmmer2.py
3. Collect insertions
4. Gather evidence sequences and paired sequencing mates to identify retrotransposons
5. Manually annotate novel insertion sites
6. Measure Source Virus Coverage

## Step 1: Generate bam files

You need to generate bam files from your sequencing data. 
* trim sequencing adapters using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
* align to the [GRCm39](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/) reference genome using [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
* index the bam files using samtools 

There are plenty of tutorials on the internet indicating how to do this and you can have a look at the file [scripts/make_bams.sh](scripts/make_bams.sh) for an example.

## Step 2: Extract insertions isa_hmmer2.py
```{bash}
python peartree/main.py --step discover --bam sample_name.bam --out output_folder/sample_name.txt.gz --threads 6
```
Requirements (install with pip):
* see [env.txt](env.txt), use python [venv](https://docs.python.org/3/library/venv.html) to build your environment. I recommend python3.13, although 3.10+ is expected to work.
* HMMER, install from http://hmmer.org
* iap.hmm file, the path is set in the script on line 491. the model has to be "pressed" using hmmpress prior to usage `hmpress iap.hmm`

Parameter | Description
--- | ---
bam | path to bam file of an experiment, aligned, clipped and sorted by coordinate
output | path to output file of an experiment, tsv file
threads | number of threads to use, I recommend 6 for best performance

## Step 3: Collect insertions

`scripts/03_collect_insertions.Rmd`

It is best to run this script step-by-step in RStudio. First, download the .tsv files and save them in a folder `breakpoints`. The script then creates a file `isa_hmmer2_out.xlsx` containing all novel insertions that you need for the next step. The script also generates a file "isa.audit.xlsx" that is an aggregated view of novel insertions that can be used for manual checking (see Step 5)

## Step 4: Gather evidence sequences and paired sequencing mates to identify retrotransposons

```{bash}
python partridge/collect_evidence.py $SAMPLE_NAME
```
Note: You have to adjust the path to the excel file generated in Step 3 in the script, it is hardcoded. Also the input file and output file roots are hard-coded and might have to be adapted. I ran this step using a singularity container, the singularity definition file and the slurm submission scripts are also in the partride folder. But this script can also be run outside of a container using your venv. 

This script generates non-standard fasta files that can be used to manually check novel insertions, see next step.

## Step 5: Manually annotate novel insertion sites

Use these rules to manually annotate novel insertion sites
1.	A least one supporting read with >20bp of LTR coverage for each end of the LTR is present
2.	The insertion is not within an obviously repetitive region
3.	No reads have discordant starting points of the LTR
4.	There is a TSD of more than two bp and less than 30bp
5.	Both LTR ends map to the same LTR in the reference genome (>93% sequence identity)- and this LTR is part of an LTR pair (2-8kb apart)

Some information needed for this annotation can be found in the non-standard fasta files:

Practically, start with the file "isa.audit.xlsx" and generate columns to indicate passing of filtering and also a column indicating if the source IAP could be uniquely identified.

### Example of such a novel insertion verified in the evidence file
```{evidence.fa}
@chr16:93993165-93993171 (-)
>RIGHT_CLIP
TGTTATTCGACGCGTTCTCACGACCGGCCAGGAA
TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACC
TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACC
TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACC
GAGACTGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACC
TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACCAGAATCTTCTGCGACAAAGCTTTATTCTTACATCT
>LEFT_CLIP
GGCTCATGCGCAGATTATTTGTTTACCAACTTAGAACACAGGATGTCAGCGCCATCTTGTGACGGCGAATGTGGGGGCGGCTTCCCACA
GGCTCATGCGCAGATTATTTGTTTACCAACTTAGAACACAGGATGTCAGCGCCATCTTGTGACGGCGAATGTGGGGGCGGCTTCCCACA
>LEFT_MATE
CTGTGTTACGGGAACCTTATAACCTTGATTCGCAGTTCTGGTTCTGGAATGAGGTATCCCTCCTGCGCCAGTCCGGAGTTTTTTCTCGTCCCGGATTTTCTCGTCCCGGGTTTCGGCACCAATTGTTATTCGACGCGTTCTCACGACCGGC
TTCTGAGTTGTCCTTGGCATGCGGTCAAGATGTCACTTTGTATTGCATTTTGAAATGTCTGATGTGGGAAGCCGCCCCCACATTCGCCGTCACAAGATGGCGCTGACATCCTGTGTTCTAAGTTGGTAAACAAATAATCTGCGCATGAGCC
>RIGHT_MATE
```

`>RIGHT_CLIP` lists all sequences clipped to the right (in +-strand direction) whereas `>LEFT_CLIP` does the same for clips to the left.
In this example, 6 identical sequences starting with `TGTTATT...` can be found clipped to the right (one has some additional bases clipped by the aligner) and left-clipped sequences ending with `..CGGCTTCCCACA`.
To check which IAP this sequences belong to, paste the longest RIGHT_CLIP immeditately followed by "LEFT_CLIP" into the [UCSC BLAT](https://genome.ucsc.edu/cgi-bin/hgBlat) tool and select the mm39 reference genome.
In this example, we paste `TGTTATTCGACGCGTTCTCACGACCGGCCAGGAAGAACACCACAGACCAGAATCTTCTGCGACAAAGCTTTATTCTTACATCT GGCTCATGCGCAGATTATTTGTTTACCAACTTAGAACACAGGATGTCAGCGCCATCTTGTGACGGCGAATGTGGGGGCGGCTTCCCACA`. The "mate" sequences can be used if the alignment is not unique (not necessary in this example)
This outputs
```{blat output}
ACTIONS                 QUERY   SCORE START   END QSIZE IDENTITY  CHROM  STRAND  START       END   SPAN
------------------------------------------------------------------------------------------------------------
browser new tab details YourSeq   171     1   172   172   100.0%  chr2   -    84335554  84335888    335
browser new tab details YourSeq   171     1   172   172   100.0%  chr2   -    84340431  84340765    335
browser new tab details YourSeq   171     1   172   172   100.0%  chr3   +    60397294  60397628    335
browser new tab details YourSeq   171     1   172   172   100.0%  chr3   +    60402171  60402505    335
browser new tab details YourSeq   169     1   172   172    99.5%  chr2   -   154201119 154201512    394
browser new tab details YourSeq   169     1   172   172    99.5%  chr16  -    31971790  31972157    368
browser new tab details YourSeq   169     1   172   172    99.5%  chr8   +    72982073  72982409    337
browser new tab details YourSeq   169     1   172   172    99.5%  chr8   +    72986458  72986794    337
browser new tab details YourSeq   169     1   172   172    99.5%  chr3   +   156505180 156505515    336
browser new tab details YourSeq   169     1   172   172    99.5%  chr3   +   156510063 156510398    336
```
We have four 100% hits, covering the entire sequence. Two of these hits are very close, these are both LTRs of a virus. In this case, the virus on chr3 is in the _Gm13710_ gene, whereas the IAP on chr2 is near the _Mbln1_ gene. This example was deliberately chosen, since both of these viruses are 100% sequence identical. Although it is not possible to identify which of these two viruses caused this novel insertion, so I count this very special case of IAP source virus (that jumps a lot!) as a single unique source.

## Step 6: Measure Source Virus Coverage

For this step, we only consider viruses manually validated (see step above). Make sure to annotate at least one source virus. If the source virus is not unique, annotate one hit, since it is important to calculate the coverage of the source virus to estimate VAF. The assumption is, that hybridisation capture enrichment has the same efficieny on the source and the destination virus, since both are sequence identical and the binding sequence of the hybridisation probe is thus also identical. Therefore, the ratio of coverage between destination (unknown zygosity / subclonal?) and the source (fixed in the genome and thus homozygous in inbred mouse strains) can be used to infer allelic frequency.

Use the R script `scripts/06_cluster_source_coverage.R`, which I run on the computing cluster. This file will calculate the source locus coverage. This script has to be lightly adapted to your circumstances and you especially have to adapt the list of source viruses you are interested in. It will output an excel file per bam file indicating the coverage per source locus.

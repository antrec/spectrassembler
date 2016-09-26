# spectrassembler
Doc and full code to come (in a couple of days)

## Introduction
We provide the code used in [our paper] (http://arxiv.org/abs/1609.07293) to obtain the layout of uncorrected long-reads from the overlaps in an original and efficient way. It takes as input an overlap file (from [minimap] [minimap]) and uses a spectral algorithm to find an ordering of the reads from the pairwise overlap information. This ordering should roughly match the sorting of the reads by position on the genome. This ordering can be viewed as a coarse-grained layout, from which a more refined one is computed with the overlap information.
The pipeline also includes a consensus stage inspired from [nanocorrect] [nanocorrect]) ((c) 2015 Ontario Institute for Cancer Research) which uses [poa] [poa]) to compute multiple sequence alignements in a sliding window along the genome.

### Disclaimer
This tool is experimental and does not replace mainstream assemblers. Specifically, there is no sanity check preventing the layout found to be totally wrong for a bad choice of parameters (a warning is raised when the ordering is "suspicious", though). The sensitive parameter is the threshold on the similarity matrix, --sim_thr, which can be increased if a warning appears. However, the default value (--sim_thr 850) should work for bacterial genomes and allow for the layout to be computed cheaply, and for a high-quality draft genome to be derived.

## Dependencies
* overlapper : [minimap] (https://github.com/lh3/minimap)
To install (require gcc and zlib):
```sh
git clone https://github.com/lh3/minimap && (cd minimap && make)
```

* python packages : numpy, scipy, biopython. Easy install with pip for example : ```sh pip install biopython```

* multiple sequence aligner : [poa] [poa] (for consensus generation. Not needed to compute the layout)

* (optional) mapping to reference : [bwa] [bwa]. Needed only if you have a reference genome available and wish to plot the layout found by our algorithm vs the one found by mapping against the reference genome.

## Getting started
Create a working directory
```sh
mkdir oxford-test && cd oxford-test
```
The Loman lab released a set of E. Coli Oxford Nanopore reads (see their [page] (http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/)) or download from the command line:
```sh
curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta
```
Compute alignments with [minimap] [minimap]
```sh
minimap -S oxford.fasta oxford.fasta > oxford.mini
```
(Optional : If you want to check the layout found against the reference genome, you can download a FASTA file of the reference genome of E. Coli K-12 substr. MG1655. [here] (https://www.ncbi.nlm.nih.gov/nuccore/556503834). If you saved it under the name oxford_reference.fa then you can run [bwa] [bwa] to map the reads against the reference :
```sh
bwa index oxford_reference.fa
bwa bwasw oxford_reference.fa oxford.fasta > oxford.sam
python get_position_from_sam.py oxford.sam oxford.fasta
```
A file reads_position.csv containing the position of the reads in the same order as in the fasta file is created in the working directory.
)


Now run the main program, passing the fasta file (```sh -f```), minimap overlap file (```-m```), the instruction to write input files for [poa] [poa] (used to derive the consensus sequences) with the flag (```sh -w```) , setting verbosity to high (```sh -vv```) and providing the path to the csv file containing the position of the reads found by mapping with BWA if it was computed with ```sh --ref_pos_csvf```.
```sh
python spectral_layout_from_minimap.py -f oxford.fasta -m oxford.mini -w -vv --ref_pos_csvf reads_position.csv
```



[minimap] : https://github.com/lh3/minimap
[nanocorrect] : https://github.com/jts/nanocorrect/
[poa] : https://sourceforge.net/projects/poamsa/
[bwa] : http://bio-bwa.sourceforge.net/

# spectrassembler
Created by Antoine Recanati at INRIA, Paris.

## Introduction
We provide the code used in [our paper] (http://arxiv.org/abs/1609.07293) (https://doi.org/10.1093/bioinformatics/btx370) to obtain the layout of uncorrected long-reads from the overlaps in an original and efficient way. It takes as input an overlap file (from [minimap][minimap]) and uses a spectral algorithm to find an ordering of the reads from the pairwise overlap information. This ordering should roughly match the sorting of the reads by position on the genome. This ordering can be viewed as a coarse-grained layout, from which a more refined one is computed with the overlap information.
The pipeline also includes a consensus stage inspired from [nanocorrect][nanocorrect]) ((c) 2015 Ontario Institute for Cancer Research) which uses [spoa][spoa]) to compute multiple sequence alignments in a sliding window along the genome.

This tool is experimental and does not replace mainstream assemblers. Specifically, it does not take in account the topology of the overlap graph to check if the layout is consistent. A sensitive parameter is the initial threshold on the similarity matrix. The default is to set the 40% lowest values of the similarity matrix to zero, in order to remove false overlaps. This should work fine with Oxford Nanopore data of bacterial genomes with coverage not too high. For eukaryotic genomes or Pacbio data with high coverage (>80x), it is better to set 90% of the values to zero with the option --sim_qtile 0.9 (although there is a heuristic to check the layout and increase this parameter automatically before rerunning if it looks wrong, at the expense of runtime).

## Dependencies
* overlapper : [minimap] (https://github.com/lh3/minimap)
To install (require gcc and zlib):
```sh
git clone https://github.com/lh3/minimap && (cd minimap && make)
```

* python packages : numpy, scipy, [biopython][biopython]. Easy install with pip for example : ```pip install biopython```


## Walking through it
We follow the main steps of the pipeline in the following to get started. The code can also be used in a more black-box way with the shell script spectrassembler_pipeline.sh (see next section).

Get the code and save path to it to use scripts later
```sh
git clone https://github.com/antrec/spectrassembler
cd spectrassembler && srcd=`pwd` && git submodule init && git submodule update && \
cd tools/spoa && git submodule init && git submodule update && make
cd ../../../
```
Create a working directory and download data. The Loman lab released a set of E. Coli Oxford Nanopore reads (see their [page] (http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/)) which can be downloaded from the command line:
```sh
mkdir oxford-test && cd oxford-test
curl -L -o oxford.fasta http://nanopore.s3.climb.ac.uk/MAP006-PCR-1_2D_pass.fasta
```
Compute alignments with [minimap][minimap] (you may have to specify the full path to minimap)
```sh
minimap -Sw5 -L100 -m0 oxford.fasta oxford.fasta > oxford.mini.paf
```

Now run the main program, passing the fasta file (```-f```), minimap overlap file (```-m```), setting verbosity to high (```-vvv```), creating files in a temporary directory (```-r temp```) and running on 8 threads (```--nproc 8```)
```sh
python $srcd/spectrassembler.py -f oxford.fasta -m oxford.mini.paf -r temp -vvv --nproc 8 > contigs.fasta
```
This constructs a similarity matrix from the overlaps, with a threshold on the number of matches found for an overlap (this threshold can be modified with ```--sim_qtile```, default is to remove the 40% lowest values (--sim_qtile 0.4)). The layout is computed in each of the connected component of the similarity graph, and written to a file ```cc%d.layout``` where %d is the number of the connected component. For each connected component %d, a subdirectory ```./cc%d``` is created, containing input files for [spoa][spoa] in order to compute a consensus sequences in a sliding window. The contigs are written to stdout.


## Usage
* spectrassembler.py
```
python spectrassembler.py -f reads.fasta -m overlaps.mini
[-h (--help)]
[-f : file containing reads in FASTA or FASTQ format]
[-m : overlap file from minimap in PAF format]
[-r (--root) : root directory where to write layout files (default "./")]
[--sim_qtile <float(0.4)> : Threshold on overlap score (similarity matrix preprocessing) is set as the value of quantile(sim_qtile) of the values of the similarity matrix entries. This is a
important parameter that should be increased to 0.9 for repetitive genomes (e.g., eukaryotic) or with Pacbio data and high coverage. The higher it is, the fewer chances to have a erroneous layout but the more fragmented the assembly will be.]
[--spoapath <string(tools/spoa/spoa)> : path to spoa executable. If it is set to a non existent path (e.g, --spoapath ''), the consensus will not be performed and the program will stop after the layout computation.]
[-v verbosity level (-v, -vv, -vvv or none), default -v]
[--nproc <int(1)> : number of threads]
```

[minimap]: https://github.com/lh3/minimap
[nanocorrect]: https://github.com/jts/nanocorrect/
[spoa]: https://github.com/rvaser/spoa
[bwa]: http://bio-bwa.sourceforge.net/
[biopython]: http://biopython.org/wiki/Download#Easy_Install

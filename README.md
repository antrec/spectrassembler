# spectrassembler
Created by Antoine Recanati at INRIA, Paris.

## Introduction
We provide the code used in [our paper] (http://arxiv.org/abs/1609.07293) to obtain the layout of uncorrected long-reads from the overlaps in an original and efficient way. It takes as input an overlap file (from [minimap][minimap]) and uses a spectral algorithm to find an ordering of the reads from the pairwise overlap information. This ordering should roughly match the sorting of the reads by position on the genome. This ordering can be viewed as a coarse-grained layout, from which a more refined one is computed with the overlap information.
The pipeline also includes a consensus stage inspired from [nanocorrect][nanocorrect]) ((c) 2015 Ontario Institute for Cancer Research) which uses [spoa][spoa]) to compute multiple sequence alignments in a sliding window along the genome.

This tool is experimental and does not replace mainstream assemblers. Specifically, it does not take in account the topology of the overlap graph to check if the layout is consistent. A sensitive parameter is the initial threshold on the similarity matrix. The default is to set the 40% lowest values of the similarity matrix to zero, in order to remove false overlaps. This should work fine with Oxford Nanopore data of bacterial genomes with coverage not too high. For eukaryotic genomes or Pacbio data with high coverage (>80x), it is better to set 90% of the values to zero with the option --sim_qtile 0.9 (although there is a heuristic to check the layout and increase this parameter automatically before rerunning if it looks wrong, at the expense of runtime).

## Dependencies
* overlapper : [minimap] (https://github.com/lh3/minimap)
To install (require gcc and zlib):
```sh
git clone https://github.com/lh3/minimap && (cd minimap && make)
```

* python packages : numpy, scipy, [biopython][biopython]. Easy install with pip for example : ```pip install biopython```

* multiple sequence aligner : [spoa][spoa] (for consensus generation only. Not needed to compute the layout)


## Walking through it
We follow the main steps of the pipeline in the following to get started. The code can also be used in a more black-box way with the shell script spectrassembler_pipeline.sh (see next section).

Get the code and save path to it to use scripts later
```sh
git clone https://github.com/antrec/spectrassembler
cd spectrassembler && srcd=`pwd` && cd ../
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

Now run the main program, passing the fasta file (```-f```), minimap overlap file (```-m```), the path to [spoa][spoa] (used to derive the consensus sequences), setting verbosity to high (```-vvv```) and creating files in a temporary directory (```-r temp```)
```sh
python $srcd/spectrassembler.py -f oxford.fasta -m oxford.mini.paf \
-r temp -vvv --spoa_path
```
This constructs a similarity matrix from the overlaps, with a threhsold on the number of matches found for an overlap (this threshold can be modified with ```--sim_thr```, default value is 850). The layout is computed in each of the connected component of the similarity graph, and written to a file ```cc%d.layout``` where %d is the number of the connected component. As the option ```-w``` is specified, for each connected component %d, a subdirectory ```./cc%d``` is created, containing input files for [spoa][spoa] in order to compute a consensus sequences in a sliding window. If a file containing the position of the reads was given with ```--ref_pos_csvf```, scatter plots are generated, showing the position found by our algorithm vs the position found by mapping to the reference (for each read mapped in a given connected component).

The layout being computed, you can generate consensus sequences in each connected component. First we compute a consensus per window, which can be done in parallel or on a cluster, depending on the number of nodes you have. Let us say you have 4 CPUs and wish to compute a consensus sequence for the third largest connected component, you can use the following shell script or something similar adapted to your cluster (poa must be on your path, if not, add it with ```PATH=$PATH:/path/to/poa/```).
```sh
NUM_CC=3
NUM_NODES=4
SCORE_MAT="$srcd/poa-score.mat" #score matrix for the multiple sequence alignment
cd ./cc_$NUM_CC
for file in poa_in_cc_*.fasta
do
  while [`jobs | wc -l` -ge $NUM_NODES ]
  do
    sleep 2
  done
  if ! test -f $file.clustal.out #Do not recompute the same thing twice in case you stopped a computation earlier.
  then
    poa -read_fasta $file -clustal $file.clustal.out -hb $SCORE_MAT &
  fi
done
cd ../
```
Then you can join together the windows with
```sh
python $srcd/gen_cons_from_poa.py -cc $NUM_CC --poa_mat_path $SCORE_MAT -vv
```
which will create a file ```consensus_cc_3.fasta``` in the current directory (oxford-test).

## Usage
What's in the box :
* Main python script that computes layout
```
python spectral_layout_from_minimap.py -f reads.fasta -m overlaps.mini
[-h (--help)]
[-f : file containing reads in FASTA format]
[-m : overlap file from minimap in PAF format]
[-r (--root) : root directory where to write layout files (default "./")]
[-w (--write_poa_files) : Whether to write POA input files for consensus generation or not.]
[--w_len <int(2500)> : length of consensus windows for POA]
[--w_ovl_len <int(1250)> : overlap length between two successive consensus windows]
[--sim_thr <int(850)> : threshold on overlap score (similarity matrix preprocessing). This is a
crucial parameter that should be increased for repetitive genomes (e.g., eukaryotic) or if
the results are unsatisfactory. Conversely, if the assembly is too fragmented (too many contigs),
it can be decreased. It should also be modified according to the overlapper used
(this value was chosen when using minimap).]
[--len_thr <int(3500)> : threshold on the length of the overlap (similarity matrix preprocessing)]
[--ref_pos_csvf : csv file (generated with get_position_from_sam.py)
with position of reads (in same order as in the fasta file)
obtained from BWA, in order to plot reads position found vs reference.)]
[-v verbosity level (-v, -vv or none), default none]
```

* Python script to compute consensus after ```spectral_layout_from_minimap.py``` was ran with the ```-w```option
```
python gen_cons_from_poa.py -cc 3 --poa_mat_path /path/to/poa-score.mat -vv
[-h (--help)]
[-cc (--contig) : index of contig you wish to compute consensus for]
[--poa_mat_path : path to score matrix file for alignment]
[--poa_path : path to poa executable if it is not on your path
(do not specify this option if poa is on your path)]
[-r (--root) : root directory where to write layout files (default "./")]
[--w_len <int(2500)> : length of consensus windows for POA.
! MUST BE THE SAME VALUE AS IN spectral_layout_from_minimap.py !]
[--w_ovl_len <int(1250)> : overlap length between two successive consensus windows
! MUST BE THE SAME VALUE AS IN spectral_layout_from_minimap.py !]
[-v verbosity level (-v, -vv or none), default none]
```

* Python script to get position of the reads in a csv file from .sam file after mapping
```
get_position_from_sam.py mapping.sam reads.fasta
```

* And the shell script ```spectrassembler_pipeline.sh``` which can be used to perform the full pipeline.
Modify the definitions with your own paths and file names in the beginning of the script and then run ```/bin/bash spectrassembler_pipeline.sh```.

[minimap]: https://github.com/lh3/minimap
[nanocorrect]: https://github.com/jts/nanocorrect/
[spoa]: https://github.com/rvaser/spoa
[bwa]: http://bio-bwa.sourceforge.net/
[biopython]: http://biopython.org/wiki/Download#Easy_Install

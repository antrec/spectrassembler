# spectrassembler
Doc and full code to come (in a couple of days)

## Introduction
We provide the code used in [our paper] (http://arxiv.org/abs/1609.07293) to obtain the layout of uncorrected long-reads from the overlaps in an original and efficient way. It takes as input an overlap file (from [minimap][minimap]) and uses a spectral algorithm to find an ordering of the reads from the pairwise overlap information. This ordering should roughly match the sorting of the reads by position on the genome. This ordering can be viewed as a coarse-grained layout, from which a more refined one is computed with the overlap information.
The pipeline also includes a consensus stage inspired from [nanocorrect][nanocorrect]) ((c) 2015 Ontario Institute for Cancer Research) which uses [poa][poa]) to compute multiple sequence alignements in a sliding window along the genome.

This tool is experimental and does not replace mainstream assemblers. Specifically, there is no sanity check preventing the layout found to be totally wrong for a bad choice of parameters (a warning is raised when the ordering is "suspicious", though). The sensitive parameter is the threshold on the similarity matrix, --sim_thr, which can be increased if a warning appears. However, the default value (--sim_thr 850) should work for bacterial genomes and allow for the layout to be computed cheaply, and for a high-quality draft genome to be derived.

## Dependencies
* overlapper : [minimap] (https://github.com/lh3/minimap)
To install (require gcc and zlib):
```sh
git clone https://github.com/lh3/minimap && (cd minimap && make)
```

* python packages : numpy, scipy, [biopython][biopython]. Easy install with pip for example : ```pip install biopython```

* multiple sequence aligner : [poa][poa] (for consensus generation only. Not needed to compute the layout)

* (optional) mapping to reference : [bwa][bwa]. Needed only if you have a reference genome available and wish to plot the layout found by our algorithm vs the one found by mapping against the reference genome.

## Walking through it
We follow the main steps of the pipeline in the following to get started. The code can also be used in a more black-box way with the shell script main.sh (see next section).

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
Compute alignments with [minimap][minimap]
```sh
minimap -S oxford.fasta oxford.fasta > oxford.mini
```
(Optional : If you want to check the layout found against the reference genome, you can download a FASTA file of the reference genome of E. Coli K-12 substr. MG1655. [here] (https://www.ncbi.nlm.nih.gov/nuccore/556503834). If you saved it under the name oxford_reference.fa then you can run [bwa][bwa] to map the reads against the reference :
```sh
bwa index oxford_reference.fa
bwa bwasw oxford_reference.fa oxford.fasta > oxford.sam
python $srcd/get_position_from_sam.py oxford.sam oxford.fasta #Run the python file from spectrassembler folder
```
A file reads_position.csv containing the position of the reads in the same order as in the fasta file is created in the working directory.
)


Now run the main program, passing the fasta file (```-f```), minimap overlap file (```-m```), the instruction to write input files for [poa][poa] (used to derive the consensus sequences) with the flag ```-w```, setting verbosity to high (```-vv```), providing the path to the csv file containing the position of the reads found by mapping with BWA if it was computed with ```--ref_pos_csvf``` and creating files in the current directory (```-r not specified```)
```sh
python $srcd/spectral_layout_from_minimap.py -f oxford.fasta -m oxford.mini -w -vv --ref_pos_csvf reads_position.csv
```
This constructs a similarity matrix from the overlaps, with a threhsold on the number of matches found for an overlap (this threshold can be modified with ```--sim_thr```, default value is 850). The layout is computed in each of the connected component of the similarity graph, and written to a file ```cc%d.layout``` where %d is the number of the connected component. As the option ```-w``` is specified, for each connected component %d, a subdirectory ```./cc%d``` is created, containing input files for [poa][poa] in order to compute a consensus sequences in a sliding window. If a file containing the position of the reads was given with ```--ref_pos_csvf```, scatter plots are generated, showing the position found by our algorithm vs the position found by mapping to the reference (for each read mapped in a given connected component).

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

## Manual
Main python script that computes layout
```sh
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

Python script to compute consensus after ```spectral_layout_from_minimap.py``` was ran with the ```-w```option :
```sh
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

Python script to get position of the reads in a csv file from .sam file after mapping
```sh
get_position_from_sam.py mapping.sam reads.fasta
```

[minimap]: https://github.com/lh3/minimap
[nanocorrect]: https://github.com/jts/nanocorrect/
[poa]: https://sourceforge.net/projects/poamsa/
[bwa]: http://bio-bwa.sourceforge.net/
[biopython]: http://biopython.org/wiki/Download#Easy_Install

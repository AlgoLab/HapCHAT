# HapCHAT

Adaptive haplotype assembly for efficiently leveraging high coverage
in long reads

## <a name="cite"></a> Citation ##

A description of the algorithm, as well as a detailed comparison
experiment with other haplotype assembly tools is presented in:

Stefano Beretta*, Murray Patterson*, Simone Zaccaria, Gianluca Della
Vedova and Paola Bonizzoni.  _HapCHAT: Adaptive haplotype assembly for
efficiently leveraging high coverage in long reads_.  bioRxiv 170225.
*_Joint first authors_

DOI: [10.1101/170225](https://doi.org/10.1101/170225)

To replicate the paper: TODO

## Quick Install

If you have docker installed, to install and run HapCHAT you only have to run `docker run
-v DATADIR:/data algolab/hapchat` where `DATADIR` is a directory with the input data files
`genome.fasta`, `file.bam`, `file.vcf`. If `DATADIR` is empty, then HapCHAT is run on the
example files in the `example` directory.

### Input and output files

*  `genome.fasta`

*  `file.bam`: it is the input file containing all aligned reads in BAM format

*  `file.vcf`

*  `out.realigned.phased.vcf`

*  `out.phased.vcf`

### <a name="install"></a> Installation for Experts

HapCHAT has been developed and tested on Ubuntu Linux, but should work
any system which has python(3), C++(>=11), as well other utilities
that appear on most *nix-based systems (such as bash, awk, git, cmake
and make)

Some more specific dependencies that may not be installed are
`python3-dev`, `python3-networkx` and `virtualenv`.  These can be
obtained in Ubuntu with, _e.g._, the command `apt install
python3-dev`, etc.

Then, in principle, one needs to simply execute `setup.sh`, which is
located in the same directory as this README, and then HapCHAT can be
run by executing `HapCHAT.py` (located in this directory as well)

_Note:_ that `setup.sh` simply checks out a HapCHAT-specific git
branch of WhatsHap, installs it in a virtual environment, and then
builds with cmake the C++ code located in `src/`.  While this should
work automagically for most *nix-based systems with bash, solutions
could be found for other systems by slightly modifying `setup.sh`.

The main advantage of the expert installation is that you can adapt the `Snakefile` in the
`example` directory to run on your files, using custom filenames and directory tree.


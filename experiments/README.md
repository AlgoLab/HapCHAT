# Replicating the experiments #

This directory contains the instructions for reproducing the
experiments described in our study: "HapCHAT: Adaptive haplotype
assembly for efficiently leveraging high coverage in long reads", for
which a preprint is available at
[https://doi.org/10.1101/170225](https://doi.org/10.1101/170225).

The workflow is described as a "Snakefile", which needs to be run with
[snakemake](http://snakemake.bitbucket.org).

## Overview

The workflow is designed to be run from the directory in which this
README file is located.  In principle, you only need to run
`snakemake` here and then the rest is done automatically, but all
dependencies need to be installed first.  Further things to know:

* Approximately 150 GB of hard disk space is needed (better 200 GB to
  be on the safe side).

* No data is included in this repository since all required files will
  be generated (which take up said 150 GB).

* 8 cores should work, but more is better. Presumably 8 GB RAM works,
  but more is better.  We tested on a 32-core machine with 256 GB RAM.

## Obtaining/generaing the input data

More precisely, this workflow replicates what is described in the
subsection "Expermental setup" of the above "HapCHAT" study.  The
data, described in the subsection "Data description", which is used by
the former, must be generated first.

In order to do this, go to the child directory `/data` of this current
directory and follow instructions there for generating the input VCF
and BAM files required byr this workflow.

## Installing dependencies

Install the necessary system packages that the different phasing
softwares of this workflow depend on.  Some of these dependencies have
already been installed as a prerequisite to runng the workflow for
obtaining/generating the input data -- so we assume that this workflow
has finished successfully.

We also assume you use Debian or Ubuntu.  If you do not have root
access on your system, you can try to ignore this step.  These
packages are typically already installed.

    sudo apt-get install default-jdk
    conda install -y virtualenv networkx

Python 2, and the weave module (supported by python<=2.6) are needed
to run ProbHap.

    sudo apt-get install python python-pip
    sudo pip install weave

## Run the workflow

Start the workflow with

    nice snakemake -p -j 16

Adjust the 16, which is the number of cores you want to use.

## Results

The resulting files for each tool/dataset will be reported in the
`output` directory (grouped by tool), as

* a `.phased.vcf` file: the phasing produced by the tool

* a whatshap diff file (`.diff`) on this phased VCF file containing
  the accuracy measures

* a `.time` file for the output of `/usr/bin/time` for time and peak
  memory used to produce this phased VCF

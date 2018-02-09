# Replicating the experiments #

This directory contains a workflow, described as a "Snakefile", which
needs to be run with [snakemake](http://snakemake.bitbucket.org).

## Overview

The workflow is designed to be run from the directory in which this
README file is located.  In principle, you only need to run
`snakemake` here and then the rest is done automatically, but all
dependencies need to be installed first.  Further things to know:

* About 150 GB of hard disk space is needed (better 200GB to be on the
  safe side)

* No data is included in this repository since all required files will
  be generated (which take up said 150GB)

* 8 cores should work, but more is better. Presumably 8 GB RAM works,
  but more is better.  We tested on a 32-core machine with 256 GB RAM.

## Generaing the input data

One must first go to the child directory `/data` of this current
directory and follow instructions there for generating the input BAM
and VCF for this workflow

## Installing dependencies

Install necessary system packages that the different phasing softwares
depend on. We assume you use Debian or Ubuntu. If you do not have root
access on your system, you can try to ignore this step.  These
packages are typically already installed.

    snakemake (conda install -c bioconda snakemake)
    pysam (conda install -c bioconda pysam)
    virtualenv (conda install virtualenv)
    javac (sudo apt install default-jdk)

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

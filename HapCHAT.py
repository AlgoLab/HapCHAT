#!/usr/bin/env python

description = '''

   HapCHAT: Adaptive haplotype assembly for efficiently leveraging high coverage in long reads

'''

import sys
import os
import argparse
import subprocess
import logging

#
# preamble
#----------------------------------------------------------------------

dir = os.path.dirname(os.path.realpath(__file__))
wh_dir = '{}/software/whatshap'.format(dir)

whatshap = '{}/venv/bin/whatshap'.format(wh_dir)

#
# functions
#----------------------------------------------------------------------

# auxiliary function to add arguments to the argument parser
def add_arguments(parser) :

    arg = parser.add_argument

    # positional arguments
    arg('vcf', metavar = 'VCF',
        help = 'VCF file with variants to be phased')
    arg('bam', metavar = 'BAM',
        help = 'BAM file of sequencing reads')

    # optional arguments
    arg('--reference', '-r', metavar = 'FASTA',
        help = 'reference file, fo detecting alleles in realignment mode')
    arg('--max-coverage', '-H', metavar = 'MAXCOV', default = 15, type = int,
        help = 'downsample coverage to at most MAXCOV')


# shortcut
def shell(command, workdir = None) :
    
    subprocess.run(command.split(), cwd = workdir)


# for setting up whatshap
def setup_whatshap() :

    # place the git repo
    shell('mkdir -p {}'.format(wh_dir), dir)
    shell('''

  git clone --branch hapcol-experiments
    https://bitbucket.org/whatshap/whatshap.git {}

    '''.format(wh_dir), dir)

    # setup virtualenv
    shell('virtualenv -p python3 venv', wh_dir)
    shell('venv/bin/pip3 install Cython nose tox', wh_dir)
    shell('venv/bin/pip3 install -e .', wh_dir)


# use whatshap to read in a bam file
def read_bam(vcf, bam, reference) :

    # setup the parameters
    realignment = ''
    if reference :
        realignment = '--reference {}'.format(reference)

    rawreal = 'realigned' if realignment else 'raw'
    wif = '{}.hx.{}.wif'.format(bam, rawreal)
    out = '{}.transcript'.format(wif)
    err = '{}.log'.format(wif)

    # run whatshap
    subprocess.run('''

  {} phase -o /dev/null {} --output-wif {} -H 1000 {} {}

    '''.format(whatshap, realignment, wif, vcf, bam).split(),
                   stdout = open(out,'w'),
                   stderr = open(err,'w'))

    # wif file for the next step
    return wif


#
# main
#----------------------------------------------------------------------
def main(argv = sys.argv[1:]) :

    # check for whatshap, installing if necessary
    if not os.path.exists(whatshap) :
        setup_whatshap()

    # parse arguments
    parser = argparse.ArgumentParser(prog = description.split()[0].rstrip(':'),
                                     description = description,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    add_arguments(parser)
    args = parser.parse_args(argv)

    # read vcf / bam pair, returning (location of) resulting wif
    wif_file = read_bam(args.vcf, args.bam, args.reference)


if __name__ == '__main__' :
    main()

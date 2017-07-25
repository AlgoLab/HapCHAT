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
hx_dir = '{}/software/hapchat'.format(dir)
scr_dir = '{}/scripts'.format(dir)

whatshap = '{}/venv/bin/whatshap'.format(wh_dir)
hapchat = '{}/hapchat'.format(hx_dir)

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

    # output
    arg('--output', '-o', metavar = 'OUTPUT', default = 'standard output',
        help = 'output VCF file with phase information')

    # realignment
    arg('--reference', '-r', metavar = 'FASTA',
        help = 'reference file, fo detecting alleles in realignment mode')

    # merging reads
    arg('--thr', '-t', metavar = 'THRESHOLD', default = 6, type = int,
        help = 'threshold for the merging step')
    arg('--neg_thr', '-n', metavar = 'NEG_THRESHOLD', default = 3, type = int,
        help = 'negative threshold for the merging step')
    arg('--error_rate', '-e', metavar = 'ERROR_RATE', default = 0.15, type = float,
        help = 'probability that a site is wrong')
    arg('--max_err', '-m', metavar = 'MAX_ERR', default = 0.25, type = float,
        help = 'maximum error rate for a site')

    # downsampling
    arg('--seed', '-s', metavar = 'SEED', default = 1, type = int,
        help = 'seed for psuedorandom downsampling')
    arg('--max_coverage', '-H', metavar = 'MAXCOV', default = 15, type = int,
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

# for setting up hapchat
def setup_hapchat() :

    # add the hapchat dir and build the source
    shell('mkdir -p {}'.format(hx_dir), dir)
    shell('cmake ../../src', hx_dir)
    shell('make -j 8', hx_dir)


# use whatshap to read in a bam file
def read_bam(vcf, bam, reference) :

    # setup realignment mode
    realignment = ''
    if reference :
        realignment = '--reference {}'.format(reference)

    # set up dir for intermediate output
    int_dir = '{}/.{}.hx_'.format(os.path.dirname(bam), os.path.basename(bam))
    shell('mkdir -p {}'.format(int_dir))

    # run whatshap
    rawreal = 'realigned' if realignment else 'raw'
    wif = '{}/{}.wif'.format(int_dir, rawreal)
    out = '{}.transcript'.format(wif)
    log = '{}.log'.format(wif)
    subprocess.run('''

  {} phase -o /dev/null {} --output-wif {} -H 1000 {} {}

    '''.format(whatshap, realignment, wif, vcf, bam).split(),
                   stdout = open(out,'w'),
                   stderr = open(log,'w'))

    # return wif file (location)
    return wif


# the merging reads step
def merge_reads(wif, e, m, t, n) :

    # setup the parameters
    pe = str(e).split('.')[1]
    pm = str(m).split('.')[1]
    out = '{}.merged_e{}_m{}_t{}_n{}.wif'.format(wif, pe, pm, t, n)
    graph = '{}.graph'.format(out)
    log = '{}.log'.format(out)

    # merge the reads
    subprocess.run('''

  python {}/rb-merge.py -e {} -m {} -t {} -n {} -w {} -o {} -g {}

    '''.format(scr_dir, e, m, 10**t, 10**n, wif, out, graph).split(),
                   stdout = open(log,'w'),
                   stderr = subprocess.STDOUT)

    # return merged wif file
    return out


# the random downsampling step
def downsample(wif, seed, maxcov) :

    # seeded pseudorandom shuffle of the lines of wif
    shuffle = '{}.lines.shuf{}'.format(wif, seed)
    shell('bash {}/pseudorandomshuffle.bash {} {}'.format(scr_dir, wif, seed))

    # greedily downsample wif to coverage according to shuffle
    sample = '{}.sample_s{}_m{}'.format(wif, seed, maxcov)
    log = '{}.log'.format(sample)
    subprocess.run('''

  python {}/wiftools.py -s {} {} {}

    '''.format(scr_dir, maxcov, shuffle, wif).split(),
                   stdout = open(sample,'w'),
                   stderr = open(log,'w'))

    # extract this sample
    downs = '{}.downs.s{}.m{}.wif'.format(wif, seed, maxcov)
    shell('bash {}/extractsample.bash {} {} {} {}'.format(scr_dir, wif, sample, seed, maxcov))

    # return downsampled wif file
    return downs


# run hapchat
def run_hapchat(wif) :

    # run hapchat with default parameters
    hap = '{}.hap'.format(wif)
    log = '{}.log'.format(hap)
    subprocess.run('''

  {} -i {} -o {} -A -e 0.05 -a 0.1 -b 1000 -r 0

    '''.format(hapchat, wif, hap).split(),
                   stdout = open(log,'w'),
                   stderr = subprocess.STDOUT)

    # return resulting haplotypes
    return hap

# convert hapchat output to phased vcf
def phase_vcf(hap, wif, vcf) :

    # obtain blocks of the instance
    blocks = '{}.info_/block_sites_'.format(wif)
    shell('python {}/wiftools.py -i {}'.format(scr_dir, wif))

    # phase vcf with blocks and haplotype info
    phased_vcf = '{}.vcf'.format(hap)
    log = '{}.log'.format(phased_vcf)
    subprocess.run('''

  python {}/subvcf.py -p {} {} {}

    '''.format(scr_dir, hap, blocks, vcf).split(),
                   stdout = open(phased_vcf,'w'),
                   stderr = open(log,'w'))

    # return phased vcf (location)
    return phased_vcf


#
# main
#----------------------------------------------------------------------
def main(argv = sys.argv[1:]) :

    # check for whatshap, installing if necessary
    if not os.path.exists(whatshap) :
        setup_whatshap()

    # check for hapchat, building if necessary
    if not os.path.exists(hapchat) :
        setup_hapchat()

    # parse arguments
    parser = argparse.ArgumentParser(prog = description.split()[0].rstrip(':'),
                                     description = description,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    add_arguments(parser)
    args = parser.parse_args(argv)

    # read vcf/bam
    wif = read_bam(args.vcf, args.bam, args.reference)

    # merge reads
    merged_wif = merge_reads(wif, args.error_rate, args.max_err,
                             args.thr, args.neg_thr)

    # random downsampling
    downs_wif = downsample(merged_wif, args.seed, args.max_coverage) 

    # run hapchat
    hap = run_hapchat(downs_wif)

    # obtain phased vcf
    phased_vcf = phase_vcf(hap, wif, args.vcf)

    # output phased vcf to stdout (or output file, if specified)
    output = sys.stdout
    if args.output != 'standard output' :
        output = open(args.output,'w')

    for line in open(phased_vcf,'r') :
        print(line, file = output, end = '')


if __name__ == '__main__' :
    main()

#!/usr/bin/env python3

usage = '''

usage: python subvcf.py [-v] [-H] [-p P S] [file.vcf]

perform various functions on a (-n unphased) vcf file file.vcf, such
as gather all (single nucleotide) variants (SNVs), or add phase
information according to a phasing of these SNVs.  In all cases,
output goes to stdout

  -v          obtain all SNVs of file.vcf
  -H          obtain only the heterozygous SNVs of file.vcf
  -p P S      add phase information to (unphased) file.vcf according
              to phasing P of the (SNV) sites, where S are these sites
              as they appear in the haplotype blocks (as given in
              block_sites_ of wiftools -i of the appropriate wif file)

note: that vcf file must be unphased, of a single contig only (o.w.,
error), and that only the first sample is of concern (i.e., fields 9
and 10 only) with regards to the phasing done, or in determining
whether a site is heterozygous (in the -H context)

'''

import sys

#
# load an old-skool (2 lines) haplotype (phasing)
def load_haplotypes(lines) :

    alpha = ['0', '1', 'X', '|']

    h0 = lines.readline().strip()
    h1 = lines.readline().strip()

    if not h0 : # haplotype file is empty, exit quietly
        print('haplotype file is empty, exiting ..', file = sys.stderr)
        exit(0)

    for a in h0 :
        assert a in alpha, 'unexpected character '+a+' in first haplotype'
    for a in h1 :
        assert a in alpha, 'unexpected character '+a+' in second haplotype'

    assert len(h0) == len(h1), 'haplotypes differ in length'

    return h0, h1

#
# load the (haplotype blocks of) sites to which a phasing corresponds,
# and return the associated phase set info
def load_phase_set_info(lines) :

    phase_set = {} # dic[snv] = phase set of snv
    prev_key = 0
    for line in lines :

        if line.startswith('#') : # burn through header
            continue
        
        s = line.split()[1:]
        ps = int(s[0])
        
        for site in s :
            key = int(site)

            # sanity checks to ensure sites are strictly increasing
            assert key not in phase_set
            assert key > prev_key
            prev_key = key

            phase_set[key] = ps

    return phase_set

#
# set up phasing for vcf output: associating it with phase set info
def setup_phasing(h0, h1, phase_set) :

    phasing = {}
    sites = sorted(phase_set)
    print('#site h0 h1', file = sys.stderr) # log homozgyous sites

    i = 0
    for a,b in zip(h0, h1) :

        site = sites[i]
        ps = phase_set[site]

        if a == '|' : # should the haplotypes contain |'s
            assert b == '|'

            ps_prev = phase_set[sites[i-i]]
            assert ps > ps_prev, 'haplotypes discordant with blocks of S'

            continue

        i += 1 # o.w., we're at a site in the haplotypes

        if a == 'X' or b == 'X' : # but phase cannot be decided
            continue

        # log the homozygous site, don't give an error, but treat it
        # as another site where phase cannot be decided
        if a == b  :
            print(site, a, b, file = sys.stderr)
            continue

        phasing[site] = {
            'phasing' : '{}|{}'.format(a, b),
            'ps'      : ps
        }

    assert i == len(sites), 'number of sites different than haplotype length'

    return phasing

#
# PARSER
#

variantsmode = False
heterozygousmode = False
phasingmode = False
h0,h1 = None, None
phase_set = None
phasing = None
entree = sys.stdin

a = sys.argv[1:]
i = 0
while i < len(a) : # bash-style argparse

    if a[i].startswith('-') :
        if a[i] == '-v' :
            variantsmode = True
        elif a[i] == '-H' :
            variantsmode = True
            heterozygousmode = True
        elif a[i] == '-p' :
            phasingmode = True
            i += 1
            h0,h1 = load_haplotypes(open(a[i],'r'))
            i += 1
            phase_set = load_phase_set_info(open(a[i],'r'))
            phasing = setup_phasing(h0, h1, phase_set)
        else :
            assert False, usage
    else :
        entree = open(a[i],'r')

    i += 1 # shift once by default in any case

#
# MAIN
#

unphased = ['0/0', '0/1', '1/1', '1/0', '1/2', '2/1'] # (un-) phasings to expect
heterozygous = ['0/1', '1/0', '1/2', '2/1']

contig = ''
snvs = set([])
for line in entree :

    if line.startswith('#') : # burn through header
        if phasingmode :
            sys.stdout.write(line)
        continue

    # now we're in the body ...
    s = line.split()

    if contig : # ensure that only one contig is represented
        assert contig == s[0], 'vcf file concerns more than one contig'
    contig = s[0]

    gt = s[8].split(':',1)[0]
    sample = s[9].split(':',1)
    unphasing = sample[0]
    r_sample = None
    if len(sample) > 1 :
        r_sample = sample[1]

    # ensure that position is properly annotated, unphased (and 'known')
    assert gt == 'GT', s[8] + ' is an unexpected annotation'
    assert unphasing in unphased, s[9] + ' is phased or unknown'

    # for now we consider only SNVs
    ref, alt = s[3:5]
    if len(ref) == len(alt) == 1 :

        pos = int(s[1]) # pos must be a unique integer
        assert pos not in snvs
        snvs.add(pos)

        # variants mode
        if variantsmode :
            if heterozygousmode :
                if unphasing in heterozygous : # snv is heterozygous
                    print(pos, ref, alt, sep = '\t')
            else :
                print(pos, ref, alt, sep = '\t')

        # phasing mode
        if phasingmode :
            if pos in phasing :
                tagged = '{}:PS'.format(s[8])
                phased_sample = phasing[pos]['phasing']
                if r_sample :
                    phased_sample = '{}:{}'.format(phased_sample, r_sample)
                phased = '{}:{}'.format(phased_sample, phasing[pos]['ps'])
                print(*s[:8], tagged, phased, *s[11:], sep = '\t')
            else :
                sys.stdout.write(line)

    else : # it is not an SNV
        if phasingmode :
            sys.stdout.write(line)

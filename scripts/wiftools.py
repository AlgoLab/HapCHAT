#!/usr/bin/env python3

usage = '''

usage: python wiftools.py [-i] [-z] [-r] [-b B] [-c C] [-s S P] [-v V] file.wif

given a (whatshap) wif file file.wif, perform a 'sanity check' on this
input to verify if it is indeed in wif format, and output the
respective info, should the following optional parameters be provided

  -i          output various types of info to files in a directory
              file.wif.info_
  -z          output to stdout the corresponding 'zygosity' matrix (0s
              and 1s), in .mat format
  -r          remove all homozygous sites from the input file.wif, the
              output going to stdout (not reflected in -b or -i)
  -b B        output (to stdout) the block having id B, an int: empty
              output should no block have this id, i.e., is not among
              the blocks in file.wif.info_/blocks_ with the exception
              of 0: should the id be 0, it outputs all of the blocks
              1..N into the respective files: file.wif.blks_/blk_1.wif
              .. file.wif.blks_/blk_N.wif
  -c C        given a set of sets (as C: one per line) of lines of
              file.wif, merge each set of lines together, outputting
              the result to stdout
  -s S P      given an int S and permutation P of the lines of
              file.wif, according to order P, greedily sample the
              reads down to the physical coverage S favouring reads of
              highest cumulative weight first
  -v V        obtain phase information (from the GT field of the first
              sample, i.e., field 10) from a (phased) vcf file V and
              compute the MEC score of this wrt to input file.wif

note: that matrix is unweighted at the moment (something to possibly
add in the future)

'''

import sys
import os
import itertools

#
# format a set of tokens (castable to string) for output
def format(tokens) :
    return ' '.join(str(t) for t in tokens)

#
# load the clusters, sets of lines (one set per line) of file.wif that
# should be merged together
def load_clusters(lines) :

    clusters = {}
    for line in lines :
        s = set([int(x) for x in line.split()])
        r = min(s)

        clusters[r] = s

    return clusters

#
# load a permuation (of a discrete interval)
def load_permutation(lines) :

    count = 0
    permutation = []
    for line in lines :
        permutation.append(int(line.strip()))
        count += 1

    assert len(permutation) == count, 'input not a permutation'

    return permutation

#
# load a phasing from a vcf file
def load_phasing(lines) :

    contig = ''
    phasing = {}
    unphased = []
    for line in lines :

        if line.startswith('#') : # burn through header
            continue

        # now we're in the body ...
        s = line.split()

        if contig : # ensure that only one contig is represented
            assert contig == s[0], 'vcf file concerns more than one contig'
        contig = s[0]

        # we consider only SNVs (since wif files are on SNPs only)
        if len(s[3]) == len(s[4]) == 1 :

            pos = int(s[1]) # pos must be a unique integer
            assert pos not in phasing

            # obtain the phase at this position (cis or trans)
            phase = s[9].split(':')[0]
            if phase == '0|1' :
                phasing[pos] = 0
            if phase == '1|0' :
                phasing[pos] = 1
            # o.w., consider only (unphased) hetero sites
            if phase in ['0/1', '1/0'] :
                phasing[pos] = 0 # default
                unphased.append(pos)

    # note that we ignore the phase group when reading the vcf file.
    # This is because (a) we want to be able to compare to phased vcfs
    # from other methods (that do not concern phase group), and (b) a
    # read should never span more than one phase group anyway, o.w.,
    # they'd be in the same phase group by definition (so we don't
    # check)

    if not phasing : # vcf file is empty, exit quietly
        print('vcf file is empty, exiting ..', file = sys.stderr)
        exit(0)

    return phasing, unphased

#
# PARSER
#----------------------------------------------------------------------

filename = None
info = False
matrixmode = False
zygositymatrix = False
remove = False
blocksmode = False
blockid = 0
clusteringmode = False
clusters = None
samplingmode = False
coverage = 0
permutation = None
mecmode = False
phasing = None
unphased = None

a = sys.argv[1:]
assert len(a), usage
i = 0
while i < len(a) : # bash-style argparse

    if a[i].startswith('-') :
        if a[i] == '-i' :
            info = True
        elif a[i] == '-z' :
            matrixmode = True
            zygositymatrix = True
        elif a[i] == '-r' :
            remove = True
        elif a[i] == '-b' :
            blocksmode = True
            i += 1
            blockid = int(a[i])
        elif a[i] == '-c' :
            clusteringmode = True
            i += 1
            clusters = load_clusters(open(a[i],'r'))
        elif a[i] == '-s' :
            samplingmode = True
            i += 1
            phys_cov_limit = int(a[i])
            i += 1
            permutation = load_permutation(open(a[i],'r'))
        elif a[i] == '-v' :
            mecmode = True
            i += 1
            phasing, unphased = load_phasing(open(a[i],'r'))
        else :
            assert False, usage
    else :
        filename = a[i]

    i += 1 # shift once by default in any case

#
# MAIN
#----------------------------------------------------------------------

# info to gather
reads = {} # dic[read] = length
read_ends = {} # dic[read] = first,last position of read
read_site_allele = {} # dic[read][site] = allele (only stored in matrix mode)
sites = {} # dic[site] = coverage
site_alleles = {} # dic[site] = major and minor allele
site_zygosity = {} # dic[site] = number of 0's, 1's
blocks = set([])
block_reads = {} # dic[block] = reads
block_sites = {} # dic[block] = sites
blockB = '' # block B to output, specified by -b
all_blocks = [] # all blocks, should B = 0
read_site_aq = {} # dic[read][site] = allele,qual (only stored in mecmode)
read_weight = {} # dic[read] = cumulative weight (only stored in sampling mode)
read_sites = {} # dic[read] = sites (only stored in info mode)

superreads = {} # superreads given by the clusters (if clustering)
rep = {} # cluster representative of a read in a cluster
if clusteringmode :
    for r in clusters :
        superreads[r] = {}
        
        for x in clusters[r] :
            rep[x] = r

#
# read the wif file, gathering the info
#----------------------------------------------------------------------
read = 0
block = 0
max_site = -1
startpos = -1 # current start pos
for line in open(filename, 'r') :
    read += 1

    # tokenize line, get first and last site
    tokens = line.split(' : ')[:-2]
    first = int(tokens[0].split()[0])
    last = int(tokens[-1].split()[0])
    read_ends[read] = first, last

    # reads must be sorted by start position
    assert startpos <= first, 'disorder: '+line
    startpos = first

    # determine if we start a new block
    if first > max_site :
        block += 1
        blocks.add(block)

        if not blockid :
            all_blocks.append('')

    if block not in block_reads :
        block_reads[block] = set([])
    if block not in block_sites :
        block_sites[block] = set([])
    if last > max_site :
        max_site = last

    # add read len to reads, its id to block_reads
    reads[read] = len(tokens)
    block_reads[block].add(read)
    if matrixmode :
        read_site_allele[read] = {}
    if info :
        read_sites[read] = set([])

    # set up read_sites and weight dictionary if in sampling mode
    if samplingmode :
        read_weight[read] = 0

    # add line to blockB if the id matches B
    if block == blockid :
        blockB += line

    # or if B = 0, append line to last element of all blocks
    if not blockid :
        all_blocks[block - 1] += line

    # set up read_site_aq if in mec mode
    if mecmode :
        read_site_aq[read] = {}

    # loop through the sites
    currsite = -1 # current site
    for tok in tokens :
        a,b,c,d = tok.split() # must have 4 fields
        site = int(a)
        zyg = int(c)
        qual = int(d)

        # sanity check on fields
        errsuf = 'site: '+str(site)+' in line:\n'+line
        assert c in '0 1'.split(), 'non-bit at '+errsuf
        assert qual >= 0, 'negative qual score at '+errsuf

        # sites on a read must be unique and in order
        assert currsite < site, 'disorder or duplicate '+errsuf
        currsite = site

        # add increment coverage for site
        if site not in sites :
            sites[site] = 0
        sites[site] += 1

        # add to alleles dictionary, checking for discordancy (multi-allelic)
        if site not in site_alleles :
            site_alleles[site] = ['','']
        if not site_alleles[site][zyg] :
            site_alleles[site][zyg] = b
        else :
            deg = 'minor' if zyg else 'major'
            cur = str(site_alleles[site][zyg])
            errsuf = ', current '+deg+' allele: '+cur
            errsuf += '\n\tis discordant with new allele: '+ b
            assert site_alleles[site][zyg] == b, 'at site: '+str(site)+errsuf

        # add to zygosity dictionary
        if site not in site_zygosity :
            site_zygosity[site] = [0,0]
        site_zygosity[site][zyg] += 1

        # add site id to block_sites
        block_sites[block].add(site)
        if info :
            read_sites[read].add(site)

        # add site (resp. weight) to read_sites (resp. read_weight)
        if samplingmode :
            read_weight[read] += qual

        # add read,site = allele to read_sites (if in matrix mode)
        if matrixmode :
            read_site_allele[read][site] = b

            if zygositymatrix :
                read_site_allele[read][site] = c

        # possibly add to a superread (if in clustering mode)
        if clusteringmode :
            if read in rep :
                r = rep[read]

                if site not in superreads[r] :
                    superreads[r][site] = [0,0]

                # add weight of site to representative read
                superreads[r][site][zyg] += qual

        # add site with allele and quality to current read
        if mecmode :
            read_site_aq[read][site] = zyg,qual

#
# done reading the file, output the type of info corresponding to the
# optional parameters provided (if any)
#----------------------------------------------------------------------
if blocksmode :
    if blockid : # output blockB if B was specified
        sys.stdout.write(blockB)
        exit(0)

    else :
        # make a directory to store all the blocks
        blocks_dir = filename+'.blks_'
        if not os.path.exists(blocks_dir) :
            os.makedirs(blocks_dir)

        # now write each block to this directory
        for i in range(len(all_blocks)) :

            handle = open(blocks_dir+'/blk_'+str(i+1)+'.wif','w')
            for line in all_blocks[i] :
                handle.write(line)

            handle.close()

if matrixmode : # output matrix correponding to wif file

    # output by block (to save time, since matrix is so banded)
    for block in sorted(blocks) :
        the_sites = sorted(block_sites[block])

        for read in sorted(block_reads[block]) :
            row = []
            r_site_a = sorted(read_site_allele[read])
            first = r_site_a[0]

            for site in the_sites :
                if site < first or site > r_site_a[-1] :
                    continue

                if site in r_site_a :
                    row.append(read_site_allele[read][site])
                else :
                    row.append('-')

            print(str(read)+'\t'+str(first)+'\t'+' '.join(row))

    exit(0)

if clusteringmode : # output wif file with some rows merged

    read = 0
    for line in open(filename, 'r') :
        read += 1

        if read in rep :
            if read == rep[read] : # output representative superread

                for site in sorted(superreads[read]) :
                    z = superreads[read][site]

                    if z[0] >= z[1] :
                        sys.stdout.write(format([site,site_alleles[site][0],0,z[0]-z[1],':','']))
                    elif z[1] > z[0] :
                        sys.stdout.write(format([site,site_alleles[site][1],1,z[1]-z[0],':','']))

                print('# X : X') # sentinel at end of wif line

        else :
            sys.stdout.write(line)

    exit(0)

if samplingmode : # sample reads down to a specified physical coverage

    assert len(reads) == len(permutation), 'permutation different length than file'

    # add rank (weight) to the permutation, keeping permutation order otherwise
    permutation_rank = {}
    for read in permutation :
        rank = read_weight[read]

        if rank not in permutation_rank :
            permutation_rank[rank] = []

        permutation_rank[rank].append(read)

    # tack all ranks together (keeping permutation order otherwise)
    ranked_permutation = []
    for rank in sorted(permutation_rank, reverse = True) :
        for read in permutation_rank[rank] :
            ranked_permutation.append(read)

    # index the sites for quick indexing with read ends
    sorted_sites = sorted(sites)
    site_i = {}
    i = 0
    for site in sorted_sites :
        site_i[site] = i
        i += 1

    # set up the sampling
    sample_reads = [] # reads to take
    sample_sites = {} # dic[site] = coverage of sampled reads so far
    for site in sites :
        sample_sites[site] = 0

    # greedily sample, according to given ranked permutation
    for read in ranked_permutation :
        sample_read = True # by default
        a,b = read_ends[read]
        i,j = site_i[a], site_i[b]+1

        for site in sorted_sites[i:j] :
            if sample_sites[site] == phys_cov_limit :
                sample_read = False

        if sample_read :
            sample_reads.append(read) # add read to sample

            for site in sorted_sites[i:j] : # update sample sites
                sample_sites[site] += 1
                assert sample_sites[site] <= phys_cov_limit

    # output the reads (lines) in the sample
    for read in sample_reads :
        print(read)

    exit(0)

if mecmode : # compute mec score of a phasing wrt wif and exit

    mec_score = 0

    # compute mec score of each block against phasing
    for block in blocks :

        tophase = []
        for site in block_sites[block] :
            if site in unphased :
                tophase.append(site)

        # now we need to iterate through all possible phasings of
        # these sites, and compute the best mec score for this
        # block (note that this is an exponential phasing method)
        k = len(tophase)
        best_score = sys.maxsize
        for phase in itertools.product([0,1], repeat = k) :

            # patch phasing with current phase for block
            for i in range(k) :
                site = tophase[i]
                phasing[site] = phase[i]

            # compute mec score for this particular patching
            score = 0
            for read in block_reads[block] :

                c0,c1 = 0,0
                for site in read_site_aq[read] :
                    if site in phasing :

                        a,q = read_site_aq[read][site]
                        if a != phasing[site] :
                            c0 += q
                        else :
                            c1 += q

                score += min(c0, c1)

            # keep the best phasing
            if score < best_score :
                best_score = score

        # add best phasing of this block to (global) mec score
        mec_score += score

    print('mec score:', mec_score)
    exit(0)

if remove : # output file.wif with homozygous sites removed

    reads = {} # dic[first] = set of reads with first
    for line in open(filename, 'r') :
        tokens = line.split(' : ')

        pref = ''
        for tok in tokens[:-2] :
            site = int(tok.split()[0])

            if site_zygosity[site][0] and site_zygosity[site][1] :
                pref += tok + ' : '

        if pref : # should there remain heterozygous sites
            index = int(pref.split()[0])
            read = pref + ' : '.join(tokens[-2:]).strip()

            if index not in reads :
                reads[index] = []
            reads[index].append(read)

    for index in sorted(reads) : # removal may change order properties
        for r in reads[index] :
            print(r)

if not info : # exit if -i was not specified, o.w., ...
    exit(0)

# compute the physical (and bridge) coverage at each site
site_physcov = {}
site_bridgecov = {}
for block in sorted(blocks) :
    sorted_sites = sorted(block_sites[block])
    sorted_reads = sorted(block_reads[block])

    for site in sorted_sites :

        site_physcov[site] = 0
        site_bridgecov[site] = 0
        for read in sorted_reads :
            if read_ends[read][0] <= site and read_ends[read][1] >= site :
                site_physcov[site] += 1

    for i in range(len(sorted_sites) - 1) :

        for read in sorted_reads :
            if sorted_sites[i] in read_sites[read] and sorted_sites[i+1] in read_sites[read] :
                site_bridgecov[sorted_sites[i]] += 1

# make a directory to store info files
dr = filename+'.info_'
if not os.path.exists(dr) :
    os.makedirs(dr)

# general stats
stats = open(dr+'/stats_','w')
stats.write('#general statistics\n')
stats.write('number of blocks: {}\n'.format(len(blocks)))
stats.write('number of reads: {}\n'.format(len(reads)))
stats.write('number of sites: {}\n'.format(len(sites)))

phasing_intensity = 0 # measure of how much time (mem) phasing should take
for site in sites :
    phasing_intensity += 2**site_physcov[site]
stats.write('phasing intensity: {}\n'.format(phasing_intensity))

# reads
byread = open(dr+'/reads_','w')
byread.write('#read\tlength\n')
for read in sorted(reads) :
    byread.write(str(read)+'\t'+str(reads[read])+'\n')

byread.close()

# sites, site alleles, site zygosity
bysite = open(dr+'/sites_','w')
bysitebyallele = open(dr+'/site_alleles_','w')
bysitebyzygosity = open(dr+'/site_zygosity_','w')
bysite.write('#site\tcoverage\tphyscov\tbridgecov\n')
bysitebyallele.write('#site\tmaj\tmin\n')
bysitebyzygosity.write('#site\t0\t1\n')
for site in sorted(sites) :
    bysite.write(str(site)+'\t'+str(sites[site])+'\t'+str(site_physcov[site])+'\t'+str(site_bridgecov[site])+'\n')
    bysitebyallele.write(str(site)+'\t'+(site_alleles[site][0] if site_alleles[site][0] else '_')+'\t'+(site_alleles[site][1] if site_alleles[site][1] else '_')+'\n')
    bysitebyzygosity.write(str(site)+'\t'+str(site_zygosity[site][0])+
                           '\t'+str(site_zygosity[site][1])+'\n')
bysite.close()
bysitebyallele.close()

# blocks
byblock = open(dr+'/blocks_','w')
byblockbyread = open(dr+'/block_reads_','w')
byblockbysite = open(dr+'/block_sites_','w')
byblock.write('#block\tnreads\tnsites\tmaxlen\tmaxcov\tmaxphyscov\n')
byblockbyread.write('#block\treads\n')
byblockbysite.write('#block\tsites\n')
for block in sorted(blocks) :
    r,s = str(len(block_reads[block])),str(len(block_sites[block]))
    maxlen = str(max(reads[read] for read in block_reads[block]))
    maxcov = str(max(sites[site] for site in block_sites[block]))
    maxphys = str(max(site_physcov[site] for site in block_sites[block]))
    byblock.write('\t'.join([str(block),r,s,maxlen,maxcov,maxphys])+'\n')
    byblockbyread.write(str(block)+'\t'+format(sorted(block_reads[block]))+'\n')
    byblockbysite.write(str(block)+'\t'+format(sorted(block_sites[block]))+'\n')

byblock.close()
byblockbyread.close()
byblockbysite.close()

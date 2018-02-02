#
# for running all the different haplotyping softwares
#----------------------------------------------------------------------
#
data_dir = '/data/phasing-comparison-experiments-hapcol'
hap_dir = '/home/prj_rnabwt/haplotyping'

# haplotyping softwares
_corewh_ = 'programs/core_whatshap/build/dp'
_hapchat_ = 'programs/balancing-hapcol/build/hapcol'
_hapcol_ = 'programs/HapCol/build/hapcol'
_hapcut2_ = 'programs/HapCUT2/build/HAPCUT2'
_probhap_ = 'programs/ProbHap/probhap.py'

# limits on memory usage and runtime
memlimit = 64 * 1024 * 1024 # 64GB limit (in KB)
timelimit = '24h' # 24 hour time limit

# lists of methods by name
sih_methods = ['refhap', 'fasthare']
hapcut_methods = ['hapcut', 'hapcut2'] + sih_methods + ['probhap']
hap_methods = ['core_wh', 'hapcol', 'hapchat']
all_methods = hapcut_methods + hap_methods + ['whatshap']

# auxiliary programs
time = '/usr/bin/time'
timeout = '/usr/bin/timeout'
phase = 'programs/whatshap/venv/bin/whatshap phase'
compare = 'programs/whatshap/venv/bin/whatshap compare'
hapcut2vcf = 'programs/conversion/whatshap/venv/bin/whatshap hapcut2vcf'
extract_hairs = 'programs/hapcut/extractHAIRS'
extract_hairs2 = 'programs/HapCUT2/build/extractHAIRS'

# datasets
data = ['ashk', 'sim']
platforms = ['pacbio']
individuals = ['child'] # mother, father, ..
chromosomes = [1] # 21, ..
coverages = list(range(25, 65, 5)) # 25, 30, .., 60
max_covs = list(range(15, 40, 5))
chr_covs = { 1 : coverages }

# preprocessing
modes = ['raw', 'realigned'] # realignment
hs = max_covs  # whatshap read selection
indelmodes = ['indels', 'noindels'] # indel modes for hairs methods

# merging
merge_pattern = 'merged_e{err,[0-9]+}_m{max,[0-9]+}_t{thresh,[0-9]+}_n{neg,[0-9]+}'
merge_regex = 'no_merging|merged_e[0-9]+_m[0-9]+_t[0-9]+_n[0-9]+'
error_rates = [15]
max_errs = [25]
thresholds = [6] # 17, ..
neg_threshs = [3]
mergings = ['merged_e{}_m{}_t{}_n{}'.format(err, max, thresh, neg)
	for err in error_rates
	for max in max_errs
	for thresh in thresholds
	for neg in neg_threshs] + ['no_merging']

# downsampling to a max coverage (in a random greedy way)
downs_pattern = 'downs_s{seed,[0-9]+}_m{maxcov,[0-9]+}'
downs_regex = 'no_downs|downs_s[0-9]+_m[0-9]+'
sample_pattern = 'sample_s{seed,[0-9]+}_m{maxcov,[0-9]+}'
seeds = [1] # 2, 3, .. for (pseudo-) random downsampling
downsamplings = ['downs_s{}_m{}'.format(seed, max)
	for seed in seeds
	for max in max_covs] + ['no_downs']

# epislon / alpha pairs for hapcol, and variants
ea_vals = ['05_1', '05_01', '05_001', '05_0001', '05_00001']
ea_two = ['01_1', '01_01', '01_001', '01_0001', '1_1', '1_01', '1_001', '1_0001']

# common patterns
vcf_pattern = '{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,(all|[0-9]+)}'
hairs_pattern = dataset_pattern + '.{realignment,(raw|realigned)}{indelmode,(.indels|.noindels)}'
whatshap_pattern = dataset_pattern + '.{realignment,(raw|realigned)}.h{h,([0-9]+|N)}'
post_pattern = whatshap_pattern + '.{mergebefore,(' + merge_regex + ')}.{downsample,(' + downs_regex + ')}.{mergeafter,(' + merge_regex + ')}'
full_pattern = post_pattern + '{ea,(|.[0-9]+_[0-9]+)}{balancing,(|.b([0-9]+|N)_[0-9]+)}{indelmode,(|.indels|.noindels)}'

def list_regex(a) :
	return '|'.join([x for x in a])

sih_pattern = '{method,(' + list_regex(sih_methods) + ')}'
hapcut_pattern = '{method,(' + list_regex(hapcut_methods) + ')}'
hap_pattern = '{method,(' + list_regex(hap_methods) + ')}'
methods_pattern = '{method,(' + list_regex(all_methods) + ')}'

output_pattern = methods_pattern + '/' + full_pattern

#
# useful lists and list-defining functions
#----------------------------------------------------------------------

# list of files we need to generate
datasets = ['{}.pacbio.child.chr{}.cov{}'.format(data, chromosome, coverage)
	for data in data
	for chromosome in chromosomes
	for coverage in chr_covs[chromosome]]

simulated = ['{}.pacbio.child.chr{}.cov{}'.format(data, chromosome, coverage)
	for data in ['sim']
	for chromosome in chromosomes
	for coverage in chr_covs[chromosome]]

whatshap_downsample = ['{}.{}.h{}.no_merging.no_downs.no_merging'.format(dataset, mode, h)
	for dataset in datasets
	for mode in modes
	for h in hs]

post_whatshap = ['{}.{}.hN.{}.{}.no_merging'.format(dataset, mode, merging, downsampling)
	for dataset in datasets
	for mode in modes
	for merging in mergings
	for downsampling in downsamplings]

# datasets processed by whatshap to a specified list of max cov.
def whatshap(datasets_, modes_, maxs_) :
	return ['{}.{}.h{}.no_merging.no_downs.no_merging'.format(dataset_, mode_, h_)
		for dataset_ in datasets_
		for mode_ in modes_
		for h_ in maxs_]

# partial paramerization of a merging (threshold and neg. threshold)
def merging(thrs_, negthrs_) :
	return ['merged_e{}_m{}_t{}_n{}'.format(err_, max_, thresh_, neg_)
		for err_ in error_rates
		for max_ in max_errs
		for thresh_ in thrs_
		for neg_ in negthrs_]

# downsampling to a specified list of max coverages
def downs(maxs_) :
	return ['downs_s{}_m{}'.format(seed_, max_)
		for seed_ in seeds
		for max_ in maxs_]

# datasets postprocessed to a specified list of max cov
def postproc(datasets_, modes_, thrs_, negthrs_, maxs_, rnddowns = False) :
	only_rnddowns = ['.no_merging'] if rnddowns else []
	return ['{}.{}.hN.{}.{}.no_merging'.format(dataset_, mode_, merging_, downsampling_)
		for dataset_ in datasets_
		for mode_ in modes_
		for merging_ in merging(thrs_, negthrs_) + only_rnddowns
		for downsampling_ in downs(maxs_)]

# datasets processed by a hairs method (just a shortcut, really)
def hairs(datasets_, modes_, indelmodes_) :
	return ['{}.{}.hN.no_merging.no_downs.no_merging.{}'.format(dataset_, mode_, indelmode_)
		for dataset_ in datasets_
		for mode_ in modes_
		for indelmode_ in indelmodes_]

# datasets both processed by whatshap and postprocessed to list of max cov.
def sliceof(datasets_, modes_, thrs_, negthrs_, maxs_) :
	return whatshap(datasets_, modes_, maxs_) + postproc(datasets_, modes_, thrs_, negthrs_, maxs_, True)

# define a subset of the datasets in terms of chromosomes and coverages
def datasubset(chr_covs_) :
	return ['{}.pacbio.child.chr{}.cov{}'.format(data_, chromosome_, coverage_)
		for data_ in data
		for chromosome_ in chr_covs_
		for coverage_ in chr_covs_[chromosome_]]

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('output/whatshap/{pattern}.diff',
			pattern = whatshap(datasets, ['realigned'],
				[15, 20, 25])),

		expand('output/hapcol/{pattern}.diff',
			pattern = whatshap(datasets, ['raw'],
				[15, 20, 25, 30])),

		expand('output/hapchat/{pattern}.05_01.bN_0.{ext}',
			pattern = postproc(datasets, ['realigned'], [6], [3],
				[15, 20, 25, 30, 35, 40]),
			ext = ['diff', 'inc']),

		expand('output/hapcut2/{pattern}.diff',
			pattern = hairs(simulated, modes, indelmodes)),

		expand('output/{method}/{pattern}.diff',
			method = sih_methods + ['probhap'],
			pattern = hairs(datasets, ['raw'], indelmodes))

#
# master rule
#----------------------------------------------------------------------
rule setup :
	input :
		expand('wif/{pattern}.wif.info_/block_sites_',
			pattern = whatshap_downsample + post_whatshap),

		expand('hairs/{pattern}.raw.{indelmode}.hairs',
			pattern = datasets,
			indelmode = indelmodes),

		expand('hairs2/{pattern}.{mode}.{indelmode}.hairs',
			pattern = datasets,
			mode = modes,
			indelmode = indelmodes),

		expand('vcf/{data}.child.chr{chr}.phased.vcf',
			data = data,
			chr = chromosomes)

#
# run whatshap on a bam / vcf pair
#----------------------------------------------------------------------
rule run_whatshap :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf',
		ref = 'reference/human_g1k_v37.fasta'

	params :
		realignment = lambda wildcards, input :
                        '--reference '+input.ref if wildcards.realignment == 'realigned' else '',
		h = lambda wildcards :
			'1000' if wildcards.h == 'N' else wildcards.h

	output : 'output/whatshap/' + post_pattern + '.phased.vcf'

	log :
		log = 'output/whatshap/' + post_pattern + '.log',
		time = 'output/whatshap/' + post_pattern + '.time'

	message : '''

   running whatshap on bam/vcf pair :

   {input.bam} / {input.vcf}

   selecting coverage down to {wildcards.h}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      whatshap phase -o {output} {params.realignment} -H {params.h} \
         {input.vcf} {input.bam} > {log.log} 2>&1 || true
   touch {output} '''

#
# run the core whatshap dp on a wif file
#----------------------------------------------------------------------
rule run_core_whatshap :
	input : 'wif/' + post_pattern + '.wif'
	output : 'output/core_wh/' + post_pattern + '.hap'

	log :
		log = 'output/core_wh/' + post_pattern + '.log',
		time = 'output/core_wh/' + post_pattern + '.time',
		wif = 'output/core_wh/' + post_pattern + '.wif'

	message : '''

   running core whatshap dp on :

   {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_corewh_} -h {output} -a {input} \
         > {log.wif} 2> {log.log} || true
   touch {output} '''

#
# run hapchat with increase-k and balancing on a wif file
#----------------------------------------------------------------------
rule run_hapchat :
	input : 'wif/' + post_pattern + '.wif'

	params :
		epsilon = lambda wildcards :
			'0' + wildcards.ea.split('_')[0],
		alpha = lambda wildcards :
			'0.' + wildcards.ea.split('_')[1],
		balance_cov = lambda wildcards :
			'1000' if wildcards.balancing.split('_')[0].split('b')[1] == 'N' else wildcards.balancing.split('_')[0].split('b')[1],
		balance_ratio = lambda wildcards :
			'0.' + wildcards.balancing.split('_')[1]

	output :
		hap = 'output/hapchat/' + full_pattern + '.hap',
		log = 'output/hapchat/' + full_pattern + '.log',

	log :
		time = 'output/hapchat/' + full_pattern + '.time'

	message : '''

   running hapchat on :

   {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_hapchat_} -i {input} -o {output.hap} -A \
         -e {params.epsilon} -a {params.alpha} \
         -b {params.balance_cov} -r {params.balance_ratio} \
            > {output.log} 2>&1 || true
   touch {output} '''

#
# run hapcol on a wif file (from within a script that adapts alpha)
#----------------------------------------------------------------------
rule run_hapcol :
	input :
		script = 'scripts/run.hapcol.bash',
		wif = 'wif/' + post_pattern + '.wif'

	output : 'output/hapcol/' + post_pattern + '.hap'

	log :
		log = 'output/hapcol/' + post_pattern + '.log',
		time = 'output/hapcol/' + post_pattern + '.time'

	message : '''

   running hapcol on :

   {input.wif}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      bash {input.script} {_hapcol_} {input.wif} {output} \
         > {log.log} 2>&1 || true
   touch {output} '''

#
# run hapcut2 on a hairs file
#----------------------------------------------------------------------
rule run_hapcut2 :
	input :
		hairs = 'hairs2/' + hairs_pattern + '.hairs',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output :
		'output/hapcut2/' + full_pattern + '.out'

	log :
		log = 'output/hapcut2/' + full_pattern + '.log',
		time = 'output/hapcut2/' + full_pattern + '.time'

	message : '''

   running hapcut2 on hairs file: {input.hairs}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {_hapcut2_} --fragments {input.hairs} --VCF {input.vcf} \
         --output {output} > {log.log} 2>&1 || true
   touch {output} '''

#
# run probhap on a hairs file
#----------------------------------------------------------------------
rule run_probhap :
	input :
		'hairs/' + hairs_pattern + '.hairs'

	output :
		'output/probhap/' + full_pattern + '.out'

	log :
		reads = 'output/probhap/' + full_pattern + '.reads',
		assignments = 'output/probhap/' + full_pattern + '.assignments',
		log = 'output/probhap/' + full_pattern + '.log',
		time = 'output/probhap/' + full_pattern + '.time'

	message : '''

   running probhap on hairs file: {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      python2 {_probhap_} \
         --reads {input} \
         --parsed-reads {log.reads} \
         --phase {output} \
         --assignments {log.assignments} > {log.log} 2>&1 || true
   touch {output} '''

#
# run a SIH method (refhap, fasthare, ..)
#----------------------------------------------------------------------
def sih_method(wildcards) :

	sih = 'programs/refhap/SingleIndividualHaplotyper/SIH.jar'
	main = 'mpg.molgen.sih.main.SIH'
	opt = '-a FastHare' if wildcards.method == 'fasthare' else ''

        return 'java -cp {} {} {}'.format(sih, main, opt)

rule run_sih_method :
	input :
		'hairs/' + hairs_pattern + '.hairs',

	params :
		method = sih_method

	output :
		'output/' + sih_pattern + '/' + full_pattern + '.out'

	log :
		log = 'output/{method}/' + full_pattern + '.log',
		time = 'output/{method}/' + full_pattern + '.time'

	message : '''

   running {wildcards.method} on hairs file: {input}

   with memory limit: {memlimit}K

   with time limit: {timelimit} '''

	shell : '''

   ulimit -Sv {memlimit}
   {time} -v -o {log.time} {timeout} {timelimit} \
      {params.method} {input} {output} > {log.log} 2>&1 || true
   touch {output} '''

#
# compare phased vcfs to true phasing
#----------------------------------------------------------------------
rule compare_vcfs :
	input :
		true = 'vcf/' + vcf_pattern + '.phased.vcf',
		vcf = 'output/' + output_pattern + '.phased.vcf'

	output :
		diff = 'output/' + output_pattern + '.diff',
		bed = 'output/' + output_pattern + '.bed'

	log : 'output/' + output_pattern + '.diff.log'

	message : '''

   comparing inferred phasing:

   {input.vcf}

   to true phasing: {input.true} '''

	shell : '''

   {compare} --switch-error-bed {output.bed} \
      {input.true} {input.vcf} \
         > {output.diff} 2> {log} || true
   touch {output} '''

# convert old-skool hap format to a phased vcf
rule phase_vcf :
	input :
		script = 'scripts/subvcf.py',
		hap = 'output/' + hap_pattern + '/' + full_pattern + '.hap',
		blocks = 'wif/' + post_pattern + '.wif.info_/block_sites_',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'output/' + hap_pattern + '/' + full_pattern + '.phased.vcf'

	log : 'output/' + hap_pattern + '/' + full_pattern + '.phased.vcf.log'

	message : '''

   adding phase information from {input.hap}
   to {input.vcf},

   obtaining: {output} '''

	shell : '''

   python {input.script} -p {input.hap} {input.blocks} {input.vcf} \
      > {output} 2> {log} '''

# convert hapcut (or hapcut-like) output format to (phased) vcf
rule hapcut_to_vcf :
	input :
		out = 'output/' + hapcut_pattern + '/' + full_pattern + '.out',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf'

	output : 'output/' + hapcut_pattern + '/' + full_pattern + '.phased.vcf'

	log : 'output/' + hapcut_pattern + '/' + full_pattern + '.phased.vcf.log'

	message : 'converting hapcut {input} to {output}'

	shell : '''

   {hapcut2vcf} {input.vcf} {input.out} > {output} 2> {log} '''

#
# get details on increasing k from a hapchat log file
#----------------------------------------------------------------------
rule increments :
	input :
		script = 'scripts/increments.py',
		log = 'output/hapchat/' + full_pattern + '.log'

	output : 'output/hapchat/' + full_pattern + '.inc'

	log : 'output/hapchat/' + full_pattern + '.inc.log'

	message : 'obtain details on increasing k from {input.log}'

	shell : '''

   printf "%s Cov. %s: " {wildcards.dataset} {wildcards.coverage} > {output}
   python {input.script} -r {input.log} >> {output} 2> {log} '''

#
# link to files from phasing comparison experiments directory
#----------------------------------------------------------------------
rule link_vcf :
	input : data_dir + '/vcf/' + vcf_pattern + '.{phase,(phased|unphased)}.vcf'
	output : 'vcf/' + vcf_pattern + '.{phase,(phased|unphased)}.vcf'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_bam_bai :
	input : data_dir + '/bam/' + dataset_pattern + '.{ext,(bam|bai)}'
	output : 'bam/' + dataset_pattern + '.{ext,(bam|bai)}'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

rule link_reference :
	input : data_dir + '/reference/human_g1k_v37.fasta'
	output : 'reference/human_g1k_v37.fasta'
	message : 'linking {input} to {output}'
	shell : 'ln -fsrv {input} {output}'

#
# obtain a wif file from a bam / vcf pair using whatshap
#----------------------------------------------------------------------
rule get_wif :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		bai = 'bam/' + dataset_pattern + '.bai',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf',
		ref = 'reference/human_g1k_v37.fasta'

	params :
                realignment = lambda wildcards, input :
			'--reference '+input.ref if wildcards.realignment == 'realigned' else '',
		h = lambda wildcards :
			'1000' if wildcards.h == 'N' else wildcards.h

	output : 'wif/' + whatshap_pattern + '.wif'

	log :
		transcript = 'wif/' + whatshap_pattern + '.wif.transcript',
		log = 'wif/' + whatshap_pattern + '.wif.log',
		time = 'wif/' + whatshap_pattern + '.wif.time'

	message : '''

   obtaining wif file: {output}

   from: {input.bam} / {input.vcf} pair,

   after selecting reads down to max coverage {params.h} '''

	shell : '''

   {time} -v -o {log.time} \
      {phase} -o /dev/null {params.realignment} \
         --output-wif {output} -H {params.h} \
         {input.vcf} {input.bam} > {log.transcript} 2> {log.log} '''

#
# obtain a hairs file from a bam / vcf pair using hapcut extract hairs
#----------------------------------------------------------------------
rule get_hairs :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		bai = 'bam/' + dataset_pattern + '.bai',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf',
		ref = 'reference/human_g1k_v37.fasta'

	params : lambda wildcards :
		'--indels 1' if wildcards.indelmode == '.indels' else ''

	output : 'hairs/' + hairs_pattern + '.hairs'

	log :
		log = 'hairs/' + hairs_pattern + '.log',
		time = 'hairs/' + hairs_pattern + '.time'

	message : '''

   obtaining hairs file: {output}

   from: {input.bam} / {input.vcf} pair '''

	shell : '''

   {time} -v -o {log.time} \
      {extract_hairs} --ref {input.ref} {params} \
         --VCF {input.vcf} --bam {input.bam} \
            > {output} 2> {log.log} '''

# obtain hairs using hapcut2 extract hairs
rule get_hairs2 :
	input :
		bam = 'bam/' + dataset_pattern + '.bam',
		bai = 'bam/' + dataset_pattern + '.bai',
		vcf = 'vcf/' + vcf_pattern + '.unphased.vcf',
		ref = 'reference/human_g1k_v37.fasta'

	params :
		realignment = lambda wildcards :
			'--pacbio 1' if wildcards.realignment == 'realigned' else '',
		indels = lambda wildcards :
			'--indels 1' if wildcards.indelmode == '.indels' else ''

	output : 'hairs2/' + hairs_pattern + '.hairs'

	log :
		log = 'hairs2/' + hairs_pattern + '.log',
		time = 'hairs2/' + hairs_pattern + '.time'

	message : '''

   obtaining hairs file: {output} with HapCUT2 extractHAIRS

   from: {input.bam} / {input.vcf} pair '''

	shell : '''

   {time} -v -o {log.time} \
      {extract_hairs2} --ref {input.ref} {params} \
         --VCF {input.vcf} --bam {input.bam} --out {output} \
            > {log.log} 2>&1 '''

#
# downsample a wif file to a specified max coverage
#----------------------------------------------------------------------
rule extract_sample :
	input :
		source = '{path}.wif',
		sample = '{path}.wif.' + sample_pattern

	output : '{path}.' + downs_pattern + '.wif'
	message : 'extract lines {input.sample} from {input.source}'

	shell : '''

   awk '{{printf "%.20d %s\\n", NR, $0}}' {input.source} | join - \
      <(awk '{{printf "%.20d\\n", $1}}' {input.sample} | sort) | \
         sed 's/^[0-9]* //' > {output} '''

# dummy rule to ensure the naming is consistent
rule no_downsampling :
	input : '{path}.wif'
	output : '{path}.no_downs.wif'
	message : 'perform no downsampling'
	shell : 'cp {input} {output}'

# greedily downsample wif to a coverage according to a shuffle
rule downsample :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif',
		shuf = '{path}.wif.lines.shuf{seed}'

	output : '{path}.wif.' + sample_pattern

	log :
		log = '{path}.wif.' + sample_pattern + '.log',
		time = '{path}.wif.' + sample_pattern + '.time'

	message : '''

   psuedorandom downsampling of: {input.wif}

   to maximum coverage {wildcards.maxcov}, according to:

   {input.shuf}, producing the sample:

   {output} '''

	shell : '''

   {time} -v -o {log.time} \
      python {input.script} -s {wildcards.maxcov} {input.shuf} {input.wif} \
         > {output} 2> {log.log} '''

# seeded pseudorandom shuffle of lines of a file (cf. gnu.org)
rule permute_lines :
	input : '{path}.lines'
	output : '{path}.lines.shuf{seed,[0-9]+}'
	message : 'pseudorandom shuffle of {input} with seed {wildcards.seed}'
	shell : '''

   shuf {input} --random-source=<(openssl enc -aes-256-ctr \
      -pass pass:"{wildcards.seed}" -nosalt </dev/zero 2>/dev/null) > {output} '''

# get lines (numbers) from a file
rule get_lines :
	input : '{path}'
	output : '{path}.lines'
	message : 'obtain lines (numbers) from {input}'
	shell : ''' awk '{{print NR}}' {input} > {output} '''

#
# obtain a (red-blue-) merged wif from a wif
#----------------------------------------------------------------------
rule merge_wif :
	input :
		script = 'scripts/rb-merge.py',
		wif = '{path}.wif'

	params :
		e = lambda wildcards : '0.' + wildcards.err,
		m = lambda wildcards : '0.' + wildcards.max,
		t = lambda wildcards : 10 ** int(wildcards.thresh),
		n = lambda wildcards : 10 ** int(wildcards.neg)

	output : '{path}.' + merge_pattern + '.wif'

	log :
		log = '{path}.' + merge_pattern + '.wif.log',
		time = '{path}.' + merge_pattern + '.wif.time',
		graph = '{path}.' + merge_pattern + '.wif.graph'

	message : '''

   merge reads of {input.wif},
   with parameters:

   error rate = {params.e},
   max error rate = {params.m},
   threshold = {params.t},
   neg. threshold = {params.n}

   producing: {output} '''

	shell : '''

   {time} -v -o {log.time} \
      python {input.script} -v -e {params.e} -m {params.m} \
         -t {params.t} -n {params.n} -w {input.wif} -o {output} \
         -g {log.graph} > {log.log} 2>&1 '''

# dummy rule to ensure the naming is consistent
rule no_merging :
	input : '{path}.wif'
	output : '{path}.no_merging.wif'
	message : 'perform no merging'
	shell : 'cp {input} {output}'

#
# obtain properties of a wif file
#----------------------------------------------------------------------
rule wif_info :
	input :
		script = 'scripts/wiftools.py',
		wif = '{path}.wif'

	output :
		expand('{{path}}.wif.info_/{file}',
			file = ['block_reads_', 'block_sites_',
				'site_alleles_', 'site_zygosity_',
				'blocks_', 'reads_',
				'sites_', 'stats_'])

	message : 'obtaining info for {input.wif}'
	shell : 'python {input.script} -i {input.wif}'

'''

   How to re-run this workflow
   -------------------------------------------------------------------------------

   Install dependencies:
    - see README
    - make 'shapeit' binary available

'''

import pysam
import textwrap

picard_tmpdir_switch = ''
if 'TMPDIR' in os.environ:
	picard_tmpdir_switch = 'TMP_DIR=%s' % os.environ['TMPDIR']

datasets = ['ashk', 'sim']
individuals = ['child']
chromosomes = [1]
reference = 'reference/human_g1k_v37.fasta'
role_to_sampleid = {
	'mother' : 'HG004',
	'father' : 'HG003',
	'child' : 'HG002' }

time = "/usr/bin/time -f '%M kB; real %e; user %U; sys %S'"

# Paths to scripts distributed along with the Snakefile
vcf_merge_trio = srcdir('scripts/vcf_merge_trio.pl')
genomesimulator = srcdir('scripts/genomesimulator.py')

# Tools assumed to be installed somewhere on the PATH.
gzip = 'pigz'
whatshap = 'whatshap'

# Software that must be installed manually prior to running
# the Snakefile due to licensing restrictions
shapeit = 'restricted-software/shapeit'

# common patterns
pattern = '{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'
platform_pattern = '{platform,[a-z]+}.' + pattern
dataset_pattern = '{dataset,[a-z]+}.' + platform_pattern
coverage_pattern = platform_pattern + '.cov{coverage,(all|[0-9]+)}'
full_pattern = '{dataset,[a-z]+}.' + coverage_pattern

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		expand('stats/bam/{dataset}.pacbio.child.chr{chromosome}.covall.coverage',
			dataset = datasets,
			chromosome = chromosomes),

		expand('vcf/{dataset}.child.chr{chromosome}.unphased.vcf',
			dataset = datasets,
			chromosome = chromosomes),

		expand('stats/bam/{dataset}.pacbio.child.chr{chromosome}.covall.avgrlen',
			dataset = datasets,
			chromosome = chromosomes)

#
# obtain external dependencies: software and data files
#----------------------------------------------------------------------
rule shapeit_missing :
	output : shapeit

	shell : '''

   echo
   echo "Due to licensing restrictions, you need to manually download "
   echo "shapeit and make it available to this pipeline."
   echo
   echo "Go to https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download ,"
   echo "download version v2.r837 (shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz),"
   echo "unpack the tgz file and place the bin/shapeit binary into '{shapeit}'."
   echo "Then re-run snakemake."
   echo
   exit 1 '''


rule download_ashkenazim :
	threads : 100
	output :
		protected('download/AshkenazimTrio/{file}.{ext,(bam|bam.bai|vcf.gz)}')

	shell : '''

   wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/{wildcards.file}.{wildcards.ext} '''


rule download_reference :
	output : 'reference/human_g1k_v37.fasta'

	shell : '''

   wget -O {output}.gz.incomplete ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
   mv {output}.gz.incomplete {output}.gz
   # gunzip fails with "decompression OK, trailing garbage ignored" because
   # the file is razf-compressed (gzip-compatible, but with an index at end)
   gunzip -f {output}.gz || true
   cd reference && md5sum -c MD5SUM '''


rule index_reference :
	output : '{path}.fasta.bwt'
	input : '{path}.fasta'
	shell: 'bwa index {input}'

#
# downsample BAM files, add read group information, etc.
#----------------------------------------------------------------------
def ashkenazim_trio_bam(wildcards) :

	individual = {
		'mother': dict(name='mother', id='4', na='24143'),
		'father': dict(name='father', id='3', na='24149'),
		'child': dict(name='son', id='2', na='24385'),
	}[wildcards.individual]
	individual['chromosome'] = wildcards.chromosome
	individual['ext'] = wildcards.ext

	return 'download/AshkenazimTrio/HG00{id}_NA{na}_{name}/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg00{id}_gr37_{chromosome}.{ext}'.format(**individual)


rule create_ashk_pacbio_links :
	input : ashkenazim_trio_bam

	output : 'bam/incorrect-readgroup/ashk.pacbio.{individual}.chr{chromosome,[0-9]+}.covall.{ext,(bam|bam.bai)}'

	shell: 'ln -fsrv {input} {output}'


rule calculate_coverage:
	input: 'bam/' + full_pattern + '.bam'
	output: 'stats/bam/' + full_pattern + '.coverage'
	message: 'Computing coverage for {input}'

	run: 
		bam = pysam.Samfile(input[0])
		length = None
		for e in bam.header.get('SQ'):
			if e['SN'] == wildcards.chromosome:
				length = e['LN']
		assert length != None
		shell('''

   samtools depth {input} | \
      awk '{{sum+=$3}} END {{ print sum/{length} }}' > {output} ''')


rule calculate_avg_read_length :
	input : 'bam/' + full_pattern + '.bam'
	output : 'stats/bam/' + full_pattern + '.avgrlen'
	message : 'computing average mapped read length from {input}'
	shell : '''

   samtools view {input} | \
      awk '{{ sum += length($10) ; n += 1 }} END {{ print sum/n }}' \
         > {output} '''


def downsampling_source(wildcards) :

	chr = wildcards.chromosome
	dest_cov = int(wildcards.coverage)

	source_cov = 'all'
	if chr == '1' :
		if 2 < dest_cov < 60 :
			source_cov = dest_cov + 5
	elif chr == '21' :
		if 2 < dest_cov < 50 :
			source_cov = dest_cov + 5

	filebody = '{}.pacbio.{}.chr{}.cov{}'.format(
		wildcards.dataset, wildcards.individual, chr, source_cov)

	return { 'bam'      : 'bam/{}.bam'.format(filebody),
		 'coverage' : 'stats/bam/{}.coverage'.format(filebody) }


rule downsample :
	input: unpack(downsampling_source)

	output : 
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.bam',
		bai='bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.bai'

	log : 'bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.log'

	message : 'Downsampling {input.bam} to {wildcards.coverage}x'

	run :
		input_coverage = float(open(input.coverage).readline())
		p = float(wildcards.coverage) / input_coverage
		seed = hash(output)
		shell('''

   picard DownsampleSam INPUT={input.bam} RANDOM_SEED={seed} \
      CREATE_INDEX=true OUTPUT={output.bam} PROBABILITY={p} \
      VALIDATION_STRINGENCY=SILENT > {log} 2>&1 ''')


rule add_read_groups :
	'''Add RG, prepend sample names to read names, strip unused tags'''

	input : 'bam/incorrect-readgroup/ashk.' + coverage_pattern + '.bam'

	output : 
		bam = 'bam/ashk.' + coverage_pattern + '.bam',
		bai = 'bam/ashk.' + coverage_pattern + '.bai'

	log : 'bam/ashk.' + coverage_pattern + '.readgroup.log'

	run :
		sample = role_to_sampleid[wildcards.individual]
		shell(textwrap.dedent('''

   (samtools view -h {input} | \\
   sed '/^[^@]/ s|^|{wildcards.individual}_|g' | \\
   sed -re 's_\\t(dq|dt|ip|iq|st|sq|mq):Z:[^\\t]*__g' | \\
   samtools view -u - | \\
   picard AddOrReplaceReadGroups {picard_tmpdir_switch} \
      CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT \
      I=/dev/stdin O={output.bam} ID={wildcards.individual} \
      LB=UNKNOWN PL=PACBIO PU=UNKNOWN SM={sample}) > {log} 2>&1 '''))


rule bam_to_fastq :
	input : 'bam/{file}.bam'
	output : 'fastq/{file}.fastq'
	shell : 'bedtools bamtofastq -i {input} -fq {output}'

#
# prepare the downloaded VCF files: extract a chromosome, filter, etc.
#----------------------------------------------------------------------
rule unzip_vcf :
	input : 'download/AshkenazimTrio/analysis/NIST_CallsIn2Technologies_05182015/HG{id}-multiall-fullcombine.vcf.gz'

	output : 'unphased-ashk/HG{id,[0-9]+}-multiall-fullcombine.vcf'
	shell : 'zcat {input} > {output}'


rule filter_vcfs :
	input :
		mother = 'unphased-ashk/HG004-multiall-fullcombine.vcf',
		father = 'unphased-ashk/HG003-multiall-fullcombine.vcf',
		child = 'unphased-ashk/HG002-multiall-fullcombine.vcf'

	output : 'vcf/ashk.trio.unphased.vcf'
	log : 'vcf/ashk.trio.unphased.vcf.log'
	message : 'Filtering and merging input VCFs to create {output}'
	version : 2

	shell : r"{vcf_merge_trio} {input.mother} {input.father} {input.child} | sed '/^[^#]/ s| |\t|g; s_0|1_0/1_g; s_1|0_0/1_g; s_1/0_0/1_g; s_1|1_1/1_g; s_0|0_0/0_g' > {output} 2> {log}"


rule split_vcf :
	input : 'vcf/ashk.trio.unphased.vcf'
	output : 'vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf'
	message : 'Extracting chromosome {wildcards.chromosome} from {input}'
	shell : '''awk '/^#/ || ($1 == "{wildcards.chromosome}")' {input} > {output}'''


rule unphase :
	input : '{base}.phased.vcf'
	output : '{base}.unphased.vcf'
	shell : '{whatshap} unphase {input} > {output}'

#
# compute "ground truth" phasing with SHAPEIT
#----------------------------------------------------------------------
rule download_1000GP :
	output : 'download/1000GP_Phase3/{filename}'
	shell : '''

   wget -O {output} http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/{wildcards.filename} '''


rule shapeit_check :
	input :
		legend = 'download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.legend.gz',
		genmap = 'download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt',
		refhaps = 'download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.hap.gz',
		refsamples = 'download/1000GP_Phase3/1000GP_Phase3.sample',
		vcf = 'vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf',
		shapeit = shapeit,

	output :
		'shapeit/trio.chr{chromosome,[0-9]+}.snp.strand',
		'shapeit/trio.chr{chromosome,[0-9]+}.snp.strand.exclude'

	log : 'shapeit/trio.chr{chromosome,[0-9]+}.check.log'

	shell : '''

   ({shapeit} -check -V {input.vcf} -M {input.genmap} \
      --input-ref {input.refhaps} {input.legend} {input.refsamples} \
      --output-log shapeit/trio.chr{wildcards.chromosome} || true) > {log} 2>&1 '''


rule shapeit :
	input :
		legend = 'download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.legend.gz',
		genmap = 'download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt',
		refhaps = 'download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.hap.gz',
		refsamples = 'download/1000GP_Phase3/1000GP_Phase3.sample',
		vcf = 'vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf',
		exclude = 'shapeit/trio.chr{chromosome,[0-9]+}.snp.strand.exclude',
		shapeit = shapeit,

	output :
		'shapeit/trio.chr{chromosome,[0-9]+}.haps',
		'shapeit/trio.chr{chromosome,[0-9]+}.sample' 

	log : 'shapeit/trio.chr{chromosome,[0-9]+}.run.log'
	message : 'Running SHAPEIT on chromosome {wildcards.chromosome}'

	shell : '''

   {shapeit} -V {input.vcf} --exclude-snp {input.exclude} -M {input.genmap} \
      --input-ref {input.refhaps} {input.legend} {input.refsamples} \
      -O shapeit/trio.chr{wildcards.chromosome} > {log} 2>&1 '''


rule shapeit_convert_to_vcf :
	input :
		haps = 'shapeit/trio.chr{chromosome,[0-9]+}.haps',
		shapeit = shapeit,

	output :
		vcf = 'shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf'

	log : 'shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf.log'
	message : '''

   Converting SHAPEIT output for chromosome {wildcards.chromosome} to VCF '''

	shell : '''

   {shapeit} -convert --input-haps shapeit/trio.chr{wildcards.chromosome} \
      --output-vcf {output.vcf} > {log} 2>&1 '''


rule split_shapeit_vcf :
	input : 'shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf'
	output : 'vcf/ashk.' + pattern + '.phased.vcf'
	log: 'vcf/ashk.' + pattern + '.phased.vcf.log'
	run :
		sample = role_to_sampleid[wildcards.individual]
		shell('vcf-subset -c {sample} {input} > {output} 2> {log}')

#
# pipeline for simulating reads off of "ground truth" phasing
#----------------------------------------------------------------------
rule copy_sim :
	input : 'vcf/ashk.' + pattern + '.phased.vcf'
	output : 'vcf/sim.' + pattern + '.phased.vcf'
	shell : 'ln -fsrv {input} {output}'


rule sim_fastas :
	input :
		ref = reference,
		vcf = 'vcf/sim.' + pattern + '.phased.vcf',

	output :
		hap1 = 'sim/sim.' + pattern + '.haplotype1.true.fasta',
		hap2 = 'sim/sim.' + pattern + '.haplotype2.true.fasta',

	log : 'sim/sim.' + pattern + '.true.log'
	message: 'Creating true haplotypes {output}'

	run :
		sample = role_to_sampleid[wildcards.individual]
		shell('mkdir -p sim/tmp')
		shell('''

   {genomesimulator} -c {wildcards.chromosome} {input.vcf} \
      {input.ref} sim/tmp > {log} 2>&1 ''')
		shell('''

   mv sim/tmp/{sample}.chr{wildcards.chromosome}.1.fasta {output.hap1} ''')
		shell('''

   mv sim/tmp/{sample}.chr{wildcards.chromosome}.2.fasta {output.hap2}''')


rule simulate_pacbio_reads :
	input :
		sample = 'fastq/ashk.pacbio.' + pattern + '.cov2.fastq',
		haplotype = 'sim/sim.' + pattern + '.haplotype{hap}.true.fasta'

	output :
		fastq = 'sim/sim.' + pattern + '.haplotype{hap,[12]}.fastq.gz',
		maf = 'sim/sim.' + pattern + '.haplotype{hap,[12]}.maf.gz'

	log : 'sim/sim.' + pattern + '.haplotype{hap,[12]}.fastq.log'

	run :
		coverage = 60
		seed = abs(hash(output.fastq))
		shell('mkdir -p sim/tmp')
		shell('''

   time (pbsim --seed {seed} \
      --prefix sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap} \
      --depth {coverage} --sample-fastq {input.sample} \
      {input.haplotype}) > {log} 2>&1 ''')
		shell('''

   awk \'NR%4==1 {{printf("%s_HAP{wildcards.hap}\\n",$0)}} \
      NR%4!=1 {{print}}\' \
      sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_0001.fastq | \
      {gzip} > {output.fastq} ''')
		shell('''

   cat sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_0001.maf | \
      {gzip} > {output.maf} ''')
		shell('''

   rm -f sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_* ''')


rule bwa_mem_single_end_pacbio :
	input :
		fastq = 'sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.fastq.gz',
		ref = reference,
		indexed = reference + '.bwt'

	output : 'sim/sim.' + pattern + '.haplotype{hap,[12]}.bam'
	log : 'sim/sim.' + pattern + '.haplotype{hap,[12]}.bam.log'
	threads : 16

	run:
		sample = role_to_sampleid[wildcards.individual]
		shell('''

   {time} bwa mem -x pacbio -t {threads} {input.ref} {input.fastq} | \
   samtools view -u - | picard AddOrReplaceReadGroups \
   {picard_tmpdir_switch} VALIDATION_STRINGENCY=LENIENT \
   I=/dev/stdin O={output} ID={wildcards.individual}{wildcards.hap} \
   LB=UNKNOWN PL=PACBIO PU=UNKNOWN SM={sample} 2> {log} ''')


rule merge_hap_bams :
	input :
		bam1 = 'sim/sim.' + pattern + '.haplotype1.bam',
		bam2 = 'sim/sim.' + pattern + '.haplotype2.bam'

	output :
		bam = 'bam/sim.pacbio.' + pattern + '.covall.bam',
		bai = 'bam/sim.pacbio.' + pattern + '.covall.bai'

	log : 'bam/sim.pacbio.' + pattern + '.covall.bam.log'
	priority: 5
	message: 'Merging haplotype-specific BAM files to create {output.bam}'

	shell: '''

   {time} picard -Xmx8g MergeSamFiles {picard_tmpdir_switch} \
      VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000 \
      SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true \
      I={input.bam1} I={input.bam2} O={output.bam} >& {log} '''

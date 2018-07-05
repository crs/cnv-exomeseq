### align exome seq

THREADS=64
JAVA='~/jdk1.8.0_92/bin/java '
CNVNATOR='bin/CNVnator-master/cnvnator '
BIN_SIZE=100
rule all:
	#input: 'bam/FTD_P1_E02-35794035/FTD-P1-E02_S15_L003_L004.bam'
	#input: 'rawcnv/FTD_P1_E02-35794035/FTD-P1-E02_S15_L003_L004.rawcnv'
	#input: 'root/FTD_P1_E02-35794035/FTD-P1-E02_S15_L003_L004.root'
	#input: 'logs/FTD_P1_E02-35794035/statistics.FTD-P1-E02_S15_L003_L004.log', 'logs/FTD_P1_E02-35794035/histogram.FTD-P1-E02_S15_L003_L004.log','logs/FTD_P1_E02-35794035/partition.FTD-P1-E02_S15_L003_L004.log'
	#input: 'cnvs/FTD_P1_E02-35794035/FTD-P1-E02_S15_L003_L004.cnv'
	input: 'cnvs/FTD_P1_F04-35810021/FTD-P1-F04_S16_L003_L004.cnv'
	#input: 'aligned/FTD_P1_E02-35794035/FTD-P1-E02_S15_L003_R1_001.sai'

rule cnv_calling:
	input: log='logs/{sample}/partition.{file}_L{laneA}_L{laneB}.log', 		
	output: 'cnvs/{sample}/{file}_L{laneA}_L{laneB}.cnv'
	params: bin_size=BIN_SIZE, root='root/{sample}/{file}_L{laneA}_L{laneB}.root'
	shell: CNVNATOR + "-root {params.root} -call {params.bin_size} > {output}"
	
rule rd_signal_partitioning:
	input:  log='logs/{sample}/statistics.{file}_L{laneA}_L{laneB}.log'
	output: log='logs/{sample}/partition.{file}_L{laneA}_L{laneB}.log'
	params: bin_size=BIN_SIZE, root='root/{sample}/{file}_L{laneA}_L{laneB}.root',
	shell: CNVNATOR + "-root {params.root} -partition {params.bin_size} > {output.log}"

rule statistics:
	input: log='logs/{sample}/histogram.{file}_L{laneA}_L{laneB}.log'	
	output: log='logs/{sample}/statistics.{file}_L{laneA}_L{laneB}.log'
	params: bin_size=BIN_SIZE, root='root/{sample}/{file}_L{laneA}_L{laneB}.root'
	shell: CNVNATOR + "-root {params.root} -genome hg18 -stat {params.bin_size} > {output.log}"

rule generate_histogram:
	input: log='logs/{sample}/extraction.{file}_L{laneA}_L{laneB}.log', reference="indexed/chr"
	output: log='logs/{sample}/histogram.{file}_L{laneA}_L{laneB}.log'
	params: bin_size=BIN_SIZE, root='root/{sample}/{file}_L{laneA}_L{laneB}.root'
	shell: CNVNATOR + "-root {params.root} -his {params.bin_size} -d {input.reference} > {output.log}"
	
rule extractReadMapping:
	input: bam='bam/{sample}/{file}_L{laneA}_L{laneB}.bam'
	output: root='root/{sample}/{file}_L{laneA}_L{laneB}.root', log='logs/{sample}/extraction.{file}_L{laneA}_L{laneB}.log'
	shell: CNVNATOR + "-root {output.root} -tree {input.bam} -unique > {output.log}"
	
rule runPennCNV:
	input: bam='bam/{sample}/{file}_L{laneA}_L{laneB}.bam', reference='indexed/hg18.fa'
	output: 'rawcnv/{sample}/{file}_L{laneA}_L{laneB}.rawcnv'
	shell: "./bin/PennCNV-Seq-master/penncnv-seq-hg18.sh ~/CNV/PennCNV-1.0.4 ~/bin/PennCNV-Seq-master/reference hg18 EUR {input.reference} {input.bam}; touch {output}"

rule convert_and_sort_sam_to_bam:
	input: sam='sam/{sample}/{file}_L{laneA}_L{laneB}.sam'
	output: bam='bam/{sample}/{file}_L{laneA}_L{laneB}.bam', bai='bam/{sample}/{file}_L{laneA}_L{laneB}.bai'
	params: bam='bam/{sample}/{file}_L{laneA}_L{laneB}'
	shell: 'samtools view -bS {input.sam} | samtools sort - {params.bam}'
	
#rule convert_sam_to_bam:
#	input: sam='sam/{sample}/{file}_L{laneA}_L{laneB}.sam'
#	output: bam='bam/{sample}/{file}_L{laneA}_L{laneB}.bam'
#	shell: JAVA + """-Xmx4g -Djava.io.tmpdir=/tmp \
#					-jar bin/picard.jar SortSam \
#					SO=coordinate \
#					INPUT={input.sam} \
#					OUTPUT={output.bam} \
#					VALIDATION_STRINGENCY=LENIENT \
#					CREATE_INDEX=true"""
	
rule pairedEnd:
	input: fwd_sai='aligned/{sample}/{file}_L{laneA}_L{laneB}_R1.sai', rev_sai='aligned/{sample}/{file}_L{laneA}_L{laneB}_R2.sai', fwd_fastq='merged/{sample}/{file}_L{laneA}_L{laneB}_R1.fastq', rev_fastq='merged/{sample}/{file}_L{laneA}_L{laneB}_R2.fastq'
	output: sam='sam/{sample}/{file}_L{laneA}_L{laneB}.sam'
	shell: """bwa sampe -f {output.sam} -r '@RG\tID:{wildcards.sample}_{wildcards.file}_L{wildcards.laneA}_L{wildcards.laneB}\tLB:{wildcards.sample}_{wildcards.file}_L{wildcards.laneA}_L{wildcards.laneB}\tSM:{wildcards.sample}_{wildcards.file}_L{wildcards.laneA}_L{wildcards.laneB}\tPL:ILLUMINA' indexed/hg18 {input.fwd_sai} {input.rev_sai} {input.fwd_fastq} {input.rev_fastq}"""
	
rule alignSamples:
	"""align fastq samples to indexed reference"""
	input: fastq='merged/{sample}/{file}_L{laneA}_L{laneB}_R{R}.fastq', index='indexed/hg18.ann'
	params: index='indexed/hg18'
	output: sai='aligned/{sample}/{file}_L{laneA}_L{laneB}_R{R}.sai'
	threads: THREADS
	shell: 'bwa aln -t {threads} -f {output.sai} {params.index} {input.fastq}'

rule merge:
	input: laneA='fastq/{sample}/{file}_L{laneA}_R{R}_001.fastq', laneB='fastq/{sample}/{file}_L{laneB}_R{R}_001.fastq'
	output: 'merged/{sample}/{file}_L{laneA}_L{laneB}_R{R}.fastq'
	shell: 'cat {input.laneA} {input.laneB} > {output}'
	
rule index:
	"""index reference (takes a while) """
	input: 'indexed/hg18.fa'
	output: 'indexed/hg18.pac', 'indexed/hg18.amb', 'indexed/hg18.ann'
	shell: 'bwa index -a bwtsw -p indexed/hg19 indexed/hg19.fa'

rule makeFa:
	"""merge reference (http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/chromFa.zip) to single fasta """
	input: expand('chromFa/chr{chr}.fa', chr=list(range(1,22)) + ['X', 'Y','M'])
	output: 'indexed/hg18.fa'
	shell: 'cat {input} > {output}'
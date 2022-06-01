// wgs_to_bam_test_script

ref = file(params.ref)

/*
Create 'read_pairs' channel that emits for each read pair a
tuple containing 3 elements: pair_id, R1, R2
*/

Channel
    .fromFilePairs(params.reads)
    .ifEmpty{error "Cannot find any reads matching: ${params.reads}"}
    .set{read_pairs_ch}

//trimmomatic read trimming
process trimreads {
	tag "trim ${pair_id}"	
	
	publishDir "${params.outdir}/$pair_id",
		saveAs: {filename ->
			if (filename.indexOf("_paired.fq.gz") > 0) "trimmed_pairs/$filename"
			else if (filename.indexOf("_unpaired.fq.gz") > 0) "unpaired/$filename"
			else filename
	}
		
	input:
	set pair_id, file(reads) from read_pairs_ch

	output:
	set pair_id, file("trimmed_${pair_id}_R{1,2}_paired.fq.gz") into (trimmed_reads_ch, trimmed_reads2_ch)
	file("trimmed_${pair_id}_R{1,2}_unpaired.fq.gz")

	conda 'bioconda::trimmomatic'

	script:
	"""
	#!/usr/bin/env bash

	trimmomatic PE ${reads[0]} ${reads[1]} \
	"trimmed_${pair_id}_R1_paired.fq.gz" "trimmed_${pair_id}_R1_unpaired.fq.gz" \
	"trimmed_${pair_id}_R2_paired.fq.gz" "trimmed_${pair_id}_R2_unpaired.fq.gz" \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:3 SLIDINGWINDOW:5:20 -threads 12
	"""
}

//fastqc on each trimmed read pair
process fastqc {
    tag "FASTQC on ${pair_id}"

    publishDir "${params.outdir}/$pair_id"

    input:
    set pair_id, file(reads) from trimmed_reads_ch

    output:
    set pair_id, file("fastqc_${pair_id}") into fastqc_ch

    conda 'bioconda::fastqc'
   
    script:
    """
    #!/usr/bin/env bash

    mkdir fastqc_${pair_id}
    fastqc -o fastqc_${pair_id} -q ${reads}
    """  
}  

//multiqc report
process multiqc {
	tag "multiqc on all trimmed_fastqs"

    publishDir params.outdir, mode:'copy'
       
    input:
    file ('*') from fastqc_ch.collect()
    
    output:
    file('multiqc_report.html')  
    
    conda 'bioconda::multiqc'

    script:
    """
	#!/usr/bin/env bash
 	
	multiqc .
    """
}

/*
// bwa alignment
process bwa_align {
	tag "align ${pair_id}"

	publishDir "${params.outdir}/$pair_id"

	input:
	set pair_id, file(reads) from trimmed_reads2_ch

	output:
	set pair_id, "${pair_id}.bam" into bam_ch
	file "${pair_id}.bam.bai"

	conda 'bioconda::bwa bioconda::samtools'

	time '2h'
	cpus 8
	penv 'smp' 
	memory '16 GB'

	script:

	"""
	#!/usr/bin/env bash

	# alignment
	bwa mem -t 10 -M -R "@RG\tID:"${pair_id}"\tLB:"${pair_id}"\tPL:illumina\tSM:"${pair_id}"\tPU:"${pair_id}""
	${ref} ${reads} | \
		samtools view -b -o temp.bam

	# sorting reads
	samtools sort -@ 8 -o ${pair_id}.bam temp.bam

	# index bam
	samtools index ${pair_id}.bam

	# remove intermediary bams
	rm temp.bam
	"""
}
*/

workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

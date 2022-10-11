// kniare part I of WGS processing for variant calling

ref = file(params.ref)
refbed = file(params.refbed)
refdir = file(params.refdir)

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
	set pair_id, file("trimmed_${pair_id}_R{1,2}_paired.fq.gz") into trimmed_reads_ch, trimmed_reads2_ch
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

    publishDir "${params.outdir}/${pair_id}/fastqc"

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
    file ("fastqc_*") from fastqc_ch.collect()
    
    output:
    file('multiqc_report.html')  
    
    conda 'bioconda::multiqc'

    script:
    """
    #!/usr/bin/env bash
    
    multiqc .
    
    """
}


// bwa alignment
process bwa_align {
	
	tag "align ${pair_id}"

	publishDir "${params.outdir}/$pair_id"

	input:
	set pair_id, file(reads) from trimmed_reads2_ch

	output:
	set pair_id, file("${pair_id}.sorted.dup.bam") into bam_ch
	//file("${pair_id}.sorted.bam.bai")
	file("${pair_id}.sam")
	file("${pair_id}.sorted.bam")
	file("${pair_id}.bam")
	file("${pair_id}.clean.bam")

	time '1h'
	cpus 8
	penv 'smp' 
	memory '16 GB'

	script:

	"""
	#!/usr/bin/env bash

	# load modules
	module load CBI bwa gatk/4.2.2.0

	# alignment and populate read group header
	bwa mem -t 8 -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tSM:${pair_id}\\tPU:${pair_id}" $ref ${reads} > ${pair_id}.sam

	# sam file sorting
	gatk --java-options "-Xmx40g -Xms40g" SamFormatConverter -R $ref -I ${pair_id}.sam -O ${pair_id}.bam
    gatk --java-options "-Xmx40g -Xms40g" CleanSam -R $ref -I ${pair_id}.bam -O ${pair_id}.clean.bam
    gatk --java-options "-Xmx40g -Xms40g" SortSam -R $ref -I ${pair_id}.clean.bam -O ${pair_id}.sorted.bam -SO coordinate --CREATE_INDEX true
    gatk --java-options "-Xmx40g -Xms40g" MarkDuplicatesSpark -R $ref -I ${pair_id}.sorted.bam -O ${pair_id}.sorted.dup.bam

    # remove intermediary bams
    # rm ${pair_id}.sam ${pair_id}.bam ${pair_id}.clean.bam
	"""
}


// samtools sorting Pf and human reads
process sort_pf_human {
	
	tag "sort PfHs ${pair_id}"

	publishDir "${params.outdir}/$pair_id"

	input:
	set pair_id, file(bams) from bam_ch

	output:
	set pair_id, file("${pair_id}.sorted.dup.pf.bam") into pf_bam_ch1, pf_bam_ch2, pf_bam_ch3
	set pair_id, file("${pair_id}.sorted.dup.hs.bam") into hs_bam_ch
	//file("${pair_id}.sorted.dup.pf.bam.csi")


	time '1h'
	cpus 8
	penv 'smp' 
	memory '16 GB'

	script:

	"""
	#!/usr/bin/env bash

	# load modules
	module load CBI samtools gatk/4.2.2.0 

	# sorting of Pf and Hs aligned reads

	samtools view -b -h ${pair_id}.sorted.dup.bam -T $ref -L $refbed/Pf3D7_core.bed > ${pair_id}.sorted.dup.pf.bam
	samtools view -b -h ${pair_id}.sorted.dup.bam -T $ref -L $refbed/human.bed > ${pair_id}.sorted.dup.hs.bam
	
	# samtools index -bc ${pair_id}.sorted.dup.pf.bam
	
	# rm ${pair_id}.sorted.dup.bam
	"""	
}


// distribution of Pf read depth by chromosome
process pf_read_depth {
	
	tag "read depth Pf chroms ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir/by_chrom"

	input:
	set pair_id, file(bam) from pf_bam_ch1

	output:
	file("ReadCoverage_final_${pair_id}.tsv")

	script:
	"""
	#!/usr/bin/env bash

	# load modules
	module load CBI samtools gatk/4.2.2.0

	samtools index -bc $bam

	for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	    do
	       gatk --java-options "-Xmx80g -Xms80g" DepthOfCoverage -R "$refdir/Pf3D7.fasta" -O chr"\$i" -L Pf3D7_"\$i"_v3 --omit-locus-table true -I $bam
	       
	       awk -F"," -v OFS="\t" '{ print \$0, \$(NF+1) = '"chr\$i"' }' chr"\$i".sample_summary > chr"\$i".sample2_summary
	    done

	cat *.sample2_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15, chromosome' > ReadCoverage_final_${pair_id}.tsv

	"""
}

/*
// pf read depth summary
process pf_read_depth_summary {
	
	tag "read depth summary ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"

	input:
	set pair_id, file(chrom_summary) from read_cover_ch.collect()

	output:
	file("ReadCoverage_final_${pair_id}.tsv")

	script:
	"""
	#!/usr/bin/env bash

	cat $chrom_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15, chromosome' > ReadCoverage_final_${pair_id}.tsv
	"""
}
*/

// insert size calculation
process insert_sizes {
	
	tag "insert sizes ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"

	input:
	set pair_id, file(pf_bam) from pf_bam_ch2

	output:
	file("${pair_id}.insert.txt")
	file("${pair_id}.insert2.txt") into inserts_ch
	
	script:

	"""
	#!/usr/bin/env bash

	# load modules
	module load CBI gatk/4.2.2.0

	gatk CollectInsertSizeMetrics -I $pf_bam -O ${pair_id}.insert.txt -H ${pair_id}_histo.pdf -M 0.05
	awk 'FNR>=8 && FNR<=8 {print \$1,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$NF="${pair_id}"}' ${pair_id}.insert.txt > ${pair_id}.insert2.txt
	
	#rm ${pair_id}.insert.txt 
	"""
}

// insert summary
process insert_summary {

	tag "insert summary"

	publishDir params.outdir, mode:'copy'

	input:
	file('*.insert2.txt') from inserts_ch.collect()

	output:
	file('InsertSizes_Final.txt')

	"""
	#!/usr/bin/env bash

	cat *.insert2.txt > InsertSizes_Final.txt
	"""
}

// Pf bam statistics by sample
process pf_bam_stat_per_sample {
	
	tag "Pf bam stat ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"
	
	input: 
	set pair_id, file(sorted_pf_bam) from pf_bam_ch3

	output:
	file("${pair_id}_bamstat_pf_final.tsv") into pf_bamstat_ch
	file("${pair_id}_bamstat_pf.tsv")

	conda 'bioconda::samtools bioconda::datamash'

	script:
	"""
	#!/usr/bin/env bash

    samtools stats $sorted_pf_bam | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_pf.tsv
    
    datamash transpose < ${pair_id}_bamstat_pf.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_pf_final.tsv

    #rm ${pair_id}_bamstat_pf.tsv
    """
}

// Pf bam statistic summary
process pf_stat_summary {
	
	tag "Pf stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	file(bamstat_pf) from pf_bamstat_ch.collect()

	output:
	file('Bam_stat_summary_pf_final.tsv') into pf_summary_ch

	script:
	"""
	#!/usr/bin/env bash

	cat $bamstat_pf | sed '1irow_total_reads_pf,	filtered_reads_pf,	sequences_pf,	is_sorted_pf,	1st_fragments_pf,	last_fragments_pf,	reads_mapped_pf,	reads_mapped_and_paired_pf,	reads_unmapped_pf,	reads_properly_paired_pf,	reads_paired_pf,	reads_duplicated_pf,	reads_MQ0_pf,	reads_QC_failed_pf,	non_primary_alignments_pf,	supplementary_alignments_pf,	total_length_pf,	total_first_fragment_length_pf,	total_last_fragment_length_pf,	bases_mapped_pf,	bases_mapped_(cigar)_pf,	bases_trimmed_pf,	bases_duplicated_pf,	mismatches_pf,	error_rate_pf,	average_length_pf,	average_first_fragment_length_pf,	average_last_fragment_length_pf,	maximum_length_pf,	maximum_first_fragment_length_pf,	maximum_last_fragment_length_pf,	average_quality_pf,	insert_size_average_pf,	insert_size_standard_deviation_pf,	inward_oriented pairs_pf,	outward_oriented_pairs_pf,	pairs_with_other_orientation_pf,	pairs_on_different_chromosomes_pf,	percentage_of_properly_paired_reads_(%)_pf, 	sample_name' > Bam_stat_summary_pf_final.tsv

	#rm *_bamstat_pf_final.tsv
	"""
}

// Hs bam statistics by sample
process hs_bam_stat_per_sample {

	tag "Hs bam stat ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"
	
	input:
	set pair_id, file(sorted_hs_bam) from hs_bam_ch

	output:
	file("${pair_id}_bamstat_hs_final.tsv") into hs_bamstat_ch
	file("${pair_id}_bamstat_hs.tsv")

	conda 'bioconda::samtools bioconda::datamash'

	script:
	"""
	#!/usr/bin/env bash

    samtools stats $sorted_hs_bam | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_hs.tsv

    datamash transpose < ${pair_id}_bamstat_hs.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_hs_final.tsv

    #rm ${pair_id}_bamstat_hs.tsv
    """
}


// Hs bam statistic summary
process hs_stat_summary {
	tag "Hs stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	file(bamstat_hs) from hs_bamstat_ch.collect()

	output:
	file('Bam_stat_summary_hs_final.tsv') into hs_summary_ch

	script:
	"""
	#!/usr/bin/env bash

	cat $bamstat_hs | sed '1irow_total_reads_hs,	filtered_reads_hs,	sequences_hs,	is_sorted_hs,	1st_fragments_hs,	last_fragments_hs,	reads_mapped_hs,	reads_mapped_and_paired_hs,	reads_unmapped_hs,	reads_properly_paired_hs,	reads_paired_hs,	reads_duplicated_hs,	reads_MQ0_hs,	reads_QC_failed_hs,	non_primary_alignments_hs,	supplementary_alignments_hs,	total_length_hs,	total_first_fragment_length_hs,	total_last_fragment_length_hs,	bases_mapped_hs,	bases_mapped_(cigar)_hs,	bases_trimmed_hs,	bases_duplicated_hs,	mismatches_hs,	error_rate_hs,	average_length_hs,	average_first_fragment_length_hs,	average_last_fragment_length_hs,	maximum_length_hs,	maximum_first_fragment_length_hs,	maximum_last_fragment_length_hs,	average_quality_hs,	insert_size_average_hs,	insert_size_standard_deviation_hs,	inward_oriented pairs_hs,	outward_oriented_pairs_hs,	pairs_with_other_orientation_hs,	pairs_on_different_chromosomes_hs,	percentage_of_properly_paired_reads_(%)_hs, 	sample_name' > Bam_stat_summary_hs_final.tsv

	#rm *_bamstat_hs_final.tsv
	"""
}


// Pf:Hs read ratio calculation
process pf_hs_ratio_calc {
	
	tag "Pf:Hs ratio"

	publishDir params.outdir, mode:'copy'
	
	input: 
	file(bamsum_pf) from pf_summary_ch
	file(bamsum_hs) from hs_summary_ch

	output:
	file('ratios_hs_pf_reads.tsv')

	script:
	"""
	#!/usr/bin/env bash

	pr -m -t -s\\ $bamsum_pf $bamsum_hs | gawk '{print \$7,\$46}' | awk '!/reads_mapped/ {print \$0}' | awk OFS="\t" '{print \$1, \$2, \$2/\$1}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > ratios_hs_pf_reads.tsv
	"""
}


// Rmd run quality report generation
process run_report {
 
	tag "Run quality report"

	publishDir params.outdir, mode:'copy'

	input:
	file('Bam_stat_summary_pf_final.tsv') from pf_summary_ch.collect()
	file('Bam_stat_summary_hs_final.tsv') from hs_summary_ch.collect()

	output:
	file('run_quality_report.html')

	script:
	"""
	#!/usr/bin/env bash

	Rscript -e "rmarkdown::render('run_quality_report.Rmd')"
	"""
}


workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the these reports in your browser --> \nmultiqc = $params.outdir/multiqc_report.html\nrun quality report = $params.outdir/reportrmd.html\nnextflow summary = $params.outdir/report.html": "Oops .. something went wrong" )
}

params.refgenome = "/home/mahe/eesha/chr7.fa"
params.outdir = "/home/mahe/eesha/results/"
params.inputdir = "/home/mahe/eesha/"
params.sampleid = "Sample-4BC-SDNA"
params.vcf_file = "/home/mahe/eesha/known_sites.vcf"

log.info"""\
refgenome: ${params.refgenome}
outdir: ${params.outdir}
inputdir: ${params.inputdir}
sampleid: ${params.sampleid}
vcf_file: ${params.vcf_file}
"""

process indexing {
    container 'biocontainers/bwa:v0.7.17_cv1'

    publishDir params.outdir, mode: "copy"   

    input:
    path refgenome

    output:
    path 'bwa_indexes'

    script:
    """
    mkdir -p bwa_indexes
    bwa index ${refgenome} -p ${refgenome}
    mv ${refgenome.baseName}* ./bwa_indexes
    """
}

process samtools_index {
    container 'mgibio/samtools:1.15.1-buster'

    input:
    path refgenome

    output:
    path 'bwa_indexes'

    script:
    """
    mkdir -p bwa_indexes
    samtools faidx ${refgenome} -o ./bwa_indexes/${refgenome.baseName}.fai
    samtools dict ${refgenome} --output ./bwa_indexes/${refgenome.baseName}.dict
    """
}

process bwa_alignment {
    container 'biocontainers/bwa:v0.7.17_cv1'

    publishDir params.outdir, mode: "copy"

    input:
    path bwa_indexes
    path inputdir
    val sampleid

    output:
    path 'alignment_results'

    script:
    """
    mkdir alignment_results
    bwa mem -M -R "@RG\\tID:${sampleid}\\tSM:${sampleid}" $bwa_indexes/*.fa ${inputdir}/${sampleid}_S1_L001_R1_001.fastq.gz ${inputdir}/${sampleid}_S1_L001_R2_001.fastq.gz > alignment_${sampleid}.sam
    mv alignment_${sampleid}.sam ./alignment_results
    """
}

process samtools_sort_index {
    container 'mgibio/samtools:1.15.1-buster'

    input:
    path alignment_results
    val sampleid

    output:
    path 'alignment_results'

    script:
    """
    samtools view -b -S $alignment_results/alignment_${sampleid}.sam | samtools sort -o sorted_${sampleid}.bam
    samtools index sorted_${sampleid}.bam
    mv sorted_${sampleid}* ./alignment_results
    """
}

process markdups_gatk {

    container 'broadinstitute/gatk:4.5.0.0'

    publishDir params.outdir, mode: "copy"

    input:
    path alignment_results

    output:
    path 'dedup_result'

    script:
    """
    mkdir dedup_result
    gatk MarkDuplicates -I $alignment_results/sorted_${params.sampleid}.bam -O sorted_dedup_${params.sampleid}.bam -M sorted_dedup_metrics_${params.sampleid}.txt
    mv sorted_dedup_* ./dedup_result
    """
}

process markdups_samtools_index {

    container 'mgibio/samtools:1.15.1-buster'

    publishDir params.outdir, mode: "copy"

    input:
    path dedup_result

    output:
    path 'dedup_result'

    script:
    """
    samtools index $dedup_result/sorted_dedup_${params.sampleid}.bam
    mv sorted_dedup_* ./dedup_result
    """
}


process readgroups {

  container 'broadinstitute/gatk:4.5.0.0'

  publishDir params.outdir, mode: "copy"

  input:
    path dedup_result
    path vcf_file

  output:
    path 'dedup_result'

  script:
    """
    gatk AddOrReplaceReadGroups -I $dedup_result/sorted_dedup_${params.sampleid}.bam -O rg_sorted_dedup_${params.sampleid}.bam -LB readgroup -PL illumina -PU ESWJFHU537GDFJK -SM Sample-SRA
    mv rg_sorted_dedup* ./dedup_result
    gatk IndexFeatureFile --input ${params.vcf_file}
    mv known_sites* ./dedup_result
    """
}

process baserecal {

  container 'broadinstitute/gatk:4.5.0.0'

  publishDir params.outdir, mode: "copy"

  input:
    path bwa_indexes
    path dedup_result
    path vcf_file

  output:
    path 'recal_report'
  script:
    """
    mkdir recal_report
    gatk BaseRecalibrator -R $bwa_indexes/chr7.fa -I $dedup_result/rg_sorted_dedup_${params.sampleid}.bam --known-sites ${params.vcf_file} -O recal_report_${params.sampleid}.table
    mv recal_report_* ./recal_report
    """
}

process apply_bqsr {

  container 'broadinstitute/gatk:4.5.0.0'

  publishDir params.outdir, mode:"copy"

  input:
    path bwa_indexes
    path dedup_result
    path recal_report

  output:
    path 'recal_report'

  script:
    """
    gatk ApplyBQSR -R $bwa_indexes/chr7.fa -I $dedup_result/rg_sorted_dedup_${params.sampleid}.bam --bqsr-recal-file $recal_report/recal_report_${params.sampleid}.table -O recalibrated_${params.sampleid}.bam
    mv recalibrated_* ./recal_report
    """
}

process haplo {

  container 'broadinstitute/gatk:4.5.0.0'

  publishDir params.outdir, mode:"copy"

  input:
    path bwa_indexes
    path recal_report
  
  output:
    path 'variant_file'

  script:
    """
    mkdir variant_file
    gatk HaplotypeCaller -R $bwa_indexes/chr7.fa -I $recal_report/recalibrated_${params.sampleid}.bam -O variants_${params.sampleid}.vcf
    mv variants_* ./variant_file
    """
}

workflow {
  bwa_index_files = indexing(params.refgenome)
  samtools_indexes = samtools_index(params.refgenome)
  sorted_bam_file = bwa_alignment(bwa_index_files, params.inputdir, params.sampleid)
  sorted_bam_with_index = samtools_sort_index(sorted_bam_file, params.sampleid)
  dedup_file = markdups_gatk(sorted_bam_with_index)
  dedup_file_with_index = markdups_samtools_index(dedup_file)
  rg_sorted_file = readgroups(dedup_file, params.vcf_file)
  recalibration_report = baserecal(bwa_index_files, rg_sorted_file, params.vcf_file)
  recalibrated_bam_file = apply_bqsr(bwa_index_files, rg_sorted_file, recalibration_report)
  haplo(bwa_index_files, recalibrated_bam_file)
}

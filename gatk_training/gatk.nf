params.refgenome = "/home/eesha/my_pipelines/gatk_nxf_training/chr7.fa"
params.outdir = "/home/eesha/my_pipelines/gatk_results/"
params.inputdir = "/home/eesha/my_pipelines/gatk_nxf_training/"
params.sampleid = "Sample-4BC-SDNA"
params.vcf_file = "/home/eesha/my_pipelines/gatk_nxf_training/known_sites.vcf"

log.info"""\
refgenome: ${params.refgenome}
outdir: ${params.outdir}
inputdir: ${params.inputdir}
sampleid: ${params.sampleigd}
vcf_file: ${params.vcf_file}
"""
 
process indexing {

  publishDir params.outdir, mode: "copy"   

  input:
    path refgenome

  output:
    path 'bwa_indexes'

  script:
    """
    mkdir bwa_indexes
    bwa index ${refgenome} -p ${refgenome}
    samtools faidx ${refgenome}
    samtools dict ${refgenome} --output ${refgenome.baseName}.dict
    mv ${refgenome.baseName}* ./bwa_indexes
    """
}

process alignment {

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
    bwa mem -M -R "@RG\\tID:${sampleid}_S1_L001\\tSM:${sampleid}_S1_L001" $bwa_indexes/*.fa ${inputdir}/${sampleid}_S1_L001_R1_001.fastq.gz ${inputdir}/${sampleid}_S1_L001_R2_001.fastq.gz | samtools sort -o sorted_${sampleid}.bam
    samtools index sorted_${params.sampleid}.bam
    mv sorted_${params.sampleid}* ./alignment_results
    """
}

process markdups {
  publishDir params.outdir, mode: "copy"

  input:
    path alignment_results

  output:
    path 'dedup_result'
  
  script:
    """
    mkdir dedup_result
    gatk MarkDuplicates -I $alignment_results/sorted_${params.sampleid}.bam -O sorted_dedup_${params.sampleid}.bam -M sorted_dedup_metrics_${params.sampleid}.txt
    samtools index sorted_dedup_${params.sampleid}.bam
    mv sorted_dedup_* ./dedup_result
    """
}

process readgroups {
  publishDir params.outdir, mode: "copy"

  input:
    path dedup_result
    path vcf_file

  output:
    path 'dedup_result'

  script:
    """
    gatk AddOrReplaceReadGroups -I $dedup_result/sorted_dedup_${params.sampleid}.bam -O rg_sorted_dedup_${params.sampleid}.bam -LB readgroup -PL illumina -PU ESWJFHU537GDFJK -SM Sample-4BC
    mv rg_sorted_dedup* ./dedup_result
    gatk IndexFeatureFile --input ${params.vcf_file}
    mv known_sites* ./dedup_result
    """
}

process baserecal {
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
    tree bwa_indexes > abc.txt
    gatk BaseRecalibrator -R $bwa_indexes/chr7.fa -I $dedup_result/rg_sorted_dedup_${params.sampleid}.bam --known-sites ${params.vcf_file} -O recal_report_${params.sampleid}.table
    mv recal_report_* ./recal_report
    """
}

process apply_bqsr {

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
  sorted_bam_file = alignment(bwa_index_files, params.inputdir, params.sampleid)
  dedup_file = markdups(sorted_bam_file)
  rg_sorted_file = readgroups(dedup_file, params.vcf_file)
  recalibration_report = baserecal(bwa_index_files, rg_sorted_file, params.vcf_file)
  recalibrated_bam_file = apply_bqsr(bwa_index_files, rg_sorted_file, recalibration_report)
  haplo(bwa_index_files, recalibrated_bam_file)
}
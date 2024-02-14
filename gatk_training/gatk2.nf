params.refgenome = "/home/eesha/my_pipelines/gatk_nxf_training/chr7.fa"
params.outdir = "/home/eesha/my_pipelines/gatk_results/"
params.inputdir = "/home/eesha/my_pipelines/gatk_nxf_training/"
params.sampleid = "Sample-4BC-SDNA"
params.vcf_file = "/home/eesha/my_pipelines/gatk_nxf_training/known_sites.vcf"

log.info"""\
refgenome: ${params.refgenome}
outdir: ${params.outdir}
inputdir: ${params.inputdir}
sampleid: ${params.sampleid}
vcf_file: ${params.vcf_file}
"""
 
process indexing {

  publishDir params.outdir, mode: "copy"   

  input:
    path refgenome

  output:
    path 'bwa_index'

  script:
    """
    mkdir bwa_index
    bwa index ${refgenome} -p ${refgenome}
    samtools faidx ${refgenome}
    samtools dict ${refgenome} -o ${refgenome}.dict
    mv ${refgenome}* ./bwa_index
    """
}

process alignment {

  publishDir params.outdir, mode: "copy"

  input:
    path bwa_index
    path inputdir
    val sampleid

  output:
    path 'align_results'

  script:
    """
    mkdir align_results
    bwa mem -M -R "@RG\\tID:${sampleid}_S1_L001\\tSM:${sampleid}_S1_L001" $bwa_index/chr7.fa ${inputdir}/${sampleid}_S1_L001_R1_001.fastq.gz ${inputdir}/${sampleid}_S1_L001_R2_001.fastq.gz | samtools sort -o sorted_${sampleid}.bam
    mv sorted_${sampleid}* ./align_results
    """
  
}

process indexing_bam {
  publishDir params.outdir, mode: "copy"

  input:
    path align_results

  script:
    """
    samtools index $align_results/sorted_${params.sampleid}.bam
    """
}

process markdups {
  publishDir params.outdir, mode: "copy"

  input:
    path align_results

  output:
    path 'dedup_results'
  
  script:
    """
    mkdir dedup_results
    gatk MarkDuplicates -I $align_results/sorted_${params.sampleid}.bam -O sorted_dedup_${params.sampleid}.bam -M sorted_dedup_metrics_${params.sampleid}.txt
    mv sorted_dedup_* ./dedup_results
    """
}

process indexing_dedup_bam {
  publishDir params.outdir, mode: "copy"

  input:
    path dedup_results

  script:
    """
    samtools index $dedup_results/sorted_dedup_${params.sampleid}.bam
    """
}

process readgroups {
  publishDir params.outdir, mode: "copy"

  input:
    path dedup_results

  output:
    path 'dedup2'

  script:
    """
    mkdir dedup2
    gatk AddOrReplaceReadGroups -I $dedup_results/sorted_dedup_${params.sampleid}.bam -O output_sorted_dedup_${params.sampleid}.bam -LB readgroup -PL illumina -PU ESWJFHU537GDFJK -SM Sample-4BC
    mv output_sorted_dedup* ./dedup2 
    """
}

process index_feature_file {
  publishDir params.outdir, mode: "copy"

  input:
    path vcf_file
  
  script:
    """
    gatk IndexFeatureFile --input ${vcf_file}
    """
}

process baserecal {
  publishDir params.outdir, mode: "copy"

  input:
    path dedup2
    path bwa_index
  output:
    path 'recal_report'
  script:
    """
    mkdir recal_report
    gatk BaseRecalibrator -R $bwa_index/chr7.fa -I $dedup2/output_sorted_dedup_${params.sampleid}.bam --known-sites ${params.vcf_file} -O recal_report_${params.sampleid}.table
    mv recal_report_* ./recal_report
    """
}

workflow {
  var1 = indexing(params.refgenome)
  var2 = alignment(var1, params.inputdir, params.sampleid)
  indexing_bam(var2)
  var3 = markdups(var2)
  indexing_dedup_bam(var3)
  var4 = readgroups(var3)
  index_feature_file(params.vcf_file)
  baserecal(var4, var1)
}
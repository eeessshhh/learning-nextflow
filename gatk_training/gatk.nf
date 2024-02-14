params.inputdir = "/home/eesha/my_pipelines/gatk_nxf_training/"
params.refgenome = "/home/eesha/my_pipelines/gatk_nxf_training/chr7.fa"
params.sampleid = "Sample-4BC-SDNA"
params.outdir = "/home/eesha/my_pipelines/gatk_results/"
params.vcf_file = "/home/eesha/my_pipelines/gatk_nxf_training/known_sites.vcf"

process fastq_to_vcf {

  publishDir params.outdir, mode:"copy"

  script:
    """
    bwa index ${params.refgenome}
    bwa mem -M -R "@RG\\tID:${params.sampleid}_S1_L001\\tSM:${paramsinputdir}/${params.sampleid}_S1_L001" chr7.fa ${params.inputdir}/${params.sampleid}_S1_L001_R1_001.fastq.gz ${params.sampleid}_S1_L001_R2_001.fastq.gz | samtools sort -o sorted_${params.sampleid}.bam
    samtools index sorted_${params.sampleid}.bam
    gatk MarkDuplicates -I sorted_${params.sampleid}.bam -O sorted_dedup_${params.sampleid}.bam -M metrics_${params.sampleid}.txt

    samtools index sorted_dedup_${params.sampleid}.bam
    samtools faidx chr7.fa
    samtools dict chr7.fa -o chr7.dict
    
    gatk AddOrReplaceReadGroups -I sorted_dedup_${params.sampleid}.bam -O output_sorted_dedup_${params.sampleid}.bam -LB readgroup -PL illumina -PU ESWJFHU537GDFJK -SM Sample-4BC
    
    gatk IndexFeatureFile --input known_sites.vcf
    
    gatk BaseRecalibrator -R chr7.fa -I output_sorted_dedup_${params.sampleid}.bam --known-sites known_sites.vcf -O recalibration_report_${params.sampleid}.table
    
    gatk ApplyBQSR -R chr7.fa -I output_sorted_dedup_${params.sampleid}.bam --bqsr-recal-file recalibration_report_${params.sampleid}.table -O recalibrated_${params.sampleid}.bam
    
    #gatk HaplotypeCaller -R chr7.fa -I recalibrated_${params.sampleid}.bam -O variants_${params.sampleid}.vcf
    
    """
}

workflow {
    fastq_to_vcf()
}
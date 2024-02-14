params.refgenome = "/home/eesha/my_pipelines/gatk_nxf_training/chr7.fa"
params.outdir = "/home/eesha/my_pipelines/gatk_results/"
params.inputdir = "/home/eesha/my_pipelines/gatk_nxf_training/"
params.sampleid = "Sample-4BC-SDNA"

log.info"""\
refgenome: ${params.refgenome}
outdir: ${params.outdir}
inputdir: ${params.inputdir}
sampleid: ${params.sampleid}
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

workflow {
  var1 = indexing(params.refgenome)
  var1.view()
  alignment(var1, params.inputdir, params.sampleid)
}
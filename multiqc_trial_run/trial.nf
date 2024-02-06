params.reads = "/home/eesha/my_pipelines/fastq_input_files/"
params.outdir = "/home/eesha/my_pipelines/results"

log.info """

reads: ${params.reads}

"""

process fastqc {

  publishDir params.outdir, mode:"copy"

  input:
    path reads
  output:
    path 'result'
  script:
    """
    mkdir result
    fastqc -o result ${reads}/*.fq
    """
}


process multiqc {
  publishDir params.outdir, mode: "copy"   

  input:
    path 'result'
  
  output:
    path 'multiqc_report'

  script:
    """
    multiqc result -o multiqc_report .
    """
}


workflow {
        ch_reads = Channel.fromPath("${params.reads}")
        ch_reads.view()
	fastqc_results = fastqc(ch_reads)
        multiqc(fastqc_results)

}
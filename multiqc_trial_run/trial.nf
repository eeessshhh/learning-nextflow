params.reads = "/home/eesha/my_pipelines/fastq_input_files/"
params.transcriptome = "/home/eesha/my_pipelines/transcriptome_file/transcriptome.fa"
params.outdir = "/home/eesha/my_pipelines/results"

log.info """\

reads:         ${params.reads}
transcriptome: ${params.transcriptome}

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


process index_transcriptome {

  publishDir params.outdir, mode: "copy"   


  output:
    path 'index_file'

  script:
    """
    echo "this works"
    salmon index -t ${transcriptome} -i index_file
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
        index_transcriptome(params.transcriptome)
        
}

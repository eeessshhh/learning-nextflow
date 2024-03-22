params.sample_name = "WLVAP-F1-D-CE-S1"
params.bam_file = "home/eesha/msi_pipeline/sorted_Sample-4BC-SDNA.bam"
params.outdir = "home/eesha/msi_pipeline/results"
params.microlist = "home/eesha/msi_pipeline/micro.list"

log.info"""\
sample_name: ${params.sample_name}
outdir: ${params.outdir}
bam_file: ${params.bam_file}
microlist: ${params.microlist}
"""

process msi_sensor {

  publishDir params.outdir, mode: "copy"   

  input:
    val sample_name
    path microlist
    path bam_file

  output:
    path 'result'

  script:
    """
    mkdir result
    msi -b 30 -d home/eesha/msi_pipeline/micro.list -t home/eesha/msi_pipeline/sorted_Sample-4BC-SDNA.bam -o ./${sample_name}_msi
    mv *_msi .result
    """
}

process python_script {
    script:
    """
    python3 MSI-script.py
    """
}

workflow {
    msi_sensor(params.sample_name)
    python_script()
}
# Stage 1: BWA
FROM docker.io/biocontainers/bwa:v0.7.17_cv1 AS bwa_builder

# Stage 2: Samtools
FROM mgibio/samtools:1.15.1-buster AS samtools_builder

# Stage 3: GATK
FROM broadinstitute/gatk:4.5.0.0 AS gatk_builder

# Final stage
FROM ubuntu:latest

# Copy BWA binary
COPY --from=bwa_builder /opt/conda/bin/bwa /usr/local/bin/

# Copy Samtools binary
COPY --from=samtools_builder /usr/local/bin/samtools /usr/local/bin/

# Copy GATK binary
COPY --from=gatk_builder /usr/local/bin/gatk* /usr/local/bin/

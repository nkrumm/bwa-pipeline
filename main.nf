fastq_pair_ch = Channel.fromFilePairs(params.input_folder + '*{1,2}.fastq.gz', flat: true)

process trim_reads {
    container 'nkrumm/atropos-paper:latest'

    input: 
      set pair_name, 
          file(fastq1), 
          file(fastq2) from fastq_pair_ch
    output:
      set pair_name, 
          file("*.trimmed.1.fastq.gz"), 
          file("*.trimmed.2.fastq.gz") into trimmed_pair_ch

    cpus 2
    memory 2.GB
    script:
    """
    atropos \
      -T ${task.cpus} \
      -a CTGTCTCTTATACACATCT \
      -A CTGTCTCTTATACACATCT \
      --cut 18 \
      -U 18 \
      --trim-n \
      -pe1 ${fastq1} \
      -pe2 ${fastq2} \
      -o ${pair_name}.trimmed.1.fastq.gz \
      -p ${pair_name}.trimmed.2.fastq.gz
    """
}

trimmed_pair_ch.into { fastqc_ch; bwa_ch }


process fastqc {
    container 'quay.io/biocontainers/fastqc:0.11.8--1'
    echo true
    memory '1 GB'
    input:
    set pair_name, file(fastq1), file(fastq2) from fastqc_ch

    output:
    file "fastqc_${pair_name}_output" into multiqc_ch

    script:
    """
    mkdir fastqc_${pair_name}_output
    fastqc -o fastqc_${pair_name}_output $fastq1 $fastq2
    """
}


process multiqc {
   container 'quay.io/biocontainers/multiqc:1.7--py_4'
   memory '4 GB'
   input:
   file('*') from multiqc_ch.collect()

   output:
   file "multiqc_report.html"

   publishDir "s3://uwlm-personal/nkrumm/GLT_v3_run/BAM/"

   script:
   """
   multiqc .
   """
}

process bwa_and_sort {
  container 'nkrumm/alignment'
  echo true
  input:
    set pair_name, file(fastq1), file(fastq2) from trimmed_pair_ch
    file(reference_fasta) from reference_fasta
    file("*") from reference_index.collect()
  output:
    set val(pair_name), file('*.sorted.bam') into mapped_bam_ch
    set val(pair_name), file('*.sorted.bam.bai') into mapped_bai_ch

  cpus 4
  memory { 8.GB * task.attempt }
  errorStrategy 'retry' 
  
  publishDir "s3://uwlm-personal/nkrumm/GLT_v3_run/BAM/"
  
  script:
  """
  seqtk mergepe ${fastq1} ${fastq2} \
  | bwa mem -p -t${task.cpus}  -R'@RG\\tID:${pair_name}\\tSM:${pair_name}' ${reference_fasta} - 2> log.txt \
  | samtools sort -t@${task.cpus} -m4G - -o ${pair_name}.sorted.bam
  
  samtools index ${pair_name}.sorted.bam
  """
}
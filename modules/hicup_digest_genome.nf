process DIGEST_GENOME {

  tag "${fasta}"

  input:
  path genome_fasta
  val re_pattern

  output:
  path "Digest*.txt", emit: genome_digest

  script:
  """
  hicup_digester --genome !{genomeName} --re1 !{re_pattern} !{genome_fasta}
  """
}

process DIGEST_GENOME {

  tag "${fasta}"

  input:
  file(fasta) from fastaForHicupDigest

  output:
  file("Digest*.txt") into hicupDigestIndex

  shell:
  """
  echo !{task.memory}
  echo !{task.cpus}
  hicup_digester --genome !{genomeName} --re1 !{params.re} !{fasta}
  """
}

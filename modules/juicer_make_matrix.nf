process MAKE_CONTACT_MATRIX {

  tag { name }

  publishDir path: "${params.outputDir}/${name}/matrices/",
             mode: 'copy',
             overwrite: 'true',
             pattern: "${name}/*.hic",
             saveAs: { filename ->
                           if (filename.endsWith(".hic")) file(filename).getName()
                     }
  when:
  !params.skip_juicer

  input:
  tuple val(name), file(pairs), file(pairsIndex) from resultsPairixJuicer
  file(chromSizeFile) from chromSizeChannelJuicer.collect()

  output:
  file("${name}/${name}.hic") into resultsJuicer

  shell:
  juicerGenomes = ['hg18', 'hg19', 'hg38', 'dMel',
                   'mm9', 'mm10', 'anasPlat1', 'bTaurus3',
                   'canFam3', 'equCab2', 'galGal4', 'Pf3D7',
                   'sacCer3', 'sCerS288c', 'susScr3', 'TAIR10']
  genome = juicerGenomes.contains(params.genome) ? params.genome : chromSizeFile
  juicerPath = "${NXF_HOME}/assets/pavrilab/hicer-nf/bin"
  '''
  mkdir -p !{name}

  java -Xmx!{task.memory.toGiga()}G -jar !{juicerPath}/juicer_tools_1.22.01.jar pre \
       -r !{resolutions} \
       -k KR,GW_KR \
       --threads !{task.cpus} \
       !{pairs} \
       !{name}/!{name}.hic \
       !{genome}
  '''
}

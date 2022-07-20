process XML_TO_TSV {

      tag "xml2tsv"

      input:
      path chromSizeXML

      output:
      path "chromSizes.tsv", emit: chromSizeTsv

      script:
      """
      xml2tsv.py ${chromSizeXML} chromSizes.tsv
      """
}

package PLS

object geneName {//extends Enumeration {
//  type Level = Value
//  val pmb,entrez,uniport, hugo_symbol,ensembl = Value
  val biomartField = Array("pmb","uniprot_swissprot","uniprot_genename","uniparc","external_gene_name","unigene","uniprot_sptrembl","entrezgene","ucsc","hgnc_id","hgnc_symbol","ensembl_gene_id","refseq_mrna")
  val mygeneinfoField = Array("entrezgene","hgnc","symbol","pathway")
  val pathwayList = Array("biocarta", "humancyc", "kegg", "netpath", "pharmgkb", "pid", "reactome","smpdb", "wikipathways")
  def findMatch(s:String,nameList:Array[String]):Array[String] = {
    val maxl:Int = nameList.map(_.intersect(s).length).max
    val slen:Int = s.length
    if (maxl < slen/2){
      throw new IllegalArgumentException("didn't find the name in filterList")
    }
    nameList.filter(_.intersect(s).length == maxl)
    //    if (maxl >= slen/2){ return rs }
    //      else {throw new IllegalArgumentException("didn't find the name in filterList")}
    //    if (rs.length == 1 & maxl >= slen/2) return rs(0)
    //    else if (rs.length == 0 | maxl < slen/2) throw new IllegalArgumentException("didn't find the name in filterList")
    //    else {
    //      val prin = rs.mkString(",")
    //      throw new IllegalArgumentException(s"find multiple matches:( $prin ) in the filterList")
    //    }
  }
}

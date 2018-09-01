package PLS
import java.io._
object gPms {
  val sepr = java.io.File.separator
  val cdr = new File(".").getAbsolutePath().split(sepr)
  val cdrr = cdr.slice(0,cdr.length-2)//.mkString(sepr)
  val homerr = cdr.slice(0,3).mkString(sepr) +sepr
  var tp:String = Array.concat(cdrr,Array("temp","PLStemp/")).mkString(sepr)
  var op:String = Array.concat(cdrr,Array("resources","PLSresources")).mkString(sepr)+sepr
  var rp:String = Array.concat(cdrr,Array("results","PLSresults")).mkString(sepr)+sepr
  var gp:String = Array.concat(cdrr,Array("resources","genomeRef")).mkString(sepr)+sepr
  var egf:String = "hgnc_ensg.txt"
  var ef:String = "gbmExp.txt"
  var eaf:String = "GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
  var pf:String = "GBMLGG.rppa.txt"
  var df:String = "gbmsnp6.txt"
  var af:String = "snp6annoNew.txt"
  var off:String = "snp6corr.txt"
  var cf:String = "gbmsnp6chrn.txt"
  var of:String = "snp6pval.txt"
  var rf:String = "geneSnpRange.txt"
  var ggf:String = "GOGeneList.txt"
  var gnf:String = "ensgGene.txt"
  var glf:String= "geneLoc.txt"
  var slf:String= "snpLoc.txt"
  var gdf:String = "GOidDescript.txt"
  var rsf:String = "GBMsnp6Rs.txt"
}

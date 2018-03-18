package PLS

import java.io.{FileWriter, PrintWriter}

import scala.io.Source._
import org.scalatest.FlatSpec
import PLS._
import PLS.vegas2._
import PLS.gPms._
import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import org.apache.commons.math3.distribution._

import scala.sys.process._
import scala.math.Numeric//{FDistribution, NormalDistribution}
class calculationTest extends FlatSpec {
 // val tp:String = "/Users/qunliu/Nutstore/workspace/temp/PLStemp/"
 // val rp:String = "/Users/qunliu/Nutstore/workspace/resources/PLSresources/"
 // val sp:String = "/Users/qunliu/Nutstore/workspace/results/PLSresults/"
 // val egf:String = "hgnc_ensg.txt"
 // val ef:String = "GBMLGG.transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.data.txt"
 // val pf:String = "GBMLGG.rppa.txt"
 // val df:String = "gbmsnp6.txt"
 // val af:String = "snp6annoNew.txt"
 // val off:String = "snp6corr.txt"
 // val cf:String = "gbmsnp6chrn.txt"
 // val of:String = "snp6pval.txt"
 // val rf:String = "geneSnpRange.txt"

  val go = scala.io.Source.fromFile(gPms.rp+gPms.ggf).getLines.map(_.split("\t")).filter(_.length >1).
    map{i => (i(0), i.tail.toSet)}.toArray

  val godes = scala.io.Source.fromFile(gPms.rp+gPms.gdf).getLines.map(_.split("\t",2)).filter(_.length >1).map(i =>i(0) -> i(1)).toMap
  val gog = go.flatMap(_._2).toSet
  val ezmap = scala.io.Source.fromFile(gPms.rp+gPms.gnf).getLines.map(_.split(" ")).map(i =>i(1) -> i(0)).toMap
  val alln = ezmap.keys.toSet
  val rs = scala.io.Source.fromFile(gPms.rp+"GBMsnp6RsN.txt").getLines.map(_.split("\t")).toArray
  val alg = rs.filter(i => i(9).toInt - i(8).toInt >1).map(i => ezmap.getOrElse(i(0),"")).toSet.intersect(gog)
  val sigg = rs.filter(i => i(9).toInt - i(8).toInt >1).filter(i => i(20).toFloat <0.05).map(i => ezmap.getOrElse(i(0),"")).toSet.intersect(gog)
  val goset = go.map(i => (i._1,i._2.intersect(alg))).filter(_._2.size > 1)
  val golen = goset.map(_._2.size)
  val gglen = goset.map(_._2.intersect(sigg).size)
  val len = golen.length
  val siglen = sigg.size
  val alllen = alg.size
  val rss = Array(0 until len:_*).map(i =>(goset(i)._1,godes.getOrElse(goset(i)._1,""), calculation.fisherTest(gglen(i),golen(i),siglen,alllen))).sortBy(_._3)
  val sigrs = rss.filter(_._3<0.05).map(_._1)
  val kapall = calculation.kappaArray(goset).filter(_._3>0.6)

  val kap = calculation.kappaArray(goset.filter(i=>sigrs.contains(i._1))).filter(_._3>0.6)
  var ii = 0
  var tes = rss
  var glen = tes.length
  while (ii< glen){
    val tss = tes(ii)._1
    val kapp = kapall.filter(i => i._1 == tss | i._2 == tss).flatMap(i => Array(i._1,i._2)).toSet - tss
    tes = tes.filter(i => !kapp.contains(i._1))
    ii += 1
    glen = tes.length
  }
  val flwriter = new PrintWriter(new FileWriter("goResult.txt"))
  tes.map(i => i._1 +":"+ i._2 +"\t"+ i._3.formatted("%.9f")).foreach(flwriter.println)
  writer.close
  val H:Array[Float] = Array(0.01f, 0.03f, 0.05f)
  //val writer =
  def simugenNo(glists:Array[String]) = {
    import java.io.{FileWriter, PrintWriter}
    val glist = glists.slice(0, 4)
    val writer = new PrintWriter(new FileWriter("goR"+glist(3)+".txt"))
    //println("processing No."+g)
        vegas2.simuFgene(glist)
    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
    //val sl = glist(4)
    //    writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("calculation ssstarting"))

      for (h <- H) {
        var i = 0
        while (i < rl) {
          var j = 0
          while(j < 100) {
            //val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
            val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno(h, i, false)) :+ h :+ i).mkString("\t")
            writer.println(rs) //foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
            j += 1
          }
          i += 1
        }
      }
    writer.close()
    }




  val xx = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).take(1).toArray.flatten
  val yy = scala.io.Source.fromFile(gPms.op+gPms.pf).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
  val ee = scala.io.Source.fromFile(gPms.op+gPms.ef).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
  val mcol =fileOper.intersectCols(xx,yy,ee)

  val pakt = fileOper.toArrays(gPms.op+gPms.pf).filter(_ (0).contains("AKT")).toArray.apply(1)
  val yyy = mcol._2.map(pakt(_).toFloat)
  val xxx = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).drop(1)
  val pval = xxx.map(i => i(0) -> calculation.linearPval(i,yyy,mcol._1)).toMap
  //pval.values.toArray.sorted

  val egmap = fileOper.toArrays(gPms.op+gPms.egf).map(i => i(0) -> i(1)).toMap
  def getXYZ(df:String,ef:String,pf:String,rf:String,gen:String,pheno:String): Unit ={
    val x = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).take(1).toArray.flatten
    val y = scala.io.Source.fromFile(gPms.op+gPms.pf).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val e = scala.io.Source.fromFile(gPms.op+gPms.ef).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val egmap = fileOper.toArrays(gPms.op+gPms.egf).map(i => i(0) -> i(1)).toMap
    val mcol = fileOper.intersectCols(x,y,e)
    val pakt = fileOper.toArrays(gPms.op+gPms.pf).filter(_ (0).contains(pheno)).toArray.apply(1)
    val yy = mcol._2.map(pakt(_).toFloat)
    val gxx = fileOper.toArrays(gPms.op+gPms.ef).filter(_ (0).contains(gen)).toArray.apply(0)
    val gg = mcol._3.map(gxx(_).toFloat)
    val eg = egmap.getOrElse(gen,"non")
    val range = fileOper.toArrays(gPms.op+gPms.rf).filter(_(0) == eg).toArray.flatten
    val ss = fileOper.toArrays(gPms.op+gPms.af).filter(i => i(1) == range(1) & range(2).toInt <= i(2).toInt & range(3).toInt >= i(2).toInt).map(_(0)).toSet
    val xxg = fileOper.toArrays(gPms.rp+gPms.df).filter(i => ss.contains(i(0))).toArray
    val xx = xxg.map(i => mcol._1.map(i(_).toFloat))
    return(xx,yy,gg)
  }
  val c1 = toPedP !
  val snpCC = scala.io.Source.fromFile(gPms.tp+"ex_CEU.out.gen").getLines.map(_.split(" "))
  //val sampl = scala.io.Source.fromFile("/Users/qunliu/hapgen2_macosx_intel/ex_CEU.out.sample").getLines.map(_.split(" "))
  //val b = Array(0,1,2)
  //val X = snpCas.map(_.drop(5).map(_.toInt).sliding(3,3).map(_.zip(b).map{case(a,b) => a*b}.reduce(_+_)).toArray)
  val Xx = snpCC.map(_.drop(5).map(_.toFloat).sliding(3,3).map(vegas2.snpVal(_).toFloat).toArray).toArray
  val X = utils.Array2DM(Xx,false)

  val R = org.ddahl.rscala.RClient()
  R.x = X(::,150).toArray
  R.y = X(::,0).toArray
  val rrs = R.evalD0("cor(x,y)")
  val h = 0.05f
  val Y = vegas2.setPheno()(X)
  val n = Y.rows
  val k = 3
  val num = 20
  val tenFold = plsCalc.kfoldInx(n, 10, true)
  val YpredTenFold = plsCalc.plsCV(X, Y, k, tenFold)
  val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, Y, k))
  val rssfold = 0.6f * sqrt(plsCalc.rss(Y, YpredTenFold))
  val t1 = System.nanoTime()
  val gdoffold = Array(0 until k: _*).map(i => plsCalc.gdof(X, Y, rssfold(i), i + 1, 1000, 1))
  val t2 = System.nanoTime()
  val lapse1 = (t2 - t1) / 1e9d
  val t3 = System.nanoTime()
  val gdoffold0 = Array(0 until k: _*).map(i => plsCalc.gdof(X, Y, rssfold(i), i + 1, 1006, 8))
  val t4 = System.nanoTime()
  val lapse2 = (t2 - t1) / 1e9d

  val X1 = X(0 until 500, 0 until 20)
  val pcs1 = princomp(convert(X1,Double)).scores
  val mav1 = breeze.stats.meanAndVariance(pcs1(::, 0))
  var pc11 = convert((pcs1(::, 0) - mav1.mean) / sqrt(mav1.variance),Float)
  val theta1 = getTheta(pc11)
  val theta2 = getTheta(pc11)
  val Y1 = (sqrt(h) * pc11 + sqrt(1 - h) * theta1).toDenseMatrix.t
  val Y2 = (sqrt(h) * Y1.toDenseVector + sqrt(1 - h) * theta2).toDenseMatrix.t
  val Y3 = (sqrt(h) * pc11 + sqrt(1 - h) * theta2).toDenseMatrix.t

  val n1 = Y1.rows
  val tenFold1 = plsCalc.kfoldInx(n1, 10, true)
  val YpredTenFold1 = plsCalc.plsCV(X1, Y1, k, tenFold1)
  val Yhat1 = plsCalc.predict(X1, plsCalc.plsTrain(X1, Y1, k))
  val rssfold1 = 0.6f * sqrt(plsCalc.rss(Y1, YpredTenFold1))
  val t11 = System.nanoTime()
  val gdoffold01 = Array(0 until k: _*).map(i => plsCalc.gdof(X1, Y1, rssfold1(i), i + 1, 1000, 1))
  val t12 = System.nanoTime()
  val lapse11 = (t12 - t11) / 1e9d
  plsCalc.gdofPlsPval(X1, Y1, k)
  val YY = DenseMatrix.horzcat(Y1, Y2)
  plsCalc.gdofPlsPval(X1, YY, k)
  plsCalc.gdofAll(X1, YY, rssfold1(0))

  //example.bed,.fam,.bim

  // will return 0 if succeed or 1 if failed
  
  //val eg = eig(X.t * X)
  //val a = eg.eigenvalues
  //val b = eg.eigenvectors
  //val writer = new PrintWriter(new FileWriter(off))
  var rrss = scala.collection.mutable.ArrayBuffer.empty[Array[String]]
  for(i <- 1 to 22){
    val fil = cf.replace("chrn","chr"+i)
    val xc = fromFile(tp+fil).getLines.map(_.split("\t"))
    val rss = calculation.corVal(xc,mcol._1)
    rrss ++= rss
  }
  val cormap = rrss.map(i => i(0) -> i.drop(1).mkString("\t")).toMap
  val writer = new PrintWriter(new FileWriter(rp+off))
  scala.io.Source.fromFile(rp+af).getLines.map(i => i+"\t"+pval.getOrElse(i.split("\t",2).apply(0),1.0)+"\t"+cormap.getOrElse(i.split("\t",2).apply(0),"null")).foreach(writer.println)
  writer.close

  //2018-1-7 compile the p value from vegas2 and pls with h = 0.02 and 1000 samples
  if(false) {
    import PLS._
    import java.io.{FileWriter, PrintWriter}
    val svd = fileOper.toArrays(rp+"GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs0 = svd.filter(_ (4).toInt > 50).sortBy(i => i(6).toFloat / i(5).toFloat)
    val rs1 = svd.sortBy(i => i(7).toFloat / i(5).toFloat).reverse
    val glist = rs1(0).slice(0, 4)
    vegas2.simuFgene(glist)

    val H = Array(0.01f,0.015f,0.02f)
    for (h <- H) {
      val writer = new PrintWriter(new FileWriter("rs" + fileOper.timeForFile + ".txt"))
      var i = 0
      while (i < 100) {
        writer.println(vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)).mkString("\t"))
        i += 1
      }
      writer.close()
    }

  }
//
  val x = Array(20,5,10,15)
  "kappaTest" should "get correct result" in {
    assert(abs(calculation.kappaTest(x) - 0.4) < 1e-6)
  }
  "kappaTest" should "handle array and seperate elements" in {
    assert(abs(calculation.kappaTest(x) -calculation.kappaTest(x(0),x(1),x(2),x(3))) < 1e-6)
  }
  "fisherTest" should "get correct result" in {
    assert(abs(calculation.fisherTest(x) - 0.008579) < 1e-2)
  }
}
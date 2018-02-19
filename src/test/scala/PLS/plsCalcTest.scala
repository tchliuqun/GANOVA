package PLS

import java.io.{FileWriter, PrintWriter}

import breeze.linalg._
import org.scalatest._
import PLS._
import breeze.numerics.sqrt
import breeze.stats.regression.leastSquares

class plsCalcTest extends FlatSpec {

  val testData = scala.io.Source.fromFile("LifeCycleSavings.txt").getLines.map(_.split("\t")).toArray
  val colName = testData(0)
  val tdata = testData.drop(1).map(_.drop(1).map(_.toFloat))
  val ncol = tdata(0).length
  val nrow = tdata.length
  val td = new DenseMatrix[Float](ncol,nrow,tdata.flatten).t
  val y = td(::,0)
  val x1 = td(::,1).toDenseMatrix.t
  val x = td(::,1 until ncol)

  val k = 4
  val snpt = new snpCalc
  val chrs = snpt.sfile.replace(".txt","chr16.txt")
  val chro = snpt.rfile.replace(".txt","chr16.txt")
  val chrd = snpt.dfile.replace(".csv","chr16.txt")
  val chrg = snpt.gfile.replace(".txt","chr16.txt")
  val pakt = fileOper.toArrays(snpt.pfile).filter(_(0).contains("AKT")).toArray
  val pn = fileOper.toArrays(snpt.pfile).next.map(_.slice(0,15))
  val dn = fileOper.toArrays(snpt.dfile,",").next.drop(2)
  val mcol = fileOper.intersectCols(dn,pn)
  val Y = new DenseMatrix(236,1,mcol._2.map(pakt(1)).map(_.toFloat))
  val n = Y.rows
  val looInx = Array(0 until n :_*).map(Seq(_))
  val tenFold = plsCalc.kfoldInx(n,10,true)
  val testof = "test16.txt"
  val phlpp2N = "ENSG00000040199"
  val genR = fileOper.FindIndexinRange(chrg,chrs,testof)
  val phlpp2 = genR.map(_.split("\t")).filter(_(0) == phlpp2N).flatten

  val snpd = fileOper.toArrays(chrd).toArray
  val snpa = fileOper.toArrays(chrs).toArray
  val tst = snpd.map(_(0)).deep == snpa.map(_(0)).deep
  val dm1 = snpd.map(_.drop(1).map(_.toFloat)).map(mcol._1.map(_))
  val dm = utils.Array2DM(dm1,false)

  val X = dm(::,phlpp2(4).toInt until phlpp2(5).toInt)
  val YpredTenFold = plsCalc.plsCV(X,Y,k,tenFold)
  val YpredLoo = plsCalc.plsCV(X,Y,k,looInx)
  val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,k))
  val pdofloo = plsCalc.pdof(Y,Yhat,YpredLoo)
  val pdoffold = plsCalc.pdof(Y,Yhat,YpredTenFold)
  val rsscv = 0.6*sqrt(plsCalc.rss(Y,YpredLoo))
  val gdofloo = Array(0 until k:_*).map(i => plsCalc.gdof(X,Y,rsscv(i),1000,i+1)).map(i => new DenseMatrix[Float](1,1,Array(i)))
  val gdofPval = Array(0 until k:_*).map(i => plsCalc.anova(Y,Yhat(::,i).toDenseMatrix.t,gdofloo(i))(0))
  val pdofPval =Array(0 until k:_*).map(i => plsCalc.anova(Y,Yhat(::,i).toDenseMatrix.t,new DenseMatrix[Float](1,1,Array(pdofloo(i))))(0))
  val tanova = Array(1 until k:_*).map(i => plsCalc.twoAnova(Y,Yhat(::,i-1).toDenseMatrix.t,Yhat(::,i).toDenseMatrix.t,gdofloo(i-1),gdofloo(i)))

  if(false) {
    import breeze.linalg._
    import PLS._
    val chrnum = "16"
    val snpt = new snpCalc

    val pakt = fileOper.toArrays(snpt.pfile).filter(_ (0).contains("AKT")).toArray.apply(1)
    val ys = 1.0 +: pakt.drop(1).map(_.toFloat)
    snpt.updateY(ys)
    snpt.updateFile(chrnum)
    val X = snpt.getXs
    if (false) {
      val phlpp2N = "ENSG00000040199"
      val chrs = snpt.sfile.replace(".txt", "chr16.txt")
      val chro = snpt.ofile.replace(".txt", "chr16.txt")
      val chrd = snpt.dfile.replace(".csv", "chr16.txt")
      val chrg = snpt.gfile.replace(".txt", "chr16.txt")
      val testof = "test16.txt"
      val genR = fileOper.FindIndexinRange(snpt.gfile, snpt.sfile, testof)
      val phlpp2 = genR.map(_.split("\t")).filter(_ (0) == phlpp2N).flatten.drop(1).map(_.toInt)
      val srtp = phlpp2(3)
      val endp = phlpp2(4)
    }
    val colname = Array("ucscName","chr","start","end","nsnp")
    val permn = Array(1 to snpt.k:_*).map("perm"+_)
    val pdofloo = Array(1 to snpt.k:_*).map("pdofLoo"+_)
    val pdoflooPval = Array(1 to snpt.k:_*).map("pdofLooPval"+_)
    val pdoften = Array(1 to snpt.k:_*).map("pdofTenf"+_)
    val pdoftenPval = Array(1 to snpt.k:_*).map("pdofTenfPval"+_)
    val gdofloo = Array(1 to snpt.k:_*).map("gdofLoo"+_)
    val gdoflooPval = Array(1 to snpt.k:_*).map("gdofLooPval"+_)
    val gdoften = Array(1 to snpt.k:_*).map("gdofTenf"+_)
    val gdoftenPval = Array(1 to snpt.k:_*).map("gdofTenfPval"+_)
    val colnames = colname ++ permn ++ gdofloo++gdoflooPval++gdoften++gdoftenPval++ pdofloo ++pdoflooPval++ pdoften++pdoftenPval
    val out = new PrintWriter(new FileWriter(snpt.ofile))
    out.println(colnames.mkString("\t"))
    val k = snpt.k
    val gens = fileOper.toArrays(snpt.rfile)
    while (gens.hasNext) {
      val gen = gens.next
      val srt = gen(4).toInt
      val end = gen(5).toInt
      val num = end - srt
      if (num > 1) {
        if(num < k) snpt.k = num
        val X1 = snpt.getX(srt, end, X)
        val permPval = snpt.permT(X1)
        val dofPred = snpt.getGeneSnp(X1)
        val pval = snpt.genePval(dofPred)
        val rs1 = pval.map(i => (i._2 ++ Array.fill(k-num)(0.0),i._1 ++ Array.fill(k-num)(1.0))).flatMap(i => i._1++i._2)
        val rs2 = permPval++Array.fill(k-num)(1.0) ++ rs1
        val rs3 = num +: rs2
        val rss = gen.take(4) ++ rs3.map(_.toString)
        out.println(rss.mkString("\t"))
      }
      snpt.k = k
    }
    out.close()

  }

  def fittedLSR:(DenseMatrix[Float],DenseVector[Float]) => DenseVector[Float] = (x,y) =>{
    val nrows = x.rows
    val xx = DenseMatrix.horzcat(DenseMatrix.ones[Float](nrows,1),x)//td(::,1 until ncol))
    val fit = leastSquares(xx,y)
    fit(xx)
  }

  //val fit1 = leastSquares(x1,y)
  //val fit = leastSquares(x,y)
  val yhat = fittedLSR(x,y)
  val yhat1 = fittedLSR(x1,y)
  val dof1 = DenseMatrix.ones[Float](1,1)
  val dof2 = DenseMatrix.ones[Float](1,1) * 4
  val myAnova = plsCalc.anova(y.toDenseMatrix.t,yhat1.toDenseMatrix.t,DenseMatrix.ones[Float](1,1))
  val mAnova = plsCalc.twoAnova(y.toDenseMatrix.t,yhat1.toDenseMatrix.t,yhat.toDenseMatrix.t,dof1,dof2)
  val mnova = plsCalc.twoAnova(y.toDenseMatrix.t,yhat.toDenseMatrix.t,yhat1.toDenseMatrix.t,dof2,dof1)

  "anova" should "get same result with R" in {
    assert(myAnova(0) - 0.0008866 < 1e-6)
  }
  "anova with two models" should "get same result with R" in {
    assert(mAnova(0) -0.04177 < 1e-6)
  }
  "anova with two models" should "get same result with different order" in {
    assert(mAnova(0) -mnova(0) < 1e-6)
  }
}

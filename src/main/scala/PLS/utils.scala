package PLS

import breeze.linalg._//{*, DenseMatrix, sum}
import breeze.stats.{mean, stddev}
import breeze.math._
import scala.collection.mutable.ArrayBuffer
import scala.annotation.tailrec
/**
  * Created by davidkellis on 11/16/16.
  */

object utils {
  def currentTime =java.time.LocalDateTime.now().toString.split("\\.").apply(0).replace("T","_").replace(":","_")
  def currentTimeIn =java.time.LocalDateTime.now().toString.split("\\.").apply(0).replace("T","_")+" "
  // perform column-wise sum; given M (RxC), yields a (1xC) matrix of per-column sums
  def sumColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = sum(M(::, *)).t.toDenseMatrix

  // perform column-wise mean; given M (RxC), yields a (1xC) matrix of per-column means
  def meanColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = mean(M(::, *)).t.asDenseMatrix

  // perform column-wise standard deviation; given M (RxC), yields a (1xC) matrix of per-column standard deviations
  def stddevColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = stddev(M(::, *)).t.asDenseMatrix.mapValues(_.toFloat)

  // vertically concatenate matrix M to itself n-1 times.
  // returns a matrix with n vertically stacked copies of M.
  def vertcat(M: DenseMatrix[Float], n: Int): DenseMatrix[Float] = {
    @tailrec
    def vertcatR(A: DenseMatrix[Float], i: Int): DenseMatrix[Float] = {
      if (i <= 1)
        A
      else
        vertcatR(DenseMatrix.vertcat(M, A), i - 1)
    }
    vertcatR(M, n)
  }
  def horzcat(M: DenseMatrix[Float], n: Int): DenseMatrix[Float] = {
    @tailrec
    def horzcatR(A: DenseMatrix[Float], i: Int): DenseMatrix[Float] = {
      if (i <= 1)
        A
      else
        horzcatR(DenseMatrix.horzcat(M, A), i - 1)
    }
    horzcatR(M, n)
  }
  def cvSplit[T](data:Iterable[T],count: Int): Stream[(Iterable[T], Iterable[T])] = {
    val size = data.size
    def piece(i: Int) = i * size / count
    Stream.range(0, size - 1) map { i =>
      val (prefix, rest) = data.splitAt(piece(i))
      val (test, postfix) = rest.splitAt(piece(i + 1) - piece(i))
      val train = prefix ++ postfix
      (test, train)
    }
  }
  def getGOMatrix(rsf:String = gPms.rp + "GBMsnp6Rs_2018_05_14_15_14_49_569.txt") = {
    var idg = scala.io.Source.fromFile(gPms.op + "goa_human.gaf").getLines.drop(23).map(_.split("\t")).map(i => (i(1), i(4))).toArray.groupBy(_._1).map(i => (i._1, i._2.map(_._2)))

    var idm = scala.io.Source.fromFile(gPms.op + "HUMAN_9606_idmapping_selected.tab").getLines.map(_.split("\t")).map(i => (i.filter(_.contains("ENSG")), i(0))).filter(i => i._1.length > 0 & i._2.length > 0).map(i => (i._1.apply(0), idg.getOrElse(i._2, Array()))).filter(_._2.length > 0).toArray
    idg = null
    var gom = scala.io.Source.fromFile(gPms.op + "goose").getLines.map(_.split("\t")).dropWhile(_ (0) == "all").map(i => Array(4, 0).map(i(_))).toArray.groupBy(_ (0)).map(i => (i._1, i._2.map(_ (1))))
    var genego = idm.map(i => (i._1, i._2.flatMap(i => gom.getOrElse(i, Array())).toSet))
    idm = null
    gom = null
    var goset = genego.map(_._2).flatten.toSet.toArray.sorted
    var gomxp = genego.map(i => (i._1, goset.map(ii => if (i._2.contains(ii)) 1f else 0f)))

    genego = null
    //goset = null
    var rs = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).toArray
    var rsset = rs.map(_ (3)).toSet.intersect(gomxp.map(_._1).toSet)
    var rss0 = rs.filter(i => rsset.contains(i(3)))
    var rsorder  = rs.map(i => if(i(5).toInt >2)i(9) else if(i(5).toInt ==2) i(8) else i(7)).map(_.toFloat)
    //val rss = Array( 0 until rss0.length :_*).map( i => rss0(i).slice(0,6) :+ rsorder(i).toString).sortBy(_(6).toFloat)
    rs = null
    var rss = Array( 0 until rss0.length :_*).map( i => rss0(i).slice(0,6) :+ rsorder(i).toString).filter(!_(6).toFloat.isNaN).sortBy(_(6).toFloat)
    val chisq = DenseVector(calculation.reverseChisq(rss.map(_(6).toFloat),1))


    rss0 = null
    rsorder = null
    var gomx = gomxp.filter(i => rsset.contains(i._1))
    rsset = null
    gomxp = null
    var indm = gomx.map(_._1).zipWithIndex.toMap
    var goinx = rss.map(i => indm(i(3))).map(gomx(_))//.filter(i => sum(i._2) > 15 & sum(i._2)<500)
    gomx = null
    rss = null
    var gm = new DenseMatrix(goinx(0)._2.length, goinx.length, goinx.flatMap(_._2)).t
    goinx = null
    val goidx = sum(gm(::,*)).t.toArray.zipWithIndex.filter(i => i._1 > 15 & i._1 < 500).map(_._2)
    var goano = scala.io.Source.fromFile(gPms.op + "goose").getLines.map(_.split("\t").take(3).mkString("\t")).toArray.distinct.map(_.split("\t")).map( i => (i(0) -> i)).toMap
    val gmm = gm(::,goidx.toIndexedSeq).toDenseMatrix

    gm = null
    val coln = goidx.map(goset(_)).map(goano(_))
    //goano = null

    //val golen = sum(gm(::,*))
    //val res = calculation.gseaP(gm,chisq)



    (chisq,gmm,coln)
  }
  def fisherPathway(rsf:String = gPms.rp + "GBMsnp6Rs_2018_05_14_15_14_49_569.txt",ngene:Int = 150,gonum:Int = 300,kshld:Float = 0.75f) = {
    val (rs,gs,ano) = utils.getGOMatrix(rsf)
    val indm = ano.map(_(0)).zipWithIndex.toMap
    val rss = calculation.fisherMatrix(gs,ano,ngene)
    val idx = rss.slice(0,gonum).map(_._1.apply(0)).map(indm(_)).toIndexedSeq
    val mxk = gs(::,idx).toDenseMatrix
    val didx = kappaMatrix(mxk,kshld)
    for ( i <- Array(0 until gonum :_*); if !didx.contains(i)) yield rss(i)
  }
  def cvSplitInx[T](size:Int,count: Int,random:Boolean = false): Array[Seq[Int]] = {
    import scala.util.Random
    val piece = size / count
    val indx = if(random) Random.shuffle(Seq(0 until size :_*)) else Seq(0 until size :_*)
    Array(0 until count:_*).map{i => indx.slice(i*piece, (i+1) * piece)}
  }
  def Array2DM(arr:Array[Array[Float]],transpose:Boolean = true):DenseMatrix[Float] = {
    val ncol = arr.length
    val nrow = arr(0).length
    if (transpose) (new DenseMatrix(nrow,ncol,arr.flatten)).t else new DenseMatrix(nrow,ncol,arr.flatten)
  }
  def dataSplit(arr:Array[Array[Float]],k:Int = 10):Array[(DenseMatrix[Float],DenseMatrix[Float])] = {
    val stream = cvSplit(arr,k).take(k).toArray
    stream.map{case(x,y) => (Array2DM(x.toArray),Array2DM(y.toArray))}
  }
  def subColMeans(X:DenseMatrix[Float]):DenseMatrix[Float] = {
    X - utils.vertcat(utils.meanColumns(X), X.rows)
  }
  def kappaMatrix(ds:DenseMatrix[Float],shld:Float = 0.75f) = {
    //val indm = rss._3.map(_(0)).zipWithIndex.toMap
    //rs.slice(0,300).map(_._1.apply(0)).map(indm(_)).toIndexedSeq
    val num = ds.cols - 1
    var rs = ArrayBuffer[Int]()
    var kp = 0f
    var i = num
    var j = i-1
    while (i >= 0 ){
      kp = 0f
      j = i -1
      val d1 = ds(::,i)
      while(kp < shld & j >= 0 ){
        kp = calculation.getkappa(d1,ds(::,j))
        if(kp >= shld) {
          rs :+= i
        }
        j -= 1
      }
      i -= 1
    }
    rs.toArray.reverse
    //ds.delete(res100.toIndexedSeq,Axis._1)
    //for ( i <- Array(0 until 300 :_*); if !rs.contains(i)) yield ds(i)
  }
  def getTimeForFile = {
    val today = java.time.LocalDate.now.toString.split("-")
    val totime = java.time.LocalTime.now.toString.replace(".",":").split(":")
    //val totime = java.time.LocalTime.now.toString.split("\\.")(0).split(":")
    (today++totime).mkString("_")
  }
  def colProduct(x:DenseMatrix[Float],y:DenseMatrix[Float]):DenseMatrix[Float] = {
    val x1 = x *:* y
    sum(x1(::,*)).t.toDenseMatrix
  }
}


package PLS

import breeze.linalg._//{*, DenseMatrix, sum}
import breeze.stats.{mean, stddev}
import breeze.math._

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
  def getTimeForFile = {
    val today = java.time.LocalDate.now.toString.split("-")
    val totime = java.time.LocalTime.now.toString.split("\\.")(0).split(":")
    (today++totime).mkString("_")
  }
  def colProduct(x:DenseMatrix[Float],y:DenseMatrix[Float]):DenseMatrix[Float] = {
    val x1 = x *:* y
    sum(x1(::,*)).t.toDenseMatrix
  }
}


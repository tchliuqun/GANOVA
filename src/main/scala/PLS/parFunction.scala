package PLS


//import org.scalameter._
import myParallel._
import scala.reflect.ClassTag

//import scala.collection.parallel

object parFunction {
  var threshold:Int = 5000
  def time[F](f: => F) = {
    val t0 = System.nanoTime
    val ans = f
    printf("Elapsed: %.3f\n",1e-9*(System.nanoTime-t0))
    ans
  }

  def lots[F](n: Int, f: => F): F = if (n <= 1) f else { f; lots(n-1,f) }
  def mapsegSeq[@specialized(Double, Int, Float, Long) A,@specialized(Double, Int, Float, Long) B](inp: Array[A], left: Int, right: Int, f : A => B,
                                                                                                   out: Array[B]) = {
    // Writes to out(i) for left <= i <= right-1
    var i= left
    while (i < right) {
      out(i)= f(inp(i))
      i= i+1
    }
  }
  def mapsegPar[@specialized(Double, Int, Float, Long) A,@specialized(Double, Int, Float, Long) B](inp: Array[A], left: Int, right: Int, f : A => B,
                                                                                                   out: Array[B]): Unit = {
    // Writes to out(i) for left <= i <= right-1
    if (right - left < threshold)
      mapsegSeq(inp, left, right, f, out)
    else {
      val mid = left + (right - left)/2
      parallel(mapsegPar(inp, left, mid, f, out),
        mapsegPar(inp, mid, right, f, out))
    }
  }

  def mapSeq[@specialized(Double, Int, Float, Long) A,@specialized(Double, Int, Float, Long) B: ClassTag](inp:Array[A])(f:A =>B):Array[B] = {
    val right = inp.length
    var out = new Array[B](right)
    mapsegSeq(inp,0,right,f,out)
    out
  }
  def reduceSeg[@specialized(Double, Int, Float, Long) A](inp: Array[A], left: Int, right: Int, f: (A,A) => A): A = {
    if (right - left < threshold) {
      var res= inp(left); var i= left+1
      while (i < right) { res= f(res, inp(i)); i= i+1 }
      res
    } else {
      val mid = left + (right - left)/2
      val (a1,a2) = parallel(reduceSeg(inp, left, mid, f),
        reduceSeg(inp, mid, right, f))
      f(a1,a2) }
  }
  def mapPar[@specialized(Double, Int, Float, Long) A,@specialized(Double, Int, Float, Long) B: ClassTag](inp:Array[A])(f:A =>B):Array[B] = {
    val right = inp.length
    var out = new Array[B](right)
    mapsegPar(inp,0,right,f,out)
    out
  }
  def reducePar[@specialized(Double, Int, Float, Long) A](inp: Array[A])(f: (A,A) => A): A = {
    reduceSeg(inp, 0, inp.length, f)
  }

  def mapReduce[@specialized(Double, Int, Float, Long) A:ClassTag](inp: Array[A])(f: A => A)(fr:(A,A) => A):Array[A] = {
    val right = inp.length
    var res= f(inp(0))
    var out = new Array[A](right)
    out(0) = res
    var i= 1
    while (i < right) {
      res = fr(res, f(inp(i)))
      out(i) = res
      i= i+1
    }
    out
  }
  def fplus:(Float,Float) => Float = (i,j) => i+j
  def mapPlusReduce(inp: Array[Float])(f: Float => Float):Array[Float] = {
    val right = inp.length
    var res= f(inp(0))
    var out = new Array[Float](right)
    out(0) = res
    var i= 1
    while (i < right) {
      res += f(inp(i))
      out(i) = res
      i= i+1
    }
    out
  }
}


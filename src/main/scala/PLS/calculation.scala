package PLS

import java.io.{FileWriter, PrintWriter}

import breeze.linalg._
//import scala.math._
//import scala.math._
import breeze.stats._
import breeze.math._
import breeze.numerics._
import org.apache.commons.math3.stat.inference.TestUtils
import org.apache.commons.math3.distribution._
import org.apache.commons.math3.stat.regression._
import PLS.parFunction._
import scala.collection.mutable.ArrayBuffer
import scala.io.Source.fromFile
import scala.math.Numeric

//import scala.annotation.tailrec
//import org.scalameter._
object calculation {
  def mySum(inp:Array[Float],len:Int):Float = {
    var i = 0
    var res = 0f
    while (i<len){
      res += inp(i)
      i += 1
    }
    res
  }

  def Variance(xs: Array[Float]):Float = {

    val counts = xs.length
    val accu = mySum(xs,counts)

    val mea = accu / counts
    var i = 0
    var variance = 0f
    while ( i< counts){
      variance += (xs(i) - mea)*(xs(i) - mea)
      i += 1
    }
    variance/counts
  }
  def pearsonCorrGeneric[@specialized(Float,Double) T](as: Array[T], bs: Array[T])(implicit n: Numeric[T]): Float = {
    import n._

    val countsA = as.length
    val countsB = bs.length
    if (countsA != countsB)
      sys.error("you fucked up")
    var i = 0
    var accuA = 0f
    var accuB = 0f
    while (i < countsA){
      accuA += as(i).toFloat
      accuB += bs(i).toFloat
      i += 1
    }

    val meaA = accuA / countsA
    val meaB = accuB / countsB
    var variA = 0f
    var variB = 0f
    var ncorrelat = 0f


    while (i< countsA){
      ncorrelat = times(as(i),bs(i)).toFloat + ncorrelat
      variA = (as(i).toFloat - meaA)*(as(i).toFloat - meaA) + variA
      variB = (bs(i).toFloat - meaB)*(bs(i).toFloat - meaB) + variB
      i += 1
    }

    (ncorrelat - countsA.toFloat * meaA * meaB)/scala.math.sqrt(variA * variB).toFloat

  }

  def pearsonCorr(as: Array[Float], bs: Array[Float]): Float = {
    val countsA = as.length
    val countsB = bs.length
    if (countsA != countsB)
      sys.error("you fucked up")
    val accuA = mySum(as,countsA)
    val accuB = mySum(bs,countsA)
    val meaA = accuA / countsA
    val meaB = accuB / countsB
    var i  = 0
    var variA = 0f
    var variB = 0f
    var ncorrelat = 0f

    while (i < countsA){
      ncorrelat += as(i)*bs(i)
      variA += (as(i) - meaA)*(as(i) - meaA)
      variB += (bs(i) - meaB)*(bs(i) - meaB)
      i += 1
    }

    (ncorrelat - countsA * meaA * meaB)/scala.math.sqrt(variA*variB).toFloat

  }

  def filterDuplicateColumn(dm:DenseMatrix[Float]):DenseMatrix[Float] ={
    val covi = breeze.linalg.cov(standardize.scaleColumns(dm).mapValues(_.toDouble))
    val inx:Seq[Int] = covi.findAll(_ > 0.9999).filter(i => i._2 > i._1).map(_._2).distinct
    dm.delete(inx,Axis._1)
  }

  def standardization(X:DenseMatrix[Float]) = {
    val (m,n) = dim(X)
    var rs = DenseMatrix.zeros[Float](m,n)
    var i = 0
    while(i< n){
      val mv = meanAndVariance(X(::,i))
      rs(::,i) := (X(::,i) - mv.mean.toFloat)/sqrt(mv.variance).toFloat
      i += 1
    }
    rs
  }
  def permPvalue[@specialized(Float,Double) T](x:T,dv: Array[T])(implicit n: Numeric[T]): Float = {
    import n._
    if (gteq(x, zero) && dv.count(gteq(_, zero)) > 0) dv.count { i => gteq(i, x) }.toFloat / dv.count(gteq(_, zero)).toFloat
    else if (lteq(x, zero) && dv.count(lteq(_, zero)) > 0) dv.count { i => lteq(i, x) }.toFloat / dv.count(lteq(_, zero)).toFloat
    else 1f
  }

  def permPvalue(dv:DenseVector[Float]):Float = {
    val ts = dv(0)
    val ct = sum(dv.map(i => if (i < ts) 1 else 0))
    ct.toFloat / (dv.size -1)
  }


  def csum(end:Int,start:Int):Int = {
    var msum = 0
    var i = start + 1
    while (i <= end){
      msum = msum + i
      i += 1
    }
    msum
  }

  def cprod(end:BigInt,start:BigInt):BigInt = {
    var prod = BigInt(1)
    var i = start + 1
    while (i <= end){
      prod = prod * i
      i += 1
    }
    prod
  }

  def eigVal(X:DenseMatrix[Float],n:Int = 0) = {
    val m = X.cols
    val X1 = convert(princomp(convert(X,Double)).propvar,Float)
    if (n == 0)X1
    else if(n <= m) X1.apply(0 until n)
    else DenseVector.vertcat(X1,DenseVector.zeros[Float](n-m))
  }

  def runeig(X:DenseMatrix[Float],Y:DenseMatrix[Float],pm:Any) = {
    val n = pm.asInstanceOf[Int]
    val m = X.cols
    val eig = eigVal(X,n)
    m.toString + "\t" + eig.toArray.map(_.toString).mkString("\t")
  }

  def cca(x:DenseMatrix[Float],y:DenseMatrix[Float]):DenseVector[Float] = {
    val yn= y.rows
    val xn= x.rows
    require(xn == yn)
    val qr.QR(xQ, xR) = qr(x)
    val qr.QR(yQ, yR) = qr(y)
    val xk = rank(x.mapValues(_.toDouble))
    val yk = rank(y.mapValues(_.toDouble))
    val e = DenseMatrix.eye[Float](yk)
    val ee = DenseMatrix.zeros[Float](yn-yk, yk)
    val d = DenseMatrix.vertcat(e, ee)
    val d1 = yQ * d
    val res2 = xQ.t * d1
    val res4 = res2(0 until xk, ::)
    val svd.SVD(u, s, v) = svd(res4)
    return(s)
  }

  def HotellingLawleyTrace(rho:DenseVector[Float],p:Int,q:Int):DenseVector[Float] = {
    val minpq = scala.math.min(p,q)
    val rhosq = rho *rho
    val rh = rhosq/(1.0f-rhosq)
    val hotell = DenseVector.zeros[Float](minpq)
    for (i <- 0 until minpq){
      hotell(i) = sum(rh(i until minpq))
    }
    return(hotell)
  }

  def HotellingTest(rho:DenseVector[Float],N:Int,p:Int,q:Int) = {
    val minpq = scala.math.min(p,q)
    val hotell = HotellingLawleyTrace(rho,p,q)
    val Hotelling = DenseVector.zeros[Float](minpq)
    val df1 = DenseVector.zeros[Float](minpq)
    val df2 = DenseVector.zeros[Float](minpq)
    val pVal = DenseVector.zeros[Float](minpq)
    for (i <- 1 to minpq){
      val k = i - 1
      df1(k) = (p-k)*(q-k)
      df2(k) = minpq*(N-2-p-q+2 * k) +2
      Hotelling(k) = df2(k)* hotell(k)/minpq/df1(k)
      val fdist: FDistribution = new FDistribution(df1(k), df2(k))
      pVal(k) = 1.0f - fdist.cumulativeProbability(Hotelling(k)).toFloat
    }
    (hotell,Hotelling,pVal,df1,df2)
  }

  def Brown(pval:Array[Float],cor:Array[Float])={
    val N = pval.length
    val T0 = pval.filter(i=> !i.isNaN & i != 0f).map(log10(_)/log10(Math.E)).sum
    val tc = cor.map(i => 2 * i * (3.25 +0.75 * i)).sum
    val thetaS = 4*N + tc
    val c = thetaS / (4*N)
    val dof = 2 * N/c
    val T = -2 * T0/c
    (T,dof)
  }

  def brownTest(pval:Array[Float],cor:Array[Float]):Float = {
    if (pval.filter(i => i.isNaN | i.isInfinite).length > 0) {
      return 1f
    }else {
      if(pval.filter(_ == 0).length >0){
      for (i <- 0 until pval.length) if (pval(i) == 0f) pval(i) = 1e-16f
     // return 0f
    }
    val (ts,dof) = Brown(pval,cor)
    return (1 - new ChiSquaredDistribution(dof).cumulativeProbability(ts)).toFloat
  }

  }
  def permY(Y: DenseVector[Float],times:Int):DenseMatrix[Float] = {
    val n = Y.length
    //val y = Y.toDenseVector
    val inx = Seq(0 until n: _*)
    val res = DenseMatrix.zeros[Float](n,times+1)
    res(::,0) := Y
    if(times > 0){
      var i:Int = 1
      while(i <= times) {
        res(::,i) := Y(scala.util.Random.shuffle(inx))
        i += 1
      }
    }
    res
  }

  def wilksLambda(rho:DenseVector[Float],N:Int,p:Int,q:Int)= {
    val minpq = scala.math.min(p,q)
    val raoF = DenseVector.zeros[Float](minpq)
    val wilks = DenseVector.zeros[Float](minpq)
    val df1 = DenseVector.zeros[Float](minpq)
    val df2 = DenseVector.zeros[Float](minpq)
    val pval = DenseVector.zeros[Float](minpq)
    val wilksAll = 1.0f - rho *:* rho
    val de = N-1.0f - 0.5f*(p + q+1.0f)
    //var ii = 0
    for (i <- 0 until minpq){
      wilks(i) = wilksAll(i until minpq).reduce(_ * _)
      val k = i.toFloat
      df1(i) = (p-k) *(q-k)
      val nu= if(p-k == 1 | q-k == 1) 1.0f else Math.sqrt((df1(i)*df1(i) -4)/((p-k)*(p-k) + (q-k)*(q-k)-5.0f)).toFloat
      val nu1 = 1.0f/nu
      df2(i) = de * nu - 0.5f * df1(i) + 1.0f
      val w = Math.pow(wilks(i),nu1).toFloat
      raoF(i) = (1.0f-w)/w * df2(i)/df1(i)
      val fdist: FDistribution = new FDistribution(df1(i), df2(i))
      pval(i) = 1.0f - fdist.cumulativeProbability(raoF(i)).toFloat
      //ii += 1
    }
    pval
    //val t = if(p *q == 2.0) 1.0 else Math.sqrt((p*p * q*q - 4.0)/(p*p+q*q - 5.0))
    //val df1 = p*q
    //val df2 = w * t - 0.5*p*q + 1.0
    //val wt = Math.pow(wilks,1.0/t)
    //val f = ((1-wt)/wt) * df2/df1
    //val fdist: FDistribution = new FDistribution(df1, df2)
    //val pVal = 1.0 - fdist.cumulativeProbability(f)
    //pVal
  }


  def chisqTest(a:Int,b:Int,c:Int,d:Int):Float = {
    val inp:Array[Array[Long]] = Array(Array(a.toLong,b.toLong),Array(c.toLong,d.toLong))
    TestUtils.chiSquareTest(inp).toFloat
  }
  def linearPval(x:Array[Double],y:Array[Double]):Float = {
//    val x = if (mcol.length > 0)mcol.map(pd(_).toFloat)else pd.map(_.toFloat)
    if (calculation.Variance(x.map(_.toFloat)) < 0.01f){
      return 1.0f
    }else {
      val n = y.length
      val xy = Array(0 until n: _*).map(i => Array(x(i), y(i)))
      val regression = new SimpleRegression()
      regression.addData(xy)
      return(regression.getSignificance.toFloat)
    }
  }
  def linearPval(pd:Array[String],y:Array[Float],mcol:Array[Int] = Array[Int]()):Float = {
    val x = if (mcol.length > 0)mcol.map(pd(_).toDouble)else pd.map(_.toDouble)
    return linearPval(x,y.map(_.toDouble))
  }
  def pcr(X:DenseMatrix[Float],Y:DenseVector[Float],k:Int) = {
    val rss = princomp(convert(X,Double)).scores
    Array( 0 until k :_*).map(i => calculation.linearPval(rss(::,i).toArray,Y.toArray.map(_.toDouble)))
  }
//  def filterDuplicateColumn(dm:DenseMatrix[Double]):DenseMatrix[Double] ={
//    val covi = breeze.linalg.cov(Standardize.scaleColumns(dm))
//    val inx:Seq[Int] = covi.findAll(_ > 0.9999).filter(i => i._2 > i._1).map(_._2).distinct
//    dm.delete(inx,Axis._1)
//  }
  def sgsea(gs:DenseMatrix[Float],rs:DenseVector[Float]) = {
    val psum = gs.rows
    val gcol = gs.cols
    val asum = sum(gs(::,*)).t
    val nsum = psum.toFloat - asum
    val rst = DenseVector.zeros[Float](gcol)
    val nrst = DenseVector.zeros[Float](gcol)
    val rss = breeze.numerics.abs(rs)
    val gss = gs(::,*) *:* rss
    val ngs = 1f - gs
    val posd = sum(gss(::,*)).t
    val res =  gss(*,::) /:/ posd - ngs(*,::) /:/ nsum
    val results = accumulate(res(::,*))
    val mins = min(results(::,*))
    val maxs = max(results(::,*))
    for ( i <- 0 until gcol) {
      val es = if (abs(mins(i)) > maxs(i)) mins(i) else maxs(i)
      rst(i) = es
      nrst(i) = nes(es,results(::,i).toArray)
    }
    (rst,nrst,results)
  }
  def gseaP(gs:DenseMatrix[Float],rs:DenseVector[Float],perm:Int = 1000,cores:Int = 1) = {
    val rsp = calculation.permY(rs,perm)
    val rsd = DenseMatrix.zeros[Float](perm,gs.cols)
    var i = 0
//    while(i< perm){
//      rsd(i,::) := calculation.sgsea(gs,rsp(::,i).toDenseVector)._2.t
//      i += 1
//    }
    parFunction.threshold = perm/cores
    parFunction.mapPar(Array(1 to perm :_*)){ i =>
//      rsd(i,::) := DenseVector(Array.fill(gs.cols)(0.1f)).t;0
     rsd(i,::) := calculation.sgsea(gs,rsp(::,i).toDenseVector)._2.t;0
    }
//    Array(1 to perm :_*).par.foreach{ i =>
//      rsd(i,::) := DenseVector(Array.fill(gs.cols)(0.1f)).t
////     rsd(i,::) := calculation.sgsea(gs,rsp(::,i).toDenseVector)._2.t;0
//    }
    val rss = calculation.sgsea(gs,rs)
    val pv = DenseVector.fill(gs.cols,0f)
    i = 0
    while (i< gs.cols){
      val rsv = rss._2.apply(i)
      if (rsv < 0){
        val sm = rsd(::,i) <:< rsv //._2.apply(i)
        pv(i) = sum(sm.map(i => if (i)1f else 0f))/perm
      } else{
        val sm = rsd(::,i) >:> rsv //._2.apply(i)
        pv(i) = sum(sm.map(i => if (i)1f else 0f))/perm
      }
      i += 1
    }
    (rss._1,rss._2,pv)
  }


  def gsea(gss: Array[Float], rss: Array[Float]):(Float,Array[Float]) = {
    //    val rs = rss.map(abs(_))
    val len = rss.length
    var rs = new Array[Float](len)
    var posd = 0f
    var i:Int = 0
    var ngss = 0f
    var gsssum = 0f
    while( i < len) {
      val temp = if (rss(i) < 0) -rss(i) else rss(i)
      rs(i) = temp
      posd = posd+ gss(i) * temp
      gsssum = gsssum + gss(i)
      ngss = ngss + 1f
      i += 1
    }
    val negd = ngss - gsssum

    var j = 0
    var res = new Array[Float](len)
    var minresult = 0f
    var maxresult = 0f
    var result = 0f
    while( j < len)  {
      result = result + gss(j) * rs(j) /posd  - (1f - gss(j))/negd
      res(j) = result
      minresult = if (minresult < result) minresult else result
      maxresult = if (maxresult > result) maxresult else result
      j += 1
    }
    val reslt = if(scala.math.abs(minresult) > maxresult) minresult else maxresult
    (reslt,res)
  }
  def gseaPerm(gss: Array[Float], rss: Array[Float]):Float = {
    //    import n._
    //val rs = rss.map(abs(_))
    val len = rss.length
    var rs = new Array[Float](len)
    var posd = 0f
    var i:Int = 0
    var ngss = 0f
    var gsssum = 0f
    while( i < len) {
      val temp = if (rss(i) < 0) -rss(i) else rss(i)

      rs(i) = temp
      posd = posd+ gss(i) * temp
      gsssum = gsssum + gss(i)
      ngss = ngss + 1f
      i += 1

    }
    val negd = ngss - gsssum

    //@tailrec def gseainner(gss:Array[Float],rss:Array[Float],posdiv:Float,negdiv:Float,sigresult:Float,minresult: Float,maxresult :Float,result: List[Float]):(Float,Float,List[Float]) ={
    i = 0
    var minresult = 0f
    var maxresult = 0f
    var result = 0f
    while( i < len)  {
      result = result + gss(i) * rs(i) /posd  - (1f - gss(i))/negd
      minresult = if (minresult < result) minresult else result
      maxresult = if (maxresult > result) maxresult else result
      i += 1
    }
    val reslt = if(scala.math.abs(minresult) > maxresult) minresult else maxresult
    reslt
  }
  def nes[@specialized(Float) T](x:T,dv: Array[T])(implicit n: Numeric[T]): Float = {
    import n._
    //    val dvar = dv.toArray
    if (gteq(x, zero)) x.toFloat * dv.count(gteq(_, zero)).toFloat / dv.filter(gteq(_, zero)).sum.toFloat
    else
      x.toFloat * dv.count(lteq(_, zero)).toFloat / dv.filter(lteq(_, zero)).sum.toFloat
  }

//  def nes(x:Float,dv: DenseVector[Float]): Float = {
//    //import n._
//    //    val dvar = dv.toArray
//    num = dv.size
//    if (x > 0f) x * dv.filter(_ > 0).length.toFloat / dv.filter(gteq(_, zero)).sum.toFloat
//    else
//      x.toFloat * dv.count(lteq(_, zero)).toFloat / dv.filter(lteq(_, zero)).sum.toFloat
//  }
  def linearPvall(pf:String = gPms.op+gPms.pf,df:String =gPms.op+gPms.df,af:String =gPms.op+gPms.af,of:String = gPms.op+"snp6pval.txt") = {
    val x = scala.io.Source.fromFile(df).getLines.map(_.split(",")).take(1).toArray.flatten
    val y = scala.io.Source.fromFile(pf).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val mcol = fileOper.intersectCols(x,y)
    val writer = new PrintWriter(new FileWriter(of))

    val pakt = fileOper.toArrays(pf).filter(_ (0).contains("AKT")).toArray.apply(1)
    val yy = mcol._2.map(pakt(_).toFloat)
    val xx = scala.io.Source.fromFile(df).getLines.map(_.split(",")).drop(1)
    val pval = xx.map(i => i(0) -> linearPval(i,yy,mcol._1)).toMap
    scala.io.Source.fromFile(af).getLines.map(i => i+"\t"+pval.getOrElse(i.split("\t",2).apply(0),1.0)).foreach(writer.println)
    writer.close()
  }
  def corVal(xc:Iterator[Array[String]] = fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")),mcoll:Array[Int] = null,num:Int = 20)={
    //val xc = scala.io.Source.fromFile(df).getLines.map(_.split("\t"))
    var snp =scala.collection.mutable.ArrayBuffer.empty[Array[Float]]
    var rs = scala.collection.mutable.ArrayBuffer.empty[Array[String]]
    var nam = Array[String]()
    var i = 0
    var ii = 0
    while( i <= num){
      val xcc = xc.next
      nam :+= xcc(0)
      val xxc = if(mcoll != null) mcoll.map(xcc.apply(_).toFloat) else xcc.drop(1).map(_.toFloat)
      snp += xxc
      i += 1
    }
    while (ii < i) {
      val corr = nam(0) +: snp.map(i => calculation.pearsonCorr(snp(0), i)).drop(1).toArray.map(_.toString)
      rs += corr
      //writer.println(corr)
      ii += 1

      if (xc.hasNext) {
        i += 1
        val xcc= xc.next
        nam :+= xcc(0)
        val xxc = if(mcoll != null) mcoll.map(xcc.apply(_).toFloat) else xcc.drop(1).map(_.toFloat)
        snp += xxc
      }
      snp = snp.drop(1)
      nam = nam.drop(1)
    }
    rs.toArray
    //writer.close
  }
  def fisherMatrix(gs:DenseMatrix[Float],anno:Array[Array[String]],num:Int) = {
    val t4 = gs.rows
    val m = gs.cols
    val tgs = gs(0 until num,::)
    val t1 = sum(tgs(::,*)).t.toArray.map(_.toInt)
    val t2 = sum(gs(::,*)).t.toArray.map(_.toInt)
    val no = anno(0).length - 1
    //Array(0 until m :_*).map(i => anno(i) +: calculation.fisherTest(t1(i),t2(i),num,t4).toString).sortBy(_(no).toFloat)
    Array(0 until m :_*).map(i => (anno(i),t1(i),t2(i), calculation.fisherTest(t1(i),t2(i),num,t4))).sortBy(_._4)
  }

  def BenjaminiHochbergFDR(p:Array[Float]):Array[Float] = {
    val po = p.sorted
    val n = p.length
    val rs = Array.fill(n)(0f)
    rs(n-1) = po(n-1)
    var i = n-2
    var un = po(n-1)
    var left = 0f
    var right = 0f
    var bm = i+1
    while( i > -1){
      un = po(i)
      bm = i+1
      left = rs(bm)
      right = (n/bm)*un
      rs(i) = math.min(left,right)
      i -= 1
    }
    return(rs)
  }

  def ogsea(gs:DenseVector[Float],rs:DenseMatrix[Float]) = {
    val psum = rs.rows
    val gcol = rs.cols
    val asum = sum(gs)
    val nsum = psum.toFloat - asum
    val rst = Array.fill(gcol)(0f)
    val rss = breeze.numerics.abs(rs)
    val gss = rss(::,*) *:* gs
    val ngs = 1f - gs
    val posd = sum(gss(::,*)).t
    val res1 =  gss(*,::) /:/ posd
    val res = res1(::,*) - ngs / nsum

    var ii = 0
    while (ii < gcol){
      val res2 = res(::,ii).toArray
      val rss2 = rs(::,ii).toArray.zipWithIndex.sortBy(_._1).reverse.map(i => res2(i._2)).scanLeft(0f)(_+_).drop(1)
      val mins = min(rss2)
      val maxs = max(rss2)
      rst(ii) = if (abs(mins) > maxs) mins else maxs
      ii += 1
    }
    val results = accumulate(res(::,0))
    //val mins = min(results(::,*))
    //val maxs = max(results(::,*))
    //for ( i <- 0 until gcol) {
    //  rst(i) = if (abs(mins(i)) > maxs(i)) mins(i) else maxs(i)
    //}
    val es = rst(0)
    val nes = if(es > 0) {
      val ls = rst.drop(1).filter(! _.isNaN).filter(_ > 0)
      ls.length.toFloat * es/ls.sum
    }else {
      val ls = rst.drop(1).filter(! _.isNaN).filter(_ < 0)
      ls.length.toFloat * es/ls.sum
    }
    val pval = if(es>0) rst.filter(_ > es).length.toFloat / (gcol-1).toFloat else rst.filter(_ < es).length.toFloat / (gcol-1).toFloat
    es +: nes +: pval +: asum.toFloat +: results.toArray
  }

  def chisqTestwithCorrect(a1:Int,b1:Int,c1:Int,d1:Int):Float = {
    val a = BigDecimal(a1)
    val c = BigDecimal(c1)
    val b = BigDecimal(b1)
    val d = BigDecimal(d1)
    val n = (a+b+c+d).toDouble
    val nom = BigDecimal(math.pow(math.max(math.abs(a1*d1-b1*c1)-n/2,0),2)) * BigDecimal(n)
    val denom = if ((a+b)*(a+c)*(b+d)*(c+d) == BigDecimal(0)) BigDecimal(0.0000001) else (a+b)*(a+c)*(b+d)*(c+d)
    val chisq = (nom / denom).toFloat
    chisqPvalue(chisq,1)
  }
  def chisqPvalue(chisq:Float,k:Int):Float = {
    val pval = new ChiSquaredDistribution(k.toDouble).cumulativeProbability(chisq).toFloat
    1-pval
  }

  def reverseChisq(pval:Array[Float],dof:Int = 1):Array[Float] = {
    val pd = new ChiSquaredDistribution(dof.toDouble)//.cumulativeProbability(chisq).toFloat
    var i = 0
    val n = pval.length
    val rs = Array.fill(n)(0f)
    while (i< n){
      rs(i) = pd.inverseCumulativeProbability(1f-pval(i)).toFloat
      i += 1
    }
    return(rs)
  }

  def fisherValue(a:Int,b:Int,c:Int,d:Int):Float = {
    val Array(af,bf,cf,df) = Array(a,b,c,d).sorted.map(BigInt(_))
    val ab = a+b
    val ac = a+c
    val bd = b+d
    val cd = c+d
    val n = ab+cd
    val Array(abf,acf,bdf,cdf) = Array(ab,ac,bd,cd).sorted

    val logab = cprod(af,0)
    val logca = cprod(abf,bf)
    val logdc = cprod(acf,cf)
    val logbd = cprod(bdf,df)
    val logn = cprod(n,cdf)

    (BigDecimal(logdc*logca*logbd) / BigDecimal(logn * logab)).toFloat
  }

  def fisherTest(a:Int,b:Int,c:Int,d:Int):Float = {
    //val row = a+b
    //val col = a+c
    //val n = a+b+c+d
    val dis = new HypergeometricDistribution(d,b,c)
    dis.upperCumulativeProbability(a).toFloat
  }
  def fisherTest(x:Seq[Int]):Float = {
    if (x.length < 4)println("The collection should have at least four elements")
    fisherTest(x(0),x(1),x(2),x(3))
  }
  def getkappa[T](x1:Set[T],x2:Set[T],n:Int) = {
    val n1 = x1.intersect(x2).size
    val n2 = x1.size - n1
    val n3 = x2.size - n1
    val n4 = n - n1 - n2 - n3
    kappaTest(n1,n2,n3,n4)
  }
  def getkappa[@specialized(Float,Int) T](x1:DenseVector[T],x2:DenseVector[T])(implicit num: Numeric[T]) = {

    import num._
    val x11 = x1.map(_.toInt)
    val x22 = x2.map(_.toInt)
    val n = x11.length
    val n1 = x11.t * x22
    val n2 = sum(x11) - n1
    val n3 = sum(x22) - n1
//    val n2 = x1.reduce(plus(_,_)) - n1
//    val n3 = x2.reduce(plus(_,_)) - n1
    val n4 = n-n1-n2-n3
    kappaTest(n1,n2,n3,n4)
  }

  def getkappa(x1:DenseVector[Float],x2:DenseVector[Float]) = {
    val n = x1.length
    val n1 = (x1.t * x2).toInt
    val n2 = sum(x1).toInt - n1
    val n3 = sum(x2).toInt - n1
    val n4 = n-n1-n2-n3
    kappaTest(n1,n2,n3,n4)
  }


  def kappaMap[T](x:Array[(String,Set[T])],g0 : (String,Set[T]),n:Int = -1,thresh:Float = 0.1f)= {
    val num = if(n<0)x.flatMap(_._2).toSet.size else n
    val len = x.length
    val rs1 = x.map(i => (g0._1,i._1,getkappa(g0._2,i._2,num))).filter(_._3 > thresh)
  }
  def kappaArray[T](x:Array[(String,Set[T])],thresh:Float = 0.1f) = {
    val n = x.flatMap(_._2).toSet.size
    val len = x.length
    var i = 0
    var rs = ArrayBuffer[(String,String,Float)]()
    while(i < len){
      val g0 = x(i)
      var j = i+1
      while( j < len){
      val g1 = x(j)
        val ka = getkappa(g0._2,g1._2,n)
        if(ka >= thresh ) {
          val rs1 = (g0._1, g1._1, getkappa(g0._2, g1._2, n))
          rs += rs1
        }
        j += 1
      }
      i += 1
    }
    rs.toArray
  }
  def kappaTest(a:Int,b:Int,c:Int,d:Int):Float = {
    val Array(af,bf,cf,df) = Array(a,b,c,d).map(BigInt(_))
    val ab = af+bf
    val ac = af+cf
    val bd = bf+df
    val cd = cf+df
    val ad = af+df
    val n = ab+cd
    val po = ad.toFloat/n.toFloat
    val pe = (ab*ac+bd*cd).toFloat/(n*n).toFloat
    val kappa = (po-pe)/(1-pe)
    return(kappa)
  }
  def kappaTest(x:Seq[Int]):Float = {
    if (x.length < 4)println("The collection should have at least four elements")
    kappaTest(x(0),x(1),x(2),x(3))
  }
}

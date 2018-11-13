package PLS

import breeze.linalg._
import breeze.linalg.eigSym.EigSym
import breeze.numerics._
import breeze.util._
//import breeze.stats.regression.leastSquares
import org.apache.commons.math3.distribution.FDistribution

object plsCalc {
  //val today = java.time.LocalDate.now.toString.split("-")
  //val totime = java.time.LocalTime.now.toString.split("\\.")(0).split(":")
  //val nowtime = (today ++ totime).mkString("_")

  //  def PlsSnp(snpData: String, snpAnno: String, Pheno: Array[Float]) = {
  //    ???
  //  }
  //
  //  def splitSnpFile(snpData: String, snpAnno: String) = {
  //    ???
  //  }
  //
  //  def getGeneList = {
  //    ???
  //  }
  //
  //  def getSnpInGene(gene: String, snpAnno: Array[Array[String]]) = {
  //    ???
  //  }

  def pls(X: DenseMatrix[Float], Y: DenseMatrix[Float], A: Int): DenseMatrix[Float] = {
    val K = X.cols // number of columns in X - number of predictor variables
    val M = Y.cols // number of columns in Y - number of response variables
    val N = X.rows // number of rows in X === number of rows in Y - number of observations of predictor and response variables
    // val W: DenseMatrix[Float] = DenseMatrix.zeros[Float](K, A) // (K x A)
    val P: DenseMatrix[Float] = DenseMatrix.zeros[Float](K, A) // (K x A)
    val Q: DenseMatrix[Float] = DenseMatrix.zeros[Float](M, A) // (M x A)
    val R: DenseMatrix[Float] = DenseMatrix.zeros[Float](K, A) // (K x A)
    var YtXXtY: DenseMatrix[Float] = DenseMatrix.zeros[Float](M, M) // (M x M) matrix
    var w_a: DenseVector[Float] = DenseVector.zeros[Float](K) // (K x 1) - column vector of W
    var p_a: DenseVector[Float] = DenseVector.zeros[Float](K) // (K x 1) - column vector of P
    var q_a: DenseVector[Float] = DenseVector.zeros[Float](M) // (M x 1) - column vector of Q
    var r_a: DenseVector[Float] = DenseVector.zeros[Float](K) // (K x 1) - column vector of R
    var t_a: DenseVector[Float] = DenseVector.zeros[Float](N) // (N x 1) - column vector of T [which is N x A]
    var tt: Float = 0.0f
    var indexOfLargestEigenvalue: Int = 0
    var beta = DenseMatrix.zeros[Float](K,M*A)

    var XY: DenseMatrix[Float] = X.t * Y // compute the covariance matrices; (K x M) matrix


    var a :Int = 0
    // A = number of PLS components to compute
    while (a < A) {
      if (M==1) {
        w_a = XY.toDenseVector
      }else {
        // cforRange(0 until reps) { i =>
        //dataSet(::, *) - input
        //val xy = XY(::,c(1,i))
        //val ytxxty = xy.t *xy
        //val EigSym(eigenvalues, eigenvectors) = eigSym(convert(ytxxty,Double))          // eigenvalues is a DenseVector[Double] and eigenvectors is a DenseMatrix[Double]
        //indexOfLargestEigenvalue = argmax(eigenvalues)                  // find index of largest eigenvalue
        //q_a = eigenvectors(::, indexOfLargestEigenvalue).mapValues(_.toFloat)
        //w_a(::,i) = xy * q_a
        //}
        YtXXtY = XY.t * XY                                              // XY.t * XY is an (M x M) matrix
        val EigSym(eigenvalues, eigenvectors) = eigSym(convert(YtXXtY,Double))          // eigenvalues is a DenseVector[Double] and eigenvectors is a DenseMatrix[Double]
        indexOfLargestEigenvalue = argmax(eigenvalues)                  // find index of largest eigenvalue
        q_a = eigenvectors(::, indexOfLargestEigenvalue).mapValues(_.toFloat)                // find the eigenvector corresponding to the largest eigenvalue; eigenvector is (M x 1)
        w_a = XY * q_a
      }// in this case, XY (K x M) is a column vector (K x 1) since M == 1
      //w_a = XY(::,*).map(x => x/sqrt(x.t * x))
      w_a = w_a /sqrt(sum(w_a *:* w_a))// sqrt(w_a.toDenseMatrix * w_a.toDenseMatrix.t).apply(0,0)//.toFloat // normalize w_a to unity - the denominator is a scalar
      r_a = w_a // loop to compute r_a
      var j: Int = 0
      while (j < a ) {
        r_a = r_a - ((P(::, j).t * w_a) * R(::, j)) // (K x 1) - ( (1 x K) * (K x 1) ) * (K x 1) === (K x 1) - scalar * (K x 1)
        j += 1
      }
      t_a = X * r_a // compute score vector
      tt = sum(t_a *:* t_a)// t_a.t * t_a // compute t't - (1 x 1) that is auto-converted to a scalar
      p_a = (X.t * t_a) / tt // X-loadings - ((K x 1)' * (K x K))' / tt === (K x 1) / tt
      q_a = (r_a.t * XY).t / tt // Y-loadings - ((K x 1)' * (K x M))' / tt === (M x 1) / tt
      XY = XY - ((p_a * q_a.t) * tt) // XtY deflation

      // update loadings and weights
      // W(::, a) := w_a
      P(::, a) := p_a
      Q(::, a) := q_a
      R(::, a) := r_a
      beta(::,a*M until (a+1)*M) := R(::,0 to a) * Q(::,0 to a).t//.toDenseVector
      a += 1
      //  beta(::, a) := R(::, 0 until a) * Q(::,0 until a).t
    }
    //val beta = R * Q.t
    // compute the regression coefficients; (K x M)
    beta
  }
  def plsP(X: DenseMatrix[Float], Y: DenseMatrix[Float], A: Int,nPerm:Int = 1000) = {
    val K = X.cols // number of columns in X - number of predictor variables
    val N = X.rows // number of rows in X === number of rows in Y - number of observations of predictor and response variables
    val L = Y.cols
    val M = L/nPerm//Y.cols // number of columns in Y - number of response variables
    var P: Array[DenseMatrix[Float]] = Array.fill(A)(DenseMatrix.zeros[Float](K, nPerm)) // (K x A)
    //val Q: DenseMatrix[Float] = DenseMatrix.zeros[Float](M, A) // (M x A)
    var Yhat: DenseMatrix[Float] = DenseMatrix.zeros[Float](N, L) // (M x A)
    var beta: DenseMatrix[Float] = DenseMatrix.zeros[Float](K, L)
    var R:Array[DenseMatrix[Float]] = Array.fill(A)(DenseMatrix.zeros[Float](K, nPerm)) // (K x A)
    var w_a: DenseMatrix[Float] = DenseMatrix.zeros[Float](K,nPerm) // (K x 1) - column vector of W
    var p_a: DenseMatrix[Float] = DenseMatrix.zeros[Float](K,nPerm) // (K x 1) - column vector of P
    var q_a: DenseVector[Float] = DenseVector.zeros[Float](L) // (M x 1) - column vector of Q
    var r_a: DenseMatrix[Float] = DenseMatrix.zeros[Float](K,nPerm) // (K x 1) - column vector of R
    var t_a: DenseMatrix[Float] = DenseMatrix.zeros[Float](N,nPerm) // (N x 1) - column vector of T [which is N x A]
    var tt: DenseVector[Float] = DenseVector.zeros[Float](nPerm)
    var XY: DenseMatrix[Float] = X.t * Y // compute the covariance matrices; (K x M) matrix
    var a :Int = 0
    // A = number of PLS components to compute
    while (a < A) {
      // w_a = XY
      if (M == 1) {
        w_a = XY(::, *).map(x => x / sqrt(sum(x *:* x)))

      } else {
        var ii = 0
        while (ii < nPerm){
          var xy = XY(::,ii*M until (ii+1)* M)
          var YtXXtY = xy.t * xy                                              // XY.t * XY is an (M x M) matrix
          var EigSym(eigenvalues, eigenvectors) = eigSym(convert(YtXXtY,Double))          // eigenvalues is a DenseVector[Double] and eigenvectors is a DenseMatrix[Double]
          var indexOfLargestEigenvalue = argmax(eigenvalues)                  // find index of largest eigenvalue
          q_a = eigenvectors(::, indexOfLargestEigenvalue).mapValues(_.toFloat)                // find the eigenvector corresponding to the largest eigenvalue; eigenvector is (M x 1)
          w_a(::,ii) := xy * q_a
          ii += 1
        }
        w_a = w_a(::, *).map(x => x / sqrt(sum(x *:* x)))
      }// loop to compute r_a

      r_a = w_a
      var j: Int = 0
      while (j < a ) {
        val md = sum(P(j) *:* w_a,Axis._0).t
        r_a = r_a - (R(j)(*,::) * md) // (K x 1) - ( (1 x K) * (K x 1) ) * (K x 1) === (K x 1) - scalar * (K x 1)
        j += 1
      }
      t_a = X * r_a // compute score vector
      tt = t_a(::,*).map(x => sum(x *:* x)).t// compute t't - (1 x 1) that is auto-converted to a scalar
      val p_aa = X.t * t_a
      p_a = p_aa(*,::) / tt // X-loadings - ((K x 1)' * (K x K))' / tt === (K x 1) / tt
      //val q_aa = sum(r_a *:* XY,Axis._0).t
      var jj = 0
      while(jj< M){
        val q_aa = sum(t_a *:* Y(::,jj until L by M),Axis._0).t
        q_a = q_aa /:/ tt
        val xya =p_a(*,::) *:* (tt *:* q_a)  // Y-loadings - ((K x 1)' * (K x M))' / tt === (M x 1) / tt
//        if (jj == M -1) {
          Yhat(::,jj until L by M) += t_a(*,::) *:* q_a
          beta(::,jj until L by M) += r_a(*,::) *:* q_a
//        }
        XY(::,jj until L by M) :-= xya // XtY deflation
        jj += 1
      }
      P(a) = p_a.copy
      R(a) = r_a.copy
      //Q(::, a) := q_a
      //val Yhata = t_a(*,::) *:* q_a
      //beta(::,a*L until (a+1)*L) := R(::,0 to a) * Q(::,0 to a).t//.toDenseVector
      a += 1
      //beta(::, a) := R(::, 0 until a) * Q(::,0 until a).t
    }
    (Yhat,beta)
  }


  def permY(Y: DenseMatrix[Float],times:Int):DenseMatrix[Float] = {
    //val n = Y.rows
    val y = Y.toDenseVector
    calculation.permY(y,times)
//    val inx = Seq(0 until n: _*)
//    val res = DenseMatrix.zeros[Float](n,times+1)
//    res(::,0) := y
//    if(times > 0){
//      var i:Int = 1
//      while(i < times) {
//        res(::,i) := y(scala.util.Random.shuffle(inx))
//        i += 1
//      }
//    }
//    res
  }


  def permT(X:DenseMatrix[Float],permY:DenseMatrix[Float],Yot:DenseMatrix[Float],k:Int = 3) = {
    var i = 0
    val perm = permY.cols
    var rss = new DenseMatrix[Float](k,perm+1)
    while (i < perm ) {
      val y1 = permY(::, i).toDenseMatrix.t
      val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, y1, k))
      rss(::,i) := plsCalc.rss(y1, Yhat)
      i += 1
    }
    val st = rss(::,0).toArray
    Array(0 until k :_*).map(i => rss(i,::).t.toArray.count(ii => ii < st(i)).toFloat / perm.toFloat)
  }
  case class Gene(entrez:String,chr:String,TSstart:Int,TSend:Int)
  case class Snip(probe:String = "None",id:String = "None",chr:String,loc:Int)
  def SnpInRange(gnfl:String,snpfl:String) = {
    val xchr = Array(1 to 22: _*).map(_.toString).toSet
    val sp = scala.io.Source.fromFile(snpfl).getLines().drop(19).
      map(_.split("\" \"").map(_.replace("\"", "")).take(4)).
      filter(i => xchr.contains(i(2)) && i(3) != "---").toArray.sortBy(i => (i(2).toInt, i(3).toInt))
    //.map(i => new Snip(i(0),i(1),i(2),i(3).toInt))
    val gn = scala.io.Source.fromFile(gnfl).getLines().drop(1).map(_.split(",").take(5)).toArray
    val snl = sp.length
    var i = 0
    var rs = Array[(String, Array[String])]()
    while (i < 23) {
      val gen = gn.filter(x => x(1) == i.toString)
      val snpls = sp.filter(x => x(2) == i.toString)
      var j = 0
      val gnl = gen.length
      while (j < gnl) {
        val gnam = gen(j)(0)
        val start = gen(j)(2).toInt
        val end = gen(j)(3).toInt
        var sls = snpls.filter(s => s(3).toInt > start & s(3).toInt < end).map(_.apply(0))
        rs = rs :+ (gnam, sls)
        j += 1
      }
      i+=1
    }
    rs
  }
  //val gnfl1 = "geneSnp6Rs.csv"
  //val snpfl1 = "snp6Annot.csv"
  //val x = SnpInRange(gnfl1,snpfl1)
  def anova(Y: DenseMatrix[Float], Ypred: DenseMatrix[Float], dof: DenseMatrix[Float]):DenseVector[Float] = {
    require(Ypred.cols == dof.cols)
    //val ny = Y.rows.toFloat
    //val YY = if (Y.cols != Ypred.cols)utils.horzcat(Y,Ypred.cols) else Y
    //val ybar = YY - Ypred
    //val nmy = utils.colProduct(ybar,ybar)
    //val my = sum(Y)/ny
    //val y2 = Ypred - my
    //val y21 = utils.colProduct(y2, y2)
    val (y21,nmy,ny) = fval(Y,Ypred)
    fdpval(y21,nmy,ny,dof)
//    val nm = ny -1.0f - dof
//    val fval = (y21 / dof) / (nmy / nm)
//    val pval = new DenseVector[Float](dof.cols)
//    for (i <- 0 until dof.cols){
//      val fdist: FDistribution = new FDistribution(dof(0,i).abs, nm(0,i).abs)
//      val pv = (1.0d - fdist.cumulativeProbability(fval(0,i).toDouble)).toFloat
//      pval(i) = if (pv == 0f) 1e-16f else pv
//    }
//
//    return(pval)
  }
  def fval(Y: DenseMatrix[Float], Ypred: DenseMatrix[Float]) = {
    val ny = Y.rows.toFloat
    val YY = if (Y.cols != Ypred.cols)utils.horzcat(Y,Ypred.cols) else Y
    val ybar = YY - Ypred
    val nmy = utils.colProduct(ybar,ybar)
    val my = sum(Y)/ny
    val y2 = Ypred - my
    val y21 = utils.colProduct(y2, y2)
    (y21,nmy,ny)
  }
  def fdpval(y21:DenseMatrix[Float],nmy:DenseMatrix[Float],ny:Float,dof:DenseMatrix[Float]):DenseVector[Float]  = {
    val nm = ny -1.0f - dof
    val fval = (y21 / dof) / (nmy / nm)
    val pval = new DenseVector[Float](dof.cols)
    for (i <- 0 until dof.cols){
      val fdist: FDistribution = new FDistribution(dof(0,i).abs, nm(0,i).abs)
      val pv = (1.0d - fdist.cumulativeProbability(fval(0,i).toDouble)).toFloat
      pval(i) = if (pv == 0f) 1e-16f else pv
    }
    return(pval)
  }
  def dofPval(Ys:DenseMatrix[Float],Ypred:DenseMatrix[Float],dof:Array[Float]) = {
    val k = dof.length
    //var rs = Array.fill(k)(0d)
    //for (i <- 0 until k){
    //  rs(i) =
    //}
    Array(0 until k:_*).map(i => plsCalc.anova(Ys(::,k % Ys.cols).toDenseMatrix.t,Ypred(::,i).toDenseMatrix.t,new DenseMatrix[Float](1,1,Array(dof(i))))(0))
  }
  def twoAnova(Y: DenseMatrix[Float],Y1: DenseMatrix[Float], Y2: DenseMatrix[Float],dof1: DenseMatrix[Float] ,dof2: DenseMatrix[Float]):DenseVector[Float] = {

    require(Y1.cols == dof1.cols)
    require(Y.cols == Y1.cols)
    require(dof2.cols == dof1.cols)
    require(Y1.cols == Y2.cols)
    val ny = Y.rows.toFloat
    val y1 = Y -Y1
    val y11 = utils.colProduct(y1,y1)
    val y2 = Y -Y2
    val y22 = utils.colProduct(y2,y2)
    val yrss = y11 - y22
    val dofd = dof2 -dof1
    val dofs = sum(dofd)
    val doff = if(dofs < 0) ny - 1.0f - dof1 else ny - 1.0f - dof2
    val rss = if(dofs < 0) y11 else y22
    val fval = (yrss / dofd) / (rss / doff)
    val pval = new DenseVector[Float](dof1.cols)
    for (i <- 0 until dof1.cols){
      val fdist: FDistribution = new FDistribution(abs(dofd(0,i)), doff(0,i))
      pval(i) = 1.0f - fdist.cumulativeProbability(fval(0,i)).toFloat
    }
    return(pval)
  }
  def kfoldInx(size:Int,count: Int,random:Boolean = false): Array[Seq[Int]] = {
    import scala.util.Random
    val piece = size / count
    val indx = if(random) Random.shuffle(Seq(0 until size :_*)) else Seq(0 until size :_*)
    Array(0 until count:_*).map{i => indx.slice(i*piece, (i+1) * piece)}
  }
  def plsCV(X: DenseMatrix[Float], Y: DenseMatrix[Float],n:Int, pn: Array[Seq[Int]]):DenseMatrix[Float] = {
    val nx = X.rows
    val mx = X.cols
    val my = Y.cols
    val lp = pn.length
    //val pn = kfoldInx(nx,k,true)
    var rs = DenseMatrix.zeros[Float](nx,n*my)
    var i : Int = 0
    while (i < lp) {
      val pni = pn(i)
      val Xtrain = X.delete(pni, Axis._0)
      val Ytrain = Y.delete(pni,Axis._0)
      val coef = plsTrain(Xtrain,Ytrain,n)
      rs(pni,::) := predict(X(pni,::).toDenseMatrix,coef)
      i += 1
    }
    return(rs)
  }

  def plsLoocv() = {
    ???
  }

  def pdof(Y:DenseMatrix[Float],Yfit:DenseMatrix[Float],Ycv:DenseMatrix[Float]):DenseVector[Float] = {
    val n = Y.rows
    val k = Yfit.cols
    val y = Y.toDenseVector
    val fit = Yfit(::,*) - y
    val cv = Ycv(::,*) - y
    var ssfit = DenseVector.zeros[Float](k)//if (k>1) diag(fit.t * fit) else DenseVector(fit.t * fit)//
    var sscv = DenseVector.zeros[Float](k)//if (k>1) diag(cv.t * cv) else DenseVector(cv.t * cv)// DenseVector.zeros[Float](k)
    for ( i <- 0 until k){
      ssfit(i) = fit(::,i).t * fit(::,i)
      sscv(i) = cv(::,i).t * cv(::,i)
    }
    val vl = sqrt(ssfit / sscv)
    //val rs =
    n.toFloat * (1.0f - vl)
    //return(rs)
  }

  def gdof(X:DenseMatrix[Float],Y:DenseMatrix[Float],theta:Float,k:Int =1,nPerm:Int = 1000,ncor:Int = 1):Float = {
    val g = breeze.stats.distributions.Gaussian(0, theta)
    val n = Y.rows
    val ron = new DenseMatrix[Float](n,nPerm,g.sample(nPerm*n).toArray.map(_.toFloat))
    val est = DenseMatrix.zeros[Float](n,nPerm)
    //val Yest = Array(1 to k :_*).map(_ => est)
    val h = DenseVector.zeros[Float](n)
    val mlen = 2 * nPerm / ncor
    def getEsts(ind:Array[Int]):Unit = {
      val nlen = ind.length
      if (nlen < mlen){
        getEst(ind)
      }else{
        val sep = nlen /2
        myParallel.parallel(getEsts(ind.slice(0,sep)),getEsts(ind.slice(sep,nlen)))
      }
    }
    def getEst(ind:Array[Int]) = {
      var ii:Int = 0
      val leng = ind.length
      while ( ii < leng){
        val i = ind(ii)
        val Yturb = DenseMatrix.horzcat(Y(::,0 until Y.cols -1),(Y(::,Y.cols -1)+ron(::,i)).toDenseMatrix.t)
        //val Yturb = Y
       // Yturb(::,Yturb.cols -1) :+= ron(::,i)
        val coef = plsCalc.plsTrain(X,Yturb,k)
        val Ye = plsCalc.predict(X,coef)
        est(::,i) := Ye(::,Ye.cols -1).toDenseVector
        ii += 1
      }
    }
    getEsts(Array(0 until nPerm :_*))
//    var jj = 0
//    while ( jj < n){
//      val ronj = ron(jj,::)
//      val xtx = ronj * ronj.t
//      val xtxD = 1.0f / xtx
//      val x = xtxD * ronj * est(jj,::).t
//      h(jj) = x(0)
//      jj += 1
//    }
    //return(sum(h))
    pgdof(ron,est)
  }
  def gdofY(Y:DenseMatrix[Float],theta:Float,nPerm:Int = 1000): DenseMatrix[Float] ={
    val g = breeze.stats.distributions.Gaussian(0, theta)
    val n = Y.rows
    val ron = new DenseMatrix[Float](n,nPerm,g.sample(nPerm*n).toArray.map(_.toFloat))
    ron(::,*) + Y(::,Y.cols -1)
  }
  def pgdof(ron:DenseMatrix[Float],est: DenseMatrix[Float]): Float ={
    val n = est.rows
    val h = DenseVector.zeros[Float](n)
    var jj = 0
    while ( jj < n){
      val ronj = ron(jj,::)
      val xtx = ronj * ronj.t
      val xtxD = 1.0f / xtx
      val x = xtxD * ronj * est(jj,::).t
      h(jj) = x(0)
      jj += 1
    }
    return(sum(h))
  }

  def kramerDof(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int = 1,R: org.ddahl.rscala.RClient= org.ddahl.rscala.RClient()):Array[Float] = {

    //val R = org.ddahl.rscala.RClient()
    R.X = JavaArrayOps.dmToArray2(X).map(_.map(_.toDouble))
    R.Y = JavaArrayOps.dmToArray2(Y).flatten.map(_.toDouble)

    val rcomm =  """
              require(plsdof)
             |pl = pls.model(X,Y,m=6,compute.DoF=TRUE)
             |pld = pl$DoF
           """.stripMargin

    R eval rcomm.replaceAll("6",k.toString)

    val rs = R.getD1("pld").map(_.toFloat)

    //R.exit()
    rs

  }
  def ngdof(X:DenseMatrix[Float], Y:DenseMatrix[Float], k:Int = 1, nPerm:Int = 1000) = {
    val n = Y.rows
    val m = Y.cols
    val Xmean = utils.meanColumns(X)
    val Xm = X - utils.vertcat(Xmean, X.rows)
    val tenFold = plsCalc.kfoldInx(n,10,true)
    val YpredTenFold = plsCalc.plsCV(Xm,Y,k,tenFold)
    val rssfold = 0.6f * sqrt(plsCalc.rss(Y,YpredTenFold))
    val gdoffold = Array(1 to k:_*).map { i => // plsCalc.gdof(X,Y,rssfold(i),i+1,1000,1))
      //val Yest = gdofY(Y, rssfold(i-1), nPerm)
      val g = breeze.stats.distributions.Gaussian(0, rssfold(i-1))
      val ron = new DenseMatrix[Float](n,nPerm,g.sample(nPerm*n).toArray.map(_.toFloat))
      val Yest = ron(::, *) + Y(::, m - 1)
      val Yp = DenseMatrix.zeros[Float](n,m * nPerm)
      (0 until m-1).foreach(i => Yp(::,i until m * nPerm by m) := tile(Y(::,i),1,nPerm))
      Yp(::,m-1 until m * nPerm by m) := Yest
      val Ymean = utils.meanColumns(Yp)
      val Ym = Yp - utils.vertcat(Ymean, n)
      val Yms = plsCalc.plsP(Xm, Ym, i)
      val Yes = Yms._1 + utils.vertcat(Ymean,n)
      val Ye = Yes(::,m-1 until m * nPerm by m)
      plsCalc.pgdof(ron, Ye)
    }
    gdoffold
  }
  def ngdofP(X:DenseMatrix[Float], Y:DenseMatrix[Float], k:Int = 1, nPerm:Int = 1000): (Array[Float],Array[Float]) ={
//    val Xmean = utils.meanColumns(Xs)
//    val X = Xs - utils.vertcat(Xmean, Xs.rows)
//    val Ymean = utils.meanColumns(Y)
//    val Ys = Y - utils.vertcat(Ymean, Y.rows)
    val Ys = calculation.standardization(Y)
    //val n = Y.rows
//    val X = calculation.standardization(Xs)
    val m = Y.cols
    val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, Ys, k))
    val YY = Ys(::, m - 1).toDenseMatrix.t
    val gdof = plsCalc.ngdof(X,Ys,k,nPerm)
   //val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, Ys, k))
    val pval = plsCalc.dofPval(YY,Yhat(::,m-1 until Yhat.cols by m),gdof)
    (gdof,pval)
  }
  def dofPvalF(Xs:DenseMatrix[Float],Ys:DenseMatrix[Float],k:Int = 1,nPerm:Int = 1000,dof:Boolean = true) = {
    //val Ys = calculation.standardization(Y)
    val X = calculation.standardization(Xs)
    val n = Ys.rows
    val m = Ys.cols
    val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, Ys, k))
    val YY = Ys(::, m - 1).toDenseMatrix.t
    //    val coef = plsTrain(X,Y,k)
    //    val Yhat = predict(X, coef)

    //val tenFold = kfoldInx(n,10,true)
    val looInx = Array(0 until n :_*).map(Seq(_))
    //val YpredTenFold = plsCV(X,Y,k,tenFold)
    val YpredLoo = plsCV(X,Ys,k,looInx)
    //val pdof= pdof(Y,Yhat(::,0).asDenseMatrix.t,YpredTenFold)
    //??
    val pdofl = if(dof) pdof(YY,Yhat(::,m-1 until Yhat.cols by m),YpredLoo).toArray else null
    //val meanY = sum(Y)/n

    val gdof = if(dof) ngdof(X,Ys,k,nPerm) else null
    val yh = Yhat(::,m-1 until Yhat.cols by m)
    //val ppval = fdpval(yup,ydn,ny,pdofl)
    //val gpval =fdpval(yup,ydn,ny,gdof)
    //fdpval(y21,nmy,ny,dof)
    //if(dof == null) {
    //  val kdof = kramerDof(X,Ys,k).drop(1).map(_.toFloat)
    //  val kpval =fdpval(yup,ydn,ny,kdof)
    //} else dof


    (YY,yh,pdofl,gdof)

  }
  def dofPvalA(X:DenseMatrix[Float],Ys:DenseMatrix[Float],k:Int = 1,nPerm:Int = 1000) = {
    val (yy,yh,pdofl,gdof) = dofPvalF(X,Ys,k,nPerm,true)
    //val gdof = ngdof(X,Ys,k,nPerm)
    val kdof = kramerDof(X,Ys,k).drop(1).map(_.toFloat)
    val kpval =dofPval(yy,yh,kdof)
    val ppval = dofPval(yy,yh,pdofl)
    val gpval =dofPval(yy,yh,gdof)
    (kdof,kpval,pdofl,ppval,gdof,gpval)
  }
  def plsPerm(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int = 1,nPerm:Int = 1000) = {
    val Xmean = utils.meanColumns(X)
    val Ymean = utils.meanColumns(Y)
    val Xm = X - utils.vertcat(Xmean, X.rows)
    val Ym = Y - utils.vertcat(Ymean, Y.rows)
    val YY = plsCalc.permY(Ym,nPerm)
    val y0 = plsCalc.plsP(Xm,YY,k,nPerm)
    val py = plsCalc.rss(YY,y0._1)
    (1.0f/nPerm) * sum(py.map(i => if(i < py(0)) 1 else 0))
  }
  def plsAdof(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float] = null)= {
    val kk = min(X.cols,k)
    val rs = plsCalc.dofPvalA(X,Y,k,1000)
    Array(Array(1 to k:_*).map(i => plsCalc.plsPerm(X,Y,i,10000)),rs._1.map(_.toFloat),rs._2,rs._3,rs._4,rs._5,rs._6).map(_.mkString("\t"))
  }

  def plsSimu(X:DenseMatrix[Float],h:Float,k:Int,dof:Array[Float] = null) = {
    val pcs1 = princomp(convert(X,Double)).scores
    val mav1 = breeze.stats.meanAndVariance(pcs1(::, 0))
    var pc11 = convert((pcs1(::, 0) - mav1.mean) / sqrt(mav1.variance),Float)
    val theta1 = vegas2.getTheta(pc11)
    //val theta2 = vegas2.getTheta(pc11)
    val Y = (sqrt(h) * pc11 + sqrt(1 - h) * theta1).toDenseMatrix.t
    val Ys = calculation.standardization(Y)
    plsAdof(X,Ys,k,dof)
  }
  case class plsResult(Xm:DenseMatrix[Float],Ym:DenseMatrix[Float],mod:DenseMatrix[Float])

  def gdofAll(X:DenseMatrix[Float],Y:DenseMatrix[Float],theta:Float,k:Int =1,nPerm:Int = 1000,ncor:Int = 1) = {
    val g = breeze.stats.distributions.Gaussian(0, theta)
    val n = Y.rows
    val m = Y.cols
    val ron = new DenseMatrix[Float](n, nPerm, g.sample(nPerm * n).toArray.map(_.toFloat))
    val est = DenseMatrix.zeros[Float](n, nPerm*m*k)
    val Yest = Array(1 to k: _*).map(_ => est)
    val h = DenseMatrix.zeros[Float](n,m*k)
    val mlen = 2 * nPerm / ncor

    def getEsts(ind: Array[Int]): Unit = {
      val nlen = ind.length
      if (nlen < mlen) {
        getEst(ind)
      } else {
        val sep = nlen / 2
        myParallel.parallel(getEsts(ind.slice(0, sep)), getEsts(ind.slice(sep, nlen)))
      }
    }
    def getEst(ind: Array[Int]) = {
      var ii: Int = 0
      val leng = ind.length
      while (ii < leng) {
        val i = ind(ii)
        val Yturb =Y(::,*) + ron(::,i)
        //val Yturb = Y
        val coef = plsCalc.plsTrain(X, Yturb, k)
        val Ye = plsCalc.predict(X, coef)
        est(::, i*m*k until (i+1)*m*k) := Ye
        ii += 1
      }
    }

    getEsts(Array(0 until nPerm: _*))
    var yinx = 0
    while (yinx < m*k) {
      var jj = 0
      while (jj < n) {
        val ronj = ron(jj, ::)
        val xtx = ronj * ronj.t
        val xtxD = 1.0f / xtx
        val x = xtxD * ronj * est(jj, yinx until est.cols by m*k).t
        h(jj, yinx) = x(0)
        jj += 1
      }
      yinx += 1
    }
    sum(h(::,*)).t.toArray
  }
  def plsTrain(X: DenseMatrix[Float], Y: DenseMatrix[Float], comp: Int = 10) = {
    val Xmean = utils.meanColumns(X)
    val Ymean = utils.meanColumns(Y)
    val Xm = X - utils.vertcat(Xmean, X.rows)
    val Ym = Y - utils.vertcat(Ymean, Y.rows)
    val mod = pls(Xm, Ym, comp)
    val rs = new plsResult(Xmean, Ymean, mod)
    rs
  }
  def predict(X:DenseMatrix[Float],rs:plsResult) ={
    var mid = (X - utils.vertcat(rs.Xm, X.rows)) * rs.mod
    val m = rs.Ym.toArray
    mid(*,::) :+= DenseVector(Array.tabulate(mid.cols/m.length)(i => m).flatten)//rs.Ym.toDenseVector
    //val m = if (n == -1) rs.Ym.toDenseVector.length - 1 else n
    //val mid3 = mid2.map(_(::,m))
    //DenseMatrix(mid3:_*).t
    mid
  }
  def getCVP(X: DenseMatrix[Float], Y: DenseMatrix[Float], testInx: Seq[Int], comp: Int = 10):DenseVector[Float]= {
    val head = testInx.head
    val last = testInx.last
    val Xtest = X(head to last, ::)
    val Xtrain = X.delete(testInx, Axis._0)
    val Ytest = Y(head to last, Y.cols -1).toDenseVector
    val Ytrain = Y.delete(testInx, Axis._0)
    val meanMod = plsTrain(Xtrain, Ytrain, comp)
    //val pre = (Xtest - utils.vertcat(meanMod.Xm, Xtest.rows)) * meanMod.mod
    val pred = predict(Xtest,meanMod)// + utils.vertcat(meanMod.Ym,pre.rows)
    val msep = breeze.numerics.pow(pred(::,*) - Ytest, 2)
    //      val msepTrain = msep.delete(testInx,Axis._0)
    sum(msep(::, *)).t
    //      (,sum(msepTrain(::,*)).t
  }
  def rss(Y:DenseMatrix[Float],Yfit:DenseMatrix[Float]):DenseVector[Float] ={
    val n = Y.rows.toFloat
    val k = Yfit.cols
    val m = Y.cols
    val fit = DenseMatrix.zeros[Float](Y.rows,k)
    //val y = Y.toDenseVector
    for ( i <- 0 until k/m){
      fit(::,i*m until (i+1)* m ) := Yfit(::,i*m until (i+1)* m ) - Y
    }
    //val fit = Yfit(::,*) - y
    var ssfit = DenseVector.zeros[Float](k)//if (k>1) diag(fit.t * fit) else DenseVector(fit.t * fit)//
    for ( i <- 0 until k){
      ssfit(i) = fit(::,i).t * fit(::,i)
    }
    ssfit/n
  }
  def plsAll(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int) = {
    val coef = plsTrain(X,Y,k)
    val Yhat = predict(X, coef)
    val n = Y.rows
    val tenFold = kfoldInx(n,10,true)
    val looInx = Array(0 until n :_*).map(Seq(_))
    val YpredTenFold = plsCV(X,Y,k,tenFold)
    val YpredLoo = plsCV(X,Y,k,looInx)
    val pdoffold = pdof(Y,Yhat(::,0).asDenseMatrix.t,YpredTenFold)
    val pdofloo = pdof(Y,Yhat(::,0).asDenseMatrix.t,YpredLoo)
    val meanY = sum(Y)/n
  }
  def gdofPlsPval(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int):(Array[Float],Array[Float]) = {

    val n = Y.rows
    val m = Y.cols
    val tenFold = plsCalc.kfoldInx(n,10,true)
    //var looInx = Array(0 until n :_*).map(Seq(_))
    val YpredTenFold = plsCalc.plsCV(X,Y,k,tenFold)
    //val YpredLoo = plsCalc.plsCV(X,Y,k,looInx)
    val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,k))
    val YY = Y(::,m -1).toDenseMatrix.t
    //val pdofloo = plsCalc.pdof(Y,Yhat,YpredLoo)
    //val pdoffold = plsCalc.pdof(Y,Yhat,YpredTenFold)
    //val rssloo = 0.6*sqrt(plsCalc.rss(Y,YpredLoo))
    //val gdofloo = Array(0 until k:_*).map(i => plsCalc.gdof(X,Y,rssloo(i),i+1,1000))
    val rssfold = 0.6f * sqrt(plsCalc.rss(Y,YpredTenFold))

    val gdoffold = Array(0 until k:_*).map(i => plsCalc.gdof(X,Y,rssfold(i),i+1,1000,1))

    val gTenPval = dofPval(YY,Yhat(::,m-1 until Yhat.cols by m),gdoffold)
    (gdoffold,gTenPval)
  }
  def gdofPlsPval(X:DenseMatrix[Float],Y:DenseVector[Float],k:Int):(Array[Float],Array[Float]) = {
    val Ys = Y.toDenseMatrix.t
    gdofPlsPval(X,Ys,k)
  }
  def gdofPlsPvalAll(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int):(Array[Float],Array[Float]) = {
    val n = Y.rows
    val m = Y.cols
    val tenFold = plsCalc.kfoldInx(n,10,true)
    val YpredTenFold = plsCalc.plsCV(X,Y,k,tenFold)
    val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,k))
    val rssfold = breeze.stats.mean(0.6f * sqrt(plsCalc.rss(Y,YpredTenFold)))
    val gdoffold = plsCalc.gdofAll(X,Y,rssfold,k,1000,1)

    val gTenPval = dofPval(Y,Yhat,gdoffold)
    (gdoffold,gTenPval)
  }
  def gdofPval2(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int) = {
    val (d,p) = gdofPlsPvalAll(X,Y,k)
    val corr = calculation.pearsonCorr(Y(::,0).toArray.map(_.toFloat),Y(::,1).toArray.map(_.toFloat))
    val up = Array(1 to k:_*).map(i => calculation.brownTest(p.slice((i-1)*2,i*2),Array(corr)))
    (d,p,up)
  }
  def ngdofPvalT(X:DenseMatrix[Float],Y:DenseMatrix[Float],nk:Int) = {
    val ny = Y.cols -1
    val inx = ny +: Seq(0 to (ny - 1):_*)
    var Yi = Y
    val k = scala.math.min(X.cols,nk)
    var (df,pl) = ngdofP(X,Yi,k)
    //var pl = Array[Float]()
    var i = 0
    while (i < ny){
      val Yii = Yi(::,inx).toDenseMatrix
      Yi = Yii
      val(dfi,pli) = plsCalc.ngdofP(X,Yi,k)
      df ++:= dfi
      pl ++:= pli
      i += 1
    }
    val corr = breeze.stats.covmat(convert(Y,Double))
    var corrr = lowerTriangular(corr)
    diag(corrr) := 0d
    var fcor = corrr.toArray.filter( _ != 0d).map(_.toFloat)
    val up = Array(0 to (k -1) :_*).map(i => calculation.brownTest(Array(i until k*(ny+1) by k :_*).map(pl(_)),fcor))
    (df,pl,up)
  }
  def ngdofPvalS(X:DenseMatrix[Float],Y:DenseMatrix[Float],nk:Any):String = {
    val nks = nk.asInstanceOf[Int]
    val rs = ngdofPvalT(X,Y,nks)
    rs._1.mkString("\t") + "\t"+rs._2.mkString("\t")+"\t"+rs._3.mkString("\t")
  }
}

package PLS

import snpCalcActor._
import myParallel.actorMessage._
import akka.actor._
import breeze.linalg._
import breeze.numerics._
import com.sun.corba.se.impl.orbutil.DenseIntMapImpl//{done, system}//{Actor, Props}
object snpCalcActor{
  val name = "snpCalcActor"
  def props(fil:snpCalcPms) = Props(classOf[snpCalcActor],fil)
  case class snpCalcPms(k:Int,n:Int,perm:Int = 0,Y:DenseMatrix[Float],
                        looInx:Array[Seq[Int]], tenFold:Array[Seq[Int]],oname:String)
  case class geneSNPResult(Ypred:DenseMatrix[Float],pdofLoo:Array[Float],gdofLoo:Array[Float],
                           pdofTenf:Array[Float],gdofTenf:Array[Float],permPval:Array[Float])
  case class Xs(gene:Array[String],X:DenseMatrix[Float])
  case class XYs(gene:Array[String],X:DenseMatrix[Float],Y:DenseMatrix[Float])
  case class Ys(Y:DenseMatrix[Float])
  case class func(f:(DenseMatrix[Float],DenseMatrix[Float],Any) => String)
  case class calcPm(pms:Any)
  case class writerName(name:String)
  def getColnames(ks:Int)= {
    val colname = Array("ucscName", "chr", "start", "end", "nsnp")
    val permn = Array(1 to ks: _*).map("perm" + _)
    val pdofloo = Array(1 to ks: _*).map("pdofLoo" + _)
    val pdoflooPval = Array(1 to ks: _*).map("pdofLooPval" + _)
    val pdoften = Array(1 to ks: _*).map("pdofTenf" + _)
    val pdoftenPval = Array(1 to ks: _*).map("pdofTenfPval" + _)
    val gdofloo = Array(1 to ks: _*).map("gdofLoo" + _)
    val gdoflooPval = Array(1 to ks: _*).map("gdofLooPval" + _)
    val gdoften = Array(1 to ks: _*).map("gdofTenf" + _)
    val gdoftenPval = Array(1 to ks: _*).map("gdofTenfPval" + _)
    colname ++ permn ++ gdofloo ++ gdoflooPval ++ gdoften ++ gdoftenPval ++ pdofloo ++ pdoflooPval ++ pdoften ++ pdoftenPval
  }
}
class snpCalcActor(pms:snpCalcPms) extends Actor{
  var k = pms.k
  //var pk = pms.k
  var n = pms.n
  var perm = pms.perm
  var Y = pms.Y//DenseMatrix.zeros[Float](n,perm)
//  var Y = permY(::,0).toDenseMatrix.t //DenseMatrix.zeros[Float](n,1)
  var cr = calculation.pearsonCorr(Y(::,0).toArray,Y(::,1).toArray)
  var looInx = pms.looInx//Array(0 until n :_*).map(Seq(_))
  var tenFold = pms.tenFold//plsCalc.kfoldInx(n,10,true)
  //var mcol = pms.mcol//:(Array[Int], Array[Int]) = (Array(0),Array(0))
  var order:Option[ActorRef] = None
    var calcPms:Any = 3//(k,n,perm,looInx,tenFold)
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/writer"))
//  def permT(X:DenseMatrix[Float]) = {
//    var i = 0
//    var nk = if(X.cols < k) X.cols else k
//    var rss = new DenseMatrix[Float](nk,perm)
//    while (i < perm ) {
//      val y1 = permY(::, i).toDenseMatrix.t
//      val Yhat = plsCalc.predict(X, plsCalc.plsTrain(X, y1, nk))
//      rss(::,i) := plsCalc.rss(y1, Yhat)
//      i += 1
//    }
//    val st = rss(::,0).toArray
//    Array(0 until nk :_*).map(i => rss(i,::).t.toArray.count(ii => ii < st(i)).toFloat / perm.toFloat)
//  }

//  def getGeneSnpSimp(X:DenseMatrix[Float])= {
//    val nk = if(X.cols < k) X.cols else k
//    val YpredTenFold = plsCalc.plsCV(X,Y,nk,tenFold)
//    val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,nk))
//    val rssfold = 0.6f * sqrt(plsCalc.rss(Y,YpredTenFold))
//    val gdoffold = Array(0 until nk:_*).map(i => plsCalc.gdof(X,Y,rssfold(i),i+1,1000))
//    val pv = dofPval(Y,Yhat,gdoffold)
//    val earray = Array.fill(nk-k)(-1.0)
//    gdoffold++earray++pv++earray
//  }
//  def getGeneSnp(X:DenseMatrix[Float]):geneSNPResult= {
//    val nk = if(X.cols < k) X.cols else k
//    val YpredTenFold = plsCalc.plsCV(X,Y,nk,tenFold)
//    val YpredLoo = plsCalc.plsCV(X,Y,nk,looInx)
//    val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,nk))
//    val pdofloo = plsCalc.pdof(Y,Yhat,YpredLoo)
//    val pdoffold = plsCalc.pdof(Y,Yhat,YpredTenFold)
//    val rssloo = 0.6f * sqrt(plsCalc.rss(Y,YpredLoo))
//    val gdofloo = Array(0 until nk:_*).map(i => plsCalc.gdof(X,Y,rssloo(i),i+1,1000))
//    val rssfold = 0.6f * sqrt(plsCalc.rss(Y,YpredTenFold))
//    val gdoffold = Array(0 until nk:_*).map(i => plsCalc.gdof(X,Y,rssfold(i),i+1,1000))
//    val permPval = permT(X)
//    //.map(i => new DenseMatrix[Float](1,1,Array(i)))
//    geneSNPResult(Yhat,pdofloo.data,gdofloo,pdoffold.data,gdoffold,permPval)
//  }
  def dofPval(Ys:DenseMatrix[Float]=Y,Ypred:DenseMatrix[Float],dof:Array[Float])={
    val nk = dof.length
    Array(0 until nk:_*).map(i => plsCalc.anova(Ys,Ypred(::,i).toDenseMatrix.t,new DenseMatrix[Float](1,1,Array(dof(i))))(0))

  }

  def genePval(grs:geneSNPResult) = {
    val nk = grs.pdofLoo.length
    val pLooPval = dofPval(Y,grs.Ypred,grs.pdofLoo)
    val pTenPval = dofPval(Y,grs.Ypred,grs.pdofTenf)
    val gLooPval = dofPval(Y,grs.Ypred,grs.gdofLoo)
    val gTenPval = dofPval(Y,grs.Ypred,grs.gdofTenf)
    val rs = Array(grs.permPval,grs.gdofLoo,gLooPval,grs.gdofTenf,gTenPval,grs.pdofLoo,pLooPval,grs.pdofTenf,pTenPval)
    rs.map(i => i ++ Array.fill(nk-k)(-1.0))

  }
  var getRes:(DenseMatrix[Float],DenseMatrix[Float],Any) => String = (X:DenseMatrix[Float],Y:DenseMatrix[Float],pm:Any) => {
    //val grs = getGeneSnp(X)
    //genePval(grs).mkString("\t")
    val pmls = pm.asInstanceOf[Int]//,Int,Int,Array[Seq[Int]],Array[Seq[Int]])]
    val m = X.cols
    val k = min(m,pmls)
    val rs = plsCalc.ngdofPvalT(X,Y,k)
    rs._2.mkString("\t")+"\t"+ rs._3.mkString("\t")


//    val rs = plsCalc.dofPvalA(X,Y,k)
//
//    val rsp = Array( 1 to k:_*).map(i => plsCalc.plsPerm(X,Y,i,10000))
//    rsp.mkString("\t")+"\t"+rs._1.mkString("\t")+"\t"+ rs._2.mkString("\t")+"\t"+ rs._3.mkString("\t")+"\t"+ rs._4.mkString("\t")+"\t"+ rs._5.mkString("\t")+"\t"+ rs._6.mkString("\t")
  }

  def receive = {
    case wrt:writerName => {
      //order = Some(sender)
      writer = Some(system.actorSelection("/user/"+wrt.name))
    }
    case func:func => {
      getRes = func.f
    }
    case pm:calcPm => {
      this.calcPms = pm.pms
    }
    case ys:Ys => {
      this.Y = ys.Y
    }
    case x:Xs => {
      val X = x.X
      if (Y.cols > 1){
        var rs = "" //x.gene.mkString("\t")
        for (i <- 0 until Y.cols){
          val y = Y(::,i).toDenseMatrix.t
          rs += "\t"+getRes(X,y,calcPms)
        }
        val kk = calcPms.asInstanceOf[Int]
        val rrs = rs.split("\t").drop(1).map(_.toFloat)
        val rss = (x.gene ++ rrs.map(_.toString) ++ Array(0 until kk :_*).map(i => calculation.
          brownTest(Array(rrs(i),rrs(i +  kk -1)),Array(cr)).toString)).mkString("\t") +"\t" + getRes(X,Y,calcPms)
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rss))
      }else {
        val prs = getRes(X,Y,calcPms)
        var rs = x.gene.mkString("\t") + "\t"+prs
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
      }
      //val prs = genePval(grs).mkString("\t")
      sender ! done(1)
    }

    case xy:XYs => {

      val X = xy.X
      val Ys = if (max(Y) == 0f & min(Y) == 0f) xy.Y else DenseMatrix.horzcat(xy.Y,Y)
      val prs = getRes(X,Ys,calcPms)
      val yrs = getRes(X,Y,calcPms)
      val Xs = DenseMatrix.horzcat(X,xy.Y)
      val xrs = getRes(Xs,Y,calcPms)
      //val prs = genePval(grs).mkString("\t")
      var rs = xy.gene.mkString("\t")+"\t" +X.cols+ "\t"+prs+"\t"+"0"+"\t"+yrs+"\t"+"0"+"\t"+xrs
      writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
      sender ! done(1)
    }

    case _ =>

  }
}

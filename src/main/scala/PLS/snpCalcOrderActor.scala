package PLS

import PLS.snpCalcActor.snpCalcPms
import PLS.snpCalcDispatchActor._
import PLS.snpCalcOrderActor._
import akka.actor._
import breeze.linalg.DenseMatrix
import myParallel.actorMessage
import myParallel.actorMessage._
import PLS.snpCalcActor._
import akka.Done

import scala.concurrent.Future
object snpCalcOrderActor{
  val name = "snpCalcOrderActor"
  def props(pms:orderPms) = Props(classOf[snpCalcOrderActor],pms)
  case class chrs(chrs:Array[String])
  case class yArray(y:Array[Array[String]])
  case class orderPms(k:Int = 5,perm:Int = 1000,amp:Int = 0,nactor:Int = 10, gfile:String=gPms.op+gPms.glf,
                           dfile:String = gPms.op+gPms.df, ifile:String=gPms.op+gPms.af,
                           pfile:String=gPms.op+gPms.pf,ofile:String = gPms.rp+gPms.rsf,efile:String = gPms.op+gPms.ef)
}

class snpCalcOrderActor(pm:orderPms) extends Actor{
  var gfile = pm.gfile//"GeneLoc.txt"
  var ifile = pm.ifile//"snp6annoNew.txt"
  var dfile = pm.dfile//"gbmsnp6.csv"
  var pfile = pm.pfile///"GBMLGG.rppa.txt"
  var ofile = pm.ofile//"GBMsnp6Rs.txt"
  var efile = pm.efile
  var process:chrs = chrs(Array(1 to 22:_*).map(_.toString))
  var dispatcher:Option[ActorRef] = None
  var count = 0
  var len:Int = 1000

  var amplength =pm.amp
  var k = pm.k
  var perm = pm.perm
  var n = 0
  var cores = Runtime.getRuntime.availableProcessors()-1
  var nActor:Int = if (cores > 50) 50 else cores//pm.nactor
  //var len = 900000
  var cnt = 0
  var diff = 0
  var mcol:(Array[Int], Array[Int],Array[Int]) = (Array(0),Array(0),Array(0))
  var dm = new DenseMatrix[Float](n,100)
  var Y = DenseMatrix.zeros[Float](n,1)
  var permY = DenseMatrix.zeros[Float](n,perm)
  var looInx = Array(0 until n :_*).map(Seq(_))
  var tenFold = plsCalc.kfoldInx(n,10,true)
  var writer:Option[ActorRef] = None
  def updateY(y:Array[Array[String]]) = {
    val nl = y.length
    val pn = fileOper.toArrays(pfile).next.map(_.slice(0,15))
    val dn = fileOper.toArrays(dfile).next
    if (efile.length >0) {
      val en = fileOper.toArrays(efile, "\t").next.map(_.slice(0,15))
      //val pakt = y
      this.mcol = fileOper.intersectCols(dn, pn, en)
    }else {
      val mcoll = fileOper.intersectCols(dn, pn)
      this.mcol = (mcoll._1, mcoll._2, Array(0))
    }
    this.n = mcol._1.length
    this.looInx = Array(0 until n :_*).map(Seq(_))
    this.tenFold = plsCalc.kfoldInx(n,10,true)
    this.Y = new DenseMatrix(n,nl,y.map(i=>mcol._2.map(i(_)).map(_.toFloat)).flatten)
    //this.permY = if(perm > 0) plsCalc.permY(Y(::,0),perm) else this.Y
  }
//  override def preStart { println("kenny: prestart") }

  def updateFile(chrn:String)= {
    val regex = "(chr[0-9]+)?\\.(txt|csv)".r
    val newn = "chr"+chrn+".txt"
    this.ifile = regex.replaceAllIn(ifile,newn)// sfileo.replace(".txt","chr"+chrn+".txt")
    //this.rfile = regex.replaceAllIn(rfile,newn)// rfileo.replace(".txt","chr"+chrn+".txt")
    this.dfile = regex.replaceAllIn(dfile,newn)// dfileo.replace(".csv","chr"+chrn+".txt")
    this.gfile = regex.replaceAllIn(gfile,newn)// gfileo.replace(".txt","chr"+chrn+".txt")
    //fileOper.FindIndexinRange(gfile,ifile,rfile)
  }

  def makeActor() = {
    val calcPm = snpCalcActor.snpCalcPms(k,n,perm,Y,looInx,tenFold,this.ofile.replaceFirst("/",""))
    Array(0 until nActor:_*).map(i => system.actorOf(snpCalcActor.props(calcPm),"calc"+i))
    val wrt = system.actorOf(myParallel.paraWriterActor.props(myParallel.paraWriterActor.fileName(this.ofile)),"writer")

    val colnames = snpCalcActor.getColnames(k).mkString("\t")
    wrt ! myParallel.paraWriterActor.WriteStr(colnames)
    //wrt ! myParallel.paraWriterActor.totalNumber(2)
    writer = Some(wrt)
  }

  def receive = {
    case y:yArray => {
      updateY(y.y)
    }

    case c:chrs =>{
      this.process = c
      this.len = c.chrs.length
      makeActor()
      updateFile(c.chrs.apply(count))

      val dispPm = snpCalcDispatchActor.dispatcherPms(amplength, nActor, mcol._1,mcol._3, this.gfile, this.dfile, this.ifile, this.ofile, this.efile)
      val dispatch = context.actorOf(snpCalcDispatchActor.props(dispPm), "dispatcher"+c.chrs.apply(count))
      dispatcher = Some(dispatch)
      count += 1
        //dispatcher.foreach(_ ! actorMessage.action)
    }

    case cfunc:func =>{
      Array(0 until nActor:_*).map("/user/calc"+_).map(system.actorSelection(_)).foreach(_ ! cfunc)
    }

    case fpm:calcPm =>{
      Array(0 until nActor:_*).map("/user/calc"+_).map(system.actorSelection(_)).foreach(_ ! fpm)
    }

    case actorMessage.action => {
      dispatcher.foreach(_ ! actorMessage.action)
    }

    case don:done => {
      if (count < len) {
        updateFile(process.chrs.apply(count))
        val dispPm = snpCalcDispatchActor.dispatcherPms(amplength, nActor, mcol._1,mcol._3, this.gfile, this.dfile, this.ifile, this.ofile, this.efile)
        val dispatch = context.actorOf(snpCalcDispatchActor.props(dispPm), "dispatcher"+process.chrs.apply(count))
        dispatcher = Some(dispatch)
        dispatcher.foreach(_ ! actorMessage.action)
        //sender !  snpCalcDispatchActor.chr(process.chrs.apply(count))
        count += 1
      } else {
       // dispatcher.foreach(_ ! PoisonPill)
        writer.foreach(_ ! actorMessage.finished)
        Array(0 until nActor:_*).map("/user/calc"+_).map(system.actorSelection(_)).foreach(_ ! PoisonPill)
        println(utils.currentTimeIn+s"calculation is done -- Order")
        self ! PoisonPill
        val done: Future[Done] = CoordinatedShutdown(system).run()

        //system.shutdown
      }
    }
    //case _ =>
  }
}

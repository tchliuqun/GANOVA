package PLS

import PLS.setDispatchActor._
import PLS.snpCalcActor._
import akka.actor._
import breeze.linalg._
import myParallel.actorMessage._
import myParallel.actorMessage
import scala.collection.mutable.ArrayBuffer
object setDispatchActor{

  val name = "snpCalcDispatchActor"
  def props(pms:dispatcherPms) = Props(classOf[setDispatchActor],pms)
  case class dispatcherPms(k:Int = 3,perm:Int = 1000,gfile:String= gPms.op+gPms.glf, sfile:String =  gPms.op+gPms.df, pfile:String= gPms.op+gPms.af,pname:Array[String],ofile:String = "setResult.txt",cores:Int = 0)

  case class dataf(file:String)
  case class phenof(file:String)
  case class genef(file:String)
  case class outf(file:String)
  case class annof(file:String)
  case class rangef(file:String)
  case class kn(num:Int)
  case class permn(num:Int)
  case class ampn(num:Int)
  case class Yarray(array:Array[Float])
}
class setDispatchActor(pm:dispatcherPms) extends Actor{
  var gfile = pm.gfile//"GeneLoc.txt"
  var sfile = pm.sfile//"snp6annoNew.txt"
  var pfile = pm.pfile//"gbmsnp6.csv"
  var ofile = gPms.rp+ pm.ofile
  var phname = pm.pname
  var goChildFile = gPms.op + "goose"
  var ncores = if (pm.cores < 1) Runtime.getRuntime.availableProcessors() - 1 else pm.cores
  var cores = if (ncores > 20 ) 20 else ncores
  var wrt:Option[ActorSelection] = Some(system.actorSelection("/user/"+pm.ofile))//"GBMsnp6Rs.txt"
  //var ampl =pm.amp
  //var nActor:Int = pm.nactor
  var len = 900000
  var cnt = 0
  var k = pm.k
  var perm = pm.perm
  var n = 0
  var looInx = Array(0 until n :_*).map(Seq(_))
  var tenFold = plsCalc.kfoldInx(n,10,true)
  //var mcol:Array[Int] = pm.mcolX
  //var mcolY:Array[Int] = pm.mcolY
  var gset = Array[Array[String]]()//
  var gexp = Array[Array[String]]()//
  //var gpheno = Array[Array[Float]]()
  var snpI = Iterator[Array[String]]()//
  var genExp = Map[String, Array[Float]]()
  var order:Option[ActorRef] = None
  var actor:Option[Array[ActorRef]] = None
  var snpd = Iterator[Array[Float]]()
  var loc = 0
  var X = ArrayBuffer[(Int,Array[Float])]()
  var XX:Array[Array[Float]] = Array()
  var Y = DenseMatrix.zeros[Float](n,1)
  var permY = DenseMatrix.zeros[Float](n,perm)

  var goano = Map[String, Array[String]]()
  var completed:Boolean = false
  var sendCont:Int = 0
  var recieveCont = 0

  var gen = Array[String]()


  def updateXs = {
    val gname = scala.io.Source.fromFile(gfile).getLines.map(_.split("\t")).next.map(_.slice(0,15))
    val pname = scala.io.Source.fromFile(pfile).getLines.map(_.split("\t")).next.map(_.slice(0,15))
    //val aname = gname.intersect(pname).toArray
    val aseq = fileOper.intersectCols(gname,pname)
    val gseq = 0 +: aseq._1
    val pseq = 0 +: aseq._2
    val parray = scala.io.Source.fromFile(pfile).getLines.map(_.split("\t")).toArray
    val y = phname.map(i => parray.filter(_(0).contains(i)).flatten)
    val gpheno = y.map(i => pseq.map(i(_)).drop(1).map(_.toFloat))

    this.gset = scala.io.Source.fromFile(sfile).getLines.map(_.split("\t")).toArray
    this.goano = scala.io.Source.fromFile(goChildFile).getLines.map(_.split("\t").take(3).mkString("\t")).toArray.distinct.map(_.split("\t")).map( i => (i(0) -> i)).toMap
    this.gexp = scala.io.Source.fromFile(gfile).getLines.map(_.split("\t")).drop(2).map(i => i(0).split("\\|").apply(1) +: i.drop(1) ).toArray.map(i => gseq.map(i(_)))
    //this.gexp = scala.io.Source.fromFile(gfile).getLines.map(_.split("\t")).toArray.map(i => gseq.map(i(_)))
    this.n = gpheno(0).length
    this.Y = new DenseMatrix(n,phname.length,gpheno.flatten)
    this.permY = if(perm > 0) plsCalc.permY(Y,perm) else this.Y
    //this.Y =

    this.cnt = 0
    this.sendCont = 0
    this.recieveCont = 0
    this.len = gset.length
    //this.X = ArrayBuffer[(Int,Array[Float])]()
  }

  def getX(counts:Int = sendCont)= {
    gexp.filter(gset(counts).toSet contains _(0)).map(_.drop(1).map(_.toFloat))
  }
  def startP= {

    while (sendCont <= cores *2) {
      //val gen = gset(cnt).apply(0)
      //gen = Array(gset(cnt).apply(0))
      //print("gene length is"+gen.length)
      XX = getX(cnt)
      gen = goano.getOrElse(gset(cnt).apply(0),Array(gset(cnt).apply(0)))
      //val ye = if(mcolY.length > 1) genExp(gen(4)) else
      if (XX.length>2) {
        gen :+= XX.length.toString
        val na = sendCont % cores
        val calcular = system.actorSelection("/user/calc" + na)
//        if(mcolY.length == 1) {
          calcular ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
//        }else {
//          calcular ! snpCalcActor.XYs(gen,utils.Array2DM(XX, false),new DenseVector(genExp(gen(4))).asDenseMatrix.t)
//        }
        sendCont += 1
      }
      cnt += 1
    }
  }
  def makeActors() = {
    val calcPm = snpCalcActor.snpCalcPms(k,n,perm,Y,looInx,tenFold,this.ofile.replaceFirst("/",""))
    Array(0 until cores:_*).map(i => system.actorOf(snpCalcActor.props(calcPm),"calc"+i))

    //val calcPm = gseaCalcActor.Pms(this.ofile.replaceFirst("/",""))
    //Array(0 until cores:_*).map(i => system.actorOf(gseaCalcActor.props(calcPm),"calc"+i))
    system.actorOf(myParallel.paraWriterActor.props(myParallel.paraWriterActor.fileName(ofile)),"setwriter")
    Array(0 until cores:_*).foreach(i =>{
      val calcular = system.actorSelection("/user/calc" + i)
      calcular ! snpCalcActor.writerName("setwriter")
      calcular ! snpCalcActor.Ys(this.Y)
      calcular ! snpCalcActor.func(plsCalc.ngdofPvalS)
    }
    )
  }

  def receive = {
    case actorMessage.action => {
      println(utils.currentTimeIn+"generate calculator")
      order = Some(sender)
//      println(utils.currentTimeIn+"processing "+ dfile)
      updateXs
      makeActors()
      startP
      println(utils.currentTimeIn+"Start dispatching")
    }
//    case f:genef => this.gfile = f.file
//    case f:annof => this.ifile = f.file
//    case f:dataf => this.dfile = f.file
//    case n:ampn => this.ampl = n.num
    case don:actorMessage.done => {
      recieveCont += 1
      if (cnt < len){

        XX = getX(cnt)
        gen = goano.getOrElse(gset(cnt).apply(0),Array(gset(cnt).apply(0),"nan","nan"))

        cnt += 1
        while (XX.length< 3 & cnt < len) {

          XX = getX(cnt)
          gen = goano.getOrElse(gset(cnt).apply(0),Array(gset(cnt).apply(0)))
          cnt += 1
        }
        if(XX.length>2) {
//          //if(mcolY.length == 1) {
          gen :+= XX.length.toString
            sender ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
//          }else {
//            sender ! snpCalcActor.XYs(gen,utils.Array2DM(XX, false),new DenseVector(genExp(gen(4))).asDenseMatrix.t)
//          }
//          //sender ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
          sendCont += 1
//          cnt += 1
        }
        if (cnt % 200 == 0) println("processing" + cnt)
      }
      else {
        if (!completed) {
          completed = true


        }
      }
      if (completed & sendCont == recieveCont){
        order.foreach(_ ! done(1))
        self ! PoisonPill
      }
    }
    case actorMessage.finished => order.foreach(_ ! done(1))
    case _ =>
  }
}

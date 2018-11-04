package PLS

import PLS.snpCalcDispatchActor._
import PLS.snpCalcActor._
import akka.actor._
import breeze.linalg._
import myParallel.actorMessage._
import myParallel.actorMessage
import scala.collection.mutable.ArrayBuffer
object snpCalcDispatchActor{

  val name = "snpCalcDispatchActor"
  def props(pms:dispatcherPms) = Props(classOf[snpCalcDispatchActor],pms)
  case class dispatcherPms(amp:Int = 50000,nactor:Int = 10, mcolX:Array[Int] = Array(0),mcolY:Array[Int] = Array(0),gfile:String= gPms.op+gPms.glf,
                      dfile:String =  gPms.op+gPms.df, ifile:String= gPms.op+gPms.af,oname:String,efile:String = gPms.op+gPms.ef)

  case class chr(chr:String)
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
class snpCalcDispatchActor(pm:dispatcherPms) extends Actor{
  var gfile = pm.gfile//"GeneLoc.txt"
  var ifile = pm.ifile//"snp6annoNew.txt"
  var dfile = pm.dfile//"gbmsnp6.csv"
  var efile = pm.efile
  var wrt:Option[ActorSelection] = Some(system.actorSelection("/user/"+pm.oname))//"GBMsnp6Rs.txt"
  var ampl =pm.amp
  var nActor:Int = pm.nactor
  var len = 900000
  var cnt = 0
  var mcol:Array[Int] = pm.mcolX
  var mcolY:Array[Int] = pm.mcolY
  var gens = Array[Array[String]]()//
  var snpI = Iterator[Array[String]]()//
  var genExp = Map[String, Array[Float]]()
  var order:Option[ActorRef] = None
  var actor:Option[Array[ActorRef]] = None
  var snpd = Iterator[Array[Float]]()
  var loc = 0
  var X = ArrayBuffer[(Int,Array[Float])]()
  var XX:Array[Array[Float]] = Array()
  var completed:Boolean = false
  var sendCont = 0
  var recieveCont = 0
  var shd = 3

  def updateXs = {
    this.snpd = if(mcol.length == 1)fileOper.toArrays(dfile).map(_.drop(1).map(_.toFloat)) else fileOper.toArrays(dfile).map(i => mcol.map(i(_)).map(_.toFloat))
    this.gens = fileOper.toArrays(gfile).toArray
    this.snpI = fileOper.toArrays(ifile)
    if(mcolY.length > 1) this.genExp = fileOper.toArrays(efile).drop(1).map(i => (i(0),mcolY.map(i(_)).map(_.toFloat))).toMap
    this.loc = snpI.next.apply(2).toInt
    this.cnt = 0
    this.sendCont = 0
    this.recieveCont = 0
    this.len = gens.length
    this.X = ArrayBuffer[(Int,Array[Float])]()
  }

  def getX(gen:Array[String])= {
    var srt = gen(1).toInt - ampl
    var end = gen(2).toInt + ampl
    while (loc < end & snpd.hasNext){
      val x0 = snpd.next
      if (breeze.stats.meanAndVariance(x0).variance > 0f) X = X :+ (loc,x0)
      if(snpI.hasNext) loc = snpI.next.apply(2).toInt
    }
    X = X.dropWhile(_._1 < srt)
    X.takeWhile(_._1 < end).map(_._2).toArray
  }
  def startP= {

    while (sendCont <= nActor *2) {
      val gen = gens(cnt)
      //print("gene length is"+gen.length)
      XX = getX(gen)
      //val ye = if(mcolY.length > 1) genExp(gen(4)) else
        if (XX.length > shd & (mcolY.length == 1 | genExp.contains(gen(4)))) {
       // if (XX.length>0 & ) {
        val na = sendCont % nActor
        val calcular = system.actorSelection("/user/calc" + na)
        if(mcolY.length == 1) {
          calcular ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
        }else {
          calcular ! snpCalcActor.XYs(gen,utils.Array2DM(XX, false),new DenseVector(genExp(gen(4))).asDenseMatrix.t)
        }
        sendCont += 1
      }
      cnt += 1
    }
  }

  def receive = {
    case actorMessage.action => {
      println(utils.currentTimeIn+"generate calculator")
      order = Some(sender)
      println(utils.currentTimeIn+"processing "+ dfile)
      updateXs
      startP
      println(utils.currentTimeIn+"Start dispatching")
    }
    case f:genef => this.gfile = f.file
    case f:annof => this.ifile = f.file
    case f:dataf => this.dfile = f.file
    case n:ampn => this.ampl = n.num
    case don:actorMessage.done => {
      recieveCont += 1
      if (cnt < len){
        var gen = gens(cnt)
        XX = getX(gen)
        while ((XX.length <= shd | (mcolY.length > 1 & !genExp.contains(gen(4)))) & cnt < (len - 1) ) {
                 // while (XX.length<1 & cnt < (len - 1) ) {
            cnt += 1
            gen = gens(cnt)
            XX = getX(gen)
        }
        if(XX.length > shd & (mcolY.length == 1 | genExp.contains(gen(4)))) {
          if(mcolY.length == 1) {
            sender ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
          }else {
            sender ! snpCalcActor.XYs(gen,utils.Array2DM(XX, false),new DenseVector(genExp(gen(4))).asDenseMatrix.t)
          }
          //sender ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
          sendCont += 1

        }
        cnt += 1
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

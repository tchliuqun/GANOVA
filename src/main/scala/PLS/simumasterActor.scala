package PLS

import akka.actor._
import simumasterActor._
import simucalculateActor._
import myParallel.actorMessage._
import SnpProcessActor._
import breeze.linalg._
import myParallel._
import myParallel.paraWriterActor._

object simumasterActor{
  val name = "simumasterActor"
  def props(fil:Pms) = Props(classOf[simumasterActor],fil)
  case class Pms(fil:String,times:Int = 100,H:Array[Float]= Array(0.02f, 0.03f, 0.04f),k:Int = 3,func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] = (X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float]) => plsCalc.ngdofP(X,Y,k)._2.map(_.toString))
  case class coreNum(num:Int)
  case class chr(chrname:Array[String])}

class simumasterActor(pms:simumasterActor.Pms) extends Actor{

  var ofile = pms.fil
  var wname = "simuWriter"
  var count = 0
  var doneNum = 0
  var times = pms.times
  var H = pms.H
  var k = pms.k
  var cores = Runtime.getRuntime.availableProcessors()-1

  var glists:Array[Array[String]] = Array(Array(""))
  var order:Option[ActorRef] = None
  var writer:Option[ActorRef] = None
  var calculaters:Array[Option[ActorRef]] = Array(None)
  var ractor = system.actorOf(Ractor.props(Ractor.Pms(3)), "ractor")

  var calculiting : (DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] =pms.func//(X:DenseMatrix[Float],Y:DenseMatrix[Float],K:Int) => plsCalc.ngdofP(X,Y,2)._2.map(_.toString)
  var tlen = 0
  def getGlist(chr:String) {
    val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs2 = svd.filter(i => i(0) == chr & i(4).toInt > 10).sortBy(_ (4).toInt)
    // val glist = rs2.filter(i => i(4).toInt > 100 & i(5).toDouble > 0.80).flatten
    val glistsInx = Range(0, 120, 3).toArray ++ Array(120 until rs2.size: _*).slice(0,100)
    //val glistsInx = Array(150 until rs2.size: _*).slice(0,70)
    //val glistsInx = Range(0, 30, 3).toArray ++ Array(30 until rs2.size: _*)
    this.glists = glistsInx.map(rs2(_))
  }
  def receive = {
    case corenums:coreNum =>{
      this.cores = corenums.num
    }
    case func:calfunc => {
      this.calculiting = func.func
    }
    case chr:simumasterActor.chr =>{
      order = Some(sender)
      val wrt = system.actorOf(paraWriterActor.props(fileName(this.ofile)), wname)
      writer = Some(wrt)
      if (cores > 50) cores = 50
      println("starting writer")
      getGlist(chr.chrname.apply(0))
      this.tlen = glists.length
      //      wrt ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting 1")
      wrt ! myParallel.paraWriterActor.totalNumber(glists.length * times * H.length)
      //      writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting"))
      println("starting processing")
      if (glists.length < cores){
        cores = glists.length
        //sleep
      }
      val pms = simucalculateActor.Pms(wname,times,H)
      var ii  = 0
      while (ii < cores){
        val actr = system.actorOf(simucalculateActor.props(pms),"calc"+ii)
        calculaters :+= Some(actr)
        system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glists(count),"user/"+wname)), glists(count)(3))
        //system.actorSelection("/user/"+glists(count)).foreach(_ ! )
        actr !  simucalculateActor.geneLists(glists(count))
        Thread.sleep(myParallel.actorMessage.fs.length+100)
        count += 1
        println("processing No." + count)
        ii += 1
      }
    }
    case chr:SnpProcessActor.chr =>{
      order = Some(sender)
      val wrt = system.actorOf(paraWriterActor.props(fileName(this.ofile)), wname)
      writer = Some(wrt)
      if (cores > 40) cores = 40

      var plsWriter = system.actorOf(plssimuWriter.props(plssimuWriter.Pms(wname)), "plswriter")
      println("starting writer")
      getGlist(chr.chrname.apply(0))
      this.tlen = glists.length
      //      wrt ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting 1")
      wrt ! myParallel.paraWriterActor.totalNumber(glists.length * times * H.length)
      //      writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting"))
      println("starting processing")
      if (glists.length < cores){
        cores = glists.length
        //sleep
      }
      val pms = simucalculateActor.Pms(wname,times,H,this.k,true,this.calculiting)
      //Array(0 until cores:_*).foreach(i =>
      var ii  = 0
        while (ii < cores){
        val actr = system.actorOf(simucalculateActor.props(pms),"calc"+ii)
        calculaters :+= Some(actr)
          system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glists(count))), glists(count)(3))
          actr !  simucalculateActor.geneList(glists(count))
        Thread.sleep(myParallel.actorMessage.fs.length+100)
        count += 1
        println("processing No." + count)
          ii += 1
      }
    }
    case gList:simucalculateActor.gList =>{
      order = Some(sender)
      val wrt = system.actorOf(paraWriterActor.props(fileName(this.ofile)), wname)
      writer = Some(wrt)
      if (cores > 40) cores = 40
      this.tlen = cores * gList.n

      println("starting writer")
      //getGlist(chr.chrname.apply(0))
      //      wrt ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting 1")
      wrt ! myParallel.paraWriterActor.totalNumber(gList.glist.apply(4).toInt * times * H.length)
      //      writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting"))
      println("starting processing")
      //if (glists.length < cores){
      //  cores = glists.length
      //}
      system.actorOf(vegas2Actor.props(vegas2Actor.Pms(gList.glist)), gList.glist.apply(3))

      val pms = simucalculateActor.Pms(wname,times,H)
      Array(0 until cores:_*).foreach(i => {
        val actr = system.actorOf(simucalculateActor.props(pms),"calc"+i)
        calculaters :+= Some(actr)
        Thread.sleep(1000)
        actr !  gList
        count += 1
        println("processing No." + count)
      })
    }

    case don:done => {
      val r = don.count
      //doneNum += 1
      if (r == 0) {
        if (count < tlen) {
          system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glists(count),"user/"+wname)), glists(count)(3))
          sender ! simucalculateActor.geneLists(glists(count))
          println("processing No." + count)
        }
        else {
          sender ! PoisonPill
        }
      }else if (r == 1){
        if (count < tlen) {
          system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glists(count),"user/"+wname)), glists(count)(3))
          sender ! simucalculateActor.geneList(glists(count))
          println("processing No." + count)
        }
        else {
          sender ! PoisonPill
        }
      }
      count += 1
      //After all jobs is done, the count should equal to number of actors add length of gene list.
      if (count == tlen + cores){
        sender ! PoisonPill
        try {
          calculaters.foreach(_.foreach(_ ! PoisonPill))
          writer.foreach(_ ! PoisonPill)
          order.foreach(_ ! "done")
        }finally {
          self ! PoisonPill
        }
        }
      }

    case _ => println("the message is not correct -simumasterActor")
    }
  override def postStop {
    calculaters.foreach(_.foreach(_ ! PoisonPill))
    writer.foreach(_ ! PoisonPill)
    order.foreach(_ ! "done")
    println(utils.currentTimeIn+s"calculating is done -simumasterActor")

  }
}

package PLS

import PLS.SnpProcessActor.chr
import akka.actor._
import myParallel.actorMessage.{done, system}
import myParallel.paraWriterActor
import myParallel.paraWriterActor._
import simumasterActor._
object simumasterActor{
  val name = "simumasterActor"
  def props(fil:Pms) = Props(classOf[simumasterActor],fil)
  case class Pms(fil:String,times:Int = 100,H:Array[Float]= Array(0.01f, 0.015f, 0.02f))
  case class coreNum(num:Int)
}

class simumasterActor(pms:Pms) extends Actor{
  var ofile = pms.fil
  var wname = "simuWriter"
  var count = 0
  var doneNum = 0
  var times = pms.times
  var H = pms.H
  var cores = Runtime.getRuntime.availableProcessors()+1
  var glists:Array[Array[String]] = Array(Array(""))
  var order:Option[ActorRef] = None
  var writer:Option[ActorRef] = None
  var calculaters:Array[Option[ActorRef]] = Array(None)
  def getGlist(chr:String) {
    val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs2 = svd.filter(i => i(0) == chr & i(4).toInt > 10).sortBy(_ (4).toInt)
    val glistsInx = Range(0, 120, 3).toArray ++ Array(120 until rs2.size: _*).slice(0,100)
    //val glistsInx = Range(0, 30, 3).toArray ++ Array(30 until rs2.size: _*)
    this.glists = glistsInx.map(rs2(_))
  }
  def receive = {
    case corenums:coreNum =>{
      this.cores = corenums.num
    }
    case chr:chr =>{
      order = Some(sender)
      val wrt = system.actorOf(paraWriterActor.props(fileName(this.ofile)), wname)
      writer = Some(wrt)
      if (cores > 50) cores = 50

      println("starting writer")
      getGlist(chr.chrname.apply(0))
//      wrt ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting 1")
      wrt ! myParallel.paraWriterActor.totalNumber(glists.length * times * H.length)
//      writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("dispatch ssstarting"))
      println("starting processing")
      if (glists.length < cores){
        cores = glists.length
      }
      val pms = simucalculateActor.Pms(wname,times,H)
      Array(0 until cores:_*).foreach(i => {
        val actr = system.actorOf(simucalculateActor.props(pms),"calc"+i)
        calculaters :+= Some(actr)
        actr !  simucalculateActor.geneList(glists(count))
        count += 1
        println("processing No." + count)
      })
    }
    case don:done => {
      //val r = don.count
      //doneNum += 1
      if (count < glists.length) {
        sender ! simucalculateActor.geneList(glists(count))
        println("processing No." + count)
      }
      else {
        sender ! PoisonPill
      }
      count += 1
      //After all jobs is done, the count should equal to number of actors add length of gene list.
      if (count == glists.length + cores){
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

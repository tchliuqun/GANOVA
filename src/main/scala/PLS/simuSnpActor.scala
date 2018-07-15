package PLS

import akka.actor._
import simuSnpActor._
import simuSnpCalActor._
import myParallel.actorMessage._
//import SnpProcessActor._
import breeze.linalg._
import myParallel._
import myParallel.paraWriterActor._

object simuSnpActor{
  val name = "simuSnpActor"
  def props(fil:Pms) = Props(classOf[simuSnpActor],fil)
  case class Pms(fil:String,times:Int = 100,H:Array[Float]= Array(0.05f),k:Int = 3,n:Int = 500,n2:Int = 40)
  case class coreNum(num:Int)
  case class ns(n1s:Array[Int])
}

class simuSnpActor(pms:simuSnpActor.Pms) extends Actor {

  var ofile = pms.fil
  var n2 = pms.n2
  var n = pms.n
  var wname = "simuWriter"
  var count = 0
  var doneNum = 0
  var times = pms.times
  var H = pms.H
  var k = pms.k
  var cores = Runtime.getRuntime.availableProcessors() - 1

  var nis:Array[Int] = Array(0)
  var order: Option[ActorRef] = None
  var writer: Option[ActorRef] = None

  var tlen = 0
  var cpms = simuSnpCalActor.Pms(wname, times, H,k,n,n2)
  def receive = {
    case corenums: coreNum => {
      this.cores = corenums.num
    }
    case ns: simuSnpActor.ns => {
      nis = ns.n1s
      order = Some(sender)
      val wrt = system.actorOf(paraWriterActor.props(fileName(this.ofile)), wname)
      writer = Some(wrt)
      if (cores > 50) cores = 50

      println("starting writer")
      this.tlen = nis.length
      wrt ! myParallel.paraWriterActor.totalNumber(tlen * times * H.length)
      println("starting processing")
      if (tlen < cores) {
        cores = tlen
      }

      var ii = 0
      while (ii < cores) {
        val actr = system.actorOf(simuSnpCalActor.props(cpms), "calc" + ii)
        //calculaters :+= Some(actr)
        actr ! simuSnpCalActor.n(nis(count))
        Thread.sleep(myParallel.actorMessage.fs.length + 100)
        count += 1
        println("processing No." + count)
        ii += 1
      }
    }
    case don:done => {
      //val r = don.count
      //doneNum += 1
      if (count < tlen) {
        val actr = system.actorOf(simuSnpCalActor.props(cpms), "calc" + count)
        //calculaters :+= Some(actr)
        actr ! simuSnpCalActor.n(nis(count))
        Thread.sleep(myParallel.actorMessage.fs.length + 100)

        println("processing No." + count)

      }
      //else {
      count += 1
      sender ! PoisonPill
      //}
      //count += 1
      //After all jobs is done, the count should equal to number of actors add length of gene list.
      if (count == tlen + cores ){
      //  sender ! PoisonPill
        try {
          //   calculaters.foreach(_.foreach(_ ! PoisonPill))
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
    //calculaters.foreach(_.foreach(_ ! PoisonPill))
    writer.foreach(_ ! PoisonPill)
    order.foreach(_ ! "done")
    println(utils.currentTimeIn+s"calculating is done -simumasterActor")

  }
}
package PLS

import scala.collection.mutable.ArrayBuffer
import scala.concurrent.duration._
import scala.util.{Failure, Success}
//import akka.actor.Status.Success
import akka.actor._
import breeze.linalg._
import myParallel.actorMessage._
import simuSnpCalActor._
import akka.pattern.ask
import scala.concurrent.Future
import myParallel.paraWriterActor._
import java.io.{FileWriter, PrintWriter}

object simuSnpCalActor{
  val name = "simucalculateActor"
  def props(pms:Pms) = Props(classOf[simuSnpCalActor],pms)
  def deffun(z:Array[Float])(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int) = plsCalc.ngdofP(X,Y,k)._2.map(_.toString)
  case class Pms(fil:String,times:Int = 1000,H:Array[Float] = Array(0.01f, 0.015f, 0.02f),k:Int = 3,n:Int,n2:Int)
  case class n(n1:Int)
  //case class
}


class simuSnpCalActor(pms:Pms) extends Actor {
  var writer: Option[ActorSelection] = Some(system.actorSelection("/user/" + pms.fil))
  var H = pms.H
  var times = pms.times
  var k = pms.k
  var n = pms.n
  var n2 = pms.n2
  var gen: String = "ENSG" + (200000 + scala.util.Random.nextInt(100000)).toString

  def receive = {
    case ns: simuSnpCalActor.n => {
      val n1 = ns.n1
      val (x, x1, glist) = vegas2.simuVegasSnp(n, n1, n2, gen)
      for (h <- H) {
        var j = 0
        while (j < times) {
          val Y = vegas2.setPheno(h, 0, false)(new DenseMatrix(n, 1, x1))
          val pval = plsCalc.ngdofP(x, Y,k)._2
          val vp = vegas2.vegasP(glist, Y)
          val prs = glist ++ (pval ++ vp :+ n1.toFloat :+ h).map(_.toString)
          writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(prs.mkString("\t")))
          j += 1
        }
      }
      sender ! done(0)
    }
  }
}


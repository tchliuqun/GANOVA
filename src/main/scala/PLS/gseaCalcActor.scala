package PLS

//import PLS.simumasterActor.done

//import akka.actor.Status.Success
import akka.actor._
import myParallel.actorMessage._
import gseaCalcActor._

import breeze.linalg._
//import scala.math._
//import scala.math._
import breeze.stats._
import breeze.math._
import breeze.numerics._
import scala.concurrent.Future
object gseaCalcActor{
  val name = "gseaCalcActor"
  def props(pms:Pms) = Props(classOf[gseaCalcActor],pms)
  case class Pms(fil:String = gPms.rp+"gseaResult.txt")
  case class goresults(data:DenseMatrix[Float])
  case class results(data:DenseVector[Float])
  case class geneset(descript:Array[String],data:DenseVector[Float])
  case class gList(glist:Array[String],n:Int = 2 )
  //case class
}


class gseaCalcActor(pms:Pms) extends Actor{
  var res:DenseMatrix[Float] = null
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/gseawriter"))
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

    val results = accumulate(res(::,*))
    val mins = min(results(::,*))
    val maxs = max(results(::,*))
    for ( i <- 0 until gcol) {
      rst(i) = if (abs(mins(i)) > maxs(i)) mins(i) else maxs(i)
    }
    val es = rst(0)
    val pval = if(es>0) rst.filter(_ > es).length.toFloat / (gcol-1).toFloat else rst.filter(_ < es).length.toFloat / (gcol-1).toFloat
    es +: pval +: results(::,0).toArray
  }

  def receive = {
    case rs:goresults =>{
      this.res = rs.data
    }
    case gs:geneset =>{
      val des = gs.descript
      val rss = ogsea(gs.data,res)
      val rsstr = des.mkString("\t")+"\t"+rss.mkString("\t")
      writer.foreach(_ !myParallel.paraWriterActor.WriteStr(rsstr))
      sender ! done
    }
    case _ => println("Someone say wrong to me!!!")
  }
}

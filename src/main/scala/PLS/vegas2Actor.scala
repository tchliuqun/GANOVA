package PLS

//import PLS.simumasterActor.done
import PLS.snpCalcActor.writerName
//import PLS.vegas2.vegasX
import akka.actor._
import breeze.linalg.DenseMatrix
import myParallel.actorMessage._
import vegas2Actor._
import simucalculateActor._

object vegas2Actor{
  val name = "vegas2Actor"
  def props(pms:Pms) = Props(classOf[vegas2Actor],pms)
  case class Pms(glist:Array[String])//,spheno:DenseMatrix[Float] => DenseMatrix[Float] = vegas2.setPheno())
  case class geneList(glist:Array[String])
  case class gList(glist:Array[String],n:Int = 2 )
  case class inp(inx:String,Y:DenseMatrix[Float])

  //case class
}


class vegas2Actor(pms:vegas2Actor.Pms) extends Actor{
  var simuwriter:Option[ActorSelection] = Some(system.actorSelection("/user/plswriter"))
  //var Y:DenseMatrix[Float] = DenseMatrix.zeros[Float](1,1)

  //var pval:Array[Float] = Array(0f)
  val glist = pms.glist
  //val spheno = pms.spheno
  //val X = vegas2.vegasX(glist)
//  def vegasY = {
//    val Y = spheno(X)
//    val pval = vegas2.vegasP(glist,Y)
//    (Y,pval)
//  }
  def receive = {
    case inp:inp =>{
      //val Y = spheno(X)
      //sender ! (X,Y)
      val pval = vegas2.vegasP(glist,inp.Y)
      simuwriter.foreach(_ ! simucalculateActor.permp(inp.inx,pval))
    }

  }
}

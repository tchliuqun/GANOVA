package PLS

import breeze.linalg._
import akka.actor._
import myParallel.actorMessage._
import Ractor._
import simucalculateActor._
object Ractor{
  val name = "vegas2Actor"
  def props(pms:Pms) = Props(classOf[Ractor],pms)
  case class Pms(k:Int)//,spheno:DenseMatrix[Float] => DenseMatrix[Float] = vegas2.setPheno())
  //case class geneList(glist:Array[String])
  //case class gList(glist:Array[String],n:Int = 2 )
  case class inp(inx:String,dt:Any)
  //case class
}

class Ractor(pms:Ractor.Pms) extends Actor{
  var simuwriter:Option[ActorSelection] = Some(system.actorSelection("/user/plswriter"))
  //var Y:DenseMatrix[Float] = DenseMatrix.zeros[Float](1,1)
  val R = org.ddahl.rscala.RClient()
  var getRs:inp => Any = (ip:inp) => {
    val (x,y,k) = ip.dt.asInstanceOf[(DenseMatrix[Float],DenseMatrix[Float],Int)]
    plsCalc.kramerDof(x,y,k,R).drop(1)
  }
  //val spheno = pms.spheno
  //  }
  def receive = {

    case inp:inp =>{
      println("rrrrrrrrrrrrrrrrrrrrr")
      val pval = getRs(inp).asInstanceOf[Array[Float]]

      println("RRRRRRRRRRRRRRRRRRRRRRR")
      simuwriter.foreach(_ ! simucalculateActor.dof(inp.inx,pval))
    }

  }
}

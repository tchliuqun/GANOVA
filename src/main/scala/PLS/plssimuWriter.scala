
package PLS

import breeze.linalg._
import akka.actor._
//import Ractor._
import myParallel.actorMessage._
import simucalculateActor._
object plssimuWriter{
  val name = "plsSimuActor"
  def props(pms:Pms) = Props(classOf[plssimuWriter],pms)
  case class Pms(fil:String)//,spheno:DenseMatrix[Float] => DenseMatrix[Float] = vegas2.setPheno())
  //case class geneList(glist:Array[String])
  //case class gList(glist:Array[String],n:Int = 2 )
  case class inp(inx:String,dt:Any)
  //case class
}

class plssimuWriter(pms:plssimuWriter.Pms) extends Actor {
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var rsm:Map[String,Array[String]] = Map()
  var rsy:Map[String,(Array[String],DenseMatrix[Float],DenseMatrix[Float],Array[Float],Array[Float],Array[Float])] = Map()

  def receive = {

    case ry:rsy =>{
      println("PPPPPPPPPPPPPPPPPPPP")
      rsy += (ry.idx-> (ry.glt, ry.yy, ry.yh, ry.pdofl, ry.gdof, ry.permp))
    }
    case df:dof => {
//      println("")
//      println("f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")

      val (glt, yy, yh, pdofl, gdof, permp) = rsy(df.idx)
      rsy -= df.idx
      val ppval = plsCalc.dofPval(yy, yh, pdofl)
      val gpval = plsCalc.dofPval(yy, yh, gdof)
      val kpval = plsCalc.dofPval(yy, yh, df.dt)
      val rss = glt ++(permp ++ df.dt ++ kpval ++ pdofl ++ ppval ++ gdof ++ gpval).map(_.toString)
      if (rsm.contains(df.idx)) {
        val vgsr = rsm(df.idx)
        val rs = (rss ++ vgsr).mkString("\t") +"\t"+ df.idx.split("_").apply(2)
        rsm -= df.idx
        println("f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
      } else {
        rsm += (df.idx -> rss)
      }
    }
    case pp:permp =>{
//      println("")
//      println("f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111")

      if (rsm.contains(pp.idx)) {
        val ppr = rsm(pp.idx)
        val rs = (ppr ++ pp.dt.map(_.toString)).mkString("\t")+"\t" + pp.idx.split("_").apply(2)
        rsm -= pp.idx
        println("f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
      } else {
        rsm += (pp.idx -> pp.dt.map(_.toString))
      }
    }

  }

}

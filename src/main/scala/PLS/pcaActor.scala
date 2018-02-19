package PLS

import PLS.pcaActor.Xs
import PLS.snpCalcActor.writerName
import breeze.linalg._
import akka.actor._
import myParallel.actorMessage._
import breeze.numerics.sqrt
object pcaActor{
  val name = "snpCalcActor"
  case class Xs(gene:Array[String],X:DenseMatrix[Float])
  case class writerName(name:String)
}
class pcaActor extends Actor{
  var n = 5
    def pcaPvar(X:DenseMatrix[Float]):Array[Float] = {
      val rss = princomp(convert(X,Double)).propvar.toArray.map(_.toFloat)
      val rs = if (rss.length >= n) rss.slice(0, n) else rss++Array.fill(n-rss.length)(0f)
      rs
    }
  var order:Option[ActorRef] = None
  var writer:Option[ActorSelection] = None

    def receive = {
      case wrt:writerName => {
        order = Some(sender)
        writer = Some(system.actorSelection("/user/"+wrt.name))
      }
      case num:Int =>{
        n = num
      }


      case x:Xs =>{

        val X = x.X
        val grs = pcaPvar(X).mkString("\t")
        var rs = x.gene.mkString("\t") + "\t"+grs
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
        sender ! done(1)

      }
      case don:done => order.foreach(_ ! finished)
      case _ =>
    }


}

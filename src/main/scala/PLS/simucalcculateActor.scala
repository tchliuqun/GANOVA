package PLS

//import PLS.simumasterActor.done
import PLS.snpCalcActor.writerName
import akka.actor._
import myParallel.actorMessage._//system
import simucalculateActor._
object simucalculateActor{
  val name = "simucalculateActor"
  def props(pms:Pms) = Props(classOf[simucalculateActor],pms)
  case class Pms(fil:String,times:Int = 100,H:Array[Float] = Array(0.01f, 0.015f, 0.02f))
  case class geneList(glist:Array[String])
  //case class
}


class simucalculateActor(pms:Pms) extends Actor{
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var H = pms.H
  var times = pms.times
  def simugenNo(glists:Array[String]) = {
    val glist = glists.slice(0, 4)
    //println("processing No."+g)
    vegas2.simuFgene(glist)
    for (h <- H) {
      var i = 0
      while (i < times) {
        val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
        i += 1
      }
    }
  }
  def receive = {
    case wrt:writerName => {
      writer = Some(system.actorSelection("/user/"+wrt.name))
    }
    case gList:geneList =>{
      simugenNo(gList.glist)
      sender ! done(0)
    }
    case don:done => sender ! done(0)
  }

}

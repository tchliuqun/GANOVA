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
  case class gList(glist:Array[String],n:Int = 2 )
  //case class
}


class simucalculateActor(pms:Pms) extends Actor{
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var H = pms.H
  var times = pms.times
  def simugenNo(glists:Array[String]) = {
    val glist = glists.slice(0, 4)
    //println("processing No."+g)
//    vegas2.simuFgene(glist)
    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
//    writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("calculation ssstarting"))

    if(rl > 0) {
      for (h <- H) {
        var i = 0
        while (i < times) {
          //val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
          val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPhenoT(h, 0,0.5f)) :+ h).mkString("\t")
            writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
            i += 1
        }
      }
    }
  }
  def simugenNo1(glists:Array[String],n:Int) = {
    import java.io.{FileWriter, PrintWriter}
    val glist = glists.slice(0, 4)
    //val writer = new PrintWriter(new FileWriter("goR"+glist(3)+".txt"))
    //println("processing No."+g)
    vegas2.simuFgene(glist)
    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
    //val sl = glist(4)
    //    writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("calculation ssstarting"))

    for (h <- H) {
      var i = 0
      while (i < rl) {
        var j = 0
        while(j < n) {
          //val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
          val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno(h, i, false)) :+ h :+ i).mkString("\t")
          writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
          //writer.println(rs) //foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
          j += 1
        }
        i += 1
      }
    }
    //writer.close()
  }

  def receive = {
    case wrt:writerName => {
      writer = Some(system.actorSelection("/user/"+wrt.name))
    }
    case gList:geneList =>{
      simugenNo(gList.glist)
      sender ! done(0)
    }
    case gList:gList =>{
      simugenNo1(gList.glist,gList.n)
      sender ! done(0)
    }
    case don:done => sender ! done(0)
  }

}

package myParallel


  import java.io.BufferedWriter

  import PLS.utils
  import akka.actor._
  //import com.typesafe.akka.actor._
  import java.io.BufferedWriter
  import java.io.FileWriter

  import myParallel.BufferWriterActor.WriteToBuffer

  /**
    * Created by liuqun on 12/10/16.
    */
  object BufferWriterActor {
    val name = "BufferedWriterActor"
    def props(fn: String) = Props(classOf[BufferWriterActor],fn)
    case class WriteToBuffer(str: String)
  }

  class BufferWriterActor(fn: String) extends Actor {
    val fw:FileWriter  = new FileWriter(fn)
    val bw:BufferedWriter  = new BufferedWriter(fw)
    var orderWorker:Option[ActorRef] = None

    def receive: Actor.Receive = {
      case WriteToBuffer(str) =>
        //      println(s"writing $str")
        bw.write(str)
        bw.newLine()
      case actorMessage.finished => {
        //orderWorker = Some(sender)
        //sender ! actorMessage.done(1)
        //bw.close()
        self ! PoisonPill
      }
      case _ => println(utils.currentTimeIn+"Someone said wrong to me. -bufferedWriterActor") // Send messages to this actor from all you threads.
    }
    override def postStop: Unit = {
      bw.close()
      println(utils.currentTimeIn+"writing is done -- BufferedWriter")
    }

}

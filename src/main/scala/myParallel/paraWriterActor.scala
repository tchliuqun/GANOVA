package myParallel


  import java.io.{BufferedWriter, File, FileWriter}

  import PLS.run.currentTime
  import PLS.utils
  import akka.actor.{Actor, ActorRef, PoisonPill, Props}
  import myParallel.paraWriterActor._
  import myParallel.actorMessage._

  /**
    * Created by liuqun on 12/11/16.
    */
  object paraWriterActor {
    val name = "paraWriterActor"
    def props(fileName:fileName) = Props(classOf[paraWriterActor], fileName)
    case class fileName(fil:String)
    case class WriteStr(str: String)
    case class totalNumber(num:Int)
    //  case class count(num:Int)
    //  case class done(count:Int)
    //  case object finished
  }
  class paraWriterActor(fileN:fileName = fileName("Buffered.txt")) extends Actor {

    var fileName = fileN.fil
    //  val system = ActorSystem("mySystem")
    //var fileName = "Buffered.txt"
    val fnL = fileName.split("\\.")//txt", utils.currentTime + ".txt")
    val profix = "."+fnL(fnL.length -1)
    val fn  = fileName.replace(profix,"_"+utils.getTimeForFile+profix)

    val file = new File(fn)
    file.setWritable(true, false)
    val fw:FileWriter  = new FileWriter(file,true)
    val bw:BufferedWriter  = new BufferedWriter(fw)
    var totalNum = 0
    var count = 0
    var ifCount = false
    var orderWorker:Option[ActorRef] = None

    //  try {
    // Share this actor across all your threads.
    //val myActor = system.actorOf(BufferWriterActor.props(fn), paraWriterActor.name)

    def receive = {
      case paraWriterActor.WriteStr(str) => {
        count += 1
        bw.write(str)
        bw.newLine()
        //myActor ! BufferWriterActor.WriteToBuffer(str)
        if (ifCount & count == totalNum) {
          sender ! done(count)
      //    myActor ! PoisonPill// actorMessage.finished
          self ! PoisonPill
        }

      } //.info(s"I was greeted by $greeter.")
      case paraWriterActor.totalNumber(i) => {
        var orderWorker:Option[ActorRef] = Some(sender)
        println(utils.currentTimeIn+s"get total record numbers $i to write")
        bw.write(utils.currentTimeIn+s"get total record numbers $i to write")
        bw.newLine()
        totalNum += i
        ifCount = true
      }
      case actorMessage.finished => {
        println(utils.currentTimeIn+s"total writed record numbers $count")

     //   orderWorker = Some(sender)
      //  myActor ! PoisonPill// actorMessage.finished
        self ! PoisonPill
      }
      case actorMessage.done => {
        //orderWorker.foreach(_!actorMessage.done)
 //       myActor ! PoisonPill// actorMessage.finished
        self ! PoisonPill
      }
      case _ => println(utils.currentTimeIn+"Someone said wrong to me. - paraWriterActor") // Send messages to this actor from all you threads.
    }
    override def postStop {
      bw.close()
      println(utils.currentTimeIn+s"writing to $fileName is done - paraWriter")

    }
    //      myActor ! riteToBuffer("The Text")
    //  }finally {


}

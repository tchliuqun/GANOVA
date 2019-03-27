package myParallel


  import java.io.{BufferedWriter, File, FileWriter}

  import PLS.run.currentTime
  import PLS.simucalculateActor.rsy
  import PLS.utils
  import akka.actor.{Actor, ActorRef, PoisonPill, Props}
  import breeze.linalg.DenseMatrix
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
    //case class rsy(idx:String,glt:Array[String],yy:DenseMatrix[Float], yh:DenseMatrix[Float],pdofl:Array[Float],gdof:Array[Float],permp:Array[Float])
    case class strBf1(idx:String = "",rs:String = "")//,s3:String = "",s4:String= "")
    case class strBf2(idx:String = "",rs:String = "")
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

    var rsm = scala.collection.mutable.Map[String,(String,String)]()
    var totalNum = 0
    var count = 0
    var ifCount = false
    var orderWorker:Option[ActorRef] = None

    val file = new File(fn)
    file.setWritable(true, false)
    //    val fw:FileWriter  =
    val bw:BufferedWriter  = new BufferedWriter(new FileWriter(file,true))

    import java.io.IOException

    Runtime.getRuntime.addShutdownHook(new Thread() {
      override def run(): Unit = {
        //System.out.println("run shutdownHook")
        // save state, resource clean,up etc.
        if (bw != null) try { // try to close the open file
          //bw.flush()
          bw.close()
          System.out.println(fn+" closed successfully")
        } catch {
          case e: IOException =>
            System.out.println("Failed to flush/close"+fn +": " + e.getMessage)
            e.printStackTrace()
        }
        //System.out.println("Shutdown hook completed...")
      }
    })
    //  try {
    // Share this actor across all your threads.
    //val myActor = system.actorOf(BufferWriterActor.props(fn), paraWriterActor.name)

    def receive = {
      case ry:strBf1 =>{
        //println("PPPPPPPPPPPPPPPPPPPP")
        if(rsm.contains(ry.idx)){
          bw.write(ry.rs + "\t"+rsm(ry.idx)._2)
          bw.newLine()
          rsm -= ry.idx
        }else{
          rsm += (ry.idx-> (ry.rs,""))
        }// += (ry.idx-> ry.rs)
      }
      case ry:strBf2 =>{
        //println("PPPPPPPPPPPPPPPPPPPP")
        if(rsm.contains(ry.idx)){
          bw.write(rsm(ry.idx)._1 + "\t" + ry.rs)
          bw.newLine()
          rsm -= ry.idx
        }else{
          rsm += (ry.idx-> ("",ry.rs))
        }// += (ry.idx-> ry.rs)
      }
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
        //println("write to "+fn)
        //bw.write(utils.currentTimeIn+s"get total record numbers $i to write")
        //bw.newLine()
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
    override def preStart(): Unit = {
      if(file.canWrite()) {
        println("writing to "+fn)
        //bw.write(utils.currentTimeIn+"start write")
        //bw.newLine()
        // write access
      } else {
        println("can't write to "+fn)
        // no write access
      }
    }
    override def postStop {
      bw.close()
      println(utils.currentTimeIn+s"writing to $fileName is done - paraWriter")

    }
    //      myActor ! riteToBuffer("The Text")
    //  }finally {


}

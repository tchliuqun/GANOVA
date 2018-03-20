package myParallel


  import java.nio.file.Paths

  import akka.actor.{ActorRef, ActorSystem}
  import akka.stream.{ActorMaterializer, IOResult}
  import akka.util.Timeout

  import scala.concurrent.duration._
  import akka.dispatch.ExecutionContexts._
  import akka.http.scaladsl.model.HttpResponse
  import akka.stream.scaladsl._

  import scala.concurrent.Future
  /**
    * Created by liuqun on 12/18/16.
    */
  object actorMessage {

    implicit val system = ActorSystem("MySystem")
    implicit val materializer = ActorMaterializer()
    implicit val ec = global
    case class host(address:String,https:Boolean = true)
    case class msgActor(actor:Option[ActorRef])
    implicit val timeout = Timeout(25000.seconds)
    val fs:FiniteDuration = (100).millis
    case class result(output:String)
    case class MysqlSetup(url:String,usr:String,pwd:String)
    case class throttle(num:Int,time:FiniteDuration)
    case class httpFail(idx:Int)
    case class httpSuccess(ind:Int,output:HttpResponse)
    case class successed(ind:Int,out:String)
    case class retryTimes(num:Int)
    case object action
    case object finished
    case object wrong
    case object busy
    case class outActor(actorName:String)
    case class done(count:Int=1)

}

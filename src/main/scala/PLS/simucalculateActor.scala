package PLS

import scala.collection.mutable.ArrayBuffer
import scala.concurrent.duration._
import scala.util.{Failure, Success}
//import akka.actor.Status.Success
import akka.actor._
import breeze.linalg._
import myParallel.actorMessage._
import simucalculateActor._
import akka.pattern.ask
import scala.concurrent.Future
import myParallel.paraWriterActor._

object simucalculateActor{
  val name = "simucalculateActor"
  def props(pms:Pms) = Props(classOf[simucalculateActor],pms)
  def deffun(z:Array[Float])(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int) = plsCalc.ngdofP(X,Y,k)._2.map(_.toString)
  case class Pms(fil:String,times:Int = 100,H:Array[Float] = Array(0.01f, 0.015f, 0.02f),k:Int = 3,rscala:Boolean = true,func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] =(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float]) =>  plsCalc.ngdofP(X,Y,k)._2.map(_.toString))
  case class geneList(glist:Array[String])
  case class gList(glist:Array[String],n:Int = 2, func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String]=(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float]) =>  plsCalc.ngdofP(X,Y,k)._2.map(_.toString))// = xx:Array[Float]=>(xyk:[(DenseMatrix[Float],DenseMatrix[Float],Int)] => plsCalc.ngdofP(xyk._1,xyk._2,xyk._3)._2.map(_.toString)))
  case class calfunc(func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String])
  //case class
}


class simucalculateActor(pms:Pms) extends Actor{
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  //var vegasA:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var H = pms.H
  var times = pms.times
  var k = pms.k
  var rscala = pms.rscala
  var calculiting : (DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] = pms.func//(X:DenseMatrix[Float],Y:DenseMatrix[Float],K:Int) => plsCalc.ngdofP(X,Y,2)._2.map(_.toString)

  def simugenNo(glists:Array[String]) = {
    val glist = glists.slice(0, 4)
    var rsm:Map[String,Array[String]] = Map()
    //val vg:ActorSelection = system.actorSelection("/user/"+glist(3))
    import java.util.concurrent.TimeUnit
    //val t = 1, TimeUnit.SECONDS)
    //val fs:FiniteDuration = (100).millis
    //val ts = vg.resolveOne(fs).value

    //if (ts.isDefined & ts.get.isFailure){
    //  system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glist)), glist(3))
    //}
    val vgs = system.actorSelection("/user/"+glist(3))
    val file = new java.io.File(gPms.tp+glist(3)+".gen")
    if(!file.exists() || file.length() == 0) vegas2.simuFgene(glist)

    //    implicit val timeout = 5000 // Timeout for the resolveOne call
//    system.actorSelection(glist(3)).resolveOne().onComplete {
//      case Success(actor) => actor ! message
//
//      case Failure(ex) =>
//        val actor = system.actorOf(Props(classOf[ActorClass]), name)
//        actor ! message
//    }
    //    vegas2.simuFgene(glist)
    //val vegas2a = system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glist)), glist(3)+utils.getTimeForFile)
    //println("processing No."+g)

    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
//    writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("calculation ssstarting"))
    val X = vegas2.vegasX(glist)
    if(rl > 0) {
      for (h <- H) {
        var i = 0
        while (i < times) {
          val Ys = vegas2.setPhenoT(h,0,0.5f)(X)
          val Y = calculation.standardization(Ys)
          val sr = i+"_"+h

          val future2: Future[(String,Array[Float])] = ask(vgs,vegas2Actor.inp(sr, Ys)).mapTo[(String,Array[Float])]

          //rsm += (sr -> plsP)
          future2 onComplete{
            case Success(f) =>{
              val rs = (glists ++f._2++ rsm(f._1) :+f._1.split("_").apply(1)).mkString("\t")
              writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
              i += 1
            }
            case Failure(t) => println("An error has occured: " + t.getMessage)
          }

          //val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
//          val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPhenoT(h, 0,0.5f)) :+ h).mkString("\t")
//            writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
//            i += 1
        }
      }
    }
  }
  def simugenNo1(glists:Array[String],n:Int) = {
    import java.io.{FileWriter, PrintWriter}
    val glist = glists.slice(0, 4)

    var rsm:Map[String,Array[String]] = Map()
    val vgs:ActorSelection = system.actorSelection("/user/"+glist(3))
    val ractor:ActorSelection = system.actorSelection("/user/ractor")
//    import java.util.concurrent.TimeUnit
//    //val t = 1, TimeUnit.SECONDS)
//    val fs:FiniteDuration = (100).millis
//    val ts = vgs.resolveOne(fs).value.get
//
//    if (ts.isFailure){
//      system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glist)), glist(3))
//      vgs = system.actorSelection("/user/"+glist(3))
//    }

//    var vgs:Option[ActorSelection] = Some(system.actorSelection("/user/"+glist(3)))
//    if (vgs.isEmpty){
//      system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glist)), glist(3))
//      vgs = Some(system.actorSelection("/user/"+glist(3)))
//    }

    val file = new java.io.File(gPms.tp+glist(3)+".gen")
    if(!file.exists() || file.length() == 0) vegas2.simuFgene(glist)
    //val vegas2a = system.actorOf(vegas2Actor.props(vegas2Actor.Pms(glist)), glist(3)+utils.getTimeForFile)
    //val writer = new PrintWriter(new FileWriter("goR"+glist(3)+".txt"))

    //println("processing No."+g)

    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
    //val sl = glist(4)
    //    writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("calculation ssstarting"))
    val X = vegas2.vegasX(glist)
    for (h <- H) {
      var i = 0
      while (i < rl) {
        var j = 0
        while(j < n) {
//          if(false) {
          val Ys = vegas2.setPheno(h, i, false)(X)
          val Y = calculation.standardization(Ys)
          val permp = Array(1 to k:_*).map(i => plsCalc.plsPerm(X,Y,i,10000))
          val sr = j + "_" + h + "_" + i
          val future2: Future[(String, Array[Float])] = ask(vgs, vegas2Actor.inp(sr, Y)).mapTo[(String, Array[Float])]

//          if(rscala) {
          val future1: Future[(String,Array[Float])] = ask(ractor,Ractor.inp(sr, (X,Y,k))).mapTo[(String,Array[Float])]
          val (yy,yh,pdofl,gdof) = plsCalc.dofPvalF(X,Ys,k,1000,true)
          //val gdof = ngdof(X,Ys,k,nPerm)
//          val kdof = kramerDof(X,Ys,k).drop(1).map(_.toFloat)
//          val kpval =dofPval(yy,yh,kdof)
          val ppval = plsCalc.dofPval(yy,yh,pdofl)
          val gpval = plsCalc.dofPval(yy,yh,gdof)
//          val (yup,ydn,ny,pdofl,gdof) = dofPvalF(X,Y,k)
//            val kpval =fdpval(yup,ydn,ny,kdof)
//            val ppval = fdpval(yup,ydn,ny,pdofl)
//            val gpval =fdpval(yup,ydn,ny,gdof)
            val future10 = future1.map{i => {
              val kpval = plsCalc.dofPval(yy,yh,i._2)
              val rss = permp ++ i._2 ++ kpval ++ pdofl ++ ppval ++ gdof ++ gpval
              rsm += (i._1 -> rss.map(_.toString))
            }
//              rsm += (i._1 -> calculiting(X,Y,k,i._2).map(_.mkString("\t")))
            }
//          }else {
//            val plsP = calculiting(X, Y, k,Array(0f))
//            rsm += (sr -> plsP)
//          }

            for {
              f2 <- future2
              f1 <- future10
            }yield{
//                case Success(f) => {
                  val rs = (glists ++ f2._2.map(_.toString) ++ rsm(f2._1) :+ f2._1.split("_").apply(1)).mkString("\t")
                  writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
                  rsm -= f2._1
                  j += 1
                }
//                case Failure(t) => println("An error has occured: " + t.getMessage)
//              }
              //          }
              //val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t")
              if (false) {
                val rs = (glists ++ vegas2.vegas(glist, 3, vegas2.setPheno(h, i, false)) :+ h :+ i).mkString("\t")
                writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))

                //writer.println(rs) //foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
                j += 1
              }
            }
        i += 1
        }

      }
    }
    //writer.close()
//}

  def receive = {
//    case wrt:writerName => {
//      writer = Some(system.actorSelection("/user/"+wrt.name))
//    }
    case gList:geneList =>{
      simugenNo(gList.glist)
      sender ! done(0)
    }
    case gList:gList =>{
      this.calculiting = gList.func
      simugenNo1(gList.glist,gList.n)
      sender ! done(0)
    }
    case func:calfunc => {
      this.calculiting = func.func
    }
    case don:done => sender ! done(0)
  }
}

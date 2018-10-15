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
import java.io.{FileWriter, PrintWriter}

object simucalculateActor{
  val name = "simucalculateActor"
  def props(pms:Pms) = Props(classOf[simucalculateActor],pms)
  def deffun(z:Array[Float])(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int) = plsCalc.ngdofP(X,Y,k)._2.map(_.toString)
  case class Pms(fil:String,times:Int = 1000,H:Array[Float] = Array(0.01f, 0.015f, 0.02f),k:Int = 3,rscala:Boolean = true,func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] =(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float]) =>  plsCalc.ngdofP(X,Y,k)._2.map(_.toString))
  case class geneList(glist:Array[String])
  case class geneLists(glist:Array[String])
  case class gList(glist:Array[String],n:Int = 2, func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String]=(X:DenseMatrix[Float],Y:DenseMatrix[Float],k:Int,dof:Array[Float]) =>  plsCalc.ngdofP(X,Y,k)._2.map(_.toString))// = xx:Array[Float]=>(xyk:[(DenseMatrix[Float],DenseMatrix[Float],Int)] => plsCalc.ngdofP(xyk._1,xyk._2,xyk._3)._2.map(_.toString)))
  case class calfunc(func:(DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String])
  case class dof(idx:String,dt:Array[Float])
  case class permp(idx:String,dt:Array[Float])
  case class rsy(idx:String,glt:Array[String],yy:DenseMatrix[Float], yh:DenseMatrix[Float],pdofl:Array[Float],gdof:Array[Float],permp:Array[Float])
  //case class
}


class simucalculateActor(pms:Pms) extends Actor{
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var simuwriter:Option[ActorSelection] = Some(system.actorSelection("/user/plswriter"))
  //var vegasA:Option[ActorSelection] = Some(system.actorSelection("/user/"+pms.fil))
  var H = pms.H
  var times = pms.times
  var k = pms.k
  var rscala = pms.rscala
  var n  = 5
  var m  = 5
  var calculiting : (DenseMatrix[Float],DenseMatrix[Float],Int,Array[Float]) => Array[String] = pms.func//(X:DenseMatrix[Float],Y:DenseMatrix[Float],K:Int) => plsCalc.ngdofP(X,Y,2)._2.map(_.toString)
  val ractor:ActorSelection = system.actorSelection("/user/ractor")
  var rsm:Map[String,Array[String]] = Map()
//  var rsy:Map[String,(Array[String],DenseMatrix[Float],DenseMatrix[Float],Array[Float],Array[Float],Array[Float])] = Map()

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
          val sr = i+"_"+h
          val future2: Future[(String,Array[Float])] = ask(vgs,vegas2Actor.inp(sr, Ys)).mapTo[(String,Array[Float])]
          val Ypm = calculation.permY(Ys.toDenseVector,1)
          val Yp = DenseMatrix.horzcat(Ypm(::,1 until Ypm.cols),Ys)

          //val Ys = calculation.standardization(Yp)
          val Y = calculation.standardization(Yp)


          val plsP = plsCalc.gdofPlsPval(X, Y, 3)._2
          rsm += (sr -> plsP.map(_.toString))
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
          if(false) {
            val Y = vegas2.setPheno(h, i, false)(X)
            val sr = j + "_" + h + "\t" + i

            val future2: Future[(String, Array[Float])] = ask(vgs, vegas2Actor.inp(sr, Y)).mapTo[(String, Array[Float])]
            val plsP = plsCalc.gdofPlsPval(X, Y, 2)._2
            rsm += (sr -> plsP.map(_.toString))
            future2 onComplete {
              case Success(f) => {
                val rs = (glists ++ f._2.map(_.toString) ++ rsm(f._1) :+ f._1.split("_").apply(1)).mkString("\t")
                writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
                rsm -= f._1
                j += 1
              }
              case Failure(t) => println("An error has occured: " + t.getMessage)
            }
          }
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

    //writer.close()
//}
  def simugenNo2(glists:Array[String]) = {
    val glist = glists.slice(0, 4)
    var rsm:Map[String,Array[String]] = Map()
    import java.util.concurrent.TimeUnit
    val vgs = system.actorSelection("/user/"+glist(3))
    val file = new java.io.File(gPms.tp+glist(3)+".gen")
    if(!file.exists() || file.length() == 0) vegas2.simuFgene(glist)

    val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
    val X = vegas2.vegasX(glist)

    if(rl > 0) {
      for (h <- H) {
        var i = 0
        while (i < times) {
          val Ys = vegas2.setPhenoT(h,0,0.5f)(X)
          val sr = i+"_"+h
          var ii = 0
          var rs  = Array("")
          while (ii < n){
            val Ypm = calculation.permY(Ys.toDenseVector,ii)
            val Yp = DenseMatrix.horzcat(Ypm(::,1 until Ypm.cols),Ys)
            val Y = calculation.standardization(Yp)
            val plsP = plsCalc.gdofPlsPval(X, Y, 3)._2
            rs ++= plsP.map(_.toString)
            ii += 1
          }
          var iii = 0
          while (iii < m){
            val Ypm = vegas2.setPhenoT(h,0,0.5f)(X)
            var Yp = DenseMatrix.horzcat(Ypm,Ys)
            while(Yp.cols < m -1){
              val Ypm = vegas2.setPhenoT(h,0,0.5f)(X)
              Yp = DenseMatrix.horzcat(Ypm,Yp)
            }
            val Y = calculation.standardization(Yp)
            val plsP = plsCalc.gdofPlsPval(X, Y, 3)._2
            rs ++= plsP.map(_.toString)
            iii += 1
          }

          writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs.mkString("\t")))
          i += 1
//          }

        }
      }
    }
  }
  def receive = {
//    case wrt:writerName => {
//      writer = Some(system.actorSelection("/user/"+wrt.name))
//    }
    case gList:geneList => {
      val glist = gList.glist.slice(0, 4)
      val file = new java.io.File(gPms.tp + glist(3) + ".gen")
      val vgs:ActorSelection = system.actorSelection("/user/"+glist(3))

      if (!file.exists() || file.length() == 0) vegas2.simuFgene(glist)
      //val rl = scala.io.Source.fromFile(gPms.tp+glist(3)+"_rsid.txt").getLines.toArray.length
      val X = vegas2.vegasX(glist)
      for (h <- H) {
        var j = 0
        while (j < times) {
          val Y = vegas2.setPhenoT(h, 0, 0.5f)(X)

          val Ys = calculation.standardization(Y)
          val sr = glist(3) + "_" + j + "_" + h
          //vgs ! vegas2Actor.inp(sr, Y)
          val future2: Future[(String, Array[Float])] = ask(vgs, vegas2Actor.inp(sr, Y)).mapTo[(String, Array[Float])]

          val (yy, yh, pdofl, gdof) = plsCalc.dofPvalF(X, Ys, k, 1000, true)
          val permp = Array(1 to k: _*).map(i => plsCalc.plsPerm(X, Ys, i, 10000))
          //rsy += (sr-> (glist, yy, yh, pdofl, gdof, permp))
          future2 onComplete {
            case Success(f) => {
              simuwriter.foreach(_ ! simucalculateActor.permp(f._1, f._2))
//              val rs = (glists ++ f._2.map(_.toString) ++ rsm(f._1) :+ f._1.split("_").apply(1)).mkString("\t")
//              writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
//              rsm -= f._1
//              j += 1
            }
            case Failure(t) => println("An error has occured: " + t.getMessage)
          }

          simuwriter.foreach(_ ! rsy(sr,glist, yy, yh, pdofl, gdof, permp))

          ractor ! Ractor.inp(sr, (X, Y, k))
          j += 1
        }
      }
      sender ! done(1)
    }
    case gList:geneLists =>{
      simugenNo2(gList.glist)
      sender ! done(0)
    }
    case gList:gList =>{
      this.calculiting = gList.func
      simugenNo1(gList.glist,gList.n)
      sender ! done(0)
    }

//    case df:dof => {
//      println("")
//      println("f0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")
//
//      val (glt, yy, yh, pdofl, gdof, permp) = rsy(df.idx)
//      rsy -= df.idx
//      val ppval = plsCalc.dofPval(yy, yh, pdofl)
//      val gpval = plsCalc.dofPval(yy, yh, gdof)
//      val kpval = plsCalc.dofPval(yy, yh, df.dt)
//      val rss = glt ++(permp ++ df.dt ++ kpval ++ pdofl ++ ppval ++ gdof ++ gpval).map(_.toString)
//      if (rsm.contains(df.idx)) {
//        val vgsr = rsm(df.idx)
//        val rs = (rss ++ vgsr).mkString("\t") +"\t"+ df.idx.split("_").apply(2)
//        rsm -= df.idx
//        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
//      } else {
//        rsm += (df.idx -> rss)
//      }
//    }
//    case pp:permp =>{
//      println("")
//      println("f111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111")
//
//      if (rsm.contains(pp.idx)) {
//        val ppr = rsm(pp.idx)
//        val rs = (ppr ++ pp.dt.map(_.toString)).mkString("\t")+"\t" + pp.idx.split("_").apply(2)
//        writer.foreach(_ ! myParallel.paraWriterActor.WriteStr(rs))
//      } else {
//        rsm += (pp.idx -> pp.dt.map(_.toString))
//      }
//    }
  case func:calfunc => {
      this.calculiting = func.func
    }
    case don:done => sender ! done(0)
  }
}

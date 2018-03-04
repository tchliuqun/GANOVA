package PLS

import java.io._

import akka.actor.{ActorSelection, PoisonPill}

import scala.util.{Failure, Success}
import scala.collection.parallel._
//import collection.parallel.ForkJoinTasks.defaultForkJoinPool._
//import PLS.gPms.rp
import PLS.snpCalcActor.calcPm
import myParallel.actorMessage._
import PLS.snpCalcOrderActor._
import akka.pattern.ask
import akka.util.Timeout
import scala.concurrent._//{Await, ExecutionContext, Future}
import scala.concurrent.duration._


//import scala.concurrent.Future
import scala.util.Try

object run extends App {
  val currentTime = java.time.LocalDateTime.now().toString.split("\\.").apply(0)
  val logfile = gPms.rp + "liuTest.jfr"
  val newlog = logfile.replaceAll(".jfr", currentTime + ".jfr")

  def mv(oldName: String, newName: String) =
    Try(new File(oldName).renameTo(new File(newName))).getOrElse(false)

  if (new java.io.File(logfile).exists) {
    mv(logfile, newlog)
  }
  if (false) {
    val ch = Array(1 to 22: _*).map(_.toString)
    val orderpms = orderPms()
    val srt = system.actorOf(snpCalcOrderActor.props(orderpms), "srt")

    srt ! snpCalcOrderActor.chrs(ch)
    srt ! snpCalcActor.func(calculation.runeig)
    srt ! snpCalcActor.calcPm(5)
    srt ! action
  }


//  def withParallelism[A](n : Int)(block : => A) : A = {

//    val defaultParLevel = getParallelism
//    collection.parallel.ForkJoinTasks.defaultForkJoinPool.setParallelism(n)
//    val ret = block
//    setParallelism(defaultParLevel)
//    ret
//  }
//
//  withParallelism(2) {
//    (1 to 100).par.map(_ * 2)
//  }

  if (false) {
    val mb = 1024 * 1024
    val runtime = Runtime.getRuntime
    val cores = runtime.availableProcessors()
    runtime.freeMemory / mb
    val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs0 = svd.filter(_ (4).toInt > 50).sortBy(i => i(6).toFloat / i(5).toFloat)
    val rs1 = svd.sortBy(i => i(7).toFloat / i(5).toFloat).reverse
    val rs2 = svd.filter(i => i(0) == "15" & i(4).toInt > 10).sortBy(_ (4).toInt)
    val glistsInx = Range(0, 60, 3).toArray ++ Array(60 until rs2.size: _*)
    val glists = glistsInx.map(rs2(_))
    //val writer = new PrintWriter(new FileWriter("rs" + fileOper.timeForFile + ".txt"))
    val H = Array(0.01f, 0.015f, 0.02f)
    val writer = new PrintWriter(new FileWriter("rs" + fileOper.timeForFile + ".txt"))

    var srt = cores
    var proc = 0
    var end = 0

    def simugenNo(g: Int) = {
      val glist = glists(g).slice(0, 4)
      println("processing No." + g)
      vegas2.simuFgene(glist)
      for (h <- H) {
        var i = 0
        while (i < 100) {
          writer.println((glists(g) ++ vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)) :+ h).mkString("\t"))
          i += 1
        }

      }
    }

//    Array(0 until cores: _*).foreach(i => {
//      val actr = system.actorOf(snpCalcActor.props(calcPm), "calc" + i)
//      actr ! writer
//      actr ! glists(i)
//    })
    while (proc < cores) {
      //val gen = gens(cnt)
      //XX = getX(gen)
      //if (XX.length>0) {
      //  val na = sendCont % nActor
      val calcular = system.actorSelection("/user/calc" + proc)
      //calcular ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
      proc += 1
    }
    //cnt += 1
    //  val f= Future{
    //
    //  }
    val pc = mutable.ParArray(glists.indices: _*)
    pc.tasksupport = new ForkJoinTaskSupport(new java.util.concurrent.ForkJoinPool(cores - 1))
    pc.foreach(simugenNo)
    writer.close()
  }
  //if (false) {
    val orderpms = simumasterActor.Pms(gPms.rp + "simuRs.txt", 100, Array(0.01f, 0.03f, 0.05f))
    val srt = system.actorOf(simumasterActor.props(orderpms), "srt")
    println("start")
    //implicit val timeout = Timeout(999 hours)
    srt ! SnpProcessActor.chr(Array("15"))
  //}
if (false) {
  val filn = myParallel.paraWriterActor.fileName(gPms.rp +"tests.txt")
  val testactor = system.actorOf(myParallel.paraWriterActor.props(filn),"testa")
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+"testa"))
  testactor ! myParallel.paraWriterActor.WriteStr("test1")
  writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("test2"))
  testactor ! done

}
  //  val future:Future[String] = ask(srt, SnpProcessActor.chr(Array("15"))).mapTo[String]
//  //val result: String = future.get()
//  future.onComplete {
//    case Success(value) => {
//      if (value == "done"){
//        srt ! PoisonPill
//      println("job done")
//    }else{
//        println("Something wrong")
//      }
//    }
//    case Failure(e) => e.printStackTrace
//  }
  //val result = Await.result(future, timeout.duration).asInstanceOf[String]
}

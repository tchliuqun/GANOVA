package PLS

import java.io._

import akka.actor.{ActorSelection, PoisonPill}

import scala.util.{Failure, Random, Success}
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
  // pakt and gene expression as phenotype data
//  if (false) {
    val xx = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).take(1).toArray.flatten
    val yy = scala.io.Source.fromFile(gPms.op+gPms.pf).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val ee = scala.io.Source.fromFile(gPms.op+gPms.ef).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val mcol =fileOper.intersectCols(xx,yy,ee)
    val pakt = fileOper.toArrays(gPms.op+gPms.pf).filter(_ (0).contains("AKT")).toArray

    val ch = Array(1 to 22: _*).map(_.toString)
    val orderpms = snpCalcOrderActor.orderPms(k = 3,nactor = 7)//,efile = "")
    val srt = system.actorOf(snpCalcOrderActor.props(orderpms), "srt")
    srt ! snpCalcOrderActor.yArray(pakt(2))
    srt ! snpCalcOrderActor.chrs(ch)


    //srt ! snpCalcActor.func(calculation.runeig)
    srt ! snpCalcActor.calcPm(3)
    srt ! action
//  }

  // 2018-5-14 MGMT status and/or gene expression as phenotype data
    if (false) {
  val pfl = "gbm_mgmt_stp27.txt"
  //val xx = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).take(1).toArray.flatten
  //val yy = scala.io.Source.fromFile(gPms.op+pfl).getLines.map(_.split("\t")).take(1).toArray.flatten//.map(_.slice(0,15))
  //val ee = scala.io.Source.fromFile(gPms.op+gPms.ef).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
  //val mcol =fileOper.intersectCols(xx,yy,ee)
  val mgmt = fileOper.toArrays(gPms.op+pfl).filter(_ (0).contains("mgmt")).toArray.flatten

  val ch = Array(1 to 22: _*).map(_.toString)
  val orderpms = snpCalcOrderActor.orderPms(k = 3,nactor = 7,pfile = gPms.op+pfl)//,efile = "")
  val srt = system.actorOf(snpCalcOrderActor.props(orderpms), "srt")
  srt ! snpCalcOrderActor.yArray(mgmt)
  srt ! snpCalcOrderActor.chrs(ch)
  srt ! snpCalcActor.calcPm(3)
  srt ! action
    }

  // comparing VEGAS2  and new GANOVA using simulating SNP within genes in chr15
  // results are before "simuRs_2018_04_02_07_01_24_582.txt"
  if (false) {
    val orderpms = simumasterActor.Pms(gPms.rp + "simuRs.txt", 100, Array( 0.03f, 0.05f))
    val srt = system.actorOf(simumasterActor.props(orderpms), "srt")
    println("start")
    //implicit val timeout = Timeout(999 hours)
    val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs2 = svd.filter(i => i(0) == "15" & i(4).toInt > 10).sortBy(_ (4).toInt)
    val glist = rs2.filter(i => i(4).toInt > 100 & i(5).toDouble > 0.80).flatten
    srt ! simucalculateActor.gList(glist, 20)

    if (false) {
      srt ! SnpProcessActor.chr(Array("15"))
    }
  }

  // testing parallel with future class
  if(false){
  val f = Future {
    Thread.sleep(Random.nextInt(500))
    if (Random.nextInt(500) > 250) throw new Exception("Yikes!") else 42
    //42
  }
  Thread.sleep(1000)
  println("before onComplete")

  f onComplete {
    case Success(result) => println(s"Success: $result")
    case Failure(t) => println(s"Exception: ${t.getMessage}")
  }

  // do the rest of your work
  println("A ..."); Thread.sleep(100)
  println("B ..."); Thread.sleep(100)
  println("C ..."); Thread.sleep(100)
  println("D ..."); Thread.sleep(100)
  println("E ..."); Thread.sleep(100)
  println("F ..."); Thread.sleep(100)
  Thread.sleep(2000)

  }


  // comparing VEGAS2 with old GANOVA using simulating SNP in chr 15
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
    val wrter = new PrintWriter(new FileWriter("rs" + fileOper.timeForFile + ".txt"))

//    var srt = cores
    var proc = 0
    var end = 0
    val g = 0
    val filnn = myParallel.paraWriterActor.fileName(gPms.rp +"trs.txt")
    val testwrtor = system.actorOf(myParallel.paraWriterActor.props(filnn),"testa")
    def simugenNo(g: Int) = {
      val glist = glists(g).slice(0, 4)
      println("processing No." + g)
      vegas2.simuFgene(glist)
      for (h <- H) {
        var i = 0
        while (i < 100) {
          val rs =(glists(g) ++ vegas2.vegas(glist, 3, vegas2.setPheno(h, 0)) :+ h).mkString("\t")
          wrter.println(rs)
          testwrtor ! myParallel.paraWriterActor.WriteStr(rs)
          i += 1
        }

      }

    }
    val pc = mutable.ParArray(glists.indices: _*)
    pc.foreach( i => simugenNo(i))

  wrter.close()
    testwrtor ! done

//    Array(0 until cores: _*).foreach(i => {
//      val actr = system.actorOf(snpCalcActor.props(calcPm), "calc" + i)
//      actr ! wrter
//      actr ! glists(i)
//    })
//    while (proc < cores) {
//      //val gen = gens(cnt)
//      //XX = getX(gen)
//      //if (XX.length>0) {
//      //  val na = sendCont % nActor
//      val calcular = system.actorSelection("/user/calc" + proc)
//      //calcular ! snpCalcActor.Xs(gen, utils.Array2DM(XX, false))
//      proc += 1
//    }
//    //cnt += 1
//    //  val f= Future{
//    //
//    //  }

//    pc.tasksupport = new ForkJoinTaskSupport(new java.util.concurrent.ForkJoinPool(cores - 1))
//    pc.foreach(simugenNo)
//  }
//  if (false) {

  }

  // testin parallel writing
if (false) {
  val filn = myParallel.paraWriterActor.fileName(gPms.rp +"tests.txt")
  val filn1 = myParallel.paraWriterActor.fileName(gPms.rp +"tests1.txt")
  val testactor = system.actorOf(myParallel.paraWriterActor.props(filn),"testa")
  val testactor1 = system.actorOf(myParallel.paraWriterActor.props(filn1),"testa1")
  var writer:Option[ActorSelection] = Some(system.actorSelection("/user/"+"testa"))
  var writer1:Option[ActorSelection] = Some(system.actorSelection("/user/"+"testa1"))
  testactor1 ! myParallel.paraWriterActor.WriteStr("test11")
  testactor ! myParallel.paraWriterActor.WriteStr("test1")
  writer1.foreach(_ ! myParallel.paraWriterActor.WriteStr("test22"))
  writer.foreach(_ ! myParallel.paraWriterActor.WriteStr("test2"))
  testactor ! done


//  val filn1 = myParallel.paraWriterActor.fileName(gPms.rp +"tests1.txt")
//  val testactor1 = system.actorOf(myParallel.paraWriterActor.props(filn),"testa1")
//  var writer1:Option[ActorSelection] = Some(system.actorSelection("/user/"+"testa1"))
//  testactor1 ! myParallel.paraWriterActor.WriteStr("test11")
//  writer1.foreach(_ ! myParallel.paraWriterActor.WriteStr("test22"))
  //testactor1 ! done


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

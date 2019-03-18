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
import breeze.linalg._

//import scala.concurrent.Future
import scala.util.Try

object run extends App {
  val funs = Set("snpPls","gsea","snpGsea")
  val currentTime = java.time.LocalDateTime.now().toString.split("\\.").apply(0)
  val logfile = gPms.rp + "liuTest.jfr"
  val newlog = logfile.replaceAll(".jfr", currentTime + ".jfr")


  def mv(oldName: String, newName: String) =
    Try(new File(oldName).renameTo(new File(newName))).getOrElse(false)

  //def main(args: Array[String]) = {
  if (false) {
    if (!funs.contains(args(0))) {
      println("first args should be one of: " + funs.mkString(" "))
      sys.exit(1)
    } else {
      val fun = args(0)


      //      val outF=args(1)
      //      val iters=args(2).toInt
      //      val thin=args(3).toInt
      ////      val out = genIters(State(0.0,0.0),iters,thin)
      //      val s = new java.io.FileWriter(outF)

      fun match {
        case "gsea" => {

          val rsf = gPms.homerr + args(1)
          val gsf = gPms.homerr + args(2)
          val gaf = gPms.homerr + args(3)
          val ncol = args(4).toInt
          val rsord = args(5).toInt
          val pv = args(6).toBoolean
          val gos = args(2) contains "goa_human.gaf"
          val rsorder = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).map(_ (rsord).toFloat).toArray
          val orderpms = gseaDispatchActor.pathwayPms(rsFile = rsf, ggFile = gsf, gcFile = gaf, go = gos, rsord = rsorder, pval = pv, namcol = ncol)
          //,efile = "")
          val srt = system.actorOf(gseaDispatchActor.props(orderpms), "srt")
          srt ! action
        }
        case "snpGsea" => {
          val rsf = gPms.homerr + args(1)
          //val gsf = gPms.homerr + args(2)
          val ncol = args(2).toInt
          val rsfirst = args(3).toBoolean
          val pv = args(4).toBoolean
          val rs = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).toArray

          val rsorder = if (rsfirst) {
            rs.map(i => if (i.length == 15) i(9) else if (i.length == 12) i(8) else if (i.length == 11) i(8) else if (i.length == 9) i(7) else i(6)).map(_.toFloat)
          } else {
            rs.map(i => if (i.length == 15) i.slice(9, 12).min else if (i.length == 12) i.slice(8, 10).min else if (i.length == 11) i.slice(8, 11).min else if (i.length == 9 & i(5).toDouble < 1d) i.slice(7, 9).min else if (i.length == 9 & i(5).toDouble < 1d) i(7) else i(6)).map(_.toFloat)
          }
          val orderpms = gseaDispatchActor.pathwayPms(rsFile = rsf, rsord = rsorder, pval = pv, namcol = ncol)
          //,efile = "")
          val srt = system.actorOf(gseaDispatchActor.props(orderpms), "srt")
          srt ! action
        }
      }
      //      s.write("x , y\n")
      //      out map { it => s.write(it.toString) }
      //      s.close
    }

  }
  if (new java.io.File(logfile).exists) {
    mv(logfile, newlog)
  }
  // pakt and gene expression as phenotype data

  if (false) {
    val xx = scala.io.Source.fromFile(gPms.op+gPms.df).getLines.map(_.split("\t")).take(1).toArray.flatten
    val yy = scala.io.Source.fromFile(gPms.op+gPms.pf).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val ee = scala.io.Source.fromFile(gPms.op+gPms.ef).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
    val mcol =fileOper.intersectCols(xx,yy,ee)
    val pakt = fileOper.toArrays(gPms.op+gPms.pf).filter(i => i(0).contains( "Akt_p")).toArray//| i(0).contains("EGFR_pY1068")).toArray
//val pakt = fileOper.toArrays(gPms.op+gPms.pf).filter(i => i(0).contains("Akt")).toArray.slice(1,3)
    val ch = Array(1 to 22: _*).map(_.toString)
    val orderpms = snpCalcOrderActor.orderPms(k = 3,efile = "")
    val srt = system.actorOf(snpCalcOrderActor.props(orderpms), "srt")
    srt ! snpCalcOrderActor.yArray(pakt)
    srt ! snpCalcOrderActor.chrs(ch)
    //srt ! snpCalcActor.func(calculation.runeig)
    srt ! snpCalcActor.calcPm(3)
    srt ! action
}
//2018-8-26 gbm and lgg: snp ,exp and pakt
  // first gbmlgg exp and pakt
  // second gbmlgg and pakt
  // third gbmlgg and pakt adjusted by exp
  // fourth gbm and pakt
  // fifth gbm exp and pakt
  // 2018-9-17pakt ,three dof
  if (false) {
    //val out = new PrintWriter(new FileWriter(gPms.op+"tcga_gbmlgg_rnaseq.txt"))
    // val expp = scala.io.Source.fromFile(gPms.op+"GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt").getLines.map(_.split("\t"))
    //out.println(expp.next.mkString("\t"))
    //val ds = expp.next
    //for (lin <- expp) out.println((lin(0).split("\\|").apply(1) +: lin.drop(1).map(i =>(log10(i.toDouble+1.0)/log10(2d)).toString)).mkString("\t"))
    // out.close()
    val dff = gPms.op+"tcga_gbmlgg_snp.txt"
//  val dff = gPms.op + gPms.df
//  val pff = gPms.op+ "tcga_gbmlgg_rppa_rnaseq.txt"
  val pff = gPms.op+ "LGG.rppa.txt"
    val eff = gPms.op+gPms.ef
//     val eff = ""
//    val eff = gPms.op+"tcga_gbmlgg_rnaseq.txt"


//    val xx = scala.io.Source.fromFile(dff).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
//    val yy = scala.io.Source.fromFile(pff).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
//    val ee = scala.io.Source.fromFile(eff).getLines.map(_.split("\t")).take(1).toArray.flatten.map(_.slice(0,15))
//    val mcol =fileOper.intersectCols(xx,yy,ee)
    val pakt = fileOper.toArrays(pff).filter(_ (0).contains("AKT")).toArray

    val ch = Array(1 to 22: _*).map(_.toString)
    val orderpms = snpCalcOrderActor.orderPms(k = 3,nactor = 7,dfile = dff,pfile = pff,efile = "")
    val srt = system.actorOf(snpCalcOrderActor.props(orderpms), "srt")
    srt ! snpCalcOrderActor.yArray(pakt.slice(1,3))
    srt ! snpCalcOrderActor.chrs(ch)


    //srt ! snpCalcActor.func(calculation.runeig)
    srt ! snpCalcActor.calcPm(3)
    srt ! action
  }
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
  srt ! snpCalcOrderActor.yArray(Array(mgmt))
  srt ! snpCalcOrderActor.chrs(ch)
  srt ! snpCalcActor.calcPm(3)
  srt ! action
    }
  // 2018-5-20 gsea analysis for the results
    if (false) {
  //2018-8-1
  val rsf =  gPms.rp + "GBMsnp6Rs_2018_05_14_15_14_49_569.txt"
  //2018-8-2
  //  val rsf = gPms.rp + "GBMsnp6Rs_2018_05_14_16_38_13_481.txt"
  //
    var rs = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).toArray
  //gseaResult_2018_08_02_11_15_52_376.txt, gseaResult_2018_08_03_00_25_16_486.txt,gseaResult_2018_08_03_11_04_17_209.txt,gseaResult_2018_08_04_20_45_11_609.txt

  //var rsorder  = rs.map(i => if(i.length == 11 )i(8) else if(i.length == 9 ) i(7) else i(6)).map(_.toFloat)
  // GBMsnp6Rs_2018_05_14_15_14_49_569.txt: gseaResult_2018_08_02_14_28_41_340.txt,gseaResult_2018_08_02_23_26_53_653.txt,gseaResult_2018_08_03_13_27_40_727.txt,gseaResult_2018_08_04_16_36_33_145.txt,gseaResult_2018_08_04_19_12_25_443.txt
  //var rsorder = rs.map(i => if(i.length == 15 )i(9) else if(i.length == 12 ) i(8) else i(7)).map(_.toFloat)
  //gseaResult_2018_08_02_12_14_30_216.txt,gseaResult_2018_08_03_09_03_28_657.txt,gseaResult_2018_08_04_22_56_13_394.txt
  //var rsorder = rs.map(i => if(i.length == 11 )i.slice(8,11).min else if(i.length == 9 ) i.slice(7,9).min else i(6)).map(_.toFloat)
  // GBMsnp6Rs_2018_05_14_15_14_49_569.txt:gseaResult_2018_08_02_15_26_26_047.txt,gseaResult_2018_08_02_16_29_50_746.txt,gseaResult_2018_08_02_21_41_00_747.txt,gseaResult_2018_08_03_14_42_00_044.txt,gseaResult_2018_08_04_17_21_00_824.txt
  var rsorder = rs.map(i => if(i.length == 15 )i.slice(9,12).min else if(i.length == 12 ) i.slice(8,10).min else i(7)).map(_.toFloat)
  val orderpms = gseaDispatchActor.pathwayPms(rsFile = rsf,rsord = rsorder)//,efile = "")
    val srt = system.actorOf(gseaDispatchActor.props(orderpms), "srt")
    srt ! action
  }
  // 2018-8-9 gbm expression Vs pakt go set base pls
    if (false) {
    var sfile = gPms.op +"GOGeneList.txt"
    var gfile = gPms.op +"gbmExp.txt"
    var pfile = gPms.op +"GBMLGG.rppa.txt"
    val pname = Array("Akt_pS473","Akt_pT308")

    val orderpms = setDispatchActor.dispatcherPms(gfile = gfile, sfile =sfile,pfile = pfile,pname = pname)
    //,efile = "")
    val srt = system.actorOf(setDispatchActor.props(orderpms), "srt")
    srt ! action
  }
  // 2018-8-13 gbmlgg expression Vs pakt go set base pls
    if (false) {
  var sfile = gPms.op +"GOGeneList.txt"
  var gfile = gPms.op +"GBMLGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
  var pfile = gPms.op +"GBMLGG.rppa.txt"
  val pname = Array("Akt_pS473","Akt_pT308")

  val orderpms = setDispatchActor.dispatcherPms(gfile = gfile, sfile =sfile,pfile = pfile,pname = pname)
  //,efile = "")
  val srt = system.actorOf(setDispatchActor.props(orderpms), "srt")
  srt ! action
    }
  // comparing VEGAS2  and new GANOVA using simulating SNP within genes in chr15
  // results are before "simuRs_2018_04_02_07_01_24_582.txt"
  if (false) {
//    val t1 = System.nanoTime()
//    //plsCalc.gdofPval2(X1,Y0,2)
//    val rss = calculation.sgsea(gs,rs)
//    val t2 = System.nanoTime()
//    val lapse1 = (t2 - t1) / 1e9d
    val orderpms = simumasterActor.Pms(gPms.rp + "simuRs.txt", 100, Array(0.05f, 0.06f, 0.07f))
    val srt = system.actorOf(simumasterActor.props(orderpms), "srt")
    println("start")
    //implicit val timeout = Timeout(999 hours)
    val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    val rs2 = svd.filter(i => i(0) == "15" & i(4).toInt > 10).sortBy(_ (4).toInt)
    val glist = rs2.filter(i => i(4).toInt > 100 & i(5).toDouble > 0.80).flatten
    srt ! simucalculateActor.gList(glist, 10)

    if (false) {
      srt ! simumasterActor.chr(Array("15"))
    }
  }
  // 2018-6-11 15:20 testing different dof calculating method and compare to permutation results of PLS
  if (false) {
    val orderpms = simumasterActor.Pms(gPms.rp + "simuRs.txt", 10, 0.02f.to(0.06f,0.005f).toArray,3,plsCalc.plsAdof _ )
    val srt = system.actorOf(simumasterActor.props(orderpms), "srt")
    println("start")
    //implicit val timeout = Timeout(999 hours)
    //val svd = fileOper.toArrays(gPms.rp + "GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
    //val rs2 = svd.filter(i => i(0) == "15" & i(4).toInt > 10).sortBy(_ (4).toInt)
    //al glist = rs2.filter(i => i(4).toInt > 100 & i(5).toDouble > 0.80).flatten
    srt ! SnpProcessActor.chr(Array("15"))//simucalculateActor.gList(glist, 20)


  }
  //2018-7-11
  // 2018-10-18 simulation with different length of Y and different correlation in Ys
  // 2019-3-18 simulation PLS,PCR,VEGAS2 with 100 genes in chr15
//  if (false) {
    val orderpms = simumasterActor.Pms(gPms.rp + "simuRs.txt", 100, Array(0.01f,0.03f, 0.05f))
    val srt = system.actorOf(simumasterActor.props(orderpms), "srt")
    println("start")
    //implicit val timeout = Timeout(999 hours)

      srt ! simumasterActor.chr(Array("15"))
//  }
//2018-7-15
  if (false) {
  val orderpms = simuSnpActor.Pms(gPms.rp + "simuRs.txt", 500, Array(0.03f),3,500,50)
  val srt = system.actorOf(simuSnpActor.props(orderpms), "srt")
  val nis = (5 to 100 by 3).toArray
  println("start")
  //implicit val timeout = Timeout(999 hours)

  srt ! simuSnpActor.ns(nis)
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

  // testing parallel writing
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

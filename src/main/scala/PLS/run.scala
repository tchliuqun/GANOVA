package PLS

import java.io.{File, FileWriter, PrintWriter}

import PLS.gPms.rp
import myParallel.actorMessage._
import PLS.snpCalcOrderActor._

import scala.util.Try

object run extends App {
  val currentTime = java.time.LocalDateTime.now().toString.split("\\.").apply(0)
  val logfile = gPms.rp+"liuTest.jfr"
  val newlog = logfile.replaceAll(".jfr", currentTime + ".jfr")

  def mv(oldName: String, newName: String) =
    Try(new File(oldName).renameTo(new File(newName))).getOrElse(false)

  if (new java.io.File(logfile).exists) {
    mv(logfile, newlog)
  }
  val ch = Array(1 to 22:_*).map(_.toString)
  val orderpms = orderPms()
  val srt = system.actorOf(snpCalcOrderActor.props(orderpms),"srt")

  srt ! snpCalcOrderActor.chrs(ch)
  srt ! snpCalcActor.func(calculation.runeig)
  srt ! snpCalcActor.calcPm(5)
  srt ! action

  val svd = fileOper.toArrays(gPms.rp+"GBMsnp6Rs_2018-01-01_23.txt").drop(1).toArray
  val rs0 = svd.filter(_ (4).toInt > 50).sortBy(i => i(6).toFloat / i(5).toFloat)
  val rs1 = svd.sortBy(i => i(7).toFloat / i(5).toFloat).reverse
  val glist = rs1(0).slice(0, 4)

  vegas2.simuFgene(glist)
  val H = Array(0.01f,0.015f,0.02f)
  for (h <- H) {
    val writer = new PrintWriter(new FileWriter("rs" + fileOper.timeForFile + ".txt"))
    var i = 0
    while (i < 100) {
      writer.println(vegas2.vegas(glist, 3, vegas2.setPheno2(h, 2)).mkString("\t"))
      i += 1
    }
    writer.close()
  }
}

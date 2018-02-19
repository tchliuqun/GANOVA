package PLS

import java.io.File

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
}

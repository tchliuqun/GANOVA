package PLS


import java.io.{FileWriter, PrintWriter}

import akka.actor.{Actor, ActorRef, ActorSystem, PoisonPill, Props}
import myParallel.paraWriterActor._
import myParallel.actorMessage._
import myParallel.{actorMessage, paraWriterActor}

import scala.util.Try

/**
  * Created by liuqun on 12/11/16.
  */
object SnpProcessActor{
  val name = "SnpProcessActor"
  def props(fil:threeFile) = Props(classOf[SnpProcessActor],fil)
  case class threeFile(sfile: String,gfile:String,dfile:String)
  case class chr(chrname:Array[String])
  case class totalNumber(num:Int)
  def sepSnpDataFile(ddfile:String,ssfile:String)  = {
    val dpre = ddfile.split("\\.")
    val spre = ssfile.split("\\.")
    val dpref = dpre(0)+"Tchr"
    val spref = spre(0) + "Tchr"
    val sf = fileOper.toArrays(ssfile).toArray
    val chr = sf.map(_(1)).toSet.filter(i => Try(i.toInt).isSuccess).toArray.sortBy(_.toInt)
    val dnam = chr.map(i => dpref + i)
    val snam = chr.map(i => spref + i)
    val dactors = dnam.map(i => system.actorOf(paraWriterActor.props(fileName(i+".txt")),i))
    val sactors = snam.map(i => system.actorOf(paraWriterActor.props(fileName(i+".txt")),i))
    val smap = sf.map(i => i(0)->i(1)).toMap
    val lines = fileOper.toArrays(ddfile,",")
    var i = 0
    while (lines.hasNext){
      val line = lines.next
      val sline = sf(i)
      val chrn = smap.get(line(0))
      chrn match {
        case Some(x) => {
          val dwriter = system.actorSelection("/user/"+dpref+x)
          val swriter = system.actorSelection("/user/"+spref+x)
          swriter ! paraWriterActor.WriteStr(sline.mkString("\t"))
          dwriter ! paraWriterActor.WriteStr(line.drop(1).mkString("\t"))
        }
        case _ =>
      }
      i += 1
    }
    dactors.foreach(_ ! actorMessage.finished)
    sactors.foreach(_ ! actorMessage.finished)
  }
  def sortSNP(dataf:String,annof:String): Unit ={
    //val dpref = dataf.split("\\.")
    //val apref = annof.split("\\.")
    val dnam = dataf.replace(".csv","Tchrn.txt")
    val snam = annof.replace(".txt","Tchrn.txt")
    val sf = fileOper.toArrays(annof).toArray
    val chr = sf.map(_(1)).toSet.filter(i => Try(i.toInt).isSuccess).toArray.sortBy(_.toInt).map("chr"+_)
    for( i <- chr){
      val ddfil = fileOper.toArrays(dnam.replace("chrn",i))(2).toArray
      val ssfil = fileOper.toArrays(snam.replace("chrn",i)).toArray
      val xy = ddfil.map(_(0)).toSet.intersect(ssfil.map(_(0)).toSet)
      val ssfilt = ssfil.filter(i =>xy.contains(i(0))).sortBy(i =>i(2).toInt)
      val sortMap = ssfilt.map(_(0)).zipWithIndex.toMap
      val ddfilt = ddfil.filter(i => xy.contains(i(0))).sortBy(i => sortMap(i(0)))
      val dout = new PrintWriter(new FileWriter(dnam.replace("Tchrn",i)))
      val sout = new PrintWriter(new FileWriter(snam.replace("Tchrn",i)))
      ssfilt.foreach(i => sout.println(i.mkString("\t")))
      ddfilt.foreach(i => dout.println(i.mkString("\t")))
      dout.close()
      sout.close()
    }
  }
}
class SnpProcessActor(fil:SnpProcessActor.threeFile) extends Actor {

  var chrName:Array[String] = Array(1 to 22:_*).map(_.toString)
  var sfile:String = fil.sfile// "snp6annoNew.txt"
  var gfile:String = fil.gfile //"geneLoc.txt"
  var dfile:String = fil.dfile//"gbmsnp6.csv"
  val chrdf = "dfilechrn"
  def sepName(s:Array[String],f:String) = {
    s.map(i => f.split("\\.").apply(0)+"chr"+i)
  }
  def setWriter(chr:Array[String] = chrName) = {
    val dnam = chr.map(i => dfile.split("\\.").apply(0)+"chr"+i)
    dnam.map(i => system.actorOf(paraWriterActor.props(fileName(i+".csv")),i))
  }

  def receive = {
    case _ =>

  }


}

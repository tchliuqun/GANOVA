package PLS

import akka.actor._
import breeze.linalg._
import PLS.gseaDispatchActor._
import myParallel.actorMessage._

import scala.concurrent.Future
//import akka.Done

object gseaDispatchActor{
  val name = "pathwayDispatchActor"
  def props(pms:pathwayPms) = Props(classOf[gseaDispatchActor],pms)
  case class yArray(y:Array[String])
  case class pathwayPms(cores:Int = 0,perm:Int = 1000,namcol:Int = 3,pval:Boolean = true,rsFile:String =  gPms.rp + "GBMsnp6Rs_2018_05_14_15_14_49_569.txt",
                        rsord:Array[Float] = Array[Float](),go:Boolean= true,
                        ggFile:String = gPms.op +"goa_human.gaf",gaFile:String = gPms.op + "HUMAN_9606_idmapping_selected.tab",
                        gcFile:String = gPms.op + "goose",outFile:String = gPms.rp+"gseaResult.txt")
}

class gseaDispatchActor(pm:pathwayPms) extends Actor {
  var ncores = if (pm.cores < 1) Runtime.getRuntime.availableProcessors() - 1 else pm.cores

  var cores = if (ncores > 50 ) 50 else ncores
  var total = 0
  var counts = 0
  var donecount = 0
  var yarray:DenseMatrix[Float] = null
  var goSet:DenseMatrix[Float] = null
  var genesetAno:Array[Array[String]] = null

  var rsords = pm.rsord
  var ofile:String = pm.outFile
  var perm = pm.perm
  var pval = pm.pval
  var namcol = pm.namcol
  var gos = pm.go
  var resultFile:String = pm.rsFile
  var goGeneFile:String = pm.ggFile
  var geneAnno:String = pm.gaFile
  var goChildFile = pm.gcFile
  def getGOMatrix(rsf:String = this.resultFile,rsord:Array[Float] = this.rsords,pval:Boolean = this.pval) = {
    var idg = scala.io.Source.fromFile(goGeneFile).getLines.drop(23).map(_.split("\t")).map(i => (i(1), i(4))).toArray.groupBy(_._1).map(i => (i._1, i._2.map(_._2).distinct))

    var idm = scala.io.Source.fromFile(geneAnno).getLines.map(_.split("\t")).map(i => (i.filter(_.contains("ENSG")), i(0))).filter(i => i._1.length > 0 & i._2.length > 0).map(i => (i._1.apply(0), idg.getOrElse(i._2, Array()))).filter(_._2.length > 0).toArray
    idg = null
    var gom = scala.io.Source.fromFile(goChildFile).getLines.map(_.split("\t")).dropWhile(_ (0) == "all").map(i => Array(4, 0).map(i(_))).toArray.groupBy(_ (0)).map(i => (i._1, i._2.map(_ (1))))
    var genego = idm.map(i => (i._1, i._2.flatMap(i => gom.getOrElse(i, Array())).toSet))
    idm = null
    gom = null
    var goset = genego.map(_._2).flatten.toSet.toArray.sorted
    var gomxp = genego.map(i => (i._1, goset.map(ii => if (i._2.contains(ii)) 1f else 0f)))

    genego = null
    //goset = null
    var rs = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).toArray
    var rsset = rs.map(_ (namcol)).toSet.intersect(gomxp.map(_._1).toSet)
    //var rss0 = rs.filter(i => rsset.contains(i(3)))
    var rsorder  = if(rsord.length != 0) rsord else rs.map(i => if(i(5).toInt >2)i(9) else if(i(5).toInt ==2) i(8) else i(7)).map(_.toFloat)
    //val rss = Array( 0 until rss0.length :_*).map( i => rss0(i).slice(0,6) :+ rsorder(i).toString).sortBy(_(6).toFloat)

    var rss = Array( 0 until rs.length :_*).map( i => Array(rs(i).apply(namcol), rsorder(i).toString)).filter(i => rsset.contains(i(0))).filter(!_(1).toFloat.isNaN).sortBy(_(1).toFloat)

    val chisq = if(pval)DenseVector(calculation.reverseChisq(rss.map(_(1).toFloat),1))else DenseVector(rss.map(_(1).toFloat))
    rs = null

    //rss0 = null
    rsorder = null
    var gomx = gomxp.filter(i => rsset.contains(i._1))
    rsset = null
    gomxp = null
    var indm = gomx.map(_._1).zipWithIndex.toMap

    var goinx = rss.map(i => indm(i(0))).map(gomx(_))//.filter(i => sum(i._2) > 15 & sum(i._2)<500)
    //.filter(i => sum(i._2) > 15 & sum(i._2)<500)
    gomx = null
    rss = null

    var gm = new DenseMatrix(goinx(0)._2.length, goinx.length, goinx.flatMap(_._2)).t
    goinx = null
    val goidx = sum(gm(::,*)).t.toArray.zipWithIndex.filter(i => i._1 > 15 & i._1 < 500).map(_._2)
    var goano = scala.io.Source.fromFile(goChildFile).getLines.map(_.split("\t").take(3).mkString("\t")).toArray.distinct.map(_.split("\t")).map( i => (i(0) -> i)).toMap
    val gmm = gm(::,goidx.toIndexedSeq).toDenseMatrix

    gm = null
    val coln = goidx.map(goset(_)).map(goano(_))
    //goano = null

    //val golen = sum(gm(::,*))
    //val res = calculation.gseaP(gm,chisq)

    (chisq,gmm,coln)
  }
  def getGSMatrix(rsf:String = this.resultFile,gsf:String = this.goGeneFile,rsord:Array[Float] = this.rsords,pval:Boolean = this.pval) = {
    var idg = fileOper.fileFlip(gsf)// scala.io.Source.fromFile(gsf).getLines.map(_.split("\t")).map(i => (i(0), i.drop(1)).toArray
    var goset = idg.map(_._2).flatten.toSet.toArray.sorted


    var rs = scala.io.Source.fromFile(rsf).getLines.drop(1).map(_.split("\t")).toArray
    var rsset = rs.map(_ (namcol)).toSet//filter(genego.contains.intersect(idg.map(_._1).toSet)
    var rsst = rsset.toArray
    var gomxp = idg.map(i => (i._1, goset.map(ii => if (i._2.contains(ii)) 1f else 0f)))

    var rsorder  = if(rsord.length != 0) rsord else rs.map(i => if(i(5).toInt >2)i(9) else if(i(5).toInt ==2) i(8) else i(7)).map(_.toFloat)

    var rss = Array( 0 until rs.length :_*).map( i => Array(rs(i).apply(namcol), rsorder(i).toString)).filter(!_(1).toFloat.isNaN).sortBy(_(1).toFloat)
    val chisq = if(pval)DenseVector(calculation.reverseChisq(rss.map(_(1).toFloat),1))else DenseVector(rss.map(_(1).toFloat))

    var gomx = gomxp.filter(i => rsset.contains(i._1))
    var nogen = rsset.diff(gomx.map(_._1).toSet).toArray
    val ncol = gomx(0)._2.length
    nogen.foreach(i => gomx +:= (i,Array.fill(ncol)(0f)))

    var indm = gomx.map(_._1).zipWithIndex.toMap
    var goinx = rss.map(i => indm(i(0))).map(gomx(_))//.filter(i => sum(i._2) > 15 & sum(i._2)<500)
    var gm = new DenseMatrix(goinx(0)._2.length, goinx.length, goinx.flatMap(_._2)).t
    val goidx = sum(gm(::,*)).t.toArray.zipWithIndex.filter(i => i._1 > 15 & i._1 < 500).map(_._2)
    var goano = scala.io.Source.fromFile(goChildFile).getLines.map(_.split("\t")).map( i => (i(0) -> i)).toMap
    val gmm = gm(::,goidx.toIndexedSeq).toDenseMatrix

    val coln = goidx.map(goset(_)).map(goano(_))


    (chisq,gmm,coln)
  }
  def makeActors() = {
    val calcPm = gseaCalcActor.Pms(this.ofile.replaceFirst("/",""))
    Array(0 until cores:_*).map(i => system.actorOf(gseaCalcActor.props(calcPm),"calc"+i))
    system.actorOf(myParallel.paraWriterActor.props(myParallel.paraWriterActor.fileName(this.ofile)),"gseawriter")
  }

  def receive = {
    case myParallel.actorMessage.action => {
      val goset = if(this.gos)getGOMatrix(resultFile) else getGSMatrix(resultFile)
      this.yarray = calculation.permY(goset._1, perm)
      this.goSet  = goset._2
      this.genesetAno = goset._3
      this.total = genesetAno.length
      makeActors()
      while (counts < cores *2 & counts < total){

        val calculator = system.actorSelection("/user/calc" +(counts % cores))
        val wrt = system.actorSelection("/user/gseawriter")
        if (counts == 0) wrt ! myParallel.paraWriterActor.totalNumber(total)
        if (counts < cores) calculator ! gseaCalcActor.goresults(yarray)
        calculator ! gseaCalcActor.geneset(genesetAno(counts),goSet(::,counts).toDenseVector)
        counts += 1

      }

      //updateY(y.y)
    }

    case dn:myParallel.actorMessage.done => {
      donecount += 1

      if (counts < total) {
        if (counts % 100 == 0) println("calculating No."+counts)

        sender ! gseaCalcActor.geneset(genesetAno(counts), goSet(::, counts).toDenseVector)
        counts += 1

      }else {
        //println("gsea job done")
        sender ! PoisonPill

      }
      if (donecount >= total ){
        println("gsea job done")
        val don : Future[akka.Done] = CoordinatedShutdown(system).run()
        self ! PoisonPill
      }
    }
  }
}
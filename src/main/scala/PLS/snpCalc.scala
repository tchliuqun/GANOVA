package PLS

  import java.io.{FileWriter, PrintWriter}

  import breeze.numerics.sqrt

  //import SNP.fileOper.FindIndexinRange
  import breeze.linalg._

  import scala.collection.mutable.ListBuffer
  import scala.io.Source
  import scala.collection.mutable.ArrayBuffer

  /**
    * Created by liuqun on 11/20/16.
    */
  class snpCalc{
    var gfile = gPms.op+gPms.glf
    var sfile = gPms.op+gPms.af
    var dfile = gPms.op+gPms.df
    var rfile = gPms.op+gPms.rf
    var pfile = gPms.op+gPms.pf
    val ofile = gPms.rp+gPms.rsf
    var mcol:(Array[Int], Array[Int]) = (Array(0),Array(0))
    var amplength = 50000
    var k = 5
    var n = 236
    var perm = 1000
    var Y = DenseMatrix.zeros[Float](n,1)
    var permY = DenseMatrix.zeros[Float](n,perm)
    var looInx = Array(0 until n :_*).map(Seq(_))
    var tenFold = plsCalc.kfoldInx(n,10,true)
    case class geneSNPResult(Ypred:DenseMatrix[Float],pdofLoo:Array[Float],gdofLoo:Array[Float],pdofTenf:Array[Float],gdofTenf:Array[Float])
    def updateFile(chrn:String)= {
      this.sfile = sfile.replace(".txt","chr"+chrn+".txt")
      this.rfile = rfile.replace(".txt","chr"+chrn+".txt")
      this.dfile = dfile.replace(".csv","chr"+chrn+".txt")
      this.gfile = gfile.replace(".txt","chr"+chrn+".txt")
    }
    def update(y:Array[Float],chrn:String) = {
      this.updateY(y)
      this.updateFile(chrn)
    }
    def updateY(y:Array[Float]) = {
      val pn = fileOper.toArrays(pfile).next.map(_.slice(0,15))
      val dn = fileOper.toArrays(dfile,",").next.drop(2)
      val pakt = y
      this.mcol = fileOper.intersectCols(dn,pn)
      this.n = mcol._1.length
      this.looInx = Array(0 until n :_*).map(Seq(_))
      this.tenFold = plsCalc.kfoldInx(n,10,true)
      this.Y = new DenseMatrix(n,1,mcol._2.map(y))
      this.permY = plsCalc.permY(Y,perm)
    }
    def getN = {
      mcol._1.length
    }
    def getY = {
      val pakt = fileOper.toArrays(pfile).filter(_(0).contains("AKT")).toArray
      new DenseMatrix(n,1,mcol._2.map(pakt(1)).map(_.toFloat))
    }
    def setY(Ys:DenseVector[Float])= {
      this.Y = Ys.toDenseMatrix.t
    }
    def setY(Ys:Array[Float]) = {
      this.Y = new DenseMatrix(1,Ys.length,Ys)
    }
    def getYperm = {
      plsCalc.permY(Y,perm)
    }
    def getXs = {
      val snpd = fileOper.toArrays(dfile).toArray
      val dm1 = snpd.map(_.drop(1).map(_.toFloat)).map(mcol._1.map(_))
      utils.Array2DM(dm1,false)
    }
    def getX(start:Int,end:Int,dm:DenseMatrix[Float])= {
      dm(::,start until end)
    }

    def getGeneSnp(X:DenseMatrix[Float]):geneSNPResult= {
      val YpredTenFold = plsCalc.plsCV(X,Y,k,tenFold)
      val YpredLoo = plsCalc.plsCV(X,Y,k,looInx)
      val Yhat= plsCalc.predict(X,plsCalc.plsTrain(X,Y,k))
      val pdofloo = plsCalc.pdof(Y,Yhat,YpredLoo)
      val pdoffold = plsCalc.pdof(Y,Yhat,YpredTenFold)
      val rssloo = 0.6f * sqrt(plsCalc.rss(Y,YpredLoo))
      val gdofloo = Array(0 until k:_*).map(i => plsCalc.gdof(X,Y,rssloo(i),i+1,1000))
      val rssfold = 0.6f * sqrt(plsCalc.rss(Y,YpredTenFold))
      val gdoffold = Array(0 until k:_*).map(i => plsCalc.gdof(X,Y,rssfold(i),i+1,1000))
      //.map(i => new DenseMatrix[Float](1,1,Array(i)))
      geneSNPResult(Yhat,pdofloo.data,gdofloo,pdoffold.data,gdoffold)
    }

    def genePval(grs:geneSNPResult) = {
      val pLooPval = plsCalc.dofPval(Y,grs.Ypred,grs.pdofLoo)
      val pTenPval = plsCalc.dofPval(Y,grs.Ypred,grs.pdofTenf)
      val gLooPval = plsCalc.dofPval(Y,grs.Ypred,grs.gdofLoo)
      val gTenPval = plsCalc.dofPval(Y,grs.Ypred,grs.gdofTenf)
      Array((gLooPval,grs.gdofLoo),(gTenPval,grs.gdofTenf),(pLooPval,grs.pdofLoo),(pTenPval,grs.pdofTenf))
    }

    def chrfile(str:String) = str.replaceAll("(.txt|.csv)","chrn.txt")
    //def
    def validRows(sfile:String,dfile:String) = {
      val ann = fileOper.toArrays(sfile).map(_(0)).toArray
      val snp = fileOper.toArrays(dfile,",")(2).map(_(0)).toArray
      ann.deep == snp.deep
    }
    //val x = fileOper.FindIndexinRange(gfile,sfile,ofile,0)

    //val paktFile = scala.io.Source.fromFile("./resources/aktrppawithSnp.txt").getLines().drop(1).map(_.split("\t").drop(1).map(_.toFloat))
    //val pakt473org = utils.Array2DM(Array(paktFile.next)).t
    //val pakt308 = utils.Array2DM(Array(paktFile.next)).t


   // def getresult(Y:DenseMatrix[Float] = pakt473org,adj:Boolean = false,times:Int = 0,parallel:Int = 2) = {
   //   println(new java.util.Date + "\t" + "getting SNP in the range")
   //   val out = fileOper.FindIndexinRange(ifile,amp= amplength).map(_.split("\t")).groupBy(_(1))//.par
   //   //    val tes = out.map { i => i.slice(i.length - 2, i.length).map(_.toInt) }
   //   //    val res = tes.map(i => i(1) - i(0)).groupBy(identity).mapValues(_.size).toArray.sortBy(_._1).takeWhile(_._1 < 9).map(_._2).sum
   //   //  res.takeWhile(_._1 <9).map(_._2)
   //   val Yts = if (times > 0) plsCalc.permY(Y,times) else Y
   //   //    val output2 = new ListBuffer[(String, DenseVector[Float], DenseVector[Float])]
   //   //    val output1 = new ListBuffer[String]//, DenseVector[Float])]
   //   //    val output3 = new ListBuffer[(String, Float,Float)]
   //   val result = ArrayBuffer[String]()
   //   val chrCol = 1
   //   val nameCol = 0
   //   val startCol = out("X").apply(0).length - 2
   //   val endCol =out("X").apply(0).length - 1
   //   val para = parallel -1
   //   def getCvPval(snpInx:Array[Array[Float]]):Array[String] => String = { gene =>
   //     val chr = gene(chrCol)
   //     val start = gene(startCol).toInt
   //     val end = gene(endCol).toInt
   //     var outp = new String
   //     if (end - start > ncomp) {
   //       val X = calculation.filterDuplicateColumn(utils.Array2DM(snpInx.slice(start, end), false))
   //       if (X.cols > ncomp) {
   //         if (times > 0) {
   //           val res = plsCalc.cvforPvalue(X, Yts, 10, ncomp)
   //           val pvalue = calculation.permPvalue(res)
   //           outp = gene(nameCol)+ "\t" +(end-start) + "\t"+res(0)+"\t"+ pvalue
   //           //              fileOut.println(outp) // +"\t"+outp._3.toArray.mkString("\t") )
   //           //            output1 += outp
   //         } else {
   //           if (adj) {
   //             val res = plsMain.KfoldCV(X, Yts, 10, ncomp)
   //             outp = gene(nameCol) + "\t" +(end-start) + "\t"+ res._1.toArray.mkString("\t") + "\t" + res._2.toArray.mkString("\t")
   //             //                fileOut.println(outp)
   //             //              output1 += outp

   //           } else {
   //             val res = plsMain.KfoldCVnoAdj(X, Yts, 10, ncomp)
   //             outp = gene(nameCol) + "\t" +(end-start) + "\t"+ res.toArray.mkString("\t")
   //             //                fileOut.println(outp) // +"\t"+outp._3.toArray.mkString("\t") )
   //             //              output1 += outp
   //           }
   //         }
   //       }else{

   //         if (times > 0) outp = gene(nameCol)+ "\t" +2+"\t"+ 2
   //         //          output1 += outp
   //       }
   //     }else{
   //       if (times > 0) outp = gene(nameCol)+ "\t" +2+"\t"+ 2
   //       //        output1 += outp

   //     }
   //     //        val reslt = if (times > 0) output3.toList else {
   //     //          if (adj) output2.toList else output1.toList
   //     //        }
   //     //      output1.toList.toArray
   //     outp
   //   }

   //   println(new java.util.Date + "\t" + "Starting calculation")
   //   for  ((k,v) <- out) {
   //     println(new java.util.Date + "\t" + "Starting calculation chromosome " + k)
   //     val snpInx = Source.fromFile("./resources/Tcga_lggSNPCHR" + k + "_out.txt").getLines().map(_.split("\t").drop(4).map(_.toFloat)).toArray
   //     //      var i = 0
   //     val inp = v.toArray
   //     val len = inp.length
   //     //      val parFunc = new parFunctions
   //     parFunction.threshold = if (para < 1 ) len+1 else len/para

   //     val reslt = parFunction.mapPar(inp)(getCvPval(snpInx))// if (parallel) v.toArray.par.map(getCvPval(_,snpInx)).toArray else v.toArray.map(getCvPval(_,snpInx))
   //     result ++= reslt
   //   }

   //   val fileOut = new PrintWriter(new FileWriter("./resources/SNP_result"+utils.getTimeForFile+"_"+amplength+".txt"))
   //   result.foreach(fileOut.println(_)) // +"\t"+outp._3.toArray.mkString("\t") )
   //   fileOut.close()
   //   result.toArray
   // }
   // def getResultbySeq(Y:DenseMatrix[Float] = pakt473org,adj:Boolean = false,times:Int = 0)  = {
   //   println(new java.util.Date + "\t" + "getting SNP in the range")
   //   val out = fileOper.FindIndexinRange(ifile, amp = amplength).map(_.split("\t")).toArray
   //   val Yts = if (times > 0) plsCalc.permY(Y,times) else Y
   //   val output1 = new ListBuffer[String]//, DenseVector[Float])]
   //   var chrRef = out(0)(1)
   //   var snpInx = Source.fromFile("./resources/Tcga_lggSNPCHR" + chrRef + "_out.txt").getLines().map(_.split("\t").drop(4).map(_.toFloat)).toArray
   //   var i = 0
   //   val chrCol = 1
   //   val nameCol = 0
   //   val startCol = out(0).length - 2
   //   val endCol = out(0).length - 1
   //   val fileOut = new PrintWriter(new FileWriter("./resources/SNP_result"+utils.getTimeForFile+"_"+amplength+".txt"))
   //   for (gene <- out) {
   //     val chr = gene(chrCol)
   //     val start = gene(startCol).toInt
   //     val end = gene(endCol).toInt
   //     var outp = new String
   //     if (chr != chrRef) {
   //       chrRef = chr
   //       snpInx = Source.fromFile("./resources/Tcga_lggSNPCHR" + chrRef + "_out.txt").getLines().map(_.split("\t").drop(4).map(_.toFloat)).toArray
   //     }
   //     if (end - start > ncomp) {
   //       val X = calculation.filterDuplicateColumn(utils.Array2DM(snpInx.slice(start, end),false))
   //       if (X.cols > ncomp) {
   //         if (times > 0){
   //           val res = plsCalc.cvforPvalue(X, Yts, 10, ncomp)
   //           val pvalue = calculation.permPvalue(res)
   //           outp = gene(nameCol)+ "\t" +res(0)+"\t"+ pvalue
   //           //            val outp = (gene(nameCol), res(0),pvalue)
   //           fileOut.println(outp) // +"\t"+outp._3.toArray.mkString("\t") )
   //           output1 += outp
   //         }else{
   //           if (adj){
   //             val res = plsCalc.KfoldCV(X, Yts, 10, ncomp)
   //             outp = gene(nameCol) + "\t" + res._1.toArray.mkString("\t") + "\t" + res._2.toArray.mkString("\t")
   //             //            val outp = (gene(nameCol), res._1, res._2)
   //             fileOut.println(outp)
   //             output1 += outp

   //           } else {
   //             val res = plsMain.KfoldCVnoAdj(X, Yts, 10, ncomp)
   //             outp = gene(nameCol) + "\t" + res.toArray.mkString("\t")
   //             //            val outp = (gene(nameCol), res)
   //           }
   //         }
   //       } else{
   //         outp = gene(nameCol)+ "\t" +2+"\t"+ 2
   //       }
   //     } else {
   //       outp = gene(nameCol) + "\t" + 2 + "\t" + 2
   //     }
   //     fileOut.println(outp) // +"\t"+outp._3.toArray.mkString("\t") )
   //     output1 += outp
   //     if (i % 2000 == 0) {
   //       println(new java.util.Date + "\t" + "completed "+100 * i.toFloat/out.length + "%")
   //     }
   //     i += 1
   //   }
   //   fileOut.close()
    // }
  }


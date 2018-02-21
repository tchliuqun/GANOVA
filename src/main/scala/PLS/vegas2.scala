package PLS

import java.io._

import breeze.linalg.{DenseMatrix, DenseVector, princomp}
import breeze.linalg._
import breeze.numerics._
import breeze.stats._
import scala.math.Numeric
import scala.sys.process._
object vegas2 {
  var chr = 2
  val pPath = "/Users/qunliu/VEGAS2offline"
  val dPath = "VEGAS2database/1000GEURO/0kbloc/geneset"
  val aGene = "allgene"
  val sGene = "snpgene"
  val comm = "pathtovegas/vegas2 test_vegas2input.txt -genelist genelist.txt -pop 1000GEURO -subpop EURO -genesize 0kbloc -top 100 –sex BothMnF –max 1000000 –out genebased.V2out"
  val simuSnpP = "hapgen2 -l refpath/CEU.0908/CEU.0908.chr10.legend -m  refpath/genetic_map_chr10_combined_b36.txt -h refpath/CEU.0908.chr10.hap -dl 8050000 1 1 1 -n 500 500 -int 8000000 8100000 -o ex_CEU.out"
  val simu37P = "hapgen2 -l refpath/ALL.chr15.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend -m refpath/genetic_map_chr15_combined_b37.txt -h refpath/ALL.chr15.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes -dl 8050000 1 1 1 -n 500 500 -int 8000000 8100000 -o ex_CEU.out"
  val mergeFile = "gtool -M --g ex_CEU.out.controls.gen ex_CEU.out.cases.gen --s ex_CEU.out.controls.sample ex_CEU.out.cases.sample  --og ex_CEU.out.gen  --os ex_CEU.out.sample --log ex_CEU.out.log"
  val snpG = "/Users/qunliu/VEGAS2offline/VEGAS2database/1000GEURO/0kbloc/geneset/allgene1"
  val subsetFileP = "gtool -S --g ex_CEU.out.gen  --s ex_CEU.out.sample --og ex_CEU.gen  --inclusion ex_CEU_rsid.txt"
  val ptest = "snptest -data ex_CEU.out.gen ex_CEU.out.sample -o ex.out -frequentist 1 -method score -pheno pheno"
  val modifyPheno = "awk 'NR<3{print $0}NR>2{print $1,$2,$3,$4+1}' ex_CEU.out.sample > ex_CEU.sample && mv ex_CEU.out0.sample ex_CEU.out.sample"
  var toPedP = "gtool -G --g ex_CEU.gen --s ex_CEU.sample --chr 10 --phenotype pheno --ped ex_CEU.ped --map ex_CEU.map "
  val tobed = "plink --file ex_CEU --out ex_CEU --make-bed"
  val snpTest = "plink --bfile ex_CEU --out ex_CEU --assoc --allow-no-sex"
  //list file have 4 columns: chromosome, start, end, name, same as glist file of VEGAS2
  val makeSet = "plink --file ex_CEU.out --make-set ex_CEU.glist --write-set"
  val vegas2v2 = "vegas2v2 -G -snpandp example.txt -custom fullexample -glist example.glist -out example_vegas2out"

  def writeGl(gl:Array[String]) {
    val pw = new PrintWriter(new File(pPath+"genelist.txt"))
    gl.foreach(pw.println)
    pw.close()
  }

  def setchr(c:Int) = {
    chr = c
  }
  def gFile = pPath+dPath+aGene+chr
  def sFile = pPath+dPath+sGene+chr
  def gSet(gfil:String) = scala.io.Source.fromFile(gfil).getLines().toSet
  def checkName(g:String) ={
    val gset = gSet(gFile)
    gset.contains(g)
  }
  def getSnps(g:String,c:Int) = {
    setchr(c)
    require(checkName(g))
    scala.io.Source.fromFile(sFile).getLines().map(_.split(" ")).filter(_.contains(g)).map(_.apply(1)).toArray
  }

  def setComm = {
    val pvRs = "rs.txt"
    val glfil = "genelist.txt"
    val comm1 = comm.replace("pathtovegas/",pPath).replace("test_vegas2input.txt",pvRs).replace("genelist.txt",glfil)

  }
  def snpVal(x:Array[Float],c:Float = 0.9f):Int = {
    //val b = Array(0,1,2)
    val a:Int = if(x(0)>c) 0 else if(x(1)>c) 1 else if(x(2) >c) 2 else -1
    return(a)
  }
  def getTheta[T](xv:DenseVector[T])(implicit n: Numeric[T])= {
    import org.apache.commons.math3.distribution.NormalDistribution
    import n._
    val xv1 = xv.map(_.toFloat).toArray
    val len = xv1.length
    val nom = new NormalDistribution(0, 1)
    var thetaVal = nom.sample(len).map(_.toFloat)
    while (breeze.numerics.abs(calculation.pearsonCorr(xv1, thetaVal)) > 0.001) {
      thetaVal = nom.sample(len).map(_.toFloat)
    }
    val thetaV = DenseVector[Float](thetaVal)
    val mv = meanAndVariance(thetaV)
    (thetaV - mv.mean.toFloat) / sqrt(mv.variance).toFloat
  }
  def repPheno(sampleF:String = gPms.tp+"ex_CEU.out.sample",pheno:Array[Float],outf:String) = {
    val sampl = scala.io.Source.fromFile(sampleF).getLines.map(_.split(" ")).toArray
    val phenoVal = sampl.slice(0,2).map(_(3)) ++ pheno.map(_.toString)
    val writer = new PrintWriter(new FileWriter(outf))
    sampl.indices.toArray.map{i => sampl(i).slice(0,3) :+ phenoVal(i)}.foreach( i => writer.println(i.mkString("\t")))
    writer.close()
  }
  def setPheno(h:Float = 0.05f,num:Int = 0,pca:Boolean = true)(X:DenseMatrix[Float]):DenseMatrix[Float] = {
    val X1 = convert(X,Double)
    val pcs = if (pca) princomp(X1).scores else X1
    val mav = breeze.stats.meanAndVariance(pcs(::, num))
    var pc1 = convert((pcs(::, num) - mav.mean) / sqrt(mav.variance),Float)
    val theta = getTheta(pc1)
    (sqrt(h) * pc1 + sqrt(1 - h) * theta).toDenseMatrix.t
  }
  def setPheno2(h:Float = 0.05f,num:Int = 2,pca:Boolean = true)(X:DenseMatrix[Float]):DenseMatrix[Float] = {
    val X1 = convert(X,Double)
    val pcs = if (pca) princomp(X1).scores else X1
    val mav = breeze.stats.meanAndVariance(pcs(::, 0 until num))
    var pc1 = convert((pcs(::, 0 until num) - mav.mean) / sqrt(mav.variance),Float)
    val theta = getTheta(pc1(::,0))
    (sqrt(h/num.toFloat) * pc1(::,0) + sqrt(h/num.toFloat) * pc1(::,1)+ sqrt(1 - h) * theta).toDenseMatrix.t
  }
  def simuFgene(glist:Array[String]) = {
    val chr = glist(0)
    val stt = glist(1).toInt
    val end = glist(2).toInt
    val gname = glist(3)
    makeSetF(glist)
    val snpL = fileOper.toArrays(gPms.gp+"ALL.chr"+chr+".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend"," ").drop(1).filter(i => i(1).toInt > stt & i(1).toInt < end).toArray
    val snpp = snpL(snpL.length/2+1).apply(1)
    val simu = simu37P.replace("chr10","chr"+chr).replace("8000000 8100000",stt + " "+end).replace("8050000",snpp).replace("refpath/",gPms.gp).replace("ex_CEU",gname)
    val comm1 = Process(simu,new File(gPms.tp)).!
    val comm2 = Process(mergeFile.replace("ex_CEU",gname),new File(gPms.tp)).!

    getSnp(glist)
    Process(subsetFileP.replace("ex_CEU",gname),new File(gPms.tp)).!

  }
  def makeSetF(glist:Array[String],outf:String = gPms.tp+"ex_CEU.glist") = {
    val gname = glist(3)
    val writer = new PrintWriter(new FileWriter(outf.replace("ex_CEU",gname)))
    writer.println(glist.mkString(" "))
    writer.close()
  }

  def getPvalF(fil:String = gPms.tp+"ex_CEU.out.qassoc") = {
    val qa = "[q]?assoc".r
    val outf =qa.replaceFirstIn(fil,"txt")
    val writer = new PrintWriter(new FileWriter(outf))
    fileOper.toArrays(fil, " ").drop(1).map(_.filter(_.length >0)).filter(i => i(i.length -1) != "NA").map(i => i(1) + "\t"+i(i.length-1)).foreach(writer.println(_))
    writer.close()
  }
  def getSnp(glist:Array[String],outf:String = gPms.tp+"ex_CEU_rsid.txt") = {
    val chr = glist(0)
    val stt = glist(1).toInt
    val end = glist(2).toInt
    val gname = glist(3)
    val snpL = fileOper.toArrays(gPms.op+"snp6annoNewchr"+chr+".txt").filter(i => i(2).toInt > stt & i(2).toInt < end).map(_(2)).toSet
    val snpp = fileOper.toArrays(gPms.tp+gname+".out.gen"," ").filter(i => snpL.contains(i(2))).map(_(1)).toArray
    val writer = new PrintWriter(new FileWriter(outf.replace("ex_CEU",gname)))
    snpp.foreach(writer.println)
    writer.close()
  }
  def getVegas(ph:Array[Float],chr:String) = {
    repPheno(gPms.tp+"ex_CEU.out.sample",ph,gPms.tp+"ex_CEU.sample")

  }
  def vegas(glist:Array[String],n:Int = 1,spheno:DenseMatrix[Float] => DenseMatrix[Float] = setPheno()) = {

    val gname = glist(3)
    val snpCC = scala.io.Source.fromFile(gPms.tp+gname+".gen").getLines.map(_.split(" "))
    val Xx = snpCC.map(_.drop(5).map(_.toFloat).sliding(3,3).map(snpVal(_).toFloat).toArray).toArray
    val X = utils.Array2DM(Xx,false)
    val Y = spheno(X)
    repPheno(gPms.tp+gname+".out.sample",pheno = Y.toArray,outf = gPms.tp+gname+".sample")
    val toPed = toPedP.replace("chr 10","chr "+glist(0))
    val comm3 = Process(toPed.replace("ex_CEU",gname),new File(gPms.tp)).!
    val comm4 = Process(tobed.replace("ex_CEU",gname),new File(gPms.tp)).!
    val comm5 = Process(snpTest.replace("ex_CEU",gname), new File(gPms.tp)).!

    if (Y.toArray.toSet.size > Y.rows/2){
      getPvalF(gPms.tp+gname+".qassoc")
    }else{
      getPvalF(gPms.tp+gname+".assoc")
    }
    val vegas2 = vegas2v2.replace("fullexample",gPms.tp+gname).replace("example",gname)
    val comm6 = Process(vegas2,new File(gPms.tp)).!
    val pval = fileOper.toArrays(gPms.tp+gname+"_vegas2out.out"," ").drop(1).toArray.map(i =>Array(i(7).toFloat,i(9).toFloat)).apply(0)
    pval ++ plsCalc.gdofPlsPval(X,Y,n)._2
  }


}

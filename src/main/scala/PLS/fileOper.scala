package PLS



import java.io.{FileWriter, PrintWriter}

  import myParallel._

  import scala.collection.mutable.{ArrayBuffer, ListBuffer}
  import scala.io.Source
  import scala.reflect.ClassTag
  import scala.reflect.io.File
  import scala.util.Try
  //import scala.reflect.io.File

  /**
    * Created by liuqun on 6/27/16.
    */
  object fileOper {
    def toLines: String => Iterator[String] = f => scala.io.Source.fromFile(f).getLines
    //def getArray:(String,String) => Iterator[Array[String]] = (f,s) => getLines(f).map(_.split(s))
    def toArrays(f:String, s:String = "\t")(implicit ncol:Int = -1):Iterator[Array[String]] = {
      if(ncol<0)toLines(f).map(_.split(s))else toLines(f).map(_.split(s,ncol))
    }
    def filetoMap:String => Map[String, String] = f => toArrays(f).filter(_.length == 2).map(k => k(0) -> k(1)).toMap
    def filetoTuple(file:String) = {
      val GOtuple = toLines(file).drop(3).map(_.split("\t")).toList.groupBy(_(0)).map{ case (k,v) => (k, v(0)(1),v(0)(2),v.map(_(4)))}
      GOtuple
    }

    def csvToTxt(csvf:String,dropl:Int) = {
      val txtf = csvf.replace(".csv",".txt")
      val out = new PrintWriter(new FileWriter(txtf))
      val line = toArrays(csvf,",").drop(dropl)
      line.foreach(i => out.println(i.mkString("\t")))
    }
    def fileFlip(fl:String):Array[(String,Array[String])] ={
      var idg = scala.io.Source.fromFile(fl).getLines.map(_.split("\t")).map(i => (i(0), i.drop(2))).toArray
      val gen = idg.flatMap(_._2).toSet.toArray.sorted
      val gens = gen.map(i => (i,scala.collection.mutable.ArrayBuffer.empty[String])).toMap
      idg.foreach(i => i._2.map(ii => gens(ii) += i._1))
      return(gens.toArray.map(i => (i._1,i._2.toArray)))
    }
    def select[T:ClassTag](in:Array[T],sel:Array[Int]):Array[T] = {
      val out = new ArrayBuffer[T]()
      val selLen = sel.length
      var i = 0
      while (i < selLen){
        out += in(sel(i))
      }
      out.toArray
      //sel.map(in)
    }
    def intersectCols(x:Array[String],y:Array[String]):(Array[Int],Array[Int]) = {
      val xy = x.intersect(y).toSet
      val x1 = Array.tabulate(x.length){ i => (x(i),i) }.filter(xy contains _._1).groupBy(_._1).
        map{case(k,v) => v(0)._2}.toArray.sorted
      val y1 = Array.tabulate(y.length){ i => (y(i),i) }.filter(xy contains _._1).groupBy(_._1).
        map{case(k,v) => v(0)._2}.toArray.sorted
      val nmap = x1.map(x).zipWithIndex.toMap
      val yn = y1.map(i => nmap(y(i)))
      val ynm = y1.zip(yn).sortBy(_._2).map(_._1)
      (x1,ynm)
    }
    def intersectCols(x:Array[String],y:Array[String],z:Array[String]):(Array[Int],Array[Int],Array[Int]) = {
      val xy = x.intersect(y).intersect(z).toSet
      val x1 = Array.tabulate(x.length){ i => (x(i),i) }.filter(xy contains _._1).groupBy(_._1).
        map{case(k,v) => v(0)._2}.toArray.sorted
      val y1 = Array.tabulate(y.length){ i => (y(i),i) }.filter(xy contains _._1).groupBy(_._1).
        map{case(k,v) => v(0)._2}.toArray.sorted
      val z1 = Array.tabulate(z.length){ i => (z(i),i) }.filter(xy contains _._1).groupBy(_._1).
        map{case(k,v) => v(0)._2}.toArray.sorted
      val nmap = x1.map(x).zipWithIndex.toMap
      val yn = y1.map(i => nmap(y(i)))
      val ynm = y1.zip(yn).sortBy(_._2).map(_._1)
      val zn = z1.map(i => nmap(z(i)))
      val znm = z1.zip(zn).sortBy(_._2).map(_._1)
      (x1,ynm,znm)
    }
    def sepFile(fil:String,chrCol:Int,sortCol:Int = -1) = {
      val f = fileOper.toArrays(fil).toArray.groupBy(_(chrCol))
      for ((k,v) <- f){
        val writer = new PrintWriter(new FileWriter(fil.replace(".txt","Tchr"+k+".txt")))
        val outT = v.filter(i=> Try(i(sortCol).toInt).isSuccess)
        val out = if (sortCol >0) outT.sortBy(i => i(sortCol).toInt) else outT
        out.foreach(i => writer.println(i.mkString("\t")))
        writer.close
      }
    }
    def parseAndSortFile(inFile:String, outFile:String,ratingCol:Array[Int],head:Boolean = true): Unit = {
      val writer = new PrintWriter(new FileWriter(outFile))
      var values = toArrays(inFile)//Source.fromFile(inFile).getLines().map { x => x.split("\t") }
      if (head) {
        val header = values.next
        writer.write(header.mkString("\t") + "\n")
      }
      //    val ratingColumn = header.indexOf(ratingCol)
      //    if (ratingColumn == -1) {
      //      println("could not find ratings column");
      //    } else {

      val value = values.toArray
      val sortind = if (ratingCol.length == 1) value.map(_ (ratingCol(0))).zipWithIndex.sortBy(_._1).map(_._2)
        else value.map(select(_,ratingCol)).zipWithIndex.sortBy(i => (i._1.apply(0),i._1.apply(1))).map(_._2)
      val len = sortind.length
      val outArray = sortind.map(value)// Array.fill[String](len)("")
      outArray.foreach { l => writer.write(l.mkString("\t") + "\n") }
      //    }
      writer.close()
    }
    def FileApppend(inFile:String,outFile:String,addFile:String,keyCol:String,valueCol:String): Unit= {
      val lines = Source.fromFile(inFile).getLines().toArray
      val reflines = Source.fromFile(addFile).getLines().toArray
      val refheader = reflines.head.split("\t")
      val header = lines.head.split("\t")
      val tailer = lines.tail.map{l => l.split("\t")}
      val refkeyColumn = header.indexOf(keyCol)
      val refvalueColumn = header.indexOf(valueCol)
      val keyColumn = header.indexOf(keyCol)
      val remainSnp = lines.tail.map(_(keyColumn)).toSet -- reflines.tail.map(_(refkeyColumn)).toSet
      remainSnp.toList.map(println)
      if(remainSnp.toList.length == 0){
        val refmap = reflines.map(l=> l(refkeyColumn).toString() -> l(refvalueColumn).toString()).toMap
        val writer = new PrintWriter(new FileWriter(outFile))
        if (keyColumn == -1) {
          println("could not find ratings column");
        }else {
          writer.write(header.mkString("\t")+"\t"+valueCol + "\n")
          val values = tailer.map { x  => writer.write(x.mkString("\t")+"\t"+ refmap(x(refkeyColumn)) + "\n")}
        }
      }
    }
    def caseclasstoFile(tabfile:String,myList:List[Any]) = {

      if(!File(tabfile).exists)File(tabfile).createFile()
      val out = new PrintWriter(new FileWriter(tabfile))

      try{
        myList.map(out.println(_))
      }
      finally{
        out.close
      }
    }
    def seperateAsChr(filenm:String = gPms.op+"lggsnpti.txt",map:Map[String,Array[String]],outFile:String = gPms.op+"Tcga_lgg_snpChrn.txt",col:Int = 0,filter:Float = 1.0f):Unit={
      val writeChr1 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr1")))
      val writeChr2 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr2")))
      val writeChr3 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr3")))
      val writeChr4 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr4")))
      val writeChr5 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr5")))
      val writeChr6 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr6")))
      val writeChr7 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr7")))
      val writeChr8 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr8")))
      val writeChr9 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr9")))
      val writeChr10 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr10")))
      val writeChr11 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr11")))
      val writeChr12 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr12")))
      val writeChr13 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr13")))
      val writeChr14 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr14")))
      val writeChr15 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr15")))
      val writeChr16 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr16")))
      val writeChr17 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr17")))
      val writeChr18 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr18")))
      val writeChr19 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr19")))
      val writeChr20 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr20")))
      val writeChr21 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr21")))
      val writeChr22 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr22")))
      val writeChrX = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrX")))
      val writeChrY = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrY")))
      val writeChrMT = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrMT")))

      val lines =scala.io.Source.fromFile(filenm).getLines().map(_.split("\t"))
      val name = lines.next.map(_.slice(0,15))
      val rppa = scala.io.Source.fromFile(gPms.op+"GBMLGG.rppa.txt").getLines().map(_.split("\t"))
      val rppaname = rppa.next.map(_.slice(0,15))
      val commname = name.toSet.intersect(rppaname.toSet).toArray
      val nameMap = name.zipWithIndex.toMap
      val rppaMap = rppaname.zipWithIndex.toMap
      val ind = commname.map(nameMap(_)).sorted
      val rppaind = commname.map(rppaMap(_)).sorted
      val outF = new PrintWriter(new FileWriter(gPms.op+"aktrppawithSnp.txt"))
      outF.println(select(rppaname,rppaind.+:(0)).mkString("\t"))
      rppa.filter(_(0).matches(".*Akt_p.*")).toArray.map(select(_,rppaind.+:(0)).mkString("\t")).foreach(outF.println(_))
      outF.close()

      for (line <- lines) {

        val linPre =  select(line,ind)


        val vari =  if (filter == 1.0f) 2.0f else calculation.Variance(linPre.map(_.toFloat))
        if (vari > filter) {

          val lin = linPre.mkString("\t")

          val mat = map.get(line.apply(col))
          mat match {
            case Some(x) => x(2) match {
              case "1" => writeChr1.write(x.mkString("\t") + "\t" + lin + "\n")
              case "2" => writeChr2.write(x.mkString("\t") + "\t" + lin + "\n")
              case "3" => writeChr3.write(x.mkString("\t") + "\t" + lin + "\n")
              case "4" => writeChr4.write(x.mkString("\t") + "\t" + lin + "\n")
              case "5" => writeChr5.write(x.mkString("\t") + "\t" + lin + "\n")
              case "6" => writeChr6.write(x.mkString("\t") + "\t" + lin + "\n")
              case "7" => writeChr7.write(x.mkString("\t") + "\t" + lin + "\n")
              case "8" => writeChr8.write(x.mkString("\t") + "\t" + lin + "\n")
              case "9" => writeChr9.write(x.mkString("\t") + "\t" + lin + "\n")
              case "10" => writeChr10.write(x.mkString("\t") + "\t" + lin + "\n")
              case "11" => writeChr11.write(x.mkString("\t") + "\t" + lin + "\n")
              case "12" => writeChr12.write(x.mkString("\t") + "\t" + lin + "\n")
              case "13" => writeChr13.write(x.mkString("\t") + "\t" + lin + "\n")
              case "14" => writeChr14.write(x.mkString("\t") + "\t" + lin + "\n")
              case "15" => writeChr15.write(x.mkString("\t") + "\t" + lin + "\n")
              case "16" => writeChr16.write(x.mkString("\t") + "\t" + lin + "\n")
              case "17" => writeChr17.write(x.mkString("\t") + "\t" + lin + "\n")
              case "18" => writeChr18.write(x.mkString("\t") + "\t" + lin + "\n")
              case "19" => writeChr19.write(x.mkString("\t") + "\t" + lin + "\n")
              case "20" => writeChr20.write(x.mkString("\t") + "\t" + lin + "\n")
              case "21" => writeChr21.write(x.mkString("\t") + "\t" + lin + "\n")
              case "22" => writeChr22.write(x.mkString("\t") + "\t" + lin + "\n")
              case "X" => writeChrX.write(x.mkString("\t") + "\t" + lin + "\n")
              case "Y" => writeChrY.write(x.mkString("\t") + "\t" + lin + "\n")
              case "MT" => writeChrMT.write(x.mkString("\t") + "\t" + lin + "\n")
              //        case chrMTmatch(line) => writeChrMT.write(line+"\n")
              case _ => //{
            }
            case None =>
          }
        }
      }

      writeChr1.close()
      writeChr2.close()
      writeChr3.close()
      writeChr4.close()
      writeChr5.close()
      writeChr6.close()
      writeChr7.close()
      writeChr8.close()
      writeChr9.close()
      writeChr10.close()
      writeChr11.close()
      writeChr12.close()
      writeChr13.close()
      writeChr14.close()
      writeChr15.close()
      writeChr16.close()
      writeChr17.close()
      writeChr18.close()
      writeChr19.close()
      writeChr20.close()
      writeChr21.close()
      writeChr22.close()
      writeChrX.close()
      writeChrY.close()
      writeChrMT.close()

    }
    def simpleSeperateChr(filenm:String = gPms.op+"geneLoc.txt",outFile:String = gPms.op+"geneLocChrn.txt",col:Int = 0):Unit={
      val writeChr1 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr1")))
      val writeChr2 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr2")))
      val writeChr3 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr3")))
      val writeChr4 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr4")))
      val writeChr5 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr5")))
      val writeChr6 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr6")))
      val writeChr7 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr7")))
      val writeChr8 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr8")))
      val writeChr9 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr9")))
      val writeChr10 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr10")))
      val writeChr11 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr11")))
      val writeChr12 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr12")))
      val writeChr13 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr13")))
      val writeChr14 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr14")))
      val writeChr15 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr15")))
      val writeChr16 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr16")))
      val writeChr17 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr17")))
      val writeChr18 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr18")))
      val writeChr19 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr19")))
      val writeChr20 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr20")))
      val writeChr21 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr21")))
      val writeChr22 = new PrintWriter(new FileWriter(outFile.replace("Chrn","chr22")))
      val writeChrX = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrX")))
      val writeChrY = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrY")))
      val writeChrMT = new PrintWriter(new FileWriter(outFile.replace("Chrn","chrMT")))
      val lines =scala.io.Source.fromFile(filenm).getLines().map(_.split("\t"))
      for (line <- lines) {

        line(col) match {
          case "1" => writeChr1.write(line.mkString("\t") + "\n")
          case "2" => writeChr2.write(line.mkString("\t") + "\n")
          case "3" => writeChr3.write(line.mkString("\t") + "\n")
          case "4" => writeChr4.write(line.mkString("\t") + "\n")
          case "5" => writeChr5.write(line.mkString("\t") + "\n")
          case "6" => writeChr6.write(line.mkString("\t") + "\n")
          case "7" => writeChr7.write(line.mkString("\t") + "\n")
          case "8" => writeChr8.write(line.mkString("\t") + "\n")
          case "9" => writeChr9.write(line.mkString("\t") + "\n")
          case "10" => writeChr10.write(line.mkString("\t") + "\n")
          case "11" => writeChr11.write(line.mkString("\t") + "\n")
          case "12" => writeChr12.write(line.mkString("\t") + "\n")
          case "13" => writeChr13.write(line.mkString("\t") + "\n")
          case "14" => writeChr14.write(line.mkString("\t") + "\n")
          case "15" => writeChr15.write(line.mkString("\t") + "\n")
          case "16" => writeChr16.write(line.mkString("\t") + "\n")
          case "17" => writeChr17.write(line.mkString("\t") + "\n")
          case "18" => writeChr18.write(line.mkString("\t") + "\n")
          case "19" => writeChr19.write(line.mkString("\t") + "\n")
          case "20" => writeChr20.write(line.mkString("\t") + "\n")
          case "21" => writeChr21.write(line.mkString("\t") + "\n")
          case "22" => writeChr22.write(line.mkString("\t") + "\n")
          case "X" => writeChrX.write(line.mkString("\t") + "\n")
          case "Y" => writeChrY.write(line.mkString("\t") + "\n")
          case "MT" => writeChrMT.write(line.mkString("\t") + "\n")
          //        case chrMTmatch(line) => writeChrMT.write(line+"\n")
          case _ => //{
        }
      }
        writeChr1.close()
        writeChr2.close()
        writeChr3.close()
        writeChr4.close()
        writeChr5.close()
        writeChr6.close()
        writeChr7.close()
        writeChr8.close()
        writeChr9.close()
        writeChr10.close()
        writeChr11.close()
        writeChr12.close()
        writeChr13.close()
        writeChr14.close()
        writeChr15.close()
        writeChr16.close()
        writeChr17.close()
        writeChr18.close()
        writeChr19.close()
        writeChr20.close()
        writeChr21.close()
        writeChr22.close()
        writeChrX.close()
        writeChrY.close()
        writeChrMT.close()
      }


      case class GENE(gene_id:String,stable_id:String,description:String,Chr:String,Start:Int,End:Int,Snp:List[Any])
    case class SNP(variation_id: String, variation_name: String, chromosome: String, loc: Int)

    def fileSortSplit(filen:String,head:Boolean = true) = {
      val affymap = scala.io.Source.fromFile(filen).getLines.map(_.split("\t")).toArray.map(i => i(0) -> i).toMap
      val nam = Array(1 to 22:_*).map(_.toString) ++ Array("X","Y","MT")

      fileOper.seperateAsChr(gPms.op+"lggsnpti.txt",affymap,filter = 0.003f)

      nam.foreach(i => fileOper.parseAndSortFile(gPms.op+"Tcga_lggSNPCHR"+i+".txt", gPms.op+"Tcga_lggSNPCHR"+i+"_out.txt",Array(3),false))
    }
    def FindSNPinRange(file:String,amp:Int):List[Any]= {
      val Gene = Source.fromFile(file).getLines().toArray
      val GeneStartCol = Gene.head.split("\t").indexOf("seq_region_start")
      val GeneEndCol = Gene.head.split("\t").indexOf("seq_region_end")
      val GeneDescriptCol = Gene.head.split("\t").indexOf("description")
      val GeneIdCol = Gene.head.split("\t").indexOf("gene_id")
      val GeneStableCol = Gene.head.split("\t").indexOf("stable_id")
      val GeneChrCol = Gene.head.split("\t").indexOf("chromosome")
      val GeneValues = Gene.tail.map { x => x.split("\t") }

      var chrRef = "1"
      var ChrSnp = Source.fromFile(gPms.op+"Sorted_SNPCHR" + chrRef + ".txt").getLines().toArray
      var snpInx = Source.fromFile(gPms.op+"/Tcga_lggSNPCHR" + chrRef + "_out.txt").getLines().map(_.split("\t")).map(_(3).toInt).toArray.zipWithIndex
      val SnpLocCol = ChrSnp.head.split("\t").indexOf("seq_region_start")
      val SnpIdCol = ChrSnp.head.split("\t").indexOf("variation_id")
      val SnpNameCol = ChrSnp.head.split("\t").indexOf("variation_name")
      val SnpChrCol = ChrSnp.head.split("\t").indexOf("chromosome")
      println(SnpChrCol,SnpLocCol,SnpIdCol,SnpNameCol)
      var SnpValues = ChrSnp.tail.map { x => x.split("\t") }
      var loc = 0
      var i = 0
      var geneList = new ListBuffer[GENE]
      for (l <- 0 to GeneValues.length-1){
        val chr = GeneValues(l)(GeneChrCol)
        val start = GeneValues(l)(GeneStartCol).toInt - amp
        val end = GeneValues(l)(GeneEndCol).toInt + amp
        var snps = new ListBuffer[SNP]
        if (chr != chrRef) {
          chrRef = chr
          var ChrSnp = Source.fromFile(gPms.op+"Sorted_SNPCHR" + chrRef + ".txt").getLines().toArray
          var SnpValues = ChrSnp.tail.map { x => x.split("\t") }
          i = 0
          loc = 0
        }
        while ((i < SnpValues.length) && (SnpValues(i)(SnpLocCol).toInt < end)) {
          if (SnpValues(i)(SnpLocCol).toInt < start) {
            loc = i
          }
          else {
            snps += new SNP(SnpValues(i)(SnpIdCol), SnpValues(i)(SnpNameCol),
              SnpValues(i)(SnpChrCol), SnpValues(i)(SnpLocCol).toInt)
          }
          if (i == 5)println(SnpValues(i)(SnpIdCol),SnpValues(i)(SnpNameCol),SnpValues(i)(SnpChrCol),SnpValues(i)(SnpLocCol).toInt)
          i += 1
        }
        i = loc
        if (snps.length > 0) println(snps)
        geneList += new GENE(GeneValues(l)(GeneIdCol),GeneValues(l)(GeneStableCol),
          GeneValues(l)(GeneDescriptCol),GeneValues(l)(GeneChrCol),
          GeneValues(l)(GeneStartCol).toInt,GeneValues(l)(GeneEndCol).toInt,snps.toList)
      }
      val genes = geneList.toList
      genes
      //List.range(1, 30).map(i => println(genes(i)))
    }

    case class GeneSnp(stable_id:String,Chr:String,Start:Int,End:Int,snpStart:Int,snpEnd:Int)
    def FindIndexinRange(gfile:String = gPms.op+gPms.glf,sfile:String = gPms.op+gPms.slf,ofile:String =gPms.op+gPms.rf,
                         amp:Int = 0,ampEnd:Int = -1,byend:Boolean = true)= {
      val gene = scala.io.Source.fromFile(gfile).getLines().map(_.split("\t")).drop(0).toArray

      val snpInx = scala.io.Source.fromFile(sfile).getLines().drop(0).map(_.split("\t")).toArray
      var loc1 = 0
      var i = 0
      var l = 0
      val GeneChrCol = 0
      val GeneStartCol = 1
      val GeneEndCol = 2
      val GeneNameCol = 3
      var geneList = new ArrayBuffer[GeneSnp]
      val geneL = gene.length
      val snpL = snpInx.length -1
      while (l < geneL) {
        var chr = gene(l)(0).toInt
        val starts = gene(l)(1).toInt
        val start = starts- amp
        val ends = if (byend) gene(l)(2).toInt else starts
        val end = if (ampEnd < 0 ) ends + amp else ends +ampEnd
        while (i < snpL && (snpInx(i)(1).toInt == chr && snpInx(i)(2).toInt < start)|snpInx(i)(1).toInt < chr ) {
          i += 1
        }
        loc1 = i
        while (snpInx(loc1)(2).toInt < end && snpInx(loc1)(1).toInt == chr && loc1 < snpL){
          loc1 += 1
        }
        if (loc1 == snpL) loc1 +=1
        geneList += new GeneSnp( gene(l)(GeneNameCol), gene(l)(GeneChrCol),
          gene(l)(GeneStartCol).toInt, gene(l)(GeneEndCol).toInt, i,loc1)
        l += 1
      }
      val genes = geneList.toArray
      val out = new PrintWriter(new FileWriter(ofile))
      val outputs = genes.map(_.productIterator.mkString("\t"))
      outputs foreach (out.println(_))
      out.close()
    }
    def timeForFile:String= {
      val today = java.time.LocalDate.now.toString.split("-")
      val totime = java.time.LocalTime.now.toString.split("\\.")(0).split(":")
      "_"+(today++totime).mkString("_")
    }
}

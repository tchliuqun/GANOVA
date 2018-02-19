package PLS


import scala.io.Source
import scala.util.matching.Regex
import java.io._

//import pathway.FileOper
/**
  * Created by liuqun on 6/25/16.
  */
object ensembl {

  val myquery =
    """Select b.name as chromosome, a.*
                   from `homo_sapiens_core_84_37`.`gene` a
                   inner join `homo_sapiens_core_84_37`.`seq_region` b
                   on a.seq_region_id = b.seq_region_id
                   inner join `homo_sapiens_core_84_37`.`coord_system` c
                   on b.coord_system_id = c.coord_system_id
                   where a.biotype = 'protein_coding'
                   and c.version = 'GRCh37'
                   and c.name = 'chromosome' """
  val myquery1 =
    """Select b.name as chromosome, b.coord_system_id, a.*
                   from `homo_sapiens_core_84_37`.`gene` a
                   inner join `homo_sapiens_core_84_37`.`seq_region` b
                   on a.seq_region_id = b.seq_region_id
                   where a.biotype = 'protein_coding' """
  val myquery20 =
    """
        Select coalesce(b1.name, b.name) as chromosome, b1.name as chromosomeb,b.name as location, a.*
        from `homo_sapiens_core_84_37`.`gene` a
        inner join `homo_sapiens_core_84_37`.`seq_region` b
        on a.seq_region_id = b.seq_region_id
        left join `homo_sapiens_core_84_37`.`assembly_exception` s
        on b.seq_region_id = s.seq_region_id
        left join `homo_sapiens_core_84_37`.`seq_region` b1
        on s.exc_seq_region_id = b1.seq_region_id
        inner join `homo_sapiens_core_84_37`.`coord_system` c
        on b.coord_system_id = c.coord_system_id
        where a.biotype = 'protein_coding'
        and c.version = 'GRCh37'
        and c.name = 'chromosome'
        order by chromosome, seq_region_start
    """
  val myquery2 =
    """
      SELECT DISTINCT s1.name AS patch,
       s2.name AS chromosome
       FROM `homo_sapiens_core_84_37`.`assembly_exception` AS a,
      `homo_sapiens_core_84_37`.`seq_region` AS s1,
      `homo_sapiens_core_84_37`.`seq_region` AS s2
       WHERE a.seq_region_id = s1.seq_region_id
       AND a.exc_seq_region_id = s2.seq_region_id"""
  //AND a.exc_type IN ('PATCH_FIX', 'PATCH_NOVEL')"""
  val myquery3 = "select * from `homo_sapiens_core_84_37`.`seq_region`"

  val myquery30 =
    "select * from `information_schema`.`tables` where TABLE_NAME ='variation'"// TABLE_SCHEMA = 'homo_sapiens_core_84_37'"// and "
  val myquery4 =
    "select * from `homo_sapiens_variation_84_37`.`variation_set_variation` where variation_set_id = 4"// between 27520 and 27529  "
  val myquery40 = """
                   select o.linkage_type AS evid, b.*,g.* from `homo_sapiens_core_84_37`.`xref` AS g,
                  `homo_sapiens_core_84_37`.`ontology_xref` AS o,
                  `homo_sapiens_core_84_37`.`object_xref` AS b
                  WHERE o.object_xref_id = b.object_xref_id
                  AND o.source_xref_id = g.xref_id
                  """
  val myquery5 = "select * from `homo_sapiens_core_84_37`.`assembly_exception`"
  val myquery400 = """
                   SELECT b.*,o.source_xref_id AS source from
                   `homo_sapiens_core_84_37`.`ontology_xref` AS o,
                   `homo_sapiens_core_84_37`.`object_xref` AS b
                   WHERE o.object_xref_id = b.object_xref_id
                   select * from
                  `homo_sapiens_core_84_37`.`ontology_xref` AS o,
                  `homo_sapiens_core_84_37`.`object_xref` AS b
                  WHERE o.object_xref_id = b.object_xref_id
                  """
  // where variation_set_id = 4"// between 27520 and 27529  "
  val myquery6 = "select * from `homo_sapiens_core_84_37`.`coord_system`"
  // where variation_set_id = 4"// between 27520 and 27529  "
  val myquery7 =
    """
                    Select coalesce(b1.name, c.name) as chromosome, b1.name as chromosomeb,c.name as location, a.*, b.variation_set_id as uniset_id
                   from `homo_sapiens_variation_84_37`.`variation_feature` a
                   inner join `homo_sapiens_variation_84_37`.`variation_set_variation` b
                   on a.variation_id = b.variation_id
                   inner join `homo_sapiens_variation_84_37`.`seq_region` c
                   on a.seq_region_id = c.seq_region_id
                   left join `homo_sapiens_core_84_37`.`assembly_exception` s
                   on c.seq_region_id = s.seq_region_id
                   left join `homo_sapiens_variation_84_37`.`seq_region` b1
                   on s.exc_seq_region_id = b1.seq_region_id
                   where b.variation_set_id = 4"""
  //                   and (b1.name = '1' or c.name = '1')
  //                   order by seq_region_start
  //                     """
  //var querydtable = " select * from `information_schema`.`tables` where TABLE_SCHEMA = 'homo_sapiens_variation_84_37'"// and TABLE_NAME ='variation'"
  //var querydtable = "select * from `homo_sapiens_core_84_37`.`seq_region` where seq_region_id between 27520 and 27529  "
  //  val rs = MysqlFunc.connect(MysqlAddress.ENSEMBLHG19HOST, myquery2)
  //  val chrPatch = MysqlFunc.RS2Map(rs)
  //  for (i <- 1 to 22) chrPatch += (i.toString -> i.toString)
  //  chrPatch += ("MT" -> "MT")
  //  chrPatch += ("X" -> "X")
  //  chrPatch += ("Y" -> "Y")
  //  chrPatch.foreach(println)
  //  println(chrPatch("HG1453_PATCH"))
  def getGeneTable() = {
    val mysqlFunc = new mysqlFunc
    val rs2 = mysqlFunc.connect(myquery40)
    mysqlFunc.RS2File(rs2, gPms.op+"GeneTGO.txt") //,chrPatch)
  }

  val Chr = List.range(1, 23).map(_.toString)++List("X","Y","MT")
  val inFilenames = Chr.map(n => s"SNPCHR$n.txt")
  //FileOper.parseAndSortFile("SNPCHR22.TXT","Sorted_SNPCHR22.txt","seq_region_start")
  //  val genes = FileOper.FindSNPinRange("GeneTable.txt",0)
  //  List.range(0, 10).map(i => println(genes(i)))

  //Chr.map(n => FileOper.parseAndSortFile(s"SNPCHR$n.TXT",s"Sorted_SNPCHR$n.txt","seq_region_start"))
}



//val RSStream = new Iterator[ResultSet] { def hasNext = rs.next ; def next = rs }.toStream
//RSStream.foreach()
//sqlconn.RS2File(rs2,"generow.txt",loc = -1,chrPatch)
//  val genes = sql"""SELECT * FROM gene LIMIT 5""".as[List[String]]
//  for (gene <- genes){println(gene)}

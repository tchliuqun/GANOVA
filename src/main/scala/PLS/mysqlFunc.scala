package PLS


  import java.io.{FileWriter, PrintWriter}
  import java.sql._
  import java.sql.DriverManager

  import com.mysql.jdbc.jdbc2.optional.MysqlDataSource

  import scala.language.experimental.macros
  import scala.reflect.io.File


class mysqlFunc {
    val MYSQLDRIVER = "com.mysql.jdbc.Driver"
    //var querydtable = " select * from `information_schema`.`tables` where TABLE_SCHEMA = 'homo_sapiens_core_84_37'"
    var usr = "anonymous"
    var pwd = ""
    var Host = MysqlAddress.ENSEMBLHG19HOST
    def url = "jdbc:mysql://" + Host
    def getDataSource():MysqlDataSource = {
      val dataSource = new MysqlDataSource()
      dataSource.setUser(usr)
      dataSource.setPassword(pwd)
      dataSource.setUrl(url)
      dataSource
    }


    //def getScriptEngine(): ScriptEngine = {
    //  import scala.collection.JavaConversions._

    //  val factories = javax.imageio.spi.ServiceRegistry.lookupProviders(classOf[ScriptEngineFactory])
    //  val scalaEngineFactory = factories.find(_.getEngineName == "Scala Scripting Engine")

    //  scalaEngineFactory.map(_.getScriptEngine).getOrElse(
    //    throw new AssertionError("Scala Scripting Engine not found")
    //  )
    //}

    // create a scala engine

    //val scriptEngine: ScriptEngine = getScriptEngine()
    //val myname = "name"


    //  val say = "hello";

    //  val b = scriptEngine.getBindings(ScriptContext.ENGINE_SCOPE);
    //  b.put("obj", say)

    //  val writer = new StringWriter();
    //  scriptEngine.getContext().setWriter(writer);

    //  scriptEngine.eval(code.toString(), b)
    //  writer.toString.trim()

    def connect(myquery:String):ResultSet = {
      val dataSource = getDataSource()
      val conn = dataSource.getConnection();
      val stmt:Statement = conn.createStatement();
      val rs:ResultSet = stmt.executeQuery(myquery);

      //    Class.forName(MYSQLDRIVER)
      //    val conn = DriverManager.getConnection(url, usr, pwd)
      //    var st = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE, ResultSet.CONCUR_READ_ONLY)
      //
      //    val rs = st.executeQuery(myquery)
      //    //conn.close()
      return rs
    }
    def RS2CC(rs:ResultSet) {
      import javax.script.{ScriptEngine, ScriptEngineManager}
      import scala.collection.mutable.ListBuffer


      val rsMeta = rs.getMetaData()
      val rscols = rsMeta.getColumnCount()
      //val engine: ScriptEngine = new ScriptEngineManager().getEngineByName("scala")
      //val settings = engine.asInstanceOf[scala.tools.nsc.interpreter.IMain].settings



      //settings.embeddedDefaults[SqlTest]

      //settings.usejavacp.value = true

      var code = new StringBuilder()
      //case class Gene()
      code ++= "case class Gene("


      for (i <- 1 to rscols) {
        //columnNames ++= rsMeta.getColumnName(i)
        println(rsMeta.getColumnLabel(i) + "\t" + rsMeta.getColumnClassName(i))
        code ++= rsMeta.getColumnLabel(i)
        if (rsMeta.getColumnClassName(i) == "java.lang.String") code ++= ":String"
        if (rsMeta.getColumnClassName(i) == "java.lang.Integer") code ++= ":Int"
        if (rsMeta.getColumnClassName(i) == "java.lang.Long") code ++= ":Long"
        if (rsMeta.getColumnClassName(i) == "java.lang.Boolean") code ++= ":Boolean"
        if (rsMeta.getColumnClassName(i) == "java.lang.Float") code ++= ":Float"
        if (rsMeta.getColumnClassName(i) == "java.lang.Double") code ++= ":Double"
        if (rsMeta.getColumnClassName(i) == "java.sql.Timestamp") code ++= ":java.util.Date"
        if (i < rscols) {
          code ++= ","
        } else {
          code ++= ")"
        }
      }
      println(code.toString)

      //engine.eval(code.toString)
      //engine.put("rscols:Int",rscols)
      // engine.eval("import java.sql.ResultSet")
      //case class Gene(engine.get("Gene"))
      var genes = new ListBuffer[Any]()
      while (rs.next){
        val genelist = (1 to rscols).map(n => rs.getObject(n)).toList
        genes += genelist
        //genelist.foreach(println)
        // engine.put("genelist:List[Any]",genelist)
        //engine.eval("""val gene = Gene.getClass.getMethods.find(x => x.getName == "apply" && x.isBridge).get.invoke(Gene, genelist map (_.asInstanceOf[AnyRef]): _*).asInstanceOf[Gene]""")
        //genes += engine.get("gene")//Gene.tupled((1 to rscols).map(i => rs.getObject(i)))
      }
      if (rs.isAfterLast()) rs.first()
      println(genes.toList.length)

    }

    def RS2List(rs:ResultSet):List[List[String]]= {
      import scala.collection.mutable.ListBuffer


      val rsMeta = rs.getMetaData()

      val rscols = rsMeta.getColumnCount()

      var rslist = new ListBuffer[List[String]]()
      while (rs.next){
        val resultset = (1 to rscols).map(n => rs.getString(n)).toList
        rslist += resultset
      }
      if (rs.isAfterLast()) rs.first()


      val rsList = rslist.toList
      //rsList.foreach(println)
      println(rsList.length)
      return rsList
    }
    def RS2File(rs:ResultSet,file:String,loc:Int = -1,map:scala.collection.mutable.Map[String,String] = scala.collection.mutable.Map[String,String]("a"->"a")) = {
      val rsMeta = rs.getMetaData()
      val rscols = rsMeta.getColumnCount()
      val fileOut = new PrintWriter(new FileWriter(file))
      var columnNames = new StringBuilder()
      for (i <- 1 to rscols) {
        columnNames ++= rsMeta.getColumnLabel(i)
        if (i < rscols) {
          columnNames ++= "\t"
        }
        //      else {
        //        columnNames ++= "\n"
        //      }
      }
      fileOut.println(columnNames.toString)
      var count = 1
      while (rs.next) {
        count += 1

        var genepr = new StringBuilder
        for (i <- 1 to rscols) {
          if (i == loc){
            try{
              genepr ++= map(rs.getString(i)).toString
            } catch{
              case e:NoSuchElementException => genepr ++= rs.getString(i)
            }
          } else{
            genepr ++= rs.getString(i)
          }
          if (i < rscols) genepr ++= "\t"
          //if (i == rscols) genepr ++= "\n"
        }
        //if (count <= 20) println(genepr.toString)
        fileOut.println(genepr.toString)

      }
      if (rs.isAfterLast()) rs.first()
      println(count)
      fileOut.close()
    }
    def RS2Map(rs:ResultSet): scala.collection.mutable.Map[String,String] ={
      var map = scala.collection.mutable.Map[String,String]()
      while (rs.next){
        map += (rs.getString(1) -> rs.getString(2))
        // if (rs.isLast()) rs.first()
      }
      if (rs.isAfterLast()) rs.first()
      map
    }
}

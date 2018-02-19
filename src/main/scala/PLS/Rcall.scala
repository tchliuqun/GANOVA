package PLS

object Rcall {
    val R = org.ddahl.rscala.RClient()
    var n = 1000
    val sz = 10
    lazy val Rscript =  s"""
         source("limma.R")
         x = getPermLimma(permN = $n)
         """
    //  def getEngine():Rengine = {
    //   // System.setProperty("
    //    //System.setProperty("java.library.path", "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/rJava/jri/")
    //    new Rengine(Array( "--no-save" ), false, null)
    //  }

    //  val javaVector = "c(1,2,3,4,5)"
    //  engine.eval("rVector=" + javaVector)
    //  engine.eval("meanVal=mean(rVector)")
    //  val mean = engine.eval("meanVal").asDouble()
    def get2DarrayDouble(script:String = Rscript) = {
      R eval script
      R.getD2("x")
    }
    def eval(script:String = Rscript) = {
      R eval script
    }
}

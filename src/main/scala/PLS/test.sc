import PLS.plsCalc
import breeze.linalg._
import breeze.numerics._

val y11 = new DenseMatrix(4,3,Array(1 to 12:_*))
val y2 = 3 * (1 - y11)
val y4 = diag(y11.t * y11)
//DenseVector(y1(::,1).t * y1(::,1))

val y0 = new DenseMatrix(4,3,linspace(1.0,12.0,12).toArray)
val y41 = y0(::,0)
val y42 = y0(::,1)
val y43 = y41 / y42

val y44 = sqrt(y43)
var y45 = DenseVector.zeros[Double](1)
y45(0) = y11(::,0).t * y11(::,0)

val y5 = 3.0 * (1.0 - y0)
val y6 = plsCalc.kfoldInx(100,10,true)
//val fl = scala.io.Source.fromFile("gbmsnp6.csv").getLines
//println(fl.next.mkStrings)
val fl = scala.io.Source.fromFile("gbmsnp6.csv").getLines
val snm = fl.next.split(",")
val rp = scala.io.Source.fromFile("GBMLGG.rppa.txt").getLines
val rpn = rp.next.split("\t").map(_.substring(0,15))
val xy = snm.intersect(rpn).toSet
val x1 = Array.tabulate(rpn.length){ i => (rpn(i),i) }.filter(xy contains _._1).map(_._2)
val y1 = Array.tabulate(snm.length){ i => (snm(i),i) }.filter(xy contains _._1).groupBy(_._1).map{case(k,v) => v(v.length-1)._2}
y1.map(snm)







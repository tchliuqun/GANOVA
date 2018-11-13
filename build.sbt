name := "GANOVA"

version := "0.1"

scalaVersion := "2.12.3"
import com.github.retronym.SbtOneJar._

oneJarSettings
resolvers ++= Seq(
  "Artima Maven Repository" at "http://repo.artima.com/releases",
  "Sonatype OSS Snapshots" at
    "https://oss.sonatype.org/content/repositories/releases"
)
libraryDependencies ++= Seq(

  "commons-httpclient" % "commons-httpclient" % "3.1",
  "commons-io" % "commons-io" % "2.5",
  "commons-net" % "commons-net" % "3.4",
  "com.storm-enroute" %% "scalameter" % "0.8.2",
  "mysql" % "mysql-connector-java" % "5.1.22",
  "org.scalanlp" %% "breeze" % "0.13.2",
  "org.scalanlp" %% "breeze-natives" % "0.13.2",
  "org.ddahl" % "rscala_2.12" % "2.4.0",
  "org.scalactic" %% "scalactic" % "3.0.1",
  "org.scalatest" %% "scalatest" % "3.0.1" % "test",
  "com.typesafe.akka" % "akka-actor_2.12" % "2.5.3",
  "com.typesafe.akka" %% "akka-testkit" % "2.5.3" % Test,
  "com.typesafe.akka" %% "akka-http" % "10.0.9",
  "com.typesafe.akka" %% "akka-http-testkit" % "10.0.9" % Test,
  "org.apache.httpcomponents" % "httpclient" % "4.5.2",
  "org.apache.httpcomponents" % "httpcore" % "4.4.5"
)


dependencyOverrides += "org.scala-lang" % "scala-compiler" % scalaVersion.value

val scalaMeterFramework = new TestFramework("org.scalameter.ScalaMeterFramework")
testFrameworks in ThisBuild += scalaMeterFramework
testOptions in ThisBuild += Tests.Argument(scalaMeterFramework, "-silent")
logBuffered := false
parallelExecution in Test := false

fork in run := true
//javaOptions ++= Seq(
//  "-Xmx6G",
//  "-Xms6G",
//  "-XX:+PrintGCTimeStamps",
//  "-XX:+UnlockCommercialFeatures",
//  "-XX:+FlightRecorder",
//  "-XX:+UnlockDiagnosticVMOptions",
//  "-XX:+DebugNonSafepoints",
//  "-XX:FlightRecorderOptions=defaultrecording=true,dumponexit=true,dumponexitpath=/Users/qunliu/Nutstore/workspace/PLS/resources/liuTest.jfr",
//  "-Xss2048k"
//)
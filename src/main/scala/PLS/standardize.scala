package PLS

import breeze.linalg.DenseMatrix
/**
  * Created by davidkellis on 11/16/16.
  */


object standardize {
  case class StandardizationFactors(
                                     meanOfColumns: DenseMatrix[Float],     // meanOfColumns is a row vector
                                     stdDevOfColumns: DenseMatrix[Float]    // stdDevOfColumns is a row vector
                                   )


  ////////////////////// Mean Centering //////////////////////

  // M is a (R x C) matrix
  // meanOfColumnsInM is a (1 x C) matrix
  def meanCenterColumnsWithMeans(M: DenseMatrix[Float], meanOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = M - utils.vertcat(meanOfColumnsInM, M.rows)

  // The average of each column is subtracted from all the values in corresponding column.
  // As a result of centering, the columns of the centered matrix have a mean of 0. This property does not hold in the case that a column contains a single distinct value, V - in
  // that case, the column is treated as having mean 0, and the mean of the mean-centered column would be V.
  // See Applied Predictive Modeling, page 30-31.
  // Returns a pair of matrices (meanCenteredM, meanOfColumnsInM)
  def meanCenterColumns2(M: DenseMatrix[Float]): (DenseMatrix[Float], DenseMatrix[Float]) = {
    val meanOfColumnsInM = utils.meanColumns(M)

    // We don't want the mean-centered columns to be a column of all zeros, which will will be case when a column in M is filled with a single distinct value.
    // To prevent the mean-centered columns to be a column of all zeros, when we detect that a column contains only one distinct value, then we just
    // consider the mean of that column to be zero. This screws up the property that every mean centered column has a mean of 0, but it allows us to have
    // a column consisting of a single constant value. This is especially important when a predictor matrix, X, has a column of 1s, which will be the case
    // any time we want to learn the equation of a line and we're interested in learning the y-intercept: e.g. y = x*b1 + 1*b0, where b0 is the y intercept.
    val stdDevOfColumnsInM = utils.stddevColumns(M)
    val adjustedColumnMeans = meanOfColumnsInM.mapPairs { (matrixIndex, columnMean) =>
      val (_, col) = matrixIndex
      if (stdDevOfColumnsInM(0, col) == 0.0f) {    // if std dev of column == 0.0, then all values in the column are the same, and we don't want to mean center that column
        0.0f
      } else {
        columnMean
      }
    }
    val meanCenteredM = meanCenterColumnsWithMeans(M, adjustedColumnMeans)
    (meanCenteredM, adjustedColumnMeans)
  }

  // Returns the mean-centered version of M
  def meanCenterColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = meanCenterColumns2(M)._1

  def denormalizeMeanCenteredColumns(meanCenteredM: DenseMatrix[Float], meanOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = {
    meanCenteredM + utils.vertcat(meanOfColumnsInM, meanCenteredM.rows)
  }

  ////////////////////// Scaling //////////////////////

  // M is a (R x C) matrix
  // meanOfColumnsInM is a (1 x C) matrix
  def scaleColumnsWithStdDevs(M: DenseMatrix[Float], stdDevOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = M :/ utils.vertcat(stdDevOfColumnsInM, M.rows)

  // To scale the data, every value within each specific column is divided by the column-specific standard deviation.
  // Each column of the scaled matrix has a standard deviation of 1 or 0. The only time a column will have a std. dev. of 0 is when
  // the column has only one distinct value.
  // See Applied Predictive Modeling, page 30-31.
  // Returns a pair of matrices (scaledM, stdDevOfColumnsInM)
  def scaleColumns2(M: DenseMatrix[Float]): (DenseMatrix[Float], DenseMatrix[Float]) = {
    val stdDevOfColumnsInM = utils.stddevColumns(M)
    // Any column in M that has a standard deviation of 0 is caused by either (1) M having zero rows,
    // or (2) every row in the column is the same value (i.e. the column has one distinct value)
    // If M has no rows, that doesn't cause a problem for us because an empty matrix divided by 0 results in an empty matrix.
    // If every row in the column is the same value, V, then that's a special case, because we can't divide V by zero. To handle it,
    // we treat the column with 0 standard deviation as having a standard deviation of V, so the scaled value, scaledV, is 1. In other words,
    // if the column is full of a single value, V, then we treat the column as having a standard deviation of V, and then to compute the
    // scaled value of V, we compute it with: scaledV = V / V === scaledV = 1.
    val zeroAdjustedStandardDeviations = stdDevOfColumnsInM.mapPairs {(matrixIndex, stdDev) =>
      val (_, col) = matrixIndex
      if (stdDev == 0.0f) {
        if (M(0, col) == 0.0f) {
          1.0f
        } else {
          M(0, col)
        }
      } else {
        stdDev
      }
    }

    val scaledM = scaleColumnsWithStdDevs(M, zeroAdjustedStandardDeviations)
    (scaledM, zeroAdjustedStandardDeviations)
  }

  // Returns the standard-deviation-scaled version of M
  def scaleColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = scaleColumns2(M)._1

  def denormalizeScaledColumns(scaledM: DenseMatrix[Float], stdDevOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = {
    scaledM :* utils.vertcat(stdDevOfColumnsInM, scaledM.rows)
  }


  ////////////////////// Scaling and Mean Centering //////////////////////

  def centerAndScaleColumnsWithFactors(M: DenseMatrix[Float], meanOfColumnsInM: DenseMatrix[Float], stdDevOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = {
    scaleColumnsWithStdDevs(meanCenterColumnsWithMeans(M, meanOfColumnsInM), stdDevOfColumnsInM)
  }

  // mean center and then scale each column of the given matrix, M, returning the pair (scaledMeanCenteredM, StandardizationFactors(meanOfColumnsInM, stdDevOfMeanCenteredColumnsInM))
  def centerAndScaleColumnsReturningFactors(M: DenseMatrix[Float]): (DenseMatrix[Float], StandardizationFactors) = {
    val (meanCenteredM, meanOfColumnsInM) = meanCenterColumns2(M)
    val (scaledMeanCenteredM, stdDevOfMeanCenteredColumnsInM) = scaleColumns2(meanCenteredM)
    (scaledMeanCenteredM, StandardizationFactors(meanOfColumnsInM, stdDevOfMeanCenteredColumnsInM))
  }

  // mean center and then scale each column of the given matrix, M
  def centerAndScaleColumns(M: DenseMatrix[Float]): DenseMatrix[Float] = scaleColumns(meanCenterColumns(M))

  def denormalizeCenteredAndScaledColumns(standardizedM: DenseMatrix[Float], meanOfColumnsInM: DenseMatrix[Float], stdDevOfColumnsInM: DenseMatrix[Float]): DenseMatrix[Float] = {
    denormalizeMeanCenteredColumns(denormalizeScaledColumns(standardizedM, stdDevOfColumnsInM), meanOfColumnsInM)
  }
}

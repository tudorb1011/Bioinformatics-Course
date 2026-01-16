case class Matrix(data: Array[Array[Double]]) {
  val size: Int = data.length

  def multiply(vector: Array[Double]): Array[Double] = {
    require(vector.length == size, "Vector length must match matrix size")
    data.map(row => row.zip(vector).map { case (a, b) => a * b }.sum)
  }

  override def toString: String = {
    data.map(row => row.mkString("[", " ", "]")).mkString("\n")
  }
}

class Predictor(matrix: Matrix, initialVector: Array[Double]) {
  require(matrix.size == initialVector.length, "Matrix and vector dimensions must match")

  def predictSingleStep(currentState: Array[Double]): Array[Double] = {
    matrix.multiply(currentState)
  }

  def predictNSteps(nSteps: Int = 5): List[Array[Double]] = {
    var states = List(initialVector)
    var currentState = initialVector

    for (_ <- 1 to nSteps) {
      val nextState = predictSingleStep(currentState)
      states = states :+ nextState
      currentState = nextState
    }

    states
  }

  def displayPredictions(nSteps: Int = 5): Unit = {
    val states = predictNSteps(nSteps)

    println("\nTransition Matrix:")
    println(matrix)
    println(s"\nInitial State (t=0):")
    println(states.head.mkString("[", " ", "]"))
    println()

    for (i <- 1 until states.length) {
      println(s"t=$i: ${states(i).mkString("[", " ", "]")}")
    }
  }
}

object Main {
  def main(args: Array[String]): Unit = {
    // Example 1: Weather prediction
    println("Example 1: Weather prediction")
    val transitionMatrix1 = Matrix(Array(
      Array(0.7, 0.2, 0.1),
      Array(0.3, 0.4, 0.3),
      Array(0.2, 0.3, 0.5)
    ))
    val initialState1 = Array(1.0, 0.0, 0.0)
    val predictor1 = new Predictor(transitionMatrix1, initialState1)
    predictor1.displayPredictions(5)

    println("\n")

    // Example 2: Population dynamics
    println("Example 2: Population dynamics")
    val leslieMatrix = Matrix(Array(
      Array(0.0, 1.5, 0.8),
      Array(0.6, 0.9, 0.0),
      Array(0.0, 0.3, 0.7)
    ))
    val initialPopulation = Array(100.0, 50.0, 20.0)
    val predictor2 = new Predictor(leslieMatrix, initialPopulation)
    predictor2.displayPredictions(5)

    println("\n")

    // Example 3: Linear system
    println("Example 3: Linear system")
    val simpleMatrix = Matrix(Array(
      Array(0.8, 0.2),
      Array(0.1, 0.9)
    ))
    val simpleVector = Array(10.0, 5.0)
    val predictor3 = new Predictor(simpleMatrix, simpleVector)
    predictor3.displayPredictions(5)
  }
}
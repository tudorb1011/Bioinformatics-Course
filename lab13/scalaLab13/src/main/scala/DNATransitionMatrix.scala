import scala.io.Source
import java.io.PrintWriter
import scala.collection.mutable

case class TransitionMatrix(matrix: Map[Char, Map[Char, Double]], sequence: String) {

  def toJSON: String = {
    val nucleotides = List('A', 'C', 'G', 'T')
    val matrixJson = nucleotides.map { from =>
      val transitions = nucleotides.map { to =>
        val prob = matrix.get(from).flatMap(_.get(to)).getOrElse(0.0)
        s""""$to": $prob"""
      }.mkString(", ")
      s""""$from": {$transitions}"""
    }.mkString(",\n    ")

    s"""{
  "sequence": "$sequence",
  "length": ${sequence.length},
  "transition_matrix": {
    $matrixJson
  }
}"""
  }

  def saveToFile(filename: String): Unit = {
    val writer = new PrintWriter(filename)
    try {
      writer.write(toJSON)
      println(s"Transition matrix saved to $filename")
    } finally {
      writer.close()
    }
  }
}

object DNATransitionMatrix {

  def calculateTransitionMatrix(sequence: String): TransitionMatrix = {
    val nucleotides = List('A', 'C', 'G', 'T')
    val counts = mutable.Map[Char, mutable.Map[Char, Int]]()

    // Initialize counts
    for (from <- nucleotides) {
      counts(from) = mutable.Map[Char, Int]()
      for (to <- nucleotides) {
        counts(from)(to) = 0
      }
    }

    // Count transitions
    for (i <- 0 until sequence.length - 1) {
      val from = sequence(i)
      val to = sequence(i + 1)
      if (nucleotides.contains(from) && nucleotides.contains(to)) {
        counts(from)(to) += 1
      }
    }

    // Calculate probabilities
    val matrix = counts.map { case (from, transitions) =>
      val total = transitions.values.sum.toDouble
      val probs = if (total > 0) {
        transitions.map { case (to, count) => (to, count / total) }.toMap
      } else {
        transitions.map { case (to, _) => (to, 0.0) }.toMap
      }
      (from, probs)
    }.toMap

    TransitionMatrix(matrix, sequence)
  }

  def generateRandomDNASequence(length: Int): String = {
    val nucleotides = Array('A', 'C', 'G', 'T')
    val random = new scala.util.Random()
    (1 to length).map(_ => nucleotides(random.nextInt(4))).mkString
  }

  def displayMatrix(tm: TransitionMatrix): Unit = {
    println(s"\nDNA Sequence: ${tm.sequence}")
    println(s"Length: ${tm.sequence.length}")
    println("\nTransition Matrix:")
    println("From -> To:  A      C      G      T")
    println("----------------------------------------")

    val nucleotides = List('A', 'C', 'G', 'T')
    for (from <- nucleotides) {
      val probs = nucleotides.map { to =>
        val prob = tm.matrix.get(from).flatMap(_.get(to)).getOrElse(0.0)
        f"$prob%.3f"
      }.mkString("  ")
      println(s"$from        $probs")
    }
  }

  def main(args: Array[String]): Unit = {
    // Generate a random DNA sequence of 50 letters
    val dnaSequence = generateRandomDNASequence(50)

    // Calculate transition matrix
    val transitionMatrix = calculateTransitionMatrix(dnaSequence)

    // Display the matrix
    displayMatrix(transitionMatrix)

    // Save to JSON file
    transitionMatrix.saveToFile("dna_transition_matrix.json")

    println("\nJSON output:")
    println(transitionMatrix.toJSON)
  }
}
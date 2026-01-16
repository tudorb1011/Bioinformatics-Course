import scala.io.Source
import scala.util.Random
import java.io.File

object DNASequenceGenerator {

  case class TransitionMatrix(
                               matrix: Map[String, Map[String, Double]],
                               sequence: String
                             )

  def parseJSON(filename: String): TransitionMatrix = {
    val source = Source.fromFile(filename)
    val jsonString = try source.mkString finally source.close()

    // Simple JSON parser for our specific format
    val sequenceLine = jsonString.split("\n").find(_.contains("\"sequence\"")).get
    val sequence = sequenceLine.split(":")(1).trim.stripPrefix("\"").stripSuffix("\",")

    val matrixLines = jsonString.split("\n")
      .dropWhile(!_.contains("\"transition_matrix\""))
      .drop(1)
      .takeWhile(!_.contains("}"))
      .filter(_.contains("\":"))

    val nucleotides = List("A", "C", "G", "T")
    val matrix = nucleotides.map { from =>
      val line = matrixLines.find(_.contains(s"\"$from\"")).getOrElse("")
      val transitions = nucleotides.map { to =>
        val pattern = s""""$to":\\s*([0-9.]+)""".r
        val prob = pattern.findFirstMatchIn(line).map(_.group(1).toDouble).getOrElse(0.0)
        (to, prob)
      }.toMap
      (from, transitions)
    }.toMap

    TransitionMatrix(matrix, sequence)
  }

  def generateSequence(matrix: TransitionMatrix, length: Int, startNucleotide: String = "A"): String = {
    val random = new Random()
    val sequence = new StringBuilder(startNucleotide)
    var current = startNucleotide

    for (_ <- 1 until length) {
      val transitions = matrix.matrix.getOrElse(current, Map.empty)

      if (transitions.isEmpty || transitions.values.sum == 0) {
        // If no transitions, pick random
        val nucleotides = Array("A", "C", "G", "T")
        current = nucleotides(random.nextInt(4))
      } else {
        // Use weighted random selection based on probabilities
        val rand = random.nextDouble()
        var cumulative = 0.0
        var selected = current

        for ((nucleotide, prob) <- transitions) {
          cumulative += prob
          if (rand <= cumulative && selected == current) {
            selected = nucleotide
          }
        }
        current = selected
      }
      sequence.append(current)
    }

    sequence.toString()
  }

  def main(args: Array[String]): Unit = {
    val jsonFile = "dna_transition_matrix.json"

    if (!new File(jsonFile).exists()) {
      println(s"Error: $jsonFile not found!")
      println("Please run DNATransitionMatrix.scala first to generate the transition matrix.")
      sys.exit(1)
    }

    println("Loading transition matrix from JSON...")
    val matrix = parseJSON(jsonFile)

    println(s"\nOriginal sequence (first 50 chars):")
    println(matrix.sequence.take(50))

    println("\nGenerating new DNA sequences using transition probabilities...")
    println("\nGenerated Sequence 1 (length 100):")
    val seq1 = generateSequence(matrix, 100)
    println(seq1)

    println("\nGenerated Sequence 2 (length 100):")
    val seq2 = generateSequence(matrix, 100, "G")
    println(seq2)

    println("\nGenerated Sequence 3 (length 200):")
    val seq3 = generateSequence(matrix, 200, "C")
    println(seq3)

    // Statistics
    println("\n--- Statistics ---")
    def countNucleotides(seq: String): Map[Char, Int] = {
      seq.groupBy(identity).mapValues(_.length).toMap
    }

    println("\nOriginal sequence composition:")
    countNucleotides(matrix.sequence).foreach { case (n, count) =>
      println(f"  $n: $count (${count * 100.0 / matrix.sequence.length}%.1f%%)")
    }

    println("\nGenerated sequence 1 composition:")
    countNucleotides(seq1).foreach { case (n, count) =>
      println(f"  $n: $count (${count * 100.0 / seq1.length}%.1f%%)")
    }
  }
}
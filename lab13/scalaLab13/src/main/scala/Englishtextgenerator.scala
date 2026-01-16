import scala.io.Source
import scala.util.Random
import java.io.File

object EnglishTextGenerator {

  case class WordTransitionMatrix(
                                   wordToSymbol: Map[String, String],
                                   symbolToWord: Map[String, String],
                                   transitionMatrix: Map[String, Map[String, Double]]
                                 )

  def parseJSON(filename: String): WordTransitionMatrix = {
    val source = Source.fromFile(filename)
    val jsonString = try source.mkString finally source.close()

    // Parse word mapping
    val wordMappingSection = jsonString
      .split("\"word_mapping\":")(1)
      .split("\"transition_matrix\":")(0)

    val symbolPattern = """"(W\d+)":\s*"([^"]+)"""".r
    val symbolToWord = symbolPattern.findAllMatchIn(wordMappingSection).map { m =>
      (m.group(1), m.group(2))
    }.toMap

    val wordToSymbol = symbolToWord.map(_.swap)

    // Parse transition matrix
    val matrixSection = jsonString
      .split("\"transition_matrix\":")(1)
      .split("},\n\\s*\"total_words\"")(0)

    val transitionPattern = """"(W\d+)":\s*\{([^}]+)\}""".r
    val transitionMatrix = transitionPattern.findAllMatchIn(matrixSection).map { m =>
      val fromSymbol = m.group(1)
      val transitionsStr = m.group(2)

      val probPattern = """"(W\d+)":\s*([0-9.]+)""".r
      val transitions = probPattern.findAllMatchIn(transitionsStr).map { pm =>
        (pm.group(1), pm.group(2).toDouble)
      }.toMap

      (fromSymbol, transitions)
    }.toMap

    WordTransitionMatrix(wordToSymbol, symbolToWord, transitionMatrix)
  }

  def generateText(matrix: WordTransitionMatrix, wordCount: Int, startWord: Option[String] = None): String = {
    val random = new Random()

    val firstWord = startWord.getOrElse {
      val allWords = matrix.symbolToWord.values.toSeq
      allWords(random.nextInt(allWords.size))
    }

    val words = new StringBuilder(firstWord)
    var currentSymbol = matrix.wordToSymbol.getOrElse(firstWord, {
      val allSymbols = matrix.symbolToWord.keys.toSeq
      allSymbols(random.nextInt(allSymbols.size))
    })

    for (_ <- 1 until wordCount) {
      val transitions = matrix.transitionMatrix.getOrElse(currentSymbol, Map.empty)

      if (transitions.isEmpty || transitions.values.sum == 0) {
        // Pick random word
        val allSymbols = matrix.symbolToWord.keys.toSeq
        currentSymbol = allSymbols(random.nextInt(allSymbols.size))
      } else {
        // Weighted random selection
        val rand = random.nextDouble()
        var cumulative = 0.0
        var selected = currentSymbol

        for ((symbol, prob) <- transitions) {
          cumulative += prob
          if (rand <= cumulative && selected == currentSymbol) {
            selected = symbol
          }
        }
        currentSymbol = selected
      }

      val nextWord = matrix.symbolToWord(currentSymbol)
      words.append(" ").append(nextWord)
    }

    words.toString()
  }

  def formatText(text: String, lineLength: Int = 80): String = {
    val words = text.split(" ")
    val lines = new StringBuilder()
    var currentLine = new StringBuilder()

    for (word <- words) {
      if (currentLine.length + word.length + 1 > lineLength && currentLine.nonEmpty) {
        lines.append(currentLine.toString().trim).append("\n")
        currentLine = new StringBuilder()
      }
      currentLine.append(word).append(" ")
    }

    if (currentLine.nonEmpty) {
      lines.append(currentLine.toString().trim)
    }

    // Capitalize first letter of each sentence
    val result = lines.toString()
    result.charAt(0).toUpper + result.substring(1)
  }

  def main(args: Array[String]): Unit = {
    val jsonFile = "word_transition_matrix.json"

    if (!new File(jsonFile).exists()) {
      println(s"Error: $jsonFile not found!")
      println("Please run WordTransitionMatrix.scala first to generate the transition matrix.")
      sys.exit(1)
    }

    println("Loading word transition matrix from JSON...")
    val matrix = parseJSON(jsonFile)

    println(s"Loaded ${matrix.symbolToWord.size} unique words")
    println(s"Transition matrix has ${matrix.transitionMatrix.size} entries")

    println("\n" + "=" * 80)
    println("Generated Text 1 (50 words, starting with 'the'):")
    println("=" * 80)
    val text1 = generateText(matrix, 50, Some("the"))
    println(formatText(text1))

    println("\n" + "=" * 80)
    println("Generated Text 2 (100 words, random start):")
    println("=" * 80)
    val text2 = generateText(matrix, 100)
    println(formatText(text2))

    println("\n" + "=" * 80)
    println("Generated Text 3 (150 words, starting with 'science'):")
    println("=" * 80)
    val text3 = generateText(matrix, 150, Some("science"))
    println(formatText(text3))

    println("\n" + "=" * 80)
    println("Generated Text 4 (200 words, starting with 'artificial'):")
    println("=" * 80)
    val text4 = generateText(matrix, 200, Some("artificial"))
    println(formatText(text4))

    // Word frequency analysis
    println("\n" + "=" * 80)
    println("Word Frequency in Generated Text 2:")
    println("=" * 80)
    val wordFreq = text2.toLowerCase().split("\\s+").groupBy(identity).mapValues(_.length)
    wordFreq.toSeq.sortBy(-_._2).take(10).foreach { case (word, count) =>
      println(f"  $word%-20s $count")
    }
  }
}
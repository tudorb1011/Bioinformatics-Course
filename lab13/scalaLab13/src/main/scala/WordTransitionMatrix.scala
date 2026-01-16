import scala.io.Source
import java.io.PrintWriter
import scala.collection.mutable

case class WordTransitionMatrix(
                                 wordToSymbol: Map[String, String],
                                 symbolToWord: Map[String, String],
                                 transitionMatrix: Map[String, Map[String, Double]],
                                 text: String
                               ) {

  def toJSON: String = {
    val wordMappingJson = wordToSymbol.map { case (word, symbol) =>
      s""""$symbol": "$word""""
    }.mkString(",\n    ")

    val matrixJson = transitionMatrix.map { case (fromSymbol, transitions) =>
      val transitionsJson = transitions.map { case (toSymbol, prob) =>
        s""""$toSymbol": $prob"""
      }.mkString(", ")
      s""""$fromSymbol": {$transitionsJson}"""
    }.mkString(",\n    ")

    s"""{
  "word_mapping": {
    $wordMappingJson
  },
  "transition_matrix": {
    $matrixJson
  },
  "total_words": ${wordToSymbol.size},
  "text_length": ${text.split("\\s+").length}
}"""
  }

  def saveToFile(filename: String): Unit = {
    val writer = new PrintWriter(filename)
    try {
      writer.write(toJSON)
      println(s"\nTransition matrix saved to $filename")
    } finally {
      writer.close()
    }
  }
}

object WordTransitionMatrix {

  def generateEnglishText(lines: Int): String = {
    val sentences = Array(
      "The quick brown fox jumps over the lazy dog.",
      "Artificial intelligence is transforming the world rapidly.",
      "Machine learning algorithms require large amounts of data.",
      "Natural language processing enables computers to understand text.",
      "Deep learning models have revolutionized computer vision.",
      "Science and technology advance human civilization forward.",
      "The universe contains billions of galaxies and stars.",
      "Education is the foundation of personal growth.",
      "Innovation drives progress in modern society today.",
      "Climate change affects ecosystems around the globe.",
      "Sustainable development ensures future generations prosper.",
      "Digital transformation reshapes business operations significantly.",
      "Cybersecurity protects sensitive information from threats.",
      "Blockchain technology enables secure decentralized transactions.",
      "Quantum computing promises unprecedented computational power.",
      "Renewable energy sources reduce carbon emissions effectively.",
      "Biotechnology advances medical treatments and diagnostics.",
      "Space exploration expands our understanding of cosmos.",
      "Cultural diversity enriches human experience profoundly.",
      "Communication technology connects people across continents."
    )

    val random = new scala.util.Random()
    val lines_text = (1 to lines).map { _ =>
      sentences(random.nextInt(sentences.length))
    }.mkString("\n")

    lines_text
  }

  def calculateTransitionMatrix(text: String): WordTransitionMatrix = {
    val words = text.toLowerCase()
      .replaceAll("[^a-z\\s]", "")
      .split("\\s+")
      .filter(_.nonEmpty)
      .toList

    val uniqueWords = words.distinct
    val wordToSymbol = uniqueWords.zipWithIndex.map { case (word, idx) =>
      (word, s"W$idx")
    }.toMap

    val symbolToWord = wordToSymbol.map(_.swap)

    val counts = mutable.Map[String, mutable.Map[String, Int]]()

    for (symbol <- wordToSymbol.values) {
      counts(symbol) = mutable.Map[String, Int]()
    }

    for (i <- 0 until words.length - 1) {
      val fromWord = words(i)
      val toWord = words(i + 1)
      val fromSymbol = wordToSymbol(fromWord)
      val toSymbol = wordToSymbol(toWord)

      counts(fromSymbol)(toSymbol) = counts(fromSymbol).getOrElse(toSymbol, 0) + 1
    }

    val matrix = counts.map { case (fromSymbol, transitions) =>
      val total = transitions.values.sum.toDouble
      val probs = if (total > 0) {
        transitions.map { case (toSymbol, count) => (toSymbol, count / total) }.toMap
      } else {
        Map[String, Double]()
      }
      (fromSymbol, probs)
    }.toMap

    WordTransitionMatrix(wordToSymbol, symbolToWord, matrix, text)
  }

  def displayMatrix(wtm: WordTransitionMatrix): Unit = {
    println(s"\nTotal unique words: ${wtm.wordToSymbol.size}")
    println(s"Total words in text: ${wtm.text.split("\\s+").length}")

    println("\nWord to Symbol Mapping (first 20):")
    wtm.wordToSymbol.take(20).foreach { case (word, symbol) =>
      println(s"  $symbol -> $word")
    }

    println("\nSample Transition Probabilities (first 10 words):")
    wtm.transitionMatrix.take(10).foreach { case (fromSymbol, transitions) =>
      val fromWord = wtm.symbolToWord(fromSymbol)
      println(s"\nFrom '$fromWord' ($fromSymbol):")
      transitions.toSeq.sortBy(-_._2).take(5).foreach { case (toSymbol, prob) =>
        val toWord = wtm.symbolToWord(toSymbol)
        println(f"  -> '$toWord' ($toSymbol): $prob%.3f")
      }
    }
  }

  def main(args: Array[String]): Unit = {
    println("Generating 300 lines of English text...")
    val text = generateEnglishText(300)

    println("Calculating word transition probabilities...")
    val transitionMatrix = calculateTransitionMatrix(text)

    displayMatrix(transitionMatrix)

    transitionMatrix.saveToFile("word_transition_matrix.json")

    println("\nDone! Check word_transition_matrix.json for full results.")
  }
}
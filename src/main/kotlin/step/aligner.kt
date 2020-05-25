package step

import mu.KotlinLogging
import util.*
import java.nio.file.*
import java.io.*
import util.CmdRunner
import java.util.zip.*

private val log = KotlinLogging.logger {}

data class AlignmentParameters (
    val r1: Path,
    val index: Path,
    val outputDirectory: Path,
    val r2: Path? = null,
    val libraryId: String? = null,
    val cores: Int = 1,
    val ram: Int = 16,
    val outputPrefix: String = "output",
    val indexTarPrefix: String? = null,
    val minReadLength: Int = 20
)

fun CmdRunner.isRSEMSorted(bam: Path): Boolean?
    = this.runCommand("rsem-sam-validator ${bam}")?.trim()?.split("\t")?.get(0)?.endsWith("is valid!")

fun CmdRunner.getFlagstats(inputPath: Path, outputPath: Path)
    = this.run("samtools flagstat $inputPath >  $outputPath")

fun CmdRunner.trimFASTQ(fastq: Path, outputDirectory: Path, minLength: Int): File {
    val trimmedFASTQ = createTempFile(suffix = ".gz", directory = outputDirectory.toFile())
    var index: Int = 0
    var currentLine: String = ""
    var shouldWrite: Boolean = true
    val inputStream = if (fastq.fileName.toString().endsWith(".gz")) (
        GZIPInputStream(fastq.toFile().inputStream())
    ) else fastq.toFile().inputStream()
    inputStream.bufferedReader().use { input ->
        GZIPOutputStream(trimmedFASTQ.outputStream()).bufferedWriter().use { output ->
            input.forEachLine {
                currentLine += it + "\n"
                if (index % 4 == 1)
                    shouldWrite = it.length >= minLength
                else if (index % 4 == 3) {
                    if (shouldWrite) output.write(currentLine)
                    currentLine = ""
                }
                ++index
            }
        }
    }
    return trimmedFASTQ
}

fun CmdRunner.align(parameters: AlignmentParameters) {

    // create output directory, unpack index
    val indexDir = Files.createDirectories(parameters.outputDirectory.resolve("index"))
    this.run("tar xvf${if (parameters.index.endsWith("gz")) "z" else ""} ${parameters.index} -C ${indexDir}")
    if (parameters.indexTarPrefix !== null) this.run("mv ${indexDir}/${parameters.indexTarPrefix}/* ${indexDir}")

    // drop reads which are too short
    val trimmedR1 = this.trimFASTQ(parameters.r1, parameters.outputDirectory, parameters.minReadLength)
    val trimmedR2: File? = if (parameters.r2 !== null) this.trimFASTQ(parameters.r2, parameters.outputDirectory, parameters.minReadLength) else null

    // run STAR
    this.run("""
        STAR \
            --genomeDir ${indexDir} \
            --readFilesIn ${trimmedR1} ${if (trimmedR2 !== null) trimmedR2.toPath() else ""} \
            --readFilesCommand zcat \
            --runThreadN ${parameters.cores} \
            --genomeLoad NoSharedMemory \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outSAMheaderCommentFile COfile.txt \
            --outSAMheaderHD @HD VN:1.4 SO:coordinate \
            --outSAMunmapped Within \
            --outFilterType BySJout \
            --outSAMattributes NH HI AS NM MD \
            ${ if (parameters.r2 !== null) "--outSAMstrandField intronMotif" else "" } \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM \
            --sjdbScore 1 \
            --limitBAMsortRAM ${parameters.ram}000000000 \
            --outFileNamePrefix ${parameters.outputDirectory}/star
    """)
    indexDir.toFile().deleteRecursively()

    // rename outputs
    Files.move(
        parameters.outputDirectory.resolve("starAligned.sortedByCoord.out.bam"),
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_genome.bam")
    )
    Files.move(
        parameters.outputDirectory.resolve("starLog.final.out"),
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_Log.final.out")
    )

    // sort for RSEM if necessary
    if (this.isRSEMSorted(parameters.outputDirectory.resolve("starAligned.toTranscriptome.out.bam")) == true)
        Files.move(
            parameters.outputDirectory.resolve("starAligned.toTranscriptome.out.bam"),
            parameters.outputDirectory.resolve("${parameters.outputPrefix}_anno.bam")
        )
    else
        this.run("""
            convert-sam-for-rsem \
                ${parameters.outputDirectory.resolve("starAligned.toTranscriptome.out.bam")} \
                ${parameters.outputDirectory.resolve("${parameters.outputPrefix}_anno")}
        """)

    // perform flagstat
    getFlagstats(
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_genome.bam"),
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_genome_flagstat.txt")
    )
    getFlagstats(
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_anno.bam"),
        parameters.outputDirectory.resolve("${parameters.outputPrefix}_anno_flagstat.txt")
    )
    
}

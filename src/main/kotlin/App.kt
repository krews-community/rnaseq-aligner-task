import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.*
import step.*
import util.*
import java.nio.file.*
import util.CmdRunner


fun main(args: Array<String>) = Cli().main(args)

class Cli : CliktCommand() {
    private val repFile1: Path by option("-repFile1", help = "path to fastq rep file")
            .path(exists = true).required()
    private val repFile2: Path? by option("-repFile2", help = "path to fastq rep2 file")
            .path()

    private val indexFile:Path by option("-indexFile", help = "path to index tar file")
            .path(exists = true).required()

    private val pairedEnd: Boolean by option("-pairedEnd", help = "Paired-end BAM.").flag()
    private val aligner:String? by option("-aligner",help = "Aligner Type").choice("star").default("star")
    private val libraryid:String? by option("-libraryid",help = "library identifier which will be added to bam header").default("libraryID")
    private val outputPrefix: String by option("-outputPrefix", help = "output file name prefix; defaults to 'output'").default("output")
    private val outDir by option("-outputDir", help = "path to output Directory")
        .path().required()
    private val ncpus: Int by option("-ncpus", help = "Number of cpus available.").int().default(4)
    private val ramGB: Int by option("-ramGB", help = "Amount of RAM available in GB").int().default(8)

    override fun run() {
        val cmdRunner = DefaultCmdRunner()
        cmdRunner.runTask(repFile1,repFile2,indexFile, pairedEnd,aligner,libraryid,ncpus,ramGB,outDir,outputPrefix)
    }
}

/**
 * Runs pre-processing and bwa for raw input files
 *
 * @param taFiles pooledTa Input
 * @param outDir Output Path
 */
fun CmdRunner.runTask(repFile1:Path,repFile2:Path?,indexFile:Path,pairedEnd:Boolean,aligner:String?,libraryid:String?,ncpus:Int,ramGB:Int, outDir:Path,outputPrefix:String) {

    aligner(repFile1,repFile2,indexFile, pairedEnd,aligner,libraryid,ncpus,ramGB,outDir,outputPrefix)
}
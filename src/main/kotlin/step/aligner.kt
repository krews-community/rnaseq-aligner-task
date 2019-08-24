package step
import mu.KotlinLogging
import util.*
import java.nio.file.*
import util.CmdRunner
private val log = KotlinLogging.logger {}

fun CmdRunner.get_flagstats(input_path:String, output_path:String) {
    log.info {""}
    val command = "samtools flagstat ${input_path} >  ${output_path}"
    this.run(command)

}
fun CmdRunner.aligner(repFile1:Path,repFile2:Path?,indexFile:Path,pairedEnd:Boolean,aligner:String?,libraryid:String?,ncpus:Int,ramGB:Int, outDir:Path,outputPrefix:String) {
    log.info { "Make output Diretory" }
    Files.createDirectories(outDir)
   this.run("tar xvf $indexFile -C $outDir")

    var cmd=""
    if(pairedEnd)
    {
        cmd="STAR --genomeDir ${outDir}/out \\\n" +
                "    --readFilesIn ${repFile1} ${repFile2} \\\n" +
                "    --readFilesCommand zcat \\\n" +
                "    --runThreadN ${ncpus} \\\n" +
                "    --genomeLoad NoSharedMemory \\\n" +
                "    --outFilterMultimapNmax 20 \\\n" +
                "    --alignSJoverhangMin 8 \\\n" +
                "    --alignSJDBoverhangMin 1 \\\n" +
                "    --outFilterMismatchNmax 999 \\\n" +
                "    --outFilterMismatchNoverReadLmax 0.04 \\\n" +
                "    --alignIntronMin 20 \\\n" +
                "    --alignIntronMax 1000000 \\\n" +
                "    --alignMatesGapMax 1000000 \\\n" +
                "    --outSAMheaderCommentFile COfile.txt \\\n" +
                "    --outSAMheaderHD @HD VN:1.4 SO:coordinate \\\n" +
                "    --outSAMunmapped Within \\\n" +
                "    --outFilterType BySJout \\\n" +
                "    --outSAMattributes NH HI AS NM MD \\\n" +
                "    --outSAMtype BAM SortedByCoordinate \\\n" +
                "    --quantMode TranscriptomeSAM \\\n" +
                "    --sjdbScore 1 \\\n" +
                "    --limitBAMsortRAM ${ramGB}000000000 \\\n" +
                " --outFileNamePrefix ${outDir}/star"
    }else {
        cmd="STAR --genomeDir ${outDir}/out \\\n" +
                "    --readFilesIn ${repFile1} \\\n" +
                "    --readFilesCommand zcat \\\n" +
                "    --runThreadN ${ncpus} \\\n" +
                "    --genomeLoad NoSharedMemory \\\n" +
                "    --outFilterMultimapNmax 20 \\\n" +
                "    --alignSJoverhangMin 8 \\\n" +
                "    --alignSJDBoverhangMin 1 \\\n" +
                "    --outFilterMismatchNmax 999 \\\n" +
                "    --outFilterMismatchNoverReadLmax 0.04 \\\n" +
                "    --alignIntronMin 20 \\\n" +
                "    --alignIntronMax 1000000 \\\n" +
                "    --alignMatesGapMax 1000000 \\\n" +
                "    --outSAMheaderCommentFile COfile.txt \\\n" +
                "    --outSAMheaderHD @HD VN:1.4 SO:coordinate \\\n" +
                "    --outSAMunmapped Within \\\n" +
                "    --outFilterType BySJout \\\n" +
                "    --outSAMattributes NH HI AS NM MD \\\n" +
                "    --outSAMstrandField intronMotif \\\n" +
                "    --outSAMtype BAM SortedByCoordinate  \\\n" +
                "    --quantMode TranscriptomeSAM \\\n" +
                "    --sjdbScore 1 \\\n" +
                "    --limitBAMsortRAM ${ramGB}000000000  \\\n" +
                "    --outFileNamePrefix ${outDir}/star"
    }
    this.run(cmd)

    //rename
   val f1 = outDir.resolve("starAligned.sortedByCoord.out.bam")
    val f2 = outDir.resolve("${outputPrefix}_genome.bam")
    this.run("mv ${f1} ${f2}")

    //rename
    val f3 = outDir.resolve("starLog.final.out")
    val f4 = outDir.resolve("${outputPrefix}_Log.final.out")
    this.run("mv ${f3} ${f4}")

    val rsem_check_cmd = "rsem-sam-validator ${outDir.resolve("starAligned.toTranscriptome.out.bam")}"

    val o = this.runCommand(rsem_check_cmd)
    val rsem_valid = o!!.trim().split("\t")[0].endsWith("is valid!")
    if(rsem_valid)
    {
        log.info { "Transcriptome bam is already rsem-sorted." }
        val f1 = outDir.resolve("starAligned.toTranscriptome.out.bam")
        val f2 = outDir.resolve("${outputPrefix}_anno.bam")
        this.run("mv ${f1} ${f2}")
    }else {
        log.info {"Transcriptome bam is not rsem-sorted."}
        val rsem_sort_cmd = "convert-sam-for-rsem ${outDir.resolve("starAligned.toTranscriptome.out.bam")} ${outDir.resolve("${outputPrefix}_anno")}"
        this.run(rsem_sort_cmd)
    }
    val genome_bam_path = outDir.resolve("${outputPrefix}_genome.bam")
    val anno_bam_path = outDir.resolve("${outputPrefix}_anno.bam")
    val genome_flagstat_path =  outDir.resolve("${outputPrefix}_genome_flagstat.txt")
    val anno_flagstat_path =  outDir.resolve("${outputPrefix}_anno_flagstat.txt")
    get_flagstats(genome_bam_path.toString(), genome_flagstat_path.toString())
    get_flagstats(anno_bam_path.toString(), anno_flagstat_path.toString())
}

package testutil

import java.security.MessageDigest
import java.nio.file.*
import util.CmdRunner
import util.*

val cmdRunner = TestCmdRunner()

class TestCmdRunner : CmdRunner {
    override fun run(cmd: String) = exec("docker", "exec", "rnaseq-align-base", "sh", "-c", cmd)
    override fun runCommand(cmd: String):String? = getCommandOutput("docker", "exec", "rnaseq-align-base", "sh", "-c", cmd)
}

fun copyDirectory(fromDir: Path, toDir: Path) {
    Files.walk(fromDir).forEach { fromFile ->
        if (Files.isRegularFile(fromFile)) {
            val toFile = toDir.resolve(fromDir.relativize(fromFile))
            Files.createDirectories(toFile.parent)
            Files.copy(fromFile, toFile)
        }
    }
}

fun setupTest() {
    exec(
        "docker", "run", "--name", "rnaseq-align-base", "--rm", "-i",
        "-t", "-d", "-v", "${testInputResourcesDir}:${testInputResourcesDir}",
        "-v", "${testDir}:${testDir}", "genomealmanac/rnaseq-align-base"
    )
}

fun cleanupTest() {
    testDir.toFile().deleteRecursively()
}

fun String.toMD5(): String {
    val bytes = MessageDigest.getInstance("MD5").digest(this.toByteArray())
    return bytes.toHex()    
}

fun ByteArray.toHex(): String {
    return joinToString("") { "%02x".format(it) }   
}

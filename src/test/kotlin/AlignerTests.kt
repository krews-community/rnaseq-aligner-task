import org.junit.jupiter.api.*
import step.*
import testutil.*
import java.nio.file.*
import org.assertj.core.api.Assertions.*

class AlignerTests {

    @BeforeEach fun setup() = setupTest()
    @AfterEach fun cleanup() = cleanupTest()

    @Test fun `test STAR`() {
        cmdRunner.align(
            AlignmentParameters(
                r1 = getResourcePath("test.fastq.gz"),
                index = getResourcePath("index.tar.gz"),
                outputDirectory = testDir
            )
        )

        assertThat(testDir.resolve("output_anno.bam")).exists()
        assertThat(testDir.resolve("output_anno_flagstat.txt")).exists()
        testDir.resolve("output_anno_flagstat.txt").toFile().bufferedReader().use {
            assertThat(it.readText().toMD5()).isEqualTo("87c78efa9ad0a76dec66ab5c8e8e7754")
        }

        assertThat(testDir.resolve("output_genome.bam")).exists()
        assertThat(testDir.resolve("output_genome_flagstat.txt")).exists()
        testDir.resolve("output_genome_flagstat.txt").toFile().bufferedReader().use {
            assertThat(it.readText().toMD5()).isEqualTo("a96d66f579ad1d88239764e059fc2554")
        }
        
    }

}

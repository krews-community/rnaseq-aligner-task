import org.junit.jupiter.api.*
import step.*
import testutil.*
import testutil.cmdRunner
import testutil.setupTest
import java.nio.file.*
import org.assertj.core.api.Assertions.*

class AlignerTests {
   @BeforeEach fun setup() = setupTest()
   @AfterEach fun cleanup() = cleanupTest()

    //TODO: Add Tests
     @Test fun `run pooledTa step singal file `() {

       //  cmdRunner.aligner(F1,null, CTL_INDEX, false,"star","libraryID",4,15, testOutputDir,"testaligner")
         //assertThat(testOutputDir.resolve("mergedta.pooled.tagAlign.gz")).exists()

    }

}
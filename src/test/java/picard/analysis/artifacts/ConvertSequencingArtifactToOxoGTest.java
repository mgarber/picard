package picard.analysis.artifacts;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectOxoGMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by farjoun on 10/29/17.
 */
public class ConvertSequencingArtifactToOxoGTest  {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
    private static final File REFERENCE_SEQUENCE = new File(TEST_DATA_DIR, "merger.fasta");

    @Test
    public void testEquivalence() throws IOException {
        final File input = SAM_FILE;

        final File outputFileOxoG = File.createTempFile("test", ".oxo_g_metrics", TEST_DATA_DIR);
        outputFileOxoG.deleteOnExit();

        final File outputFileArtifacts = File.createTempFile("test", "", TEST_DATA_DIR);
        outputFileArtifacts.deleteOnExit();

        final File convertedArtifacts = File.createTempFile("test", ".oxog_metrics", TEST_DATA_DIR);
        convertedArtifacts.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("OUTPUT=" + outputFileOxoG.getAbsolutePath());
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath());

        CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        Assert.assertEquals(collectOxoGMetrics.instanceMain(args.toArray(new String[args.size()])), 0,
                "CollectOxoGMetrics can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        args.clear();
        args.add( "OUTPUT=" + outputFileArtifacts.getAbsolutePath());
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath());

        CollectSequencingArtifactMetrics collectArtifactMetrics = new CollectSequencingArtifactMetrics();
        Assert.assertEquals(collectArtifactMetrics.instanceMain(args.toArray(new String[args.size()])), 0,
                "CollectSequencingArtifactMetrics can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        args.clear();
        args.add("OUTPUT_BASE=" + convertedArtifacts.getAbsolutePath());
        args.add("INPUT_BASE=" + outputFileArtifacts.getAbsolutePath());


        ConvertSequencingArtifactToOxoG convertSequencingArtifactToOxoG = new ConvertSequencingArtifactToOxoG();
        Assert.assertEquals(convertSequencingArtifactToOxoG.instanceMain(args.toArray(new String[args.size()])), 0,
                "ConvertSequencingArtifactToOxoG can't process base input" + outputFileArtifacts.getAbsolutePath() + " correctly");

       Assert.assertTrue(MetricsFile.areMetricsEqual(convertedArtifacts,outputFileOxoG), "Metrics differ");

    }


}
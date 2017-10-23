package picard.fingerprint;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.vcf.SamTestUtils;
import picard.vcf.VcfTestUtils;
import sun.nio.ch.IOUtil;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created by farjoun on 10/23/17.
 */
public class CheckFingerprintTest extends CommandLineProgramTest{

    @Override
    public String getCommandLineProgramName() {
        return this.getClass().getName();
    }

    private final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private final File HAPLOTYPE_MAP = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    private final File NA12891_r1_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
    private final File NA12891_r2_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");

    //this is a copy of a previous one, but with a different sample name
    private final File NA12891_named_NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12891_named_NA12892.over.fingerprints.r1.sam");

    private final File NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
    private final File NA12892_r2_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam");

    private File NA12891_r1, NA12891_r2, NA12891_named_NA12892_r1, NA12892_r1, NA12892_r2;

    private final int NA12891_r1_RGs = 27;
    private final int NA12891_r2_RGs = 26;
    private final int NA12892_r1_RGs = 25;
    private final int NA12892_r2_RGs = 26;

    private static File  NA12891_1_vcf;
    private static File  NA12891_2_vcf;
    private static File  NA12892_1_vcf;
    private static File  NA12891_swapped_nonref_g_vcf;
    private static File  NA12892_2_vcf;
    private static File  NA12891_named_NA12892_vcf;
    private static File  NA12891_g_vcf;
    private static File  NA12892_g_vcf;
    private static File  NA12892_and_NA123891_vcf;
    private static File  NA12892_and_NA123891_part1_vcf;
    private static File  NA12892_and_NA123891_part2_vcf;
    private static File  NA12892_and_NA123891_part3_vcf;

    final File na12891_r1 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
    final File na12891_r2 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");
    final File na12892_r1 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
    final File na12892_r2 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
    final File na12891_fp = new File(TEST_DATA_DIR, "NA12891.fp.vcf");


    private static final Map<CrosscheckMetric.DataType, List<String>> lookupMap = new HashMap<>(4);

    @BeforeClass
    public void setup() throws IOException {
        NA12891_r1 = SamTestUtils.createIndexedBam(NA12891_r1_sam, NA12891_r1_sam);
        NA12891_r2 = SamTestUtils.createIndexedBam(NA12891_r2_sam, NA12891_r2_sam);
        NA12891_named_NA12892_r1 = SamTestUtils.createIndexedBam(NA12891_named_NA12892_r1_sam, NA12891_named_NA12892_r1_sam);
        NA12892_r1 = SamTestUtils.createIndexedBam(NA12892_r1_sam, NA12892_r1_sam);
        NA12892_r2 = SamTestUtils.createIndexedBam(NA12892_r2_sam, NA12892_r2_sam);

        lookupMap.put(CrosscheckMetric.DataType.FILE, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.FILE).addAll(Arrays.asList("LEFT_FILE", "RIGHT_FILE"));

        lookupMap.put(CrosscheckMetric.DataType.SAMPLE, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(Arrays.asList("LEFT_SAMPLE", "RIGHT_SAMPLE"));
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(lookupMap.get(CrosscheckMetric.DataType.FILE));

        lookupMap.put(CrosscheckMetric.DataType.LIBRARY, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(Arrays.asList("LEFT_LIBRARY", "RIGHT_LIBRARY"));
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(lookupMap.get(CrosscheckMetric.DataType.SAMPLE));

        lookupMap.put(CrosscheckMetric.DataType.READGROUP, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(Arrays.asList("LEFT_RUN_BARCODE", "LEFT_LANE",
                "LEFT_MOLECULAR_BARCODE_SEQUENCE","RIGHT_RUN_BARCODE",
                "RIGHT_LANE", "RIGHT_MOLECULAR_BARCODE_SEQUENCE"));
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(lookupMap.get(CrosscheckMetric.DataType.LIBRARY));


        NA12891_1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.vcf"), "fingerprint");
        NA12891_2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.fp.vcf"), "fingerprint");
        NA12891_named_NA12892_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891_named_NA12892.vcf"), "fingerprint");
        NA12892_1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.vcf"), "fingerprint");
        NA12892_2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.fp.vcf"), "fingerprint");
        NA12891_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.g.vcf"), "fingerprint");
        NA12891_swapped_nonref_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891.with_swapped_NON_REF.g.vcf"), "fingerprint");

        NA12892_g_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.g.vcf"), "fingerprint");
        NA12892_and_NA123891_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892.vcf"), "fingerprint");
        NA12892_and_NA123891_part1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part1.vcf"), "fingerprint");
        NA12892_and_NA123891_part2_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part2.vcf"), "fingerprint");
        NA12892_and_NA123891_part3_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891andNA12892_part3.vcf"), "fingerprint");

    }
    @Test
    void testCheckFingerprint() throws IOException {
        final List<String> args = new ArrayList<>();
        final File outputSummary = File.createTempFile("fingerprint","summary_metrics");
        outputSummary.deleteOnExit();
        final File outputDetail = File.createTempFile("fingerprint","detail_metrics");
        outputSummary.deleteOnExit();

        args.add("INPUT="+na12891_r1.getAbsolutePath());
        args.add("G="+na12891_fp.getAbsolutePath());
        args.add("H="+HAPLOTYPE_MAP.getAbsolutePath());
        args.add("SUMMARY_OUTPUT="+outputSummary.getAbsolutePath());
        args.add("DETAIL_OUTPUT="+outputDetail.getAbsolutePath());

        runPicardCommandLine(args);
    }

}
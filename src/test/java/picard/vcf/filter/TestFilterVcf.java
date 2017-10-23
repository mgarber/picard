/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.vcf.filter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.ListMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.PrintWriter;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Tests for VCF filtration
 */
public class TestFilterVcf {
    private final File INPUT = new File("testdata/picard/vcf/filter/testFiltering.vcf");
    private final File BAD_INPUT = new File("testdata/picard/vcf/filter/testFilteringNoSeqDictionary.vcf");

    /* write content of javascript in the returned file */
	private File quickJavascriptFilter(String content) throws Exception {
		final File out = File.createTempFile("jsfilter", ".js");
		out.deleteOnExit();
		try (final PrintWriter pw = new PrintWriter(out)) {
			pw.println(content);
		}
		return out;
	}

	@Test
	public void testJavaScript() throws Exception {
        final File out = VcfTestUtils.createTemporaryIndexedFile("filterVcfTestJS.", ".vcf");
		final FilterVcf filterer = new FilterVcf();
		filterer.INPUT = INPUT;
		filterer.OUTPUT = out;
		filterer.JAVASCRIPT_FILE = quickJavascriptFilter("variant.getStart()%5 != 0");

		final int retval = filterer.doWork();
		Assert.assertEquals(retval, 0);

		//count the number of reads
		final int expectedNumber = 4;
		int count=0;
		VCFFileReader in = new VCFFileReader(filterer.OUTPUT, false);
		CloseableIterator<VariantContext> iter = in.iterator();
		while(iter.hasNext()) {
			final VariantContext ctx = iter.next();
			count += (ctx.isFiltered()?1:0);
		}
		iter.close();
		in.close();
		Assert.assertEquals(count, expectedNumber);
	}

    /** Returns a sorted copy of the supplied set, for safer comparison. */
    <T extends Comparable> SortedSet<T> sorted(Set<T> in) { return new TreeSet<T>(in); }

    /** Tests that all records get PASS set as their filter when extreme values are used for filtering. */
    @Test public void testNoFiltering() throws Exception {
        final File out = testFiltering(INPUT, ".vcf.gz", 0, 0, 0, Double.MAX_VALUE);
        final VCFFileReader in = new VCFFileReader(out, false);
        for (final VariantContext ctx : in) {
            if (!ctx.filtersWereApplied() || ctx.isFiltered()) {
                Assert.fail("Context should not have been filtered: " + ctx.toString());
            }
        }
        in.close();
    }

    /** Tests that sites with a het allele balance < 0.4 are marked as filtered out. */
    @Test public void testAbFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("tf2", "rs28566954", "rs28548431");
        final File out = testFiltering(INPUT, ".vcf.gz", 0.4, 0, 0, Double.MAX_VALUE);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(sorted(filters.keySet()), sorted(fails), "Failed sites did not match expected set of failed sites.");
    }

    /** Tests that genotypes with DP < 18 are marked as failed, but not >= 18. */
    @Test public void testDpFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs71509448", "rs71628926", "rs13302979", "rs2710876");
        final File out = testFiltering(INPUT, ".vcf.gz", 0, 18, 0, Double.MAX_VALUE);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(sorted(filters.keySet()), sorted(fails), "Failed sites did not match expected set of failed sites.");
    }

    /** Tests that genotypes with DP < 18 are marked as failed, but not >= 18. */
    @Test public void testDpFilteringToVcf() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs71509448", "rs71628926", "rs13302979", "rs2710876");
        final File out = testFiltering(INPUT, ".vcf", 0, 18, 0, Double.MAX_VALUE);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(sorted(filters.keySet()), sorted(fails), "Failed sites did not match expected set of failed sites.");
    }

    /** Tests that genotypes with low GQ are filtered appropriately. */
    @Test public void testGqFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs71509448"); // SNP with GQ=21; lowest GQ in file

        {
            final File out = testFiltering(INPUT, ".vcf.gz", 0, 0, 20, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(filters.size(), 0, "Should not have filtered sites: " + filters);
        }
        {
            final File out = testFiltering(INPUT, ".vcf.gz", 0, 0, 21, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(filters.size(), 0, "Should not have filtered sites: " + filters);
        }
        {
            final File out = testFiltering(INPUT, ".vcf.gz", 0, 0, 22, Double.MAX_VALUE);
            final ListMap<String, String> filters = slurpFilters(out);
            Assert.assertEquals(sorted(filters.keySet()), sorted(fails), "Failed sites did not match expected set of failed sites.");
        }
    }

    /** Tests that genotypes with DP < 18 are marked as failed, but not >= 18. */
    @Test public void testFsFiltering() throws Exception {
        final Set<String> fails = CollectionUtil.makeSet("rs13303033", "rs28548431", "rs2799066");
        final File out = testFiltering(INPUT, ".vcf.gz", 0, 0, 0, 5.0d);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(sorted(filters.keySet()), sorted(fails), "Failed sites did not match expected set of failed sites.");
    }

    @Test public void testCombinedFiltering() throws Exception {
        final TreeSet<String> fails = new TreeSet<String>(CollectionUtil.makeSet("rs13302979", "rs13303033", "rs2710876" , "rs2799066" , "rs28548431", "rs28566954", "rs71509448", "rs71628926", "tf2"));
        final File out = testFiltering(INPUT, ".vcf.gz", 0.4, 18, 22, 5.0d);
        final ListMap<String,String> filters = slurpFilters(out);
        Assert.assertEquals(new TreeSet<String>(filters.keySet()), fails, "Failed sites did not match expected set of failed sites.");
    }

    /** Utility method that takes a a VCF and a set of parameters and filters the VCF. */
    private File testFiltering(final File vcf, final String outputExtension, final double minAb, final int minDp, final int minGq, final double maxFs) throws Exception {
        final File out = VcfTestUtils.createTemporaryIndexedFile("filterVcfTest.", outputExtension);

        final FilterVcf filterer = new FilterVcf();
        filterer.CREATE_INDEX = true;
        filterer.INPUT = vcf;
        filterer.OUTPUT = out;
        filterer.MIN_AB = minAb;
        filterer.MIN_DP = minDp;
        filterer.MIN_GQ = minGq;
        filterer.MAX_FS = maxFs;

        final int retval = filterer.doWork();
        if (retval != 0) {
            throw new PicardException("Return value non-zero: " + retval);
        }

        return out;
    }

    /** Tests that attempting to write to an uncompressed vcf fails if the input has no sequence dictionary */
    @Test(expectedExceptions = PicardException.class)
    public void testFilteringToVcfWithNoSequenceDictionary() throws Exception {
        final File out = File.createTempFile("filterVcfTest.", ".vcf");
        out.deleteOnExit();

        final FilterVcf filterer = new FilterVcf();
        filterer.CREATE_INDEX = true;
        filterer.INPUT = BAD_INPUT;
        filterer.OUTPUT = out;
        filterer.MIN_AB = 0;
        filterer.MIN_DP = 18;
        filterer.MIN_GQ = 0;
        filterer.MAX_FS = Double.MAX_VALUE;

        filterer.doWork();
    }


    /** Consumes a VCF and returns a ListMap where each they keys are the IDs of filtered out sites and the values are the set of filters. */
    private ListMap<String,String> slurpFilters(final File vcf) {
        final ListMap<String,String> map = new ListMap<String, String>();
        final VCFFileReader in = new VCFFileReader(vcf, false);
        for (final VariantContext ctx : in) {
            if (ctx.isNotFiltered()) continue;
            for (final String filter : ctx.getFilters()) {
                map.add(ctx.getID(), filter);
            }
        }
        in.close();
        return map;
    }
}

/*
 * The MIT License
 *
 * Copyright (c) 2017 The University of Massachusetts Medical School
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
package picard.vcf;

import java.io.File;
import java.util.EnumSet;
import java.util.Iterator;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingDeque;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

/**
 * Takes the difference between two VCF files. Throws IllegalArgumentException if the two files are not 
 * in the input files. 
 * <p/>
 * An index file is created for the output file by default. Using an output file name with a
 * ".gz" extension will create gzip-compressed output.
 * <p/>
 * The result is the result of the first file minus the second.
 * 
 * @author Manuel Garber
 * 
 */
@CommandLineProgramProperties(
        summary = "Outputs records that are present in the first file but not the second. VCFs need to be sorted and indexed."
        		+ " An index file is created and a sequence dictionary is required by default.",
        oneLineSummary = "Takes the difference between two VCF files",
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature


public class DifferenceVcfs extends CommandLineProgram{

    @Argument(shortName= "F1", doc="First VCF or BCF file - The result will be this file minus F2", optional=false)
    public File F1;

    @Argument(shortName = "F2", doc = "Second VCF or BCF file - The result will be F1 minus this file", optional = false)
    public File F2;
    
    @Argument(shortName = "SD", doc = "The index sequence dictionary to use instead of the sequence dictionary in the input file", optional = true)
    public File SEQUENCE_DICTIONARY;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The difference between F1 and F2")
    public File OUTPUT;

    private final Log log = Log.getInstance(DifferenceVcfs.class);
    
    public static void main(final String[] argv) {
        new DifferenceVcfs().instanceMainWithExit(argv);
    }

	@Override
	protected int doWork()  {
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        SAMSequenceDictionary sequenceDictionary = null;
        final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
        if (SEQUENCE_DICTIONARY != null) {
        		SAMSequenceDictionaryExtractor.extractDictionary(SEQUENCE_DICTIONARY);
        }
        
        IOUtil.assertFileIsReadable(F1);
        final VCFFileReader F1Reader = new VCFFileReader(F1, false);
        final VCFHeader F1Header = F1Reader.getFileHeader();
        VariantContextComparator variantContextComparator = F1Header.getVCFRecordComparator();
        if(sequenceDictionary == null) {
        		sequenceDictionary = F1Header.getSequenceDictionary();
        }
        
        final VCFFileReader F2Reader = new VCFFileReader(F2, false);
        final VCFHeader F2Header = F2Reader.getFileHeader();
        
        if (!variantContextComparator.isCompatible(F2Header.getContigLines())) {
        		F1Reader.close();
        		F2Reader.close();
        		throw new IllegalArgumentException(
                        "The contig entries in F2 file " + F2.getAbsolutePath() + " are not compatible with F1 " + F1.getAbsolutePath());
        }

        VCFHeader header = new VCFHeader(F1Header);
        header.setSequenceDictionary(sequenceDictionary);
        
        CloseableIterator<VariantContext> f1Iterator = F1Reader.iterator();
        CloseableIterator<VariantContext> f2Iterator = F2Reader.iterator();
        
        
        final VariantContextWriter writer = new VariantContextWriterBuilder()
        		.setOutputFile(OUTPUT)
        		.setReferenceDictionary(sequenceDictionary)
            .setOptions(options).build();

        writer.writeHeader(header);
        /*
         * Here we assume that the files are sorted. If they are not this will fail miserably
         */
        VariantContext f2Cealing = null;
        while(f1Iterator.hasNext()) {
        		VariantContext var = f1Iterator.next();
        		if (f2Cealing == null || variantContextComparator.compare(var, f2Cealing) > 0) {
        			f2Cealing = getNextCeiling(var, f2Iterator, variantContextComparator);
        		}
        		if(f2Cealing == null || (f2Cealing != null && variantContextComparator.compare(var, f2Cealing) < 0)) {
        			writer.add(var);
        		}
        }

        f1Iterator.close();
        f2Iterator.close();
        F1Reader.close();
        F2Reader.close();
		return 0;
	}

	private VariantContext getNextCeiling(VariantContext var, Iterator<VariantContext> varIterator, VariantContextComparator comparator) {
		VariantContext result = null;
		while (varIterator.hasNext() && (result == null || comparator.compare(var, result)>0)) {
			result = varIterator.next();
		}
		
		return result != null && comparator.compare(var, result) <= 0 ? result : null;
		
	}
	

	
}

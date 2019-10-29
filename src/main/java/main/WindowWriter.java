/*
 * Copyright (C) 2014-2016 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package main;

import beagleutil.Samples;
import imp.ImputedVcfWriter;
import blbutil.BGZIPOutputStream;
import blbutil.Utilities;
import imp.ImpData;
import imp.StateProbs;
import ints.IntList;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.VcfWriter;

/**
 * <p>Class {@code WindowWriter} writes VCF and IBD output data.
 * </p>
 * <p>Instances of class {@code WindowWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WindowWriter implements Closeable {

    private final Samples samples;
    private final String outPrefix;
    private final File vcfOutFile;

    /**
     * Constructs a new {@code WindowWriter} object.
     * @param samples the sample whose data will be printed
     * @param outPrefix the output file prefix
     *
     * @throws IllegalArgumentException if {@code outPrefix.length() == 0}
     * @throws NullPointerException if
     * {@code samples == null || outPrefix == null}
     */
    public WindowWriter(Samples samples, String outPrefix) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        if (outPrefix.length()==0) {
            throw new IllegalArgumentException("outPrefix.length()==0");
        }
        this.samples = samples;
        this.outPrefix = outPrefix;
        this.vcfOutFile = new File(outPrefix + ".vcf.gz");

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter vcfOut=new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            boolean ds = true;
            boolean ap = false;
            boolean gp = true;
            boolean gl = false;
            VcfWriter.writeMetaLines(samples.ids(), Main.PROGRAM,
                    ds, ap, gp, gl, vcfOut);
        }
        try {
            try (FileOutputStream fos=new FileOutputStream(vcfOutFile)) {
                fos.write(baos.toByteArray());
            }
        } catch (IOException e) {
            fileOutputError(vcfOutFile, e);
        }
    }

    /**
     * Returns the output file prefix.
     * @return the output file prefix
     */
    public String outPrefix() {
        return outPrefix;
    }

    /**
     * Returns the samples whose data is written by {@code this}.
     * @return the samples whose data is written by {@code this}
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Prints the data in {@code alProbs} for markers
     * with index between {@code refStart} (inclusive) and
     * {@code refEnd} (exclusive) to the output
     * VCF file: {@code this.outPrefix() + ".vcf.gz"}.
     *
     * @param impData the input data for genotype imputation
     * @param stateProbs the imputed state probabilities
     * @param refStart the starting ref marker index (inclusive)
     * @param refEnd the ending ref marker index (exclusive)
     *
     * @throws IllegalArgumentException if
     * {@code stateProbs.size() != impData.nTargHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code refStart < 0 || refEnd > impData.refGT().nMarkers()}
     * @throws NullPointerException if {@code impData==null || stateProbs==null}
     * @throws NullPointerException if any element of {@code stateProbs} is
     * {@code null}
     */
    public void printImputed(ImpData impData, int refStart, int refEnd,
            AtomicReferenceArray<StateProbs> stateProbs) {
        if (stateProbs.length() != impData.nTargHaps()) {
            throw new IllegalArgumentException("inconsistent data:");
        }
        byte[][] output = IntStream.range(0, impData.nClusters())
                .parallel()
                .mapToObj(c -> toByteArray(impData, refStart, refEnd, stateProbs, c))
                .toArray(byte[][]::new);
        print(output, vcfOutFile);
    }

    private static byte[] toByteArray(ImpData impData, int refStart, int refEnd,
            AtomicReferenceArray<StateProbs> stateProbs, int m) {
        ImputedVcfWriter ivw = new ImputedVcfWriter(impData, refStart, refEnd, m);
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter out = new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            ivw.appendRecords(stateProbs, out);
        }
        return baos.toByteArray();
    }

    /**
     * Writes the data in phased genotypes for the specified markers
     * to the output VCF file: {@code this.outPrefix() + ".vcf.gz"}.
     *
     * @param phasedTarg the estimated target haplotypes
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param nThreads the number of parallel threads to use
     *
     * @throws IllegalArgumentException if
     * {@code isImputed.length != alProbs.nMarkers()}
     * @throws IllegalArgumentException if {@code phasedTarg.isPhased() == false}
     * @throws IllegalArgumentException if {@code nThreads < 1}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > phasedTarg.nMarkers() || start > end}
     */
    public void printPhased(GT phasedTarg, int start, int end, int nThreads) {
        int step = 20;
        int[] starts = starts(start, end, step);
        byte[][] output = IntStream.range(0, starts.length-1)
                .parallel()
                .mapToObj(i -> toByteArray(phasedTarg, starts[i], starts[i+1]))
                .toArray(byte[][]::new);
        print(output, vcfOutFile);
    }

    private static int[] starts(int start, int end, int step) {
        IntList starts = new IntList(2 + ((end - start) / step));
        for (int m=start ; m<end; m+=step) {
            starts.add(m);
        }
        starts.add(end);
        return starts.toArray();
    }

    private static byte[] toByteArray(GT phasedTarg, int start, int end) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter vcfOut=new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            VcfWriter.appendRecords(phasedTarg, start, end, vcfOut);
        }
        return baos.toByteArray();
    }

    private static void print(byte[][] output, File outFile)  {
        boolean append = true;
        try {
            try (OutputStream fos = new BufferedOutputStream(
                    new FileOutputStream(outFile, append))) {
                for (byte[] ba : output) {
                   fos.write(ba);
                }
            }
        } catch (IOException e) {
            fileOutputError(outFile, e);
        }
    }

    @Override
    public void close() {
        boolean append = true;
        try {
            try (FileOutputStream fos = new FileOutputStream(vcfOutFile, append);
                    BufferedOutputStream bos = new BufferedOutputStream(fos);
                    BGZIPOutputStream bgzip = new BGZIPOutputStream(bos, true)) {
                // write empty BGZIP block to bgzip by closing bgzip
            }
        } catch (IOException e) {
            Utilities.exit("Error closing file: " + vcfOutFile, e);
        }
    }

    private static void fileOutputError(File file, Exception e) {
        Utilities.exit("Error writing to file: " + file, e);
    }
}

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
package vcf;

import beagleutil.Samples;
import blbutil.Const;
import blbutil.StringUtil;
import blbutil.Utilities;
import java.io.File;
import java.util.Arrays;

/**
 * <p>Class {@code VcfRecGTParser} parses VCF records and extracts the GT format
 * field.
 * </p>
 * <p>Instances of class {@code VcfRecGTParser} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecGTParser {

    private final VcfHeader vcfHeader;
    private final String vcfRec;
    private final Marker marker;
    private final int nAlleles;
    private final int nSamples;
    private final int ninthTabPos;

    /**
     * Constructs a new {@code VcfRecGTParser} object from the specified VCF
     * record.
     * @param vcfHeader the VCF meta-information lines and header line
     * @param vcfRec the VCF record
     * @throws IllegalArgumentException if {@code vcfHeader.nSamples() == 0}
     * @throws IllegalArgumentException if a format error is detected in the
     * {@code vcfRecord}
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRec == null}
     */
    public VcfRecGTParser(VcfHeader vcfHeader, String vcfRec) {
        if (vcfHeader.nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        this.vcfHeader = vcfHeader;
        this.vcfRec = vcfRec;
        this.marker = new BasicMarker(vcfRec);
        this.nAlleles = marker.nAlleles();
        this.nSamples = vcfHeader.nSamples();
        this.ninthTabPos = ninthTabPos(vcfRec);
    }

    private static int ninthTabPos(String vcfRec) {
        int pos = -1;
        for (int j=0; j<9; ++j) {
            pos = vcfRec.indexOf(Const.tab, pos + 1);
            if (pos == -1) {
                throw new IllegalArgumentException(
                        "VCF record format error: " + vcfRec);
            }
        }
        return pos;
    }

    /**
     * Returns the VCF meta-information lines and header line for the backing
     * VCF record
     * @return the VCF meta-information lines and header line
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    /**
     * Returns the backing VCF record.
     * @return the backing VCF record
     */
    public String vcfRecord() {
        return vcfRec;
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns {@code this.marker().nAlleles()}.
     * @return the number of alleles
     */
    public int nAlleles() {
        return nAlleles;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return vcfHeader.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return nSamples;
    }

    /**
     * Stores the genotypes genotypes in the specified long arrays.  The
     * contract for this method is undefined if any bit is set in an array
     * is set when the method is invoked, or if any array is not sufficiently
     * long.
     * @param allele1 a long array in which the first allele for each
     * sample is stored
     * @param allele2 a long array set in which the second allele for each
     * sample is stored
     * @param isMissing1 a long array whose {@code k}-th bit will be set
     * if the first allele of the {@code k}-th sample is missing
     * @param isMissing2 a long array set whose {@code k}-th bit will be set
     * if the second allele of the {@code k}-th sample is missing
     * @param isPhased a along array whose {@code k}-th bit will be set
     * if the phased allele separator is present in the {@code k}-th sample
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if any parameter is {@code null}
     */
    public void storeAlleles(long[] allele1, long[] allele2,
            long[] isMissing1, long[] isMissing2, long[] isPhased) {
        int bitsPerLong = 6;
        int bitsPerAllele = bitsPerAllele(marker);
        int pos = ninthTabPos;
        int unfilt = -1;
        for (int s=0; s<nSamples; ++s) {
            if (pos == -1) {
                throwFieldCountError();
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                pos = vcfRec.indexOf(Const.tab, pos + 1);
                if (pos == -1) {
                    throwFieldCountError();
                }
            }
            int alEnd1 = alEnd1(vcfRec, pos + 1);   // checks allele separator
            int alEnd2 = alEnd2(vcfRec, alEnd1 + 1);
            int a1 = parseAllele(pos + 1, alEnd1);
            int a2 =  parseAllele(alEnd1 + 1, alEnd2);
            char sep = vcfRec.charAt(alEnd1);
            if (sep == Const.phasedSep) {
                isPhased[s >> bitsPerLong] |= (1L << s);
            }
            if (a1 == -1) {
                isMissing1[s >> bitsPerLong] |= (1L << s);
            }
            else {
                storeAllele(allele1, s, bitsPerAllele, a1);
            }
            if (a2 == -1) {
                isMissing2[s >> bitsPerLong] |= (1L << s);
            }
            else {
                storeAllele(allele2, s, bitsPerAllele, a2);
            }
            pos = vcfRec.indexOf(Const.tab, alEnd2);
        }
    }

    private static int bitsPerAllele(Marker marker) {
        int nAllelesM1 = marker.nAlleles() - 1;
        int nStorageBits = Integer.SIZE - Integer.numberOfLeadingZeros(nAllelesM1);
        return nStorageBits;
    }

    private static void storeAllele(long[] al, int sample, int bitsPerAllele,
            int allele) {
        int bitsPerLong = 6;
        int index = sample*bitsPerAllele;
        int mask = 1;
        for (int k=0; k<bitsPerAllele; ++k) {
            if ((allele & mask)==mask) {
                al[index >> bitsPerLong] |= (1L << index);
            }
            ++index;
            mask <<= 1;
        }
    }

    /**
     * Returns the list of phased alleles in the backing VCF record.
     * @return the list of phased alleles in the backing VCF record
     * @throws IllegalArgumentException if  the VCF record contains an
     * unphased or missing genotype
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     */
    public int[] phasedAlleles() {
        int[] alleles = new int[2*nSamples];
        int pos = ninthTabPos;
        int unfilt = -1;
        for (int s=0, hap=0; s<nSamples; ++s) {
            if (pos == -1) {
                throwFieldCountError();
            }
            int nextUnfiltered = vcfHeader.unfilteredSampleIndex(s);
            while (++unfilt < nextUnfiltered) {
                pos = vcfRec.indexOf(Const.tab, pos + 1);
                if (pos == -1) {
                    throwFieldCountError();
                }
            }
            int alEnd1 = alEnd1(vcfRec, pos + 1);
            int alEnd2 = alEnd2(vcfRec, alEnd1 + 1);
            int a1 = parseAllele(pos + 1, alEnd1);
            int a2 =  parseAllele(alEnd1 + 1, alEnd2);
            char sep = vcfRec.charAt(alEnd1);
            if (sep!=Const.phasedSep || a1 == -1 || a2 == -1) {
                String msg = "Unphased or missing reference genotype at marker: "
                        + marker;
                Utilities.exit(new IllegalArgumentException(msg));
            }
            alleles[hap++] = a1;
            alleles[hap++] = a2;
            pos = vcfRec.indexOf(Const.tab, alEnd2);
        }
        return alleles;
    }

    /* returns exclusive end */
    private static int alEnd1(String rec, int start) {
        if (start==rec.length()) {
            throwGTFormatError(rec, rec.length());
        }
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.unphasedSep || c == Const.phasedSep) {
                return index;
            }
            else if (c == Const.colon || c == Const.tab) {
                throwGTFormatError(rec, index+1);
            }
            ++index;
        }
        if (index==rec.length()) {
            throwGTFormatError(rec, rec.length());
        }
        return index;
    }

    /* returns exclusive end */
    private static int alEnd2(String rec, int start) {
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.colon || c == Const.tab) {
                return index;
            }
            ++index;
        }
        return index;
    }

    private int parseAllele(int start, int end) {
        if (start==end) {
            String s = "Missing sample allele: " + vcfRec;
            throw new IllegalArgumentException(s);
        }
        int al;
        if (start + 1 == end) {
            char c = vcfRec.charAt(start);
            if (c=='.') {
                return -1;
            }
            else {
                al = (c - '0');
            }
        }
        else {
            al = Integer.parseInt(vcfRec.substring(start, end));
        }
        if (al < 0 || al >= nAlleles) {
            String strAllele = vcfRec.substring(start, end);
            String s = "invalid allele [" + strAllele + "]: " + Const.nl + vcfRec;
            throw new IllegalArgumentException(s);
        }
        return al;
    }

    private static void throwGTFormatError(String rec, int index) {
        StringBuilder sb = new StringBuilder(1000);
        sb.append("ERROR: genotype is missing allele separator:");
        sb.append(Const.nl);
        sb.append(rec.substring(0, index));
        sb.append(Const.nl);
        sb.append("Exiting Program");
        sb.append(Const.nl);
        Utilities.exit(sb.toString());
    }

    private void throwFieldCountError() {
        File f = vcfHeader.file();
        String[] fields = StringUtil.getFields(vcfRec, Const.tab);
        StringBuilder sb = new StringBuilder(1000);
        sb.append("VCF header line has ");
        sb.append(vcfHeader.nHeaderFields());
        sb.append(" fields, but data line has ");
        sb.append(fields.length);
        sb.append(" fields");
        sb.append(Const.nl);
        sb.append("File source: ");
        sb.append((f!=null ? f.toString() : "stdin"));
        sb.append(Const.nl);
        sb.append(Arrays.toString(fields));
        sb.append(Const.nl);
        Utilities.exit(sb.toString());
    }

    /**
     * Returns an array of length {@code this.nAlleles()} whose
     * {@code k}-th element is the list of haplotype indices carrying
     * the {@code k}-th allele if {@code k} is a non-major allele,
     * and whose {@code k}-th element is {@code null} if {@code k} is
     * the major allele.  If there is more than one allele with maximal count,
     * the allele with maximal count having the smallest index is defined to
     * be the major allele.
     * @return the indices of the haplotypes carrying each non-major allele
     * @throws IllegalArgumentException if a format error is detected in
     * the specified VCF record or if the specified VCF header is
     * inconsistent with the specified VCF header.
     *
     * @throws NullPointerException if {@code vcfRec == null || rec == null}
     */
    public int[][] nonMajRefIndices() {
        int[] alleles = phasedAlleles();
        int[] alCnts = new int[nAlleles];
        for (int a : alleles) {
            ++alCnts[a];
        }
        int majAl = 0;
        for (int j=1; j<nAlleles; ++j) {
            if (alCnts[j]>alCnts[majAl]) {
                majAl = j;
            }
        }
        int[][] nonMajIndices = new int[nAlleles][];
        for (int al=0; al<nAlleles; ++al) {
            nonMajIndices[al] = al==majAl ? null : new int[alCnts[al]];
        }
        Arrays.fill(alCnts, 0);
        for (int j=0; j<alleles.length; ++j) {
            int al = alleles[j];
            if (al!=majAl) {
                nonMajIndices[al][alCnts[al]++] = j;
            }
        }
        return nonMajIndices;
    }
}

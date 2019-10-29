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

import blbutil.Const;
import java.io.PrintWriter;

/**
 * <p>Class {@code VcfRecBuilder} contains methods for constructing
 * and printing a VCF record in VCF 4.2 format.  The FORMAT field data
 * for each sample is added sequentially to the record via the
 * {@code addSampleData()} method.
 *
 * </p>
 * <p>Instances of class {@code VcfRecBuilder} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecBuilder {

    /**
     * The default initial size for the string buffer, which is 50
     * characters.
     */
    public static final int DEFAULT_INIT_SIZE = 50;

    private final StringBuilder sb;
    private final Marker marker;
    private final int nAlleles;

    /**
     * Constructs a new {@code VcfRecBuilder} instance with an initial
     * capacity for the specified number of samples.
     *
     * @param marker the marker
     * @param nSamples the number of samples
     * @throws IllegalArgumentException if {@code nSamples < 0}
     */
    public VcfRecBuilder(Marker marker, int nSamples) {
        if (nSamples < 0) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        StringBuilder buffer = new StringBuilder(100 + 4*nSamples);
        writeFixedFields(marker, buffer);
        this.marker = marker;
        this.nAlleles = marker.nAlleles();
        this.sb = buffer;
    }

    /**
     * Returns the current marker.  Returns {@code null} if
     * {@code this.reset()} has not been previously invoked.
     * @return the current marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Adds the FORMAT field for a sample to the VCF record for the current
     * marker.
     * @param a1 the first allele
     * @param a2 the second allele
     * @throws IllegalStateException if {@code this.marker() == null}
     * @throws IndexOutOfBoundsException if
     * {@code a1 < 0 || a1 >= this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code a2 < 0 || a2 >= this.marker().nAlleles()}
     */
    public void addSampleData(int a1, int a2) {
        if (a1 < 0 || a1 >= nAlleles) {
            throw new IndexOutOfBoundsException(String.valueOf(a1));
        }
        if (a2 < 0 || a2 >= nAlleles) {
            throw new IndexOutOfBoundsException(String.valueOf(a2));
        }
        sb.append(Const.tab);
        sb.append(a1);
        sb.append(Const.phasedSep);
        sb.append(a2);
    }

    /**
     * Prints the current VCF record for the current marker to the specified
     * {@code PrintWriter}.
     * @param out the {@code PrintWriter} to which the VCF record will be
     * printed
     * @throws NullPointerException if {@code out == null}
     */
    public void writeRec(PrintWriter out) {
        out.println(sb);
    }

    private static void writeFixedFields(Marker marker, StringBuilder sb) {
        appendMarker(marker, sb);
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);     // QUAL
        sb.append(Const.tab);
        sb.append("PASS");                      // FILTER
        sb.append(Const.tab);
        if (marker.end()==-1) {
            sb.append(Const.MISSING_DATA_CHAR); // INFO
        }
        else {
            sb.append("END=");
            sb.append(marker.end());
        }
        sb.append(Const.tab);
        sb.append("GT");                        // FORMAT
    }

    private static void appendMarker(Marker marker, StringBuilder sb) {
        sb.append(marker.chrom());
        sb.append(Const.tab);
        sb.append(marker.pos());
        int nIds = marker.nIds();
        if (nIds==0) {
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
}
        else {
            for (int j=0; j<nIds; ++j) {
                sb.append(j==0 ? Const.tab : Const.semicolon);
                sb.append(marker.id(j));
            }
        }
        int nAlleles = marker.nAlleles();
        if (nAlleles==1) {
            sb.append(Const.tab);
            sb.append(marker.allele(0));
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<nAlleles; ++j) {
                sb.append(j<2 ? Const.tab : Const.comma);
                sb.append(marker.allele(j));
            }
        }
    }
}

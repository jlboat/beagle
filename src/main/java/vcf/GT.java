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

/**
 * <p>Interface {@code GT} represents genotype data
 * for a list of markers and a list of samples.
 * </p>
 * <p>All instances of {@code GT} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GT {

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the number of haplotypes.  The returned value is equal to
     * {@code 2*this.nSamples()}.
     * @return the number of haplotypes
     */
    int nHaps();

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

   /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Returns the first allele for the specified marker and sample
     * or return -1 if the allele is missing.  The two alleles for a
     * sample are arbitrarily ordered if
     * {@code this.unphased(marker, sample) == false}.
     * @param marker the marker index
     * @param sample the sample index
     * @return the first allele for the specified marker and sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele1(int marker, int sample);

    /**
     * Returns the second allele for the specified marker and sample
     * or return -1 if the allele is missing.  The two alleles for a
     * sample are arbitrarily ordered if
     * {@code this.unphased(marker, sample) == false}.
     * @param marker the marker index
     * @param sample the sample index
     * @return the  allele for the specified marker and sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele2(int marker, int sample);

    /**
     * Returns the allele on the specified haplotype for the specified marker
     * or return -1 if the allele is missing.  The two alleles for an
     * individual are arbitrarily ordered if
     * {@code this.unphased(marker, hap/2) == false}.
     * @param marker the marker index
     * @param hap the haplotype index
     * @return the allele on the specified haplotype for the specified marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap  >= this.nHaps()}
     */
    int allele(int marker, int hap);

    /**
     * Returns {@code true} if the genotype for each marker and sample
     * is a phased, non-missing genotype, and returns {@code false} otherwise.
     * @return {@code true} if the genotype for each marker and sample
     * is a phased, non-missing genotype
     */
    boolean isPhased();

    /**
     * Returns a {@code GT} instance restricted to genotype data for
     * the specified markers.
     * @param markers the list of markers in the returned instance
     * @param indices a list of distinct marker indices (from
     * {@code this.markers())} in increasing order
     * @return a {@code GT} instance restricted to genotype data for
     * the specified markers
     *
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= this.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (this.marker(indices[j]).equals(markers.marker(j)) == false)}
     * @throws NullPointerException if {@code indices == null}
     */
    GT restrict(Markers markers, int[] indices);

    /**
     * Returns the probability of the observed data for the specified marker
     * and sample if the specified pair of unordered alleles is the true
     * genotype.
     * @param gt the genotype data
     * @param marker the marker index
     * @param sample the sample index
     * @param allele1 the first allele index
     * @param allele2 the second allele index
     * @return the probability of the observed data for the specified marker
     * and sample if the specified pair of ordered alleles is the true
     * ordered genotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || samples >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker(marker).nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker(marker).nAlleles()}
     *
     * @implSpec The default implementation returns {@code 1.0f} if the
     * corresponding genotype determined by the {@code isPhased()},
     * {@code allele1()}, and {@code allele2()} methods is consistent
     * with the specified ordered genotype, and returns {@code 0.0f} otherwise.
     */
    static float gl(GT gt, int marker, int sample, int allele1, int allele2) {
        int nAlleles = gt.marker(marker).nAlleles();
        if (allele1 < 0 || allele1 >= nAlleles)  {
            String s = "invalid alleles: (" + allele1 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        if (allele2 < 0 || allele2 >= nAlleles) {
            String s = "invalid alleles: (" + allele2 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        int a1 = gt.allele1(marker, sample);
        int a2 = gt.allele2(marker, sample);
        boolean consistent = (a1==-1 || a1==allele1) && (a2==-1 || a2==allele2);
        if (consistent==false) {
            consistent = (a1==-1 || a1==allele2) && (a2==-1 || a2==allele1);
        }
        return consistent ? 1.0f : 0.0f;
    }
}

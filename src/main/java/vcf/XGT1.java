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

// Later consider relaxing mutability, and determine whether interface
// can be replaced with implementing class

/**
 * <p>Interface {@code XGT1} (Genotype Likelihoods) represents genotype
 * likelihoods for one sample.
 * </p>
 * <p>Instances of {@code XGT1} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface XGT1 {


    /**
     * Returns the probability of the observed data for the specified marker
     * if the specified pair of ordered alleles is the true ordered genotype.
     * @param marker the marker index
     * @param allele1 the first allele index
     * @param allele2 the second allele index
     * @return the probability of the observed data for the specified marker
     * and sample if the specified pair of ordered alleles is the true
     * ordered genotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker(marker).nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker(marker).nAlleles()}
     */
    float gl(int marker, int allele1, int allele2);

    /**
     * Returns the first allele for the specified marker if the observed data
     * include a non-missing allele, and returns -1 otherwise. Alleles are
     * arbitrarily ordered if the genotype is unphased.
     * @param marker the marker index
     * @return the first allele for the specified marker if the observed data
     * include a non-missing allele, and -1 otherwise
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele1(int marker);

    /**
     * Returns the second allele for the specified marker if the observed data
     * include a non-missing allele, and returns -1 otherwise. Alleles are
     * arbitrarily ordered if the genotype is unphased.
     * @param marker the marker index
     * @return the second allele for the specified marker if the observed data
     * include a non-missing allele, and -1 otherwise
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele2(int marker);

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the specified marker.
     * @param marker the marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the sample identifier index.
     * @return the sample identifier index
     */
    int idIndex();

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    String toString();
}

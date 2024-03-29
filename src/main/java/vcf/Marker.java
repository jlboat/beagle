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

/**
 * <p>Interface {@code Marker} represents a genetic marker.
 * </p>
 * <p>All instances of class {@code Marker} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Marker extends Comparable<Marker> {

    /**
     * Returns the chromosome.
     * @return the chromosome
     */
    String chrom();

    /**
     * Returns the chromosome index.
     * @return the chromosome index
     */
    int chromIndex();

    /**
     * Returns the chromosome position coordinate.
     * @return the chromosome position coordinate
     */
    int pos();

    /**
     * Returns the number of marker identifiers.
     * @return the number of marker identifiers
     */
    int nIds();

    /**
     * Returns the specified marker identifier.
     * @param index a marker identifier index
     * @return the specified marker identifier
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nIds()}
     */
    String id(int index);

    /**
     * Returns the first marker identifier if there is at least
     * one identifier in the VCF record ID field, and returns
     * {@code this.chr() + ":" + this.pos()} otherwise.
     *
     * @return a marker identifier
     */
    String id();

    /**
     * Returns the number of alleles for the marker, including the REF
     * allele.
     * @return the number of alleles for the marker, including the REF
     * allele
     */
    int nAlleles();

    /**
     * Returns the specified allele.  The reference allele has index 0.
     * @param index an allele index
     * @return the specified allele
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nAlleles()}
     */
    String allele(int index);

    /**
     * Returns the alleles. The {@code k}-th element of the returned array
     * is equal to {@code this.allele(k)}.
     * @return the alleles
     */
    String[] alleles();

    /**
     * Returns the number of distinct genotypes, which equals
     * {@code this.nAlleles()*(1 + this.nAlleles())/2}.
     *
     * @return the number of distinct genotypes
     */
    int nGenotypes();

    /**
     * Returns the INFO END field, or -1 if there is no INFO END field.
     *
     * @return the INFO END field, or -1 if there is no INFO END field
     */
    int end();

    /**
     * Returns {@code true} if the specified object is a
     * {@code Marker} with the same chromosome,
     * position, allele lists, and INFO END field, and
     * returns {@code false} otherwise.  Equality does not
     * depend on value of the VCF record ID field.
     *
     * @param obj object to be compared with {@code this} for equality
     *
     * @return {@code true} if the specified object is a
     * {@code Marker} with the same chromosome,
     * position, and allele lists, and INFO END field
     */
    @Override
    boolean equals(Object obj);

    /**
     * <p>Returns the hash code value for this object. The hash code does not
     * depend on value of the VCF record ID field.
     * The hash code is defined by the following calculation:
     * </p>
     * <pre>
     *   int hash = 5;
     *   hash = 29 * hash + this.chromIndex();
     *   hash = 29 * hash + this.pos();
     *   for (int j=0, n=this.nAlleles(); j&lt;n; ++j) {
     *       hash = 29 * hash + alleles[j].hashCode();
     *   }
     *   hash = 29 * hash + end();
     * </pre>
     *
     * @return the hash code value for this marker
     */
    @Override
    int hashCode();

    /**
     * Compares this marker with the specified marker
     * for order, and returns a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal to,
     * or greater than the specified marker.  Comparison is
     * on chromosome index, position, allele identifier lists, and end value
     * in that order.  Allele identifier lists are compared for
     * lexicographical order, and alleles are compared using the
     * {@code String compareTo()} method.
     *
     * @param other the {@code Marker} to be compared
     * @return a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal,
     * or greater than the specified marker
     */
    @Override
    int compareTo(Marker other);

    /**
     * Returns a string equal to the first five tab-delimited fields
     * of a VCF record corresponding to this marker.
     *
     * @return a string equal to the first five tab-delimited fields
     * of a VCF record corresponding to this marker
     */
    @Override
    String toString();
}

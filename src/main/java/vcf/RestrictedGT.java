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
 * <p>Class {@code RestrictedGT} is a wrapper for a {@code GT}
 * instance that restricts the data to a subset of the VCF records.
 * </p>
 * <p>Instances of class {@code RestrictGT} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RestrictedGT implements GT {

    private final GT gt;
    private final Markers restrictedMarkers;
    private final int[] inclusionMap;

    /**
     * Constructs a new {@code NoPhaseGTWrapper} instance from the specified
     * data
     * @param gt the genotypes to be wrapped
     * @param markers the list of markers in the returned instance
     * @param indices the mapping of marker indices from {@code markers}
     * to {@code gt.markers()}
     * @throws IndexOutOfBoundsException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (indices[j] < 0 || indices[j] >= gt.nMarkers())}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (1 <= j && j < indices.length)} such that
     * {@code (indices[j] <= indice[j - 1])}
     * @throws IllegalArgumentException if there exists {@code j} such that
     * {@code (0 <= j && j < indices.length)} such that
     * {@code (gt.marker(indices[j]).equals(markers.marker(j)) == false)}
     * @throws NullPointerException if
     * {@code (gt == null || markers == null || include == null)}
     */
    public RestrictedGT(GT gt, Markers markers, int[] indices) {
        for (int j=0; j<indices.length; ++j) {
            if (j>0 && indices[j] <= indices[j-1]) {
                throw new IllegalArgumentException(String.valueOf(indices[j]));
            }
            if (gt.marker(indices[j]).equals(markers.marker(j))==false) {
                throw new IllegalArgumentException(markers.marker(j).toString());
            }
        }
        this.gt = gt;
        this.restrictedMarkers = markers;
        this.inclusionMap = indices.clone();
    }

    @Override
    public boolean isPhased() {
        return gt.isPhased();
    }

    @Override
    public int allele1(int marker, int sample) {
        return gt.allele1(inclusionMap[marker], sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return gt.allele2(inclusionMap[marker], sample);
    }


    @Override
    public int allele(int marker, int hap) {
        return gt.allele(inclusionMap[marker], hap);
    }

    @Override
    public int nMarkers() {
        return restrictedMarkers.nMarkers();
    }

    @Override
    public Marker marker(int marker) {
        return restrictedMarkers.marker(marker);
    }

    @Override
    public Markers markers() {
        return restrictedMarkers;
    }

    @Override
    public int nHaps() {
        return gt.nHaps();
    }

    @Override
    public int nSamples() {
        return gt.nSamples();
    }

    @Override
    public Samples samples() {
        return gt.samples();
    }

    @Override
    public GT restrict(Markers markers, int[] indices) {
        return new RestrictedGT(this, markers, indices);
    }

    @Override
    public String toString() {
        return RestrictedGT.class.toString() + " : " + gt.toString();
    }
}

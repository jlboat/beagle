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
import ints.LongArray;
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * <p>Class {@code HapsGT} implements the {@code vcf.GT} interface by
 * storing the haplotypes as bit arrays.
 * </p>
 * <p>Instances of {@code HapsGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapsGT implements GT {

    private final Markers markers;
    private final Samples samples;
    private final LongArray[] haps;

    /**
     * Constructs a new {@code HapsGT} instance from the specified data.
     * Two haplotypes for the {@code k}-th sample are required to
     * stored in the {@code 2*k} and {@code 2*k + 1} elements of the
     * {@code haps} array.
     * @param markers the list of markers
     * @param samples the list of samples
     * @param haps the list of haplotypes
     * @throws IllegalArgumentException if
     * {@code haps.length != 2*samples.nSamples()}
     * @throws IllegalArgumentException if there exists a {@code j} such
     * that {@code (0 <= j && j < haps.length)} and
     * {@code (haps[j].size() < (markers.sumHaplotypeBits() + 63)/64)}
     * @throws NullPointerException if
     * {@code markers == null || samples == null || haps == null}
     * @throws NullPointerException if any elements of {@code haps} is
     * {@code null}
     */
    public HapsGT(Markers markers, Samples samples, LongArray[] haps) {
        if (haps.length != 2*samples.nSamples()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int minSize = (markers.sumHaplotypeBits() + 63)/64;
        for (LongArray hap : haps) {
            if (hap.size() < minSize) {
                throw new IllegalArgumentException("inconsistent data");
            }
        }
        this.markers = markers;
        this.samples = samples;
        this.haps = haps;
    }

    /**
     * Constructs a new {@code HapsGT} instance from the specified data.
     * Two haplotypes for the {@code k}-th sample are required to
     * stored in the {@code 2*k} and {@code 2*k + 1} elements of the
     * {@code haps} array.
     * @param markers the list of markers
     * @param samples the list of samples
     * @param haps the list of haplotypes
     * @throws IllegalArgumentException if
     * {@code haps.length() != 2*samples.nSamples()}
     * @throws IllegalArgumentException if there exists a {@code j} such
     * that {@code (0 <= j && j < haps.length)} and
     * {@code (haps[j].size() < (markers.sumHaplotypeBits() + 63)/64)}
     * @throws NullPointerException if
     * {@code markers == null || samples == null || haps == null}
     * @throws NullPointerException if any elements of {@code haps} is
     * {@code null}
     */
    public HapsGT(Markers markers, Samples samples,
            AtomicReferenceArray<LongArray> haps) {
        if (haps.length() != 2*samples.nSamples()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        LongArray[] hapsArray = new LongArray[haps.length()];
        int minSize = (markers.sumHaplotypeBits() + 63)/64;
        for (int h=0; h<hapsArray.length; ++h) {
            LongArray hap = haps.get(h);
            if (hap.size() < minSize) {
                throw new IllegalArgumentException("inconsistent data");
            }
            hapsArray[h] = hap;
        }
        this.markers = markers;
        this.samples = samples;
        this.haps = hapsArray;
    }

    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public int allele1(int marker, int hapPair) {
        return markers.bitsToAllele(haps[hapPair<<1], marker);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return markers.bitsToAllele(haps[(hapPair<<1)+1], marker);
    }

    @Override
    public int allele(int marker, int haplotype) {
        return markers.bitsToAllele(haps[haplotype], marker);
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public int nHaps() {
        return haps.length;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public GT restrict(Markers markers, int[] indices) {
        return new RestrictedGT(this, markers, indices);
    }

    @Override
    public Samples samples() {
        return samples;
    }
}

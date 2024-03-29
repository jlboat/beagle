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
import main.Pedigree;

/**
 * <p>Class {@code XBasicGT} represents genotype and genotype emission
 * probabilities for a set of samples optimized by sample.
 * </p>
 * Instances of class {@code XBasicGTWindow} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class XBasicGT implements GT {

    private final Samples samples;
    private final Markers markers;
    private final boolean isRefData;
    private final XBasicGT1[] gt1;

    /**
     * Constructs a {@code XBasicGT} instance from the specified data.
     *
     * @param gt the genotype likelihoods
     * @param ped the pedigrees
     *
     * @throws IllegalArgumentException if
     * {@code gl.samples().equals(ped.samples())==false}
     * @throws NullPointerException if {@code gl == null || ped == null}
     */
    public XBasicGT(GT gt, Pedigree ped) {
        if (gt.samples().equals(ped.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        int nSamples = gt.nSamples();
        this.markers = gt.markers();
        this.samples = gt.samples();
        this.isRefData = gt.isPhased();
        this.gt1 = new XBasicGT1[nSamples];
        for (int s=0; s<nSamples; ++s) {
            int father = ped.father(s);
            int mother = ped.mother(s);
            gt1[s] = new XBasicGT1(gt, s, father, mother);
        }
    }

    @Override
    public boolean isPhased() {
        return isRefData;
    }

    @Override
    public int allele1(int marker, int sample) {
        return gt1[sample].allele1(marker);
    }

    @Override
    public int allele2(int marker, int sample) {
        return gt1[sample].allele2(marker);
    }

    @Override
    public int allele(int marker, int hap) {
        int sample = hap>>1;
        if ((hap & 1) == 0) {
            return allele1(marker, sample);
        }
        else {
            return allele2(marker, sample);
        }
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Marker marker(int markerIndex) {
        return markers.marker(markerIndex);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nHaps() {
        return 2*samples.nSamples();
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public GT restrict(Markers markers, int[] indices) {
        return new RestrictedGT(this, markers, indices);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[XBasicGTWindow: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        sb.append(Const.nl);
        for (int m=0, n=nMarkers(); m<n; ++m) {
            sb.append(markers.marker(m));
            sb.append(Const.nl);
            sb.append(Const.MISSING_DATA_CHAR);     // QUAL
            sb.append(Const.tab);
            sb.append("PASS");                      // FILTER
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);     // INFO
            sb.append(Const.tab);
            sb.append("GT");                        // FORMAT
            for (XGT1 s : gt1) {
                sb.append(Const.tab);
                sb.append(s.allele1(m));
                sb.append(Const.unphasedSep);
                sb.append(s.allele2(m));
            }
        }
        sb.append(Const.nl);
        sb.append(']');
        return sb.toString();
    }
}

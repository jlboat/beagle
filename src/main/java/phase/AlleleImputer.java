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
package phase;

import blbutil.Utilities;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.GT;

/**
 * <p>Class {@code AlleleImputer} imputes alleles from marker allele
 * frequencies.
 * </p>
 * <p>Instances of class {@code AlleleImputer} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class AlleleImputer {

    private final double[][] alFreq;

    /**
     * Constructs a new {@code AlleleImputer} instance from the specified
     * data.
     * @param ref the reference genotypes, or {@code null} if there are none
     * @param targ the target genotypes
     * @param seed the seed for random number generation
     * @throws IllegalArgumentException if
     * {@code ref.markers().equals(targ.markers()) == false}
     * @throws NullPointerException if {@code targ == null}
     */
    public AlleleImputer(GT targ, GT ref, long seed) {
        if (ref!=null && ref.markers().equals(targ.markers())==false) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int maxHaps = 2000;
        Random rand = new Random(seed);
        int[] targHaps = hapIndices(targ, maxHaps, rand);
        int[] refHaps = hapIndices(ref, (maxHaps - targHaps.length), rand);
        this.alFreq = IntStream.range(0, targ.nMarkers())
                .parallel()
                .mapToObj(m -> alFreq(targ, targHaps, ref, refHaps, m))
                .toArray(double[][]::new);
    }

    private static int[] hapIndices(GT gt, int maxHaps, Random rand) {
        if (gt==null || maxHaps<=0) {
            return new int[0];
        }
        else {
            int[] hapIndices = IntStream.range(0, gt.nHaps()).toArray();
            if (maxHaps < hapIndices.length) {
                int nHaps = Math.min(gt.nHaps(), maxHaps);
                Utilities.shuffle(hapIndices, nHaps, rand);
                hapIndices = Arrays.copyOf(hapIndices, nHaps);
                Arrays.sort(hapIndices);
            }
            return hapIndices;
        }
    }

    private static double[] alFreq(GT targ, int[] targHaps,
            GT ref, int[] refHaps, int marker) {
        int nAlleles = targ.marker(marker).nAlleles();
        int[] cnts = new int[nAlleles];
        double[] freq = new double[nAlleles];
        for (int h : targHaps) {
            int allele = targ.allele(marker, h);
            if (allele>=0) {
                ++cnts[allele];
            }
        }
        for (int h : refHaps) {
            int allele = ref.allele(marker, h);
            ++cnts[allele];
        }
        int sum = 0;
        for (int cnt : cnts) {
            sum += cnt;
        }
        if (sum>0) {
            for (int al=0; al<cnts.length; ++al) {
                freq[al] = (double) cnts[al] / sum;
            }
        }
        else {
            freq[0] = 1.0;
        }
        return freq;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return alFreq.length;
    }

    /**
     * Returns an imputed allele for the specified marker.
     * @param marker a marker index
     * @param rand a random number generator
     * @return an imputed allele for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws NullPointerException if {@code rand == null}
     */
    public int imputeAllele(int marker, Random rand) {
        int nAlleles = alFreq[marker].length;
        double d = rand.nextDouble();
        double sum = 0.0;
        for (int j=0; j<nAlleles; ++j) {
            sum += alFreq[marker][j];
            if (sum >= d) {
                return j;
            }
        }
        return nAlleles - 1;
    }
}

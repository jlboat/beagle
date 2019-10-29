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

import beagleutil.Samples;
import blbutil.FloatList;
import ints.IntList;
import ints.LongArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.HapsGT;
import vcf.Markers;

/**
 * <p>Class {@code ImputableHaps} performs imputation to estimate
 * missing alleles and missing haplotype phase.
 * </p>
 * <p>Instances of {@code ImputableHaps} are thread-safe.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapImputer {

    private final Markers markers;
    private final Samples samples;
    private final int nHaps;
    private final AtomicReferenceArray<LongArray> haps;
    private final AtomicReferenceArray<List<PartlyImputedAllele>> missing;

    /**
     * Constructs a new {@code ImputedHaps} instance from the specified data.
     * @param markers the list of markers
     * @param samples the list of target samples
     * @throws NullPointerException if
     * {@code markers == null || samples == null}
     */
    public HapImputer(Markers markers, Samples samples) {
        this.markers = markers;
        this.samples = samples;
        this.nHaps = 2*samples.nSamples();
        this.haps = new AtomicReferenceArray<>(nHaps);
        this.missing = new AtomicReferenceArray<>(nHaps);
        for (int j=0; j<nHaps; ++j) {
            missing.set(j, new ArrayList<>(4));
        }
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targSamples() {
        return samples;
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Stores the specified haplotype.
     * @param hap the haplotype index
     * @param alleles the haplotype
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= 2*this.samples()}
     * @throws IllegalArgumentException if
     * {@code alleles.length != this.markers().nMarkers()}
     * @throws NullPointerException if {@code alleles == null}
     */
    public void setHap(int hap, int[] alleles) {
        haps.set(hap, markers.allelesToBits(alleles));
    }

    /**
     * Stores a partially-imputed allele.  if
     * {@code alProbs.length < this.markers().marker(marker).nAlleles()},
     * the missing allele probabilities are assumed to be 0.0.
     * @param hap the haplotype index
     * @param marker the marker index
     * @param alProbs the posterior allele probabilities computed from a
     * subset of HMM states
     * @param refHaps the list of reference haplotypes for HMM states
     * not included in the specified allele probabilities
     * @param stateProbs the list of state probabilities for each HMM state
     * in the specified list of reference haplotypes
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= 2 * this.targSamples().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.markers().nMarkers()}
     * @throws IllegalArgumentException if
     * {@code refHaps.size() != stateProbs.size()}
     * @throws NullPointerException if
     * {@code alProbs == null || refHaps == null || stateProbs == null}
     */
    public void setPartlyImputedAllele(int hap, int marker, float[] alProbs,
                IntList refHaps, FloatList stateProbs) {
        if (marker < 0 || marker >= markers.nMarkers()) {
            throw new IndexOutOfBoundsException(String.valueOf(marker));
        }
        if (refHaps.size()!=stateProbs.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        PartlyImputedAllele pia = new PartlyImputedAllele(markers, marker,
                alProbs, refHaps, stateProbs);
        missing.get(hap).add(pia);
    }

    /**
     * Returns the imputed haplotypes.
     * @return the imputed haplotypes
     * @throws NullPointerException if {@code this.setHap()} has not previously
     * been called for each haplotype {@code h} satisfying
     * {@code (0 <= h && h <= 2 * this.targSamples().nSamples())}
     */
    public GT imputedHaps() {
        GT oldGT = new HapsGT(markers, samples, haps);
        LongArray[] newHaps = IntStream.range(0, nHaps)
                .parallel()
                .mapToObj(h -> imputeHap(h, oldGT))
                .toArray(LongArray[]::new);
        return new HapsGT(markers, samples, newHaps);
    }

    private LongArray imputeHap(int h, GT oldGT) {
        LongArray oldHap = haps.get(h);
        List<PartlyImputedAllele> partlyImputedAlleles = missing.get(h);
        if (partlyImputedAlleles.isEmpty()) {
            return oldHap;
        }
        else {
            int[] alleles = markers.bitsToAlleles(oldHap);
            for (int j=0, n=partlyImputedAlleles.size(); j<n; ++j) {
                PartlyImputedAllele pia = partlyImputedAlleles.get(j);
                alleles[pia.marker()] = pia.imputeAllele(oldGT);
            }
            return markers.allelesToBits(alleles);
        }
    }

    private static class PartlyImputedAllele {

        private final int marker;
        private final float[] alProbs;
        private final int[] refHaps;
        private final float[] stateProbs;

        private PartlyImputedAllele(Markers markers, int marker,
                float[] alProbs, IntList refHaps, FloatList stateProbs) {
            int nAlleles =  markers.marker(marker).nAlleles();
            this.marker = marker;
            this.alProbs = Arrays.copyOf(alProbs, nAlleles);
            this.refHaps = refHaps.toArray();
            this.stateProbs = stateProbs.toArray();
        }

        private int marker() {
            return marker;
        }

        private int imputeAllele(GT phasedGT) {
            assert phasedGT.isPhased();
            float[] updatedProbs = alProbs.clone();
            for (int j=0; j<refHaps.length; ++j) {
                int allele = phasedGT.allele(marker, refHaps[j]);
                updatedProbs[allele] += stateProbs[j];
            }
            int imputedAllele = 0;
            for (int j=1; j<updatedProbs.length; ++j) {
                if (updatedProbs[j] > updatedProbs[imputedAllele]) {
                    imputedAllele = j;
                }
            }
            return imputedAllele;
        }
    }
}

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

import vcf.HapsGT;
import ints.IntArray;
import ints.IntList;
import ints.LongArray;
import java.util.Random;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.function.IntConsumer;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.Markers;
import vcf. RefGT;

/**
 * <p>Class {@code EstPhase} stores original input genotype data.
 * the current estimated haplotype pair for each target sample, a list of
 * missing genotypes for each target sample, and a list of remaining
 * unphased heterozygote genotypes for each target sample.
 * </p>
 * <p>Instances of class {@code EstPhase} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class EstPhase {

    private final GT targGT;
    private final AtomicReferenceArray<LongArray> haps;
    private final AtomicReferenceArray<IntArray> unphased;
    private final AtomicReferenceArray<IntArray> missing;

    /**
     * Constructs a new {@code EstPhase} instance from the specified data.
     * @param targGT the original input target genotype data
     * @param refGT the reference genotype data or {@code null} if there
     * are no reference samples
     * @param overlap the number of initial markers with phased, non-missing
     * genotypes due to overlap with the previous marker window
     * @param seed the seed for random number generation
     * @throws IllegalArgumentException if
     * {@code refGT != null && targGT.markers().equals(refGT.markers()) == false}
     * @throws IllegalArgumentException if {@code overlap < 0}
     * @throws NullPointerException if {@code targGT == null}
     */
    public EstPhase(GT targGT, RefGT refGT, int overlap, long seed) {
        if (refGT!=null && targGT.markers().equals(refGT.markers())==false) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (overlap<0) {
            throw new IllegalArgumentException(String.valueOf(overlap));
        }
        AlleleImputer imputer = new AlleleImputer(targGT, refGT, seed);
        this.targGT = targGT;
        this.haps = new AtomicReferenceArray<>(targGT.nHaps());
        this.unphased = new AtomicReferenceArray<>(targGT.nSamples());
        this.missing = new AtomicReferenceArray<>(targGT.nSamples());
        IntStream.range(0, targGT.nSamples())
                .parallel()
                .forEach(initializer(targGT, overlap, imputer, seed));
    }

    private IntConsumer initializer(GT targ, int overlap, AlleleImputer imputer,
            long seed) {
        return s -> {
            int nMarkers = targ.nMarkers();
            Random rand = new Random(seed + s);
            IntList unphList = new IntList(nMarkers>>6 + 1);
            IntList missList = new IntList(nMarkers>>10 + 1);
            int[] hap1 = new int[targ.nMarkers()];
            int[] hap2 = new int[targ.nMarkers()];
            boolean foundFirstHet = setOverlap(targ, overlap, s, hap1, hap2);
            int start = overlap;
            for (int m=start; m<nMarkers; ++m) {
                int a1 = targ.allele1(m, s);
                int a2 = targ.allele2(m, s);
                if (a1<0 || a2<0) {
                    missList.add(m);
                    a1 = imputer.imputeAllele(m, rand);
                    a2 = imputer.imputeAllele(m, rand);
                }
                else if (a1!=a2) {
                    if (rand.nextBoolean()) {
                        int tmp = a1;
                        a1 = a2;
                        a2 = tmp;
                    }
                    if (foundFirstHet) {
                        unphList.add(m);
                    }
                    else {
                        foundFirstHet = true;
                    }
                }
                hap1[m] = a1;
                hap2[m] = a2;
            }
            setHapPair(s, hap1, hap2);
            missing.set(s, IntArray.create(missList, targ.nMarkers()));
            unphased.set(s, IntArray.create(unphList, targ.nMarkers()));
        } ;
    }

    private static boolean setOverlap(GT gt, int overlap, int sample, int[] hap1,
            int[] hap2) {
        boolean foundFirstHet = false;
        for (int m=0; m<overlap; ++m) {
            hap1[m] = gt.allele1(m, sample);
            hap2[m] = gt.allele2(m, sample);
            if (hap1[m]<0 || hap2[m]<0) {
                throw new IllegalArgumentException("inconstent data");
            }
            if (hap1[m]!=hap2[m]) {
                foundFirstHet = true;
            }
        }
        return foundFirstHet;
    }

    /**
     * Returns the input target genotype data.
     * @return the input target genotype data
     */
    public GT targGT() {
        return targGT;
    }

    /**
     * Returns a list of remaining marker indices in increasing order for which
     * the specified sample has an unphased heterozygote genotype.
     * @param sample the sample index
     * @return a list of remaining marker indices in increasing order for which
     * the specified sample has an unphased heterozygote genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    public IntArray unphased(int sample) {
        return unphased.get(sample);
    }

    /**
     * Sets the list of remaining marker indices with unphased heterozygote
     * genotypes for the specified sample to the specified list
     * @param sample the sample index
     * @param newUnphased a list of remaining marker indices with unphased
     * heterozygote genotypes
     * @throws IllegalArgumentException if the specified {@code newUnphased}
     * list is not sorted in increasing order, contains a duplicate elements,
     * or is not a subset of {@code this.getUnphased(sample)}.
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.targGT().nSamples()}
     * @throws NullPointerException if {@code newUnphased == null}
     */
    public void setUnphased(int sample, IntArray newUnphased) {
        IntArray oldUnphased = unphased.get(sample);
        int oldSize = oldUnphased.size();
        int nextOldIndex = 0;
        int oldMkr = -1;
        for (int j=0, n=newUnphased.size(); j<n; ++j) {
            int newMkr = newUnphased.get(j);
            if (nextOldIndex<oldSize) {
                oldMkr = oldUnphased.get(nextOldIndex++);
            }
            while (oldMkr!=newMkr && nextOldIndex<oldSize) {
                oldMkr = oldUnphased.get(nextOldIndex++);
            }
            if (oldMkr!=newMkr) {
                throw new IllegalArgumentException(newUnphased.toString());
            }
        }
        unphased.set(sample, newUnphased);
    }

    /**
     * Returns a list of marker indices in increasing order for which
     * the specified sample has a missing genotype.
     * @param sample the sample index
     * @return a list of marker indices in increasing order for which
     * the specified sample has a missing genotype
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.targGT().nSamples()}
     */
    public IntArray missing(int sample) {
        return missing.get(sample);
    }

    /**
     * Sets the haplotype pair for the specified sample to the specified
     * haplotypes.
     * @param  sample the sample index
     * @param hap1 an array whose elements are the estimated alleles
     * carried by the sample's first haplotype.
     * @param hap2 an array whose elements are the estimated alleles
     * carried by the sample's second haplotype.
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.targGT().nSamples()}
     * @throws IllegalArgumentException if
     * {@code hap1.length != this.targGT().nMarkers()}
     * @throws IllegalArgumentException if
     * {@code  hap2.length != this.targGT().nMarkers()}
     * @throws IllegalArgumentException if
     * {@code (hap1[k] < 0 || hap1[k] >= this.targGT().marker(k).nAlleles())}
     * for any {@code k} satisfying {@code (0 <= k && k < this.targGT().nMarkers())}
     * @throws IllegalArgumentException if
     * {@code (hap2[k] < 0 || hap2[k] >= this.targGT().marker(k).nAlleles())}
     * for any {@code k} satisfying {@code (0 <= k && k < this.targGT().nMarkers())}
     * @throws NullPointerException if {@code hap1 == null || hap2 == null}
     */
    public void setHapPair(int sample, int[] hap1, int[] hap2) {
        Markers markers = targGT.markers();
        int h1 = sample<<1;
        int h2 = h1 | 0b1;
        haps.set(h1, markers.allelesToBits(hap1));
        haps.set(h2, markers.allelesToBits(hap2));
    }

    /**
     * Sets the {@code k}-th element of the specified {@code hap1} and
     * {@code hap2} arrays to the specified sample's phased genotype
     * at the {@code k}-th marker.
     * @param sample the sample index
     * @param hap1 an array whose elements will be set to the
     * current estimated alleles carried by the sample's first haplotype.
     * @param hap2 an array whose elements will be set to the
     * current estimated allele carried by the sample's second haplotype.
     * @throws IndexOutOfBoundsException if
     * {@code hap1.length != this.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2.length != this.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.targGT().nSamples()}
     * @throws NullPointerException if {@code hap1 == null || hap2 == null}
     */
    public void getHaps(int sample, int[] hap1, int[] hap2) {
        Markers markers = targGT.markers();
        int nMarkers = markers.nMarkers();
        if (hap1.length != nMarkers || hap2.length != nMarkers) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int h1 = sample<<1;
        int h2 = h1 | 0b1;
        LongArray bits1 = haps.get(h1);
        LongArray bits2 = haps.get(h2);
        for (int m=0; m<nMarkers; ++m) {
            hap1[m] = markers.bitsToAllele(bits1, m);
            hap2[m] = markers.bitsToAllele(bits2, m);
        }
    }

    /**
     * Returns the current estimated phased genotypes for the target samples.
     * @return the current estimated phased genotypes for the target samples
     */
    public HapsGT hapsGT() {
        return new HapsGT(targGT.markers(), targGT.samples(), haps);
    }
}

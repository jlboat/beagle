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

import main.CompHapSegment;
import blbutil.Utilities;
import ints.IntIntMap;
import ints.IntList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.IntStream;

/**
 * <p>Class {@code PhaseStates} has methods for constructing a Li and
 * Stephens HMM for a target haplotype or target sample.
 * </p>
 * <p>Instances of {@code PhaseStates} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PhaseStates {

    private static final int NIL = -103;
    private final PhaseIbs ibsHaps;
    private final PhaseData phaseData;
    private final int nMarkers;
    private final int maxStates;
    private final int minSteps;

    private final IntIntMap hapToEnd;
    private final PriorityQueue<CompHapSegment> q;
    private final IntList[] compositeHapToHap;
    private final IntList[] compositeHapToEnd;

    private final int[] compHapToListIndex;
    private final int[] compHapToHap;
    private final int[] compHapToEnd;

    /**
     * Constructs a new {@code PhaseStates} object from the specified data.
     * @param ibsHaps the IBS haplotype segments
     * @param maxStates the maximum number of composite reference
     * haplotypes that will be constructed
     * @throws IllegalArgumentException if {@code maxStates < 1}
     * @throws NullPointerException if {@code ibsHaps == null}
     */
    public PhaseStates(PhaseIbs ibsHaps, int maxStates) {
        if (maxStates < 1) {
            throw new IllegalArgumentException(String.valueOf(maxStates));
        }
        this.ibsHaps = ibsHaps;
        this.phaseData = ibsHaps.phaseData();
        this.nMarkers = phaseData.targGT().nMarkers();
        this.maxStates = maxStates;
        int defaultMinSteps = 200;
        float scaleFactor = phaseData.par().scaleFactor();
        this.minSteps = (int) Math.ceil(defaultMinSteps*scaleFactor);
        this.hapToEnd = new IntIntMap(maxStates);
        this.q = new PriorityQueue<>(maxStates);
        this.compositeHapToHap = IntStream.range(0, maxStates)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        this.compositeHapToEnd = IntStream.range(0, maxStates)
                .mapToObj(j -> new IntList())
                .toArray(IntList[]::new);
        this.compHapToListIndex = new int[maxStates];
        this.compHapToHap = new int[maxStates];
        this.compHapToEnd = new int[maxStates];
    }

    /**
     * Returns the maximum number of HMM states at a marker.
     * @return the maximum number of HMM states at a marker
     */
    public int maxStates() {
        return maxStates;
    }

    /**
     * Stores the Li and Stephens HMM for the specified target sample in
     * the specified arrays.  The number of allele mismatches (0 or 1)
     * between {@code hap1[m]} and {@code hap2[m]} for the {@code j}-th state
     * are stored in {@code nMismatchs1[m][j]} and
     * {@code nMismatches2[m][j]} respectively.
     *
     * @param sample the target sample index
     * @param hap1 the alleles on the first target haplotype
     * @param hap2 the alleles on the second target haplotype
     * @param missAlleles a list of arrays for storing HMM state alleles
     * at markers for which one or both target haplotypes have a missing allele
     * @param nMismatches1 a two-dimensional array in which the number
     * of allele mismatches (0 or 1) for {@code hap1} for each HMM state
     * will be stored
     * @param nMismatches2 a two-dimensional array in which the number
     * of allele mismatches (0 or 1) for {@code hap2} for each HMM state
     * will be stored
     * @return the number of state alleles at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= phaseData.targGT().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code hap1.length < this.phaseData().targGT.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2.length < this.phaseData().targGT.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code nMismatches1.length < this.phaseData().targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code nMismatches2.length < this.phaseData().targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if {@code missAlleles.length} is
     * less than the number of markers for which
     * {@code (hap1[m] < 0 || hap2[m] < 0)}
     * @throws IndexOutOfBoundsException if {@code missAlleles.get(j)}
     * is less than the number of model states for any {@code j}
     * satisfying {@code (0 <= j && j < misAlleles.size())}
     * @throws IndexOutOfBoundsException if {@code nMismatches1[m].length}
     * is less than the number of model states for any marker {@code m}
     * satisfying {@code (0 <= m && m < nMismatches1.length)}
     * @throws IndexOutOfBoundsException if {@code nMismatches2[m].length}
     * is less than the number of model states for any marker {@code m}
     * satisfying {@code (0 <= m && m < nMismatches2.length)}
     * @throws NullPointerException if any array is {@code null} or if
     * {@code missAlleles == null}
     */
    public int ibsStates(int sample, int[] hap1, int[] hap2,
            List<int[]> missAlleles, byte[][] nMismatches1, byte[][] nMismatches2) {
        int h1 = sample << 1;
        int h2 = (h1 | 0b1);
        initializeFields();
        for (int step=0, n=phaseData.codedSteps().nSteps(); step<n; ++step) {
            ibsHaps.addIbsHaps(h1, step, this);
            ibsHaps.addIbsHaps(h2, step, this);
        }
        if (q.isEmpty()) {
            fillQWithRandomHaps(h1);
        }
        int numStates = copyData(hap1, hap2, missAlleles, nMismatches1, nMismatches2);
        return numStates;
    }

    /**
     * Stores the Li and Stephens HMM for the specified target
     * haplotype in the specified arrays.  The haplotype for the
     * {@code j}-th state at the {@code m}-th marker is stored
     * in {@code haps[m][j]}.  The number of allele mismatches (0 or 1)
     * between the haplotype for the {@code j}-th state and the
     * target haplotype at the {@code m}-th marker is stored in
     * {@code nMismatches[m][j]}.
     * The number of HMM states states at each marker is returned.
     * @param targHap the haplotype index
     * @param haps the two-dimensional array in which the
     * haplotype for each HMM state will be stored
     * @param nMismatches the two-dimensional array in which the number
     * of allele mismatches (0 or 1) for each HMM state will be stored
     * @return the number of HMM states at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code targHap < 0 || targHap >= this.phaseData().targGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code haps.length < this.phaseData.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code nMismatches.length < this.phaseData.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if {@code haps[m].length}
     * is less than the number of model states for any marker {@code m}
     * satisfying {@code (0 <= m && m < haps.length)}
     * @throws IndexOutOfBoundsException if {@code nMismatches[m].length}
     * is less than the number of model states for any marker {@code m}
     * satisfying {@code (0 <= m && m < nMismatches.length)}
     * @throws NullPointerException if any array is {@code null}
     */
    public int ibsStates(int targHap, int[][] haps, byte[][] nMismatches) {
        initializeFields();
        for (int j=0, n=phaseData.codedSteps().nSteps(); j<n; ++j) {
            ibsHaps.addIbsHaps(targHap, j, this);
        }
        if (q.isEmpty()) {
            fillQWithRandomHaps(targHap);
        }
        int numStates = copyData(targHap, haps, nMismatches);
        return numStates;
    }

    private void initializeFields() {
        hapToEnd.clear();
        for (int j=0, n=q.size(); j<n; ++j) {
            compositeHapToHap[j].clear();
            compositeHapToEnd[j].clear();
        }
        q.clear();
    }

    /**
     * Add the specified segment of the specified reference haplotype to
     * the HMM state space.
     * @param hap a haplotype index
     * @param step a genomic interval
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.phaseData().targGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.phaseData().codedStates().nSteps()}
     */
    public void updateFields(int hap, int step) {
        if (hapToEnd.get(hap, NIL)==NIL) { // hap not currently in q
            updateHeadOfQ();
            if (q.size()==maxStates
                    || (q.isEmpty()==false && step - q.peek().step() > minSteps)) {
                CompHapSegment head = q.poll();
                int modEnd = phaseData.codedSteps().stepStart((head.step() + step) >>> 1);
                hapToEnd.remove(head.hap());
                compositeHapToHap[head.compHapIndex()].add(hap);      // hap of new segment
                compositeHapToEnd[head.compHapIndex()].add(modEnd);   // end of old segment
                head.updateHap(hap);
                head.updateStep(step);
                q.offer(head);
            }
            else {
                int compHapIndex = q.size();
                compositeHapToHap[compHapIndex].add(hap);            // hap of new segment
                q.offer(new CompHapSegment(hap, step, compHapIndex));
            }
        }
        hapToEnd.put(hap, step);
    }

    private void updateHeadOfQ() {
        CompHapSegment head = q.peek();
        if (head!=null) {
            int latestEnd = hapToEnd.get(head.hap(), NIL);
            while (head.step()!=latestEnd) {
                head = q.poll();
                head.updateStep(latestEnd);
                q.offer(head);
                head = q.peek();
                latestEnd = hapToEnd.get(head.hap(), NIL);
            }
        }
    }

    private int copyData(int[] hap1, int[] hap2, List<int[]> refAl,
            byte[][] nMismatches1, byte[][] nMismatches2) {
        int nCompHaps = q.size();
        initializeCopy(nCompHaps);
        int missIndex = 0;
        for (int m=0; m<nMarkers; ++m) {
            boolean isMissing = hap1[m] == -1 || hap2[m] == -1;
            int[] copiedAlleles = isMissing ? refAl.get(missIndex++) : null;
            for (int j=0; j<nCompHaps; ++j) {
                if (m==compHapToEnd[j]) {
                    ++compHapToListIndex[j];
                    compHapToHap[j] = compositeHapToHap[j].get(compHapToListIndex[j]);
                    compHapToEnd[j] = compositeHapToEnd[j].get(compHapToListIndex[j]);
                }
                int refAllele = phaseData.allele(m, compHapToHap[j]);
                if (isMissing) {
                    copiedAlleles[j] = refAllele;
                    nMismatches1[m][j] = 0;
                    nMismatches2[m][j] = 0;
                }
                else {
                    nMismatches1[m][j] = (refAllele==hap1[m] ? (byte) 0 : (byte) 1);
                    nMismatches2[m][j] = (refAllele==hap2[m] ? (byte) 0 : (byte) 1);
                }
            }
        }
        return nCompHaps;
    }

    private int copyData(int targHap, int[][] haps, byte[][] nMismatches) {
        int nCompositeHaps = q.size();
        initializeCopy(nCompositeHaps);
        for (int m=0; m<nMarkers; ++m) {
            int targAl = phaseData.allele(m, targHap);
            for (int j=0; j<nCompositeHaps; ++j) {
                if (m==compHapToEnd[j]) {
                    ++compHapToListIndex[j];
                    compHapToHap[j] = compositeHapToHap[j].get(compHapToListIndex[j]);
                    compHapToEnd[j] = compositeHapToEnd[j].get(compHapToListIndex[j]);
                }
                int refHap = compHapToHap[j];
                haps[m][j] = refHap;
                nMismatches[m][j] = phaseData.allele(m, refHap)==targAl
                        ? (byte) 0 : (byte) 1;
            }
        }
        return nCompositeHaps;
    }

    private void initializeCopy(int nSlots) {
        for (int j=0; j<nSlots; ++j) {
            compositeHapToEnd[j].add(nMarkers); // add missing end of last segment
            compHapToListIndex[j] = 0;
            compHapToHap[j] = compositeHapToHap[j].get(0);
            compHapToEnd[j] = compositeHapToEnd[j].get(0);
        }
    }

    private void fillQWithRandomHaps(int hap) {
        assert q.isEmpty();
        int nHaps = phaseData.nHaps();
        int nStates = Math.min(nHaps-2, maxStates);
        if (nStates<=0) {
            Utilities.exit("ERROR: there is only one sample");
        }
        else {
            int sample = hap>>1;
            Random rand = new Random(phaseData.seed() + hap);
            for (int i=0; i<nStates; ++i) {
                int h = rand.nextInt(nHaps);
                while ((h>>1)==sample) {
                    h = rand.nextInt(nHaps);
                }
                compositeHapToHap[q.size()].add(h);
                q.add(new CompHapSegment(h, nMarkers, i));
            }
        }
    }
}

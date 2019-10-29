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

import ints.IndexArray;
import ints.IntArray;
import ints.IntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;
import main.Par;
import vcf.Markers;

/**
 * <p>Class {@code BestMatch} has methods for finding a haplotype
 * that shares a long IBS sequence containing a specified genomic interval
 * with a target haplotype and that originates from a sample which
 * carries a rare variant in the specified genomic interval
 * that is also carried by the target sample.</p>
 *
 * <p>Instances of {@code BestMatch} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BestMatch {

    private final FixedPhaseData fpd;
    private final int nTargHaps;
    private final CodedSteps codedSteps;
    private final int nBufferSteps;

    /**
     * Constructs a new {@code BestMatch} instance from the specified data
     * @param phaseData the input data for the next phase update
     * @throws IllegalArgumentException if
     * {@code fpd.hiFreqTargGT() != codedSteps.phaseData().targGT()}
     * @throws NullPointerException if
     * {@code phaseData == null}
     */
    public BestMatch(PhaseData phaseData) {
        Par par = phaseData.par();
        float nSteps = par.buffer()/par.phase_step();
        this.fpd = phaseData.fpd();
        this.nTargHaps = phaseData.targGT().nHaps();
        this.codedSteps = phaseData.codedSteps();
        this.nBufferSteps = (int) Math.rint(nSteps);
    }

    /**
     * Returns the coded steps.
     * @return the coded steps
     */
    public CodedSteps codedSteps() {
        return codedSteps;
    }

    /**
     * Returns an array mapping each target haplotype to a matching haplotype.
     * The target and mathching haplotypes must originate from distinct samples
     * that carry a shared low frequency variant in a specified genomic interval
     * and that share a long IBS sequence in the forward direction
     * containing the specified genomic interval. The target haplotype is
     * mapped to {@code -1} if no such matching haplotype exists.
     * @param step a step index
     * @return return an array mapping each target haplotypes to a matching
     * haplotype or to {@code -1}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.codedSteps().nSteps()}
     */
    public int[] fwdMatch(int step) {
        int[] bestHap = IntStream.range(0, nTargHaps).map(i -> -1).toArray();
        int endStep = Math.min(step + nBufferSteps, codedSteps.nSteps());
        Random rand = new Random(fpd.par().seed() + step);
        List<IntList> hapLists = lowFreqHapLists(step);
        for (int j=step; j<endStep && hapLists.size()>0; ++j) {
            hapLists = nextHapLists(hapLists, codedSteps.get(j));
            updateBestHap(hapLists, bestHap, rand);
        }
        return bestHap;
    }

    /**
     * Returns an array mapping each target haplotype to a matching haplotype.
     * The target and mathching haplotypes must originate from distinct samples
     * that carry a shared low frequency variant in a specified genomic interval
     * and that share a long IBS sequence in the backward direction
     * containing the specified genomic interval. The target haplotype is
     * mapped to {@code -1} if no such matching haplotype exists.
     * @param step a step index
     * @return return an array mapping each target haplotypes to a matching
     * haplotype or to {@code -1}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.codedSteps().nSteps()}
     */
    public int[] bwdMatch(int step) {
        int[] bestHap = IntStream.range(0, nTargHaps).map(i -> -1).toArray();
        int endStep = Math.max(step - nBufferSteps, -1);
        Random rand = new Random(fpd.par().seed() + 1009 + step);
        List<IntList> hapLists = lowFreqHapLists(step);
        for (int s=step; s>endStep && hapLists.size()>0; --s) {
            hapLists = nextHapLists(hapLists, codedSteps.get(s));
            updateBestHap(hapLists, bestHap, rand);
        }
        return bestHap;
    }

    private List<IntList> lowFreqHapLists(int step) {
        IntArray hiFreqIndices = fpd.hiFreqIndices();
        int start = step==0 ? 0 : hiFreqIndices.get(codedSteps.stepStart(step));
        int end = (step+1 < codedSteps.nSteps())
                ? hiFreqIndices.get(codedSteps.stepStart(step+1))
                : fpd.targGT().nMarkers();
        return lowFreqHapLists(fpd, start, end);
    }

    private static List<IntList> lowFreqHapLists(FixedPhaseData fpd, int start,
            int end) {
        List<IntList> hapLists = new ArrayList<>();
        Markers markers = fpd.targGT().markers();
        for (int m=start; m<end; ++m) {
            int nAlleles = markers.marker(m).nAlleles();
            for (int al=0; al<nAlleles; ++al) {
                IntArray carriers = fpd.carriers(m, al);
                if (carriers.size()>1) {
                    hapLists.add(hapList(carriers));
                }
            }
        }
        return hapLists;
    }

    private static IntList hapList(IntArray carriers) {
        IntList hapList = new IntList(2*carriers.size());
        for (int j=0, n=carriers.size(); j<n; ++j) {
            int sample = carriers.get(j);
            int h1 = sample<<1;
            hapList.add(h1);
            hapList.add(h1|0b1);
        }
        return hapList;
    }

    private List<IntList> nextHapLists(List<IntList> hapLists,
            IndexArray codedStep) {
        List<IntList> nextHapLists = new ArrayList<>();
        List<IntList> children = new ArrayList<>();
        IntArray hap2Seq = codedStep.intArray();
        int[] seq2Child = new int[codedStep.valueSize()];
        for (IntList hapList : hapLists) {
            getChildren(hapList, hap2Seq, seq2Child, children);
            setNextHapLists(children, nextHapLists);
        }
        return nextHapLists;
    }

    private void getChildren(IntList hapList, IntArray hap2Seq, int[] seq2Child,
            List<IntList> childList) {
        Arrays.fill(seq2Child, -1);
        childList.clear();
        for (int j=0, n=hapList.size(); j<n; ++j) {
            int hap = hapList.get(j);
            int seq = hap2Seq.get(hap);
            if (seq2Child[seq] == -1) {
                seq2Child[seq] = childList.size();
                childList.add(new IntList(8));
            }
            childList.get(seq2Child[seq]).add(hap);
        }
    }

    private void setNextHapLists(List<IntList> childList,
            List<IntList> nextHapLists) {
        for (int j=0, n=childList.size(); j<n; ++j) {
            IntList hapList = childList.get(j);
            if (hapList.size()>=2) {
                int hap0 = hapList.get(0);
                if (hap0<nTargHaps  // haplist is sorted in increasing order
                        && (hapList.size()>2 || (hap0^hapList.get(1))!=1)) {
                    nextHapLists.add(hapList);
                }
            }
        }
    }

    private void updateBestHap(List<IntList> hapLists, int[] bestHap,
            Random rand) {
        hapLists.stream()
                .forEach(haps -> updateBestHap(haps, bestHap, rand));
    }

    private void updateBestHap(IntList haps, int[] bestHap, Random rand) {
        for (int j=0, n=haps.size(); j<n && haps.get(j)<bestHap.length; ++j) {
            int hap = haps.get(j);
            int sample = hap >> 1;
            int i = rand.nextInt(n);
            int h = haps.get(i);
            while ((h>>1)==sample) {
                if ((++i)==n) {
                    i = 0;
                }
                h = haps.get(i);
            }
            bestHap[hap] = h;
        }
    }
}

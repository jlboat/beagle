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

import blbutil.DoubleArray;
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import vcf.GT;
import vcf.MarkerMap;
import vcf.RefGT;

/**
 * <p>Class {@code Ibs2} stores IBS2 segments that any target sample shares
 * with another target or reference sample.</p>
 *
 * <p>Instances of {@code Ibs2} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Ibs2 {

    private static final int MIN_STEP_MARKERS = 100;
    private static final int MAX_STEP_MARKERS = 1500; // value for high-frequency markers
    private static final double MAX_IBD_GAP_CM = 4.0;

    private final GT targGT;
    private final RefGT refGT;
    private final DoubleArray genPos;
    private final double minCmIbs2;
    private final SampleSeg[][] sampleSegs; // [targ sample][segment]

    /**
     * Constructs a new {@code Ibs2} instance from the specified data.
     * @param targGT the target genotype data
     * @param refGT the reference genotype data or {@code null} if there are
     * no reference data
     * @param map an array with the genetic map positions of each marker
     * @param minCmIbs2 the minimum cM length of a stored IBS2 segment
     * @throws NullPointerException if {@code (targGT == null || map == null)}
     * @throws IllegalArgumentException if
     * {@code targGT.nMarkers() != map.genPos().size()}
     * @throws IllegalArgumentException if
     * {@code (refGT != null && refGT.markers().equals(targGT.markers())==false)}
     * @throws IllegalArgumentException if
     * {@code (minCmIbs2 <= 0 || Double.isFinite(minCmIbs2) == false)}
     */
    public Ibs2(GT targGT, RefGT refGT, MarkerMap map, double minCmIbs2) {
        checkArgs(targGT, refGT, map, minCmIbs2);
        this.targGT = targGT;
        this.refGT = refGT;
        this.genPos = map.genPos();
        this.minCmIbs2 = minCmIbs2;
        IntArray windowStarts = stepStarts(map, 0.5*minCmIbs2);

        int[][][] idSets = IntStream.range(0, windowStarts.size())
                .parallel()
                .mapToObj(w -> ibsSamples(targGT, refGT, windowStarts, w))
                .toArray(int[][][]::new); //[window][targ_sample][ibs2_samples]

        int nMarkersM1 = targGT.nMarkers() - 1;
        Predicate<SampleSeg> predicate =
                ss -> (genPos.get(ss.inclEnd()) - genPos.get(ss.start())) >= minCmIbs2;
        this.sampleSegs = IntStream.range(0, targGT.nSamples())
                .parallel()
                .mapToObj(s -> {
                    SampleSeg[] list = segList(s, nMarkersM1, windowStarts, idSets);
                    SampleSeg[] merged1 = merge(list, genPos, MAX_IBD_GAP_CM);
                    SampleSeg[] extended = extend(targGT, refGT, s, merged1);
                    SampleSeg[] merged2 = merge(extended, genPos, MAX_IBD_GAP_CM);
                    return Arrays.stream(merged2)
                            .filter(predicate)
                            .toArray(SampleSeg[]::new);
                })
                .toArray(SampleSeg[][]::new);
    }

    private static void checkArgs(GT targGT, RefGT refGT, MarkerMap map,
            double minIbs2Cm) {
        if (targGT.nMarkers()!=map.genPos().size()) {
            throw new IllegalArgumentException(String.valueOf(map.genPos().size()));
        }
        if (refGT!=null && refGT.markers().equals(targGT.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (minIbs2Cm <= 0 || Double.isFinite(minIbs2Cm)==false) {
            throw new IllegalArgumentException(String.valueOf(minIbs2Cm));
        }
    }

    private static IntArray stepStarts(MarkerMap map, double minCM) {
        double[] genPos = map.genPos().toArray();
        IntList indices = new IntList(genPos.length/100);

        int nextStart = 0;
        while (nextStart<genPos.length) {
            indices.add(nextStart);
            nextStart = nextStart(genPos, nextStart, minCM);
        }
        if (popLastInterval(indices, genPos, minCM)) {
            indices.pop();
        }
        return new WrappedIntArray(indices.toArray());
    }

    private static int nextStart(double[] genPos, int start, double minCM) {
        double nextGenPos = genPos[start] + minCM;
        int nextIndex = Arrays.binarySearch(genPos, start, genPos.length,
                nextGenPos);
        if (nextIndex<0) {
            nextIndex = -nextIndex-1;
        }
        int nextStart = nextIndex;

        int minNextStart = start + MIN_STEP_MARKERS;
        int maxNextStart = start + MAX_STEP_MARKERS;
        if (nextStart < minNextStart) {
            nextStart = minNextStart;
        }
        else if (nextStart > maxNextStart) {
            nextStart = maxNextStart;
        }
        return nextStart;
    }

    private static boolean popLastInterval(IntList indices, double[] pos,
            double minCM) {
        if (indices.size()==1) {
            return false;
        }
        else {
            int lastStart = indices.get(indices.size()-1);
            return (pos.length - lastStart) < (MIN_STEP_MARKERS>>1)
                    || (pos[pos.length-1] - pos[lastStart]) < (0.5*minCM);
        }
    }

    private static int[][] ibsSamples(GT targGT, RefGT refGT, IntArray wStarts,
            int w) {
        int wP1 = w + 1;
        int start = wStarts.get(w);
        int end = (wP1 < wStarts.size()) ? wStarts.get(wP1) : targGT.nMarkers();
        int nSamples = targGT.nSamples() + (refGT!=null ? refGT.nSamples() : 0);
        List<SampClust> equivLists = new ArrayList<>(1);
        equivLists.add(new SampClust(nSamples));
        for (int m=start; m<end; ++m) {
            int mm = m;
            equivLists = equivLists.stream()
                    .flatMap(ia -> partition(targGT, refGT, ia, mm))
                    .collect(Collectors.toCollection(ArrayList::new));
        }
        return results(equivLists, targGT.nSamples());
    }

    private static Stream<SampClust> partition(GT targGT, RefGT refGT,
            SampClust parent, int m) {
        // this method assumes int[] parent.samples is sorted in increasing order
        int[] a = new int[2];   // alleles
        int nAlleles = targGT.marker(m).nAlleles();
        int nTargSamples = targGT.nSamples();
        IntList[] gtToList = new IntList[nAlleles*nAlleles];
        boolean[] areHom = new boolean[gtToList.length];
        IntList gtIndices = new IntList(16);
        IntList missingTargSamples = new IntList(32);
        for (int s : parent.samples) {
            setAlleles(m, s, targGT, refGT, a);
            if (a[0]<0 || a[1]<0) {
                assert s<nTargSamples;
                missingTargSamples.add(s);
                for (int k=0, n=gtIndices.size(); k<n; ++k) {
                    gtToList[gtIndices.get(k)].add(s);
                }
            }
            else {
                int gtIndex = a[0]<=a[1] ? a[0]*nAlleles + a[1] : a[1]*nAlleles + a[0];
                if (gtToList[gtIndex]!=null) {
                    gtToList[gtIndex].add(s);
                }
                else if (s<nTargSamples || missingTargSamples.size()>0) { // assumes int[] parent.samples is increasing
                    gtIndices.add(gtIndex);
                    areHom[gtIndex] = parent.areHomozygous && (a[0]==a[1]);
                    gtToList[gtIndex] = new IntList();
                    for (int k=0, n=missingTargSamples.size(); k<n; ++k) {
                        // if s>=nTargSamples, put s in list w/ missTargSamples
                        gtToList[gtIndex].add(missingTargSamples.get(k));
                    }
                    gtToList[gtIndex].add(s);
                }
            }
        }
        return gtIndices.stream()
                .filter(gt -> gtToList[gt].size()>1)
                .mapToObj(gt -> new SampClust(gtToList[gt].toArray(), areHom[gt]));
    }

    private static void setAlleles(int m, int s, GT targGT, RefGT refGT,
            int[] alleles) {
        int nTargSamples = targGT.nSamples();
        boolean inTarg = s<nTargSamples;
        int modS = inTarg ? s : s - nTargSamples;
        alleles[0] = inTarg ? targGT.allele1(m, s) : refGT.allele1(m, modS);
        alleles[1] = inTarg ? targGT.allele2(m, s) : refGT.allele2(m, modS);
    }

    private static int[][] results(List<SampClust> ibd2Lists, int nTargSamples) {
        final int[] EMPTY_ARRAY = new int[0];
        int[][] results = IntStream.range(0, nTargSamples)
                .mapToObj(i -> EMPTY_ARRAY)
                .toArray(int[][]::new);
        for (int j=0, n=ibd2Lists.size(); j<n; ++j) {
            SampClust ibd2List = ibd2Lists.get(j);
            if (ibd2List.areHomozygous==false) {
                int[] ia = ibd2List.samples;
                assert ia.length>1;
                for (int s : ia) {
                    if (s<nTargSamples) {
                        if (results[s]==EMPTY_ARRAY) {
                            results[s] = ia;
                        }
                        else {
                            // sample can be in >1 list due to missing genotypes
                            IntStream is1 = Arrays.stream(results[s]);
                            IntStream is2 = Arrays.stream(ia);
                            results[s] = IntStream.concat(is1, is2)
                                    .sorted()
                                    .distinct()
                                    .toArray();
                        }
                    }
                    else {
                        break;  // int[] ia is an increasing list
                    }
                }
            }
        }
        return results;
    }

    private static SampleSeg[] segList(int s, int nMarkersM1,
            IntArray windowStarts, int[][][] idSets) {
        List<SampleSeg> list = new ArrayList<>();
        for (int w=0; w<idSets.length; ++w) {
            int[] ia = idSets[w][s];
            if (ia.length>0) {
                int wP1 = w + 1;
                int start = windowStarts.get(w);
                int inclEnd = wP1 < windowStarts.size()
                        ? windowStarts.get(wP1) - 1 : nMarkersM1;
                for (int s2 : ia) {
                    if (s2!=s) {
                        list.add(new SampleSeg(s2, start, inclEnd));
                    }
                }
            }
        }
        return list.toArray(new SampleSeg[0]);
    }

    private static SampleSeg[] merge(SampleSeg[] list, DoubleArray genPos,
            double maxGapCm)  {
        if (list.length<2) {
            return list;
        }
        Arrays.sort(list, SampleSeg.sampleComp());
        List<SampleSeg> merged = new ArrayList<>();
        SampleSeg prev = list[0];
        for (int j=1; j<list.length; ++j) {
            SampleSeg next = list[j];
            if (prev.sample()==next.sample()
                    && gapCm(prev, next, genPos) <= maxGapCm) {
                assert prev.inclEnd() < next.inclEnd();
                prev = new SampleSeg(prev.sample(), prev.start(), next.inclEnd());
            }
            else {
                merged.add(prev);
                prev = next;
            }
        }
        merged.add(prev);
        return merged.toArray(new SampleSeg[0]);
    }

    private static double gapCm(SampleSeg first, SampleSeg second,
            DoubleArray genPos) {
        return genPos.get(second.start()) - genPos.get(first.inclEnd());
    }

    private static SampleSeg[] extend(GT targGT, RefGT refGT, int sample,
            SampleSeg[] list) {
        return IntStream.range(0, list.length)
                .mapToObj(j -> extend(targGT, refGT, sample, list, j))
                .toArray(SampleSeg[]::new);
    }

    private static SampleSeg extend(GT targGT, RefGT refGT, int sample,
            SampleSeg[] list, int i) {
        int iP1 = i+1;
        int iM1 = i-1;
        int sample2 = list[i].sample();
        SampleSeg ss = list[i];
        SampleSeg ssM1 = (iM1>=0 && list[iM1].sample()==sample2)
                ? list[iM1] : null;
        SampleSeg ssP1 = (iP1<list.length && list[iP1].sample()==sample2)
                ? list[iP1] : null;

        int inclStart = ss.start();
        int minStart = ssM1==null ? 0 : ssM1.inclEnd();
        int exclEnd = ss.inclEnd() + 1;
        int maxExclEnd = ssP1==null ? targGT.nMarkers() : ssP1.start();
        while (inclStart>minStart
                && ibs2(targGT, refGT, inclStart-1, sample, sample2)) {
            --inclStart;
        }
        while (exclEnd<maxExclEnd
                && ibs2(targGT, refGT, exclEnd, sample, sample2)) {
            ++exclEnd;
        }
        return new SampleSeg(sample2, inclStart, exclEnd - 1);
    }

    private static boolean ibs2(GT targGT, RefGT refGT, int m, int s1, int s2) {
        int nTargSamples = targGT.nSamples();
        boolean inTarg1 = s1<nTargSamples;
        boolean inTarg2 = s2<nTargSamples;
        int modS1 = inTarg1 ? s1 : s1 - nTargSamples;
        int modS2 = inTarg2 ? s2 : s2 - nTargSamples;
        int a1 = inTarg1 ? targGT.allele1(m, s1) : refGT.allele1(m, modS1);
        int a2 = inTarg1 ? targGT.allele2(m, s1) : refGT.allele2(m, modS1);
        int b1 = inTarg2 ? targGT.allele1(m, s2) : refGT.allele1(m, modS2);
        int b2 = inTarg2 ? targGT.allele2(m, s2) : refGT.allele2(m, modS2);
        return arePhaseConsistent(a1, a2, b1, b2)
                || arePhaseConsistent(a1, a2, b2, b1);
    }

    private static boolean arePhaseConsistent(int a1, int a2, int b1, int b2) {
        return (a1<0 || b1<0 || a1==b1) && (a2<0 || b2<0 || a2==b2);
    }

    /**
     * Returns the genotype data for the target samples.
     * @return the genotype data for the target samples
     */
    public GT targGT() {
        return targGT;
    }

    /**
     * Returns the genotype data for the reference samples.  Returns
     * {@code null} if there are no reference samples.
     * @return the genotype data for the reference samples
     */
    public GT refGT() {
        return refGT;
    }

    /**
     * Returns the total number of target and reference samples.
     * @return the total number of target and reference samples
     */
    public int nSamples() {
        return targGT.nSamples() + (refGT!=null ? refGT.nSamples() : 0);
    }

    /**
     * Returns the minimum cM length of a stored IBS2 segment.
     * @return the minimum cM length of a stored IBS2 segment
     */
    public double minCmIbs2() {
        return minCmIbs2;
    }

    /**
     * Returns {@code true} if the specified samples are estimated
     * to be IBS2 at the specified marker, and the IBS2 interval
     * is at least {@code this.minCmIbs2()} cM in length, and
     * returns {@code false} otherwise. Reference sample indices are
     * must be indexed starting with {@code this.targGT().nSamples()}.
     * @param targSample a target sample index
     * @param otherSample a target or reference sample index
     * @param marker a marker index
     * @return {@code true} if the specified samples are estimated
     * to be IBD2 at the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code (targSample < 0 || targSample >= this.targGT().nSamples())}
     * @throws IndexOutOfBoundsException if
     * {@code (otherSample < 0 || otherSample >= this.nSamples())}
     * @throws IndexOutOfBoundsException if
     * {@code (marker < 0 || marker >= this.targGT().nMarkers())}
     */
    public boolean areIbs2(int targSample, int otherSample, int marker) {
        if (targSample>=sampleSegs.length) {
            throw new IndexOutOfBoundsException(String.valueOf(targSample));
        }
        if (targSample==otherSample) {
            return true;
        }
        if (sampleSegs[targSample].length>0) {
            for (SampleSeg ss : sampleSegs[targSample]) {
                if (ss.sample()==otherSample) {
                    if (ss.start() <= marker && marker <= ss.inclEnd()) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static class SampClust {

        private final int[] samples;
        private final boolean areHomozygous;

        private SampClust(int[] samples, boolean areHomozygous) {
            this.samples = samples;
            this.areHomozygous = areHomozygous;
        }

        private SampClust(int nSamples) {
            this.samples = IntStream.range(0, nSamples).toArray();
            this.areHomozygous = true;
        }
    }
}

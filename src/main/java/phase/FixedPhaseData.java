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
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.stream.IntStream;
import main.Par;
import main.Pedigree;
import vcf.Data;
import vcf.GT;
import vcf.GeneticMap;
import vcf.MarkerMap;
import vcf.RefGT;
import vcf.Markers;
import vcf.SplicedGT;

/**
 * <p>Class {@code FixedPhaseData} stores immutable data for a
 * marker window.  The definition of low-frequency markers is
 * determined by the {@code Par.rare()} method and
 * {@code FixedPhaseData.MAX_HIFREQ_PROP} field.</p>
 *
 * <p>Instances of class {@code FixedPhaseData} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FixedPhaseData {

    private static final float MAX_HIFREQ_PROP = 0.9f;
    private static final float MIN_IBS2_CM = 2.0f;

    private final Par par;
    private final float err;
    private final Pedigree ped;
    private final int window;

    private final MarkerMap map;
    private final GT targGT;
    private final RefGT refGT;
    private final int overlap;

    private final MarkerMap hiFreqMap;
    private final GT hiFreqTargGT;
    private final RefGT hiFreqRefGT;
    private final int hiFreqOverlap;
    private final Ibs2 ibs2;

    private final int nHaps;
    private final IntArray[][] carriers;

    private final IntArray hiFreqIndices;
    private final int[] prevHiFreqMarker;
    private final float[] prevWt;   // interpolation weight


    /**
     * Constructs a new {@code FixedPhaseData} instance from the
     * specified data.
     *
     * @param par the analysis parameters
     * @param genMap the genetic map
     * @param data input data for the current marker window
     * @param phasedOverlap initial phased target genotypes due to
     * overlap with the previous window or {@code null} if there are
     * no initial phased target genotypes
     *
     * @throws IllegalArgumentException if
     * {@code (phasedOverlap != null && phasedOverlap.isPhased() == false)}
     * @throws IllegalArgumentException if
     * {@code (phasedOverlap != null
     * && data.targGT().samples().equals(phasedOverlap.samples()) == false)}
     * @throws IllegalArgumentException if
     * {@code (phasedOverlap != null
     * && data.targGT().nMarkers() < phasedOverlap.nMarkers())}
     * @throws IllegalArgumentException if
     * {@code (phasedOverlap != null &&
     * phasedOverlap.marker(j).equals(data.targGT().marker(j) == false)}
     * for some {@code j} satisfying
     * {@code (0 <= j && j <= overlapHaps.nMarkers())}
     * @throws NullPointerException if
     * {@code (par == null || genMap == null || data == null)}
     */
    public FixedPhaseData(Par par, GeneticMap genMap, Data data,
            GT phasedOverlap) {
        checkData(data, phasedOverlap);
        int nTargMarkers = data.targGT().nMarkers();
        this.par = par;
        this.ped = data.ped();
        this.window = data.windowIndex();

        this.map = MarkerMap.create(genMap, data.targGT().markers());
        this.targGT = phasedOverlap==null ? data.targGT() :
                new SplicedGT(phasedOverlap, data.targGT());
        this.refGT = data.restrictRefGT();
        this.overlap = phasedOverlap==null ? 0 : phasedOverlap.nMarkers();
        this.nHaps = nHaps(data);
        this.err = par.err(nHaps);

        IntArray[][] rareCarriers = carriers(par, data);
        int[] hiFreqMkrs = hiFreqIndices(rareCarriers);
        if (hiFreqMkrs.length<2
                || hiFreqMkrs.length > MAX_HIFREQ_PROP*nTargMarkers) {
            hiFreqMkrs = IntStream.range(0, nTargMarkers).toArray();
            ignoreLowFreqCarriers(rareCarriers);
            this.carriers = rareCarriers;
            this.hiFreqMap = map;
            this.hiFreqTargGT = targGT;
            this.hiFreqRefGT = refGT;
            this.hiFreqOverlap = overlap;
            this.hiFreqIndices = new WrappedIntArray(hiFreqMkrs);
            this.prevHiFreqMarker = IntStream.range(0, targGT.nMarkers())
                    .parallel()
                    .toArray();
            float[] fa = new float[targGT.nMarkers()];
            Arrays.fill(fa, 1.0f);
            this.prevWt = fa;
        }
        else {
            Markers hiFreqMarkers = targGT.markers().restrict(hiFreqMkrs);
            this.carriers = rareCarriers;
            this.hiFreqMap = map.restrict(hiFreqMkrs);
            this.hiFreqTargGT = targGT.restrict(hiFreqMarkers, hiFreqMkrs);
            this.hiFreqRefGT = refGT!=null ? refGT.restrict(hiFreqMarkers, hiFreqMkrs) : null;
            this.hiFreqOverlap = hiFreqTargOverlap(phasedOverlap, hiFreqMkrs);
            this.hiFreqIndices = new WrappedIntArray(hiFreqMkrs);
            this.prevHiFreqMarker = prevHiFreqMarker(targGT.nMarkers(), hiFreqIndices);
            this.prevWt = prevWt(map, hiFreqIndices);
        }
        this.ibs2 = new Ibs2(hiFreqTargGT, hiFreqRefGT, hiFreqMap, MIN_IBS2_CM);
    }

    private static void ignoreLowFreqCarriers(IntArray[][] carriers) {
        for (int j=0; j<carriers.length; ++j) {
            Arrays.fill(carriers[j], Data.HIGH_FREQ_ARRAY);
        }
    }

    private static void checkData(Data data, GT phasedOverlap) {
        if (phasedOverlap!=null) {
            GT targ = data.targGT();
            if (phasedOverlap.isPhased()==false) {
                throw new IllegalArgumentException("unphased");
            }
            if (targ.samples().equals(phasedOverlap.samples())==false) {
                throw new IllegalArgumentException("inconsistent data");
            }
            if (phasedOverlap.nMarkers() > targ.nMarkers()) {
                throw new IllegalArgumentException("inconsistent data");
            }
            for (int j=0, n=phasedOverlap.nMarkers(); j<n; ++j) {
                if (phasedOverlap.marker(j).equals(targ.marker(j))==false) {
                    throw new IllegalArgumentException("inconsistent data");
                }
            }
        }
    }

    private static int nHaps(Data data) {
        int nRefHaps = data.refGT()!=null ? data.refGT().nHaps() : 0;
        return data.targGT().nHaps() + nRefHaps;
    }

    private static IntArray[][] carriers(Par par, Data data) {
        RefGT refGT = data.refGT();
        int nRefSamples = refGT!=null ? refGT.nSamples() : 0;
        int nSamples = data.targGT().nSamples() + nRefSamples;
        int maxCarriers = (int) Math.floor(nSamples*par.rare());
        return data.carriers(maxCarriers);
    }

    private static int[] hiFreqIndices(IntArray[][] carriers) {
        return IntStream.range(0, carriers.length)
                .parallel()
                .filter(m -> Arrays.stream(carriers[m])
                        .filter(ia -> ia==Data.HIGH_FREQ_ARRAY)
                        .count()>1)
                .toArray();
    }

    private static int hiFreqTargOverlap(GT phasedOverlap, int[] hiFreqMkrs) {
        if (phasedOverlap==null) {
            return 0;
        }
        int insPt = Arrays.binarySearch(hiFreqMkrs, phasedOverlap.nMarkers());
        return (insPt<0) ? (-insPt - 1) : insPt;
    }

    private int[] prevHiFreqMarker(int nMarkers, IntArray hiFreqIndices) {
        int[] mkrA = new int[nMarkers];
        int nHiFreq = hiFreqIndices.size();
        int start = hiFreqIndices.get(1);
        for (int j=2; j<nHiFreq; ++j) {
            int end = hiFreqIndices.get(j);
            Arrays.fill(mkrA, start, end, j-1);
            start = end;
        }
        Arrays.fill(mkrA, start, nMarkers, nHiFreq-1);
        return mkrA;
    }

    private static float[] prevWt(MarkerMap map, IntArray hiFreqIndices) {
        DoubleArray genPos = map.genPos();
        float[] prevWt = new float[genPos.size()];
        Arrays.fill(prevWt, 0, hiFreqIndices.get(0), 1.0f);
        int start = hiFreqIndices.get(0);
        for (int j=1, n=hiFreqIndices.size(); j<n; ++j) {
            int end = hiFreqIndices.get(j);
            double posA = genPos.get(start);
            double posB = genPos.get(end);
            double d = posB - posA;
            prevWt[start] = 1.0f;
            for (int m=start+1; m<end; ++m) {
                prevWt[m] = (float) ((posB - genPos.get(m))/d);
            }
            start = end;
        }
        Arrays.fill(prevWt, start, genPos.size(), 1.0f);
        return prevWt;
    }

    /**
     * Return the analysis parameters.
     * @return the analysis parameters
     */
    public Par par() {
        return par;
    }

    /**
     * Returns the index of the marker window.
     * @return the index of the marker window
     */
    public int window() {
        return window;
    }

    /**
     * Returns the parent-offspring relationships.
     * @return the parent-offspring relationships
     */
    public Pedigree ped() {
        return ped;
    }

    /**
     * Returns the genetic map for the markers.
     * @return the genetic map for the markers
     */
    public MarkerMap map() {
        return map;
    }

    /**
     * Returns the input phased, nonmissing reference genotypes
     * or {@code null} if there are no reference samples.
     * @return the input phased, nonmissing reference genotypes
     * or {@code null} if there are no reference samples
     */
    public RefGT refGT() {
        return refGT;
    }

    /**
     * Returns the input target genotypes.
     * @return the input target genotypes
     */
    public GT targGT() {
        return targGT;
    }

    /**
     * Returns the number of initial markers that have phased target
     * genotypes due to overlap with the previous marker window.
     * @return the number of initial markers that have phased target
     * genotypes due to overlap with the previous marker window
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the genetic map for the high-frequency markers.
     * @return the genetic map for the high-frequency markers
     */
    public MarkerMap hiFreqMap() {
        return hiFreqMap;
    }

    /**
     * Returns the phased, nonmissing reference genotypes for the high-frequency
     * markers or {@code null} if there are no reference samples.
     * @return the phased, nonmissing reference genotypes for the high-frequency
     * markers or {@code null} if there are no reference samples
     */
    public RefGT hiFreqRefGT() {
        return hiFreqRefGT;
    }

    /**
     * Returns the input target genotypes at the high-frequency markers.
     * @return the input target genotypes at the high-frequency markers
     */
    public GT hiFreqTargGT() {
        return hiFreqTargGT;
    }

    /**
     * Returns the number of initial high-frequency markers that
     * have phased target genotypes due to overlap with the previous window.
     * @return the number of initial high-frequency markers that
     * have phased target genotypes due to overlap with the previous window
     */
    public int hiFreqOverlap() {
        return hiFreqOverlap;
    }

    /**
     * Return the sum of the number of reference and target haplotypes.
     * @return the sum of the number of reference and target haplotypes
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Return the allele mismatch emission probability.
     * @return the allele mismatch emission probability
     */
    public float err() {
        return err;
    }

    /**
     * Returns the indices of the high-frequency markers.
     * @return the indices of the high-frequency markers
     */
    public IntArray hiFreqIndices() {
        return hiFreqIndices;
    }

    /**
     * Returns the IBS2 data for the high-frequency markers.
     * @return the IBS2 data for the high-frequency markers
     */
    public Ibs2 hiFreqIbs2() {
        return ibs2;
    }

    /**
     * Returns the indices of the reference and target samples for the
     * specified low-frequency allele.  The reference sample indices will be
     * shifted by the number of target samples. so that the first reference
     * sample will have an index equal to the number of target samples.
     * The returned list will be sorted in order of increasing sample index.
     * The returned array will be empty and equal to
     * {@code vcf.Data.ZERO_FREQ_ARRAY} if the allele has no carriers, and the
     * returned array will be empty and equal to
     * {@code vcf.Data.HIGH_FREQ_ARRAY} if the allele is not a low-frequency
     * allele.
     * @param marker a marker index
     * @param allele an allele index for the specified marker
     * @return the indices of the reference and target samples that the
     * specified low-frequency allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.targGT().marker(marker).nAlleles()}
     */
    public IntArray carriers(int marker, int allele) {
        return carriers[marker][allele];
    }

    /**
     * Returns {@code true} if the specified allele is a low-frequency allele,
     * and returns {@code false} otherwise.
     * @param marker a marker index
     * @param allele an allele index for the specified marker
     * @return {@code true} if the specified allele is a low-frequency allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.targGT().marker(marker).nAlleles()}
     */
    public boolean isLowFreq(int marker, int allele) {
        return carriers[marker][allele]!=Data.HIGH_FREQ_ARRAY;
    }

    /**
     * Returns the index of the closest high-frequency marker with position
     * less than or equal to the position of the specified marker, or
     * 0 if no such high-frequency marker exists.
     * @param marker a marker index
     * @return the index of the closest hi-frequency marker with position
     * less than or equal to the position of the specified marker, or
     * 0 if no such high-frequency marker exists
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.targGT().nMarkers()}
     */
    public int prevHiFreqMarker(int marker) {
        return prevHiFreqMarker[marker];
    }

    /**
     * Returns the linear interpolation weight for the high-frequency marker
     * with index {@code this.prevHiFreqMarker(marker)}.
     * @param marker a marker index
     * @return the linear interpolation weight for the high-frequency marker
     * with index {@code this.prevHiFreqMarker(marker)}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.targGT().nMarkers()}
     */
    public float prevWt(int marker) {
        return prevWt[marker];
    }
}

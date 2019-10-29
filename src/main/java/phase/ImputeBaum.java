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

import blbutil.FloatList;
import ints.IntArray;
import ints.IntList;
import java.util.Random;
import vcf.GT;
import vcf.RefGT;

/**
 * <p>Class {@code ImputeBaum} applies the forward and backward algorithms
 * for a haploid Li and Stephens hidden Markov model at high-frequency markers,
 * and imputes missing genotypes and heterozygote phase at low-frequency
 * markers.</p>
 *
 * <p>Instances of class {@code ImputeBaum} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImputeBaum {

    private final FixedPhaseData fpd;
    private final PhaseData phaseData;
    private final FwdBwd fwdBwd;
    private final int[] targHap = new int[2];
    private final int[] nStates = new int[2];
    private final int[][][] states;
    private final float[][][] probs;

    private final IntList savedStates = new IntList(8);
    private final FloatList savedProbs = new FloatList(8);

    private final GT hiFreqPhasedGT;
    private final GT unphTargGT;
    private final RefGT refGT;
    private final int nTargHaps;
    private final int nHaps;
    private final int nHiFreqMarkers;
    private final HapImputer imputableHaps;
    private final IntArray hiFreqIndices;
    private final Random rand;

    private final int[][] outPhase;

    /**
     * Creates a {@code ImputeBaum} instance from the specified data.
     *
     * @param phaseIbs the IBS haplotypes
     * @param hapImputer an object for imputing haplotypes at low-frequency
     * markers
     * @throws NullPointerException if
     * {@code phaseIbs == null || hapImputer == null}
     */
    public ImputeBaum(PhaseIbs phaseIbs, HapImputer hapImputer) {
        this.fpd = phaseIbs.phaseData().fpd();
        this.phaseData = phaseIbs.phaseData();
        this.nHiFreqMarkers = phaseData.targGT().nMarkers();
        this.fwdBwd =  new FwdBwd(phaseIbs);
        this.states = new int[2][nHiFreqMarkers][fwdBwd.maxStates()];
        this.probs = new float[2][nHiFreqMarkers][fwdBwd.maxStates()];

        this.hiFreqPhasedGT = phaseData.estPhase().hapsGT();
        this.unphTargGT = fpd.targGT();
        this.refGT = fpd.refGT();
        this.nTargHaps = fpd.targGT().nHaps();
        this.nHaps = fpd.nHaps();
        this.imputableHaps = hapImputer;
        this.hiFreqIndices = fpd.hiFreqIndices();
        this.rand = new Random(phaseData.seed());
        this.outPhase = new int[2][unphTargGT.nMarkers()];
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return hiFreqPhasedGT.nSamples();
    }

    /**
     * Estimates phased haplotypes for the specified sample.
     * @param sample a sample index
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nTargSamples()}
     */
    public void phase(int sample) {
        rand.setSeed(phaseData.seed() + sample);
        for (int i=0; i<2; ++i) {
            targHap[i] = (sample << 1) + i;
            nStates[i] = fwdBwd.run(targHap[i], states[i], probs[i]);
        }
        int start = 0;
        for (int j=0; j<nHiFreqMarkers; ++j) {
            int end = hiFreqIndices.get(j);
            imputeInterval(sample, start, end);
            outPhase[0][end] = hiFreqPhasedGT.allele1(j, sample);
            outPhase[1][end] = hiFreqPhasedGT.allele2(j, sample);
            start = end + 1;
        }
        imputeInterval(sample, start, unphTargGT.nMarkers());
        imputableHaps.setHap(targHap[0], outPhase[0]);
        imputableHaps.setHap(targHap[1], outPhase[1]);
    }

    private void imputeInterval(int sample, int start, int end) {
        for (int m=start; m<end; ++m) {
            int a1 = unphTargGT.allele1(m, sample);
            int a2 = unphTargGT.allele2(m, sample);
            if (a1>=0 && a2>=0) {
                boolean noFlip = true;
                if (a1!=a2) {
                    float[] alProbs1 = unscaledAlProbs(m, 0, a1, a2);
                    float[] alProbs2 = unscaledAlProbs(m, 1, a1, a2);
                    float p1 = alProbs1[a1]*alProbs2[a2];
                    float p2 = alProbs1[a2]*alProbs2[a1];
                    noFlip = p1>p2 || (p1==0.0f && p2==0.0f && rand.nextBoolean());
                }
                outPhase[0][m] = (noFlip) ? a1 : a2;
                outPhase[1][m] = (noFlip) ? a2 : a1;
            }
            else {
                outPhase[0][m] = imputeAllele(m, 0);
                outPhase[1][m] = imputeAllele(m, 1);
            }
        }
    }

    private float[] unscaledAlProbs(int m, int hapBit, int a1, int a2) {
        float[] alProbs = new float[unphTargGT.marker(m).nAlleles()];
        boolean rare1 = fpd.isLowFreq(m, a1);
        boolean rare2 = fpd.isLowFreq(m, a2);
        int mkrA = fpd.prevHiFreqMarker(m);
        int mkrB = Math.min(mkrA + 1, nHiFreqMarkers - 1);
        int[] statesA = states[hapBit][mkrA];
        float[] probsA = probs[hapBit][mkrA];
        float[] probsB = probs[hapBit][mkrB];
        for (int j=0, n=nStates[hapBit]; j<n; ++j) {
            int hap = statesA[j];
            int b1 = allele(m, hap);
            int b2 = allele(m, (hap ^ 0b1));
            if (b1>=0 && b2>=0) {
                float wt = fpd.prevWt(m);
                float prob = wt*probsA[j] + (1.0f - wt)*probsB[j];
                if (b1==b2) {
                    alProbs[b1] += prob;
                }
                else {
                    boolean match1 = rare1 && (a1==b1 || a1==b2);
                    boolean match2 = rare2 && (a2==b1 || a2==b2);
                    if (match1 ^ match2) {
                        if (match1) {
                            alProbs[a1] += prob;
                        }
                        else {
                            alProbs[a2] += prob;
                        }
                    }
                }
            }
        }
        return alProbs;
    }

    private int imputeAllele(int m, int hapBit) {
        savedStates.clear();
        savedProbs.clear();
        float[] alProbs = new float[unphTargGT.marker(m).nAlleles()];
        float unknownAlProb = 0.0f;
        int mkrA = fpd.prevHiFreqMarker(m);
        int mkrB = Math.min(mkrA + 1, nHiFreqMarkers - 1);
        int[] statesA = states[hapBit][mkrA];
        float[] probsA = probs[hapBit][mkrA];
        float[] probsB = probs[hapBit][mkrB];
        for (int j=0, n=nStates[hapBit]; j<n; ++j) {
            float wt = fpd.prevWt(m);
            float prob = wt*probsA[j] + (1.0f - wt)*probsB[j];
            int hap = statesA[j];
            int b1 = allele(m, hap);
            int b2 = allele(m, hap ^ 0b1);
            if (b1>=0 && b2>=0) {
                if (b1==b2) {
                    alProbs[b1] += prob;
                }
                else {
                    unknownAlProb += prob;
                    savedStates.add(hap);
                    savedProbs.add(prob);
                }
            }
        }
        int imputedAllele = maxIndex(alProbs);
        if (alProbs[imputedAllele] < unknownAlProb) {
            imputableHaps.setPartlyImputedAllele(targHap[hapBit], m, alProbs,
                    savedStates, savedProbs);
        }
        return imputedAllele;
    }

    private int allele(int marker, int hap) {
        if (hap>=nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        if (hap < nTargHaps) {
            return unphTargGT.allele(marker, hap);
        }
        else {
            return refGT.allele(marker, hap - nTargHaps);
        }
    }

    private int maxIndex(float[] fa) {
        int maxIndex = 0;
        for (int j=1; j<fa.length; ++j) {
            if (fa[j] > fa[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }
}

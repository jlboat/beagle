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

import blbutil.FloatArray;
import blbutil.FloatList;
import ints.IntArray;
import ints.IntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import vcf.Markers;

/**
 * <p>Class {@code PhaseBaum1} implements the forward and backward algorithms
 * for a haploid Li and Stephens hidden Markov model.
 * </p>
 * <p>Instances of class {@code PhaseBaum1} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseBaum1 {

    private final PhaseData phaseData;
    private final EstPhase estPhase;
    private final Markers markers;
    private final FloatArray pRecomb;
    private final List<int[]> refAl;
    private final byte[][][] alMatch;

    private final PhaseStates states;
    private final FloatList lrList;
    private final int nMarkers;

    private final int maxStates;
    private int nStates;
    private final int[] hap1;
    private final int[] hap2;
    private final float[][] fwd;
    private final float[][] bwd;
    private final List<float[]> savedBwd1;
    private final List<float[]> savedBwd2;

    private final float[] sum;

    private final byte[] missingMatch;
    private final HmmUpdater hmmUpdater;

    private boolean swapHaps = false;
    private int savedIndex = 0;
    private int missIndex = 0;

    /**
     * Creates a {@code PhaseLSBaum} instance from the specified data.
     *
     * @param phaseIbs the IBS haplotype segments

     * @throws NullPointerException if {@code phaseIBS == null}
     */
    public PhaseBaum1(PhaseIbs phaseIbs) {
        this.phaseData = phaseIbs.phaseData();
        this.estPhase = phaseData.estPhase();
        this.markers = phaseData.targGT().markers();
        this.pRecomb = phaseData.pRecomb();
        this.lrList = new FloatList(200);
        this.nMarkers = phaseData.targGT().nMarkers();

        this.maxStates = phaseData.par().phase_states();
        this.states = new PhaseStates(phaseIbs, maxStates);
        this.refAl = new ArrayList<>();
        this.alMatch = new byte[3][nMarkers][maxStates];

        this.hap1 = new int[nMarkers];
        this.hap2 = new int[nMarkers];
        this.savedBwd1 = new ArrayList<>();
        this.savedBwd2 = new ArrayList<>();

        this.fwd = new float[5][maxStates];
        this.bwd = new float[3][maxStates];
        this.sum = new float[5];

        this.hmmUpdater = new HmmUpdater(phaseData.err());
        this.missingMatch = new byte[maxStates];
        Arrays.fill(missingMatch, (byte) 0);
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargSamples() {
        return phaseData.targGT().nSamples();
    }

    /**
     * Estimates and stores the phased haplotypes for the specified sample
     * @param sample a sample index
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nTargSamples()}
     */
    public void phase(int sample) {
        long t0 = System.nanoTime();        // xxx
        savedIndex = 0;
        missIndex = 0;
        swapHaps = false;
        boolean missingGTs = estPhase.missing(sample).size()>0;
        boolean hasUnphasedHet = estPhase.unphased(sample).size()>0;
        if (missingGTs || hasUnphasedHet) {
            setHapsWithMissingAlleles(sample);
            this.nStates = states.ibsStates(sample, hap1, hap2, refAl,
                    alMatch[1], alMatch[2]);
            phaseAndImpute(estPhase, sample);
            estPhase.setHapPair(sample, hap1, hap2);
        }
    }

    private void setHapsWithMissingAlleles(int sample) {
        estPhase.getHaps(sample, hap1, hap2);
        IntArray missing = estPhase.missing(sample);
        int nMissing = missing.size();
        int nSaved = nMissing + estPhase.unphased(sample).size();
        for (int j=0; j<nMissing; ++j) {
            int m = missing.get(j);
            hap1[m] = hap2[m] = -1;
        }
        while (refAl.size()<nMissing) {
            refAl.add(new int[maxStates]);
        }
        while (savedBwd1.size()<nSaved) {
            savedBwd1.add(new float[maxStates]);
            savedBwd2.add(new float[maxStates]);
        }
    }

    private void phaseAndImpute(EstPhase estPhase, int sample) {
        lrList.clear();
        IntArray unph = estPhase.unphased(sample);
        IntArray miss = estPhase.missing(sample);
        if (unph.size()>0 || miss.size()>0) {
            bwdAlg(unph);
            fwdAlg(unph);
            if (unph.size()>0) {
                updateUnphased(sample);
            }
        }
    }

    private void fwdAlg(IntArray unph) {
        Arrays.fill(fwd[0], 0, nStates, 1.0f/nStates);
        Arrays.fill(fwd[3], 0, nStates, 1.0f/nStates);
        Arrays.fill(fwd[4], 0, nStates, 1.0f/nStates);
        sum[0] = sum[3] = sum[4] = 1.0f;
        int start = 0;
        for (int j=0, n=unph.size(); j<n; ++j) {
            int end = unph.get(j);
            setFwd(start, end);
            phaseHet();
            start = end;
        }
        setFwd(start, nMarkers);
    }

    private void setFwd(int start, int end) {
        long t0 = System.nanoTime();
        if (swapHaps) {
            swapHaps(start, end);
        }
        System.arraycopy(fwd[0], 0, fwd[1], 0, nStates);
        System.arraycopy(fwd[0], 0, fwd[2], 0, nStates);
        sum[1] = sum[2] = sum[0];
        for (int m=start; m<end; ++m)  {
            float pRec = pRecomb.get(m);
            boolean isMissing = (hap1[m]<0 || hap2[m]<0);
            byte[] alMatch0 = (isMissing || hap1[m]!=hap2[m]) ? missingMatch : alMatch[1][m];
            sum[0] = hmmUpdater.fwdUpdate(fwd[0], sum[0], pRec, alMatch0, nStates);
            sum[1] = hmmUpdater.fwdUpdate(fwd[1], sum[1], pRec, alMatch[1][m], nStates);
            sum[2] = hmmUpdater.fwdUpdate(fwd[2], sum[2], pRec, alMatch[2][m], nStates);
            sum[3] = hmmUpdater.fwdUpdate(fwd[3], sum[3], pRec, alMatch[1][m], nStates);
            sum[4] = hmmUpdater.fwdUpdate(fwd[4], sum[4], pRec, alMatch[2][m], nStates);
            if (isMissing) {
                int nAl = markers.marker(m).nAlleles();
                float[] bwd1 = savedBwd1.get(--savedIndex);
                float[] bwd2 = savedBwd2.get(savedIndex);
                int[] refAlleles = refAl.get(missIndex++);
                hap1[m] = imputeAl(nAl, nStates, refAlleles, fwd[3],
                        (swapHaps ? bwd2 : bwd1));
                hap2[m] = imputeAl(nAl, nStates, refAlleles, fwd[4],
                        (swapHaps ? bwd1 : bwd2));
            }
        }
    }

    private static int imputeAl(int nAlleles, int nStates, int[] refAl,
            float[] fwd, float[] bwd) {
        float[] alFreq = new float[nAlleles];
        for (int k=0; k<nStates; ++k) {
            alFreq[refAl[k]] += fwd[k]*bwd[k];
        }
        int maxIndex = 0;
        for (int j=1; j<alFreq.length; ++j) {
            if (alFreq[j] > alFreq[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    private void bwdAlg(IntArray unph) {
        long t0 = System.nanoTime();
        int end = nMarkers - 1;
        Arrays.fill(bwd[0], 0, nStates, 1.0f/nStates);
        if (hap1[end]<0 || hap2[end]<0) {
            System.arraycopy(bwd[0], 0, savedBwd1.get(savedIndex), 0, nStates);
            System.arraycopy(bwd[0], 0, savedBwd2.get(savedIndex++), 0, nStates);
        }
        for (int j=unph.size()-1; j>=0; --j) {
            int start = unph.get(j) - 1;
            runBwd(start, end);
            System.arraycopy(bwd[1], 0, savedBwd1.get(savedIndex), 0, nStates);
            System.arraycopy(bwd[2], 0, savedBwd2.get(savedIndex++), 0, nStates);
            end = start;
        }
        runBwd(0, end);
    }

    private void runBwd(int start, int end) {
        System.arraycopy(bwd[0], 0, bwd[1], 0, nStates);
        System.arraycopy(bwd[0], 0, bwd[2], 0, nStates);
        for (int m=end-1; m>=start; --m) {
            int mP1 = m + 1;
            float pRec = pRecomb.get(mP1);
            boolean emit1 = hap1[mP1]<0 || hap2[mP1]<0 || hap1[mP1]!=hap2[mP1];
            byte[] alMatch0 = emit1 ? missingMatch : alMatch[1][mP1];
            hmmUpdater.bwdUpdate(bwd[0], pRec, alMatch0, nStates);
            hmmUpdater.bwdUpdate(bwd[1], pRec, alMatch[1][mP1], nStates);
            hmmUpdater.bwdUpdate(bwd[2], pRec, alMatch[2][mP1], nStates);
            if (hap1[m]<0 || hap2[m]<0) {
                System.arraycopy(bwd[1], 0, savedBwd1.get(savedIndex), 0, nStates);
                System.arraycopy(bwd[2], 0, savedBwd2.get(savedIndex++), 0, nStates);
            }
        }
    }

    private void phaseHet() {
        float[] b1 = savedBwd1.get(--savedIndex);
        float[] b2 = savedBwd2.get(savedIndex);
        float p11 = 0.0f;
        float p12 = 0.0f;
        float p21 = 0.0f;
        float p22 = 0.0f;
        for (int k=0; k<nStates; ++k) {
            p11 += fwd[1][k]*b1[k];
            p12 += fwd[1][k]*b2[k];
            p21 += fwd[2][k]*b1[k];
            p22 += fwd[2][k]*b2[k];
        }
        float num = (p11*p22);
        float den = (p12*p21);
        swapHaps = num < den;
        lrList.add(swapHaps ? den/num : num/den);
    }

    private void swapHaps(int m1, int m2) {
        for (int m=m1; m<m2; ++m) {
            int tmp = hap1[m];
            hap1[m] = hap2[m];
            hap2[m] = tmp;

            byte[] tmpMatch = alMatch[1][m];
            alMatch[1][m] = alMatch[2][m];
            alMatch[2][m] = tmpMatch;
        }
    }

    private void updateUnphased(int sample) {
        double leaveUnphasedProp = phaseData.leaveUnphasedProp(sample);
        if (leaveUnphasedProp < 1.0) {
            IntArray prevUnph = estPhase.unphased(sample);
            IntList nextUnph = new IntList();
            float threshold = threshold(lrList, leaveUnphasedProp);
            for (int j=0, n=prevUnph.size(); j<n; ++j) {
                if (lrList.get(j) < threshold) {
                    nextUnph.add(prevUnph.get(j));
                }
            }
            estPhase.setUnphased(sample, IntArray.create(nextUnph, nMarkers));
        }
    }

    private static float threshold(FloatList lrList, double propToLeaveUnphased) {
        float[] lra = lrList.toArray();
        Arrays.sort(lra);
        int rank = (int) Math.floor(propToLeaveUnphased*lra.length + 0.5);
        return lra[rank<lra.length ? rank : (lra.length-1)];
    }
}

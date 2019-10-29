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
import java.util.Arrays;

/**
 * <p>Class {@code FwdBwd} implements computes HMM reference haplotypes
 * and state probabilities.</p>
 *
 * <p>Instances of class {@code FwdBwd} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FwdBwd {

    private final PhaseIbs phaseIbs;
    private final int maxStates;
    private final PhaseStates states;
    private final FloatArray pRecomb;
    private final byte[][] nMismatches;
    private final float[] bwd;
    private final float[] emProb;

    /**
     * Creates a {@code FwdBwd} instance from the specified data.
     * @param phaseIbs the IBS haplotypes
     * @throws NullPointerException if {@code phaseIbs == null}
     */
    public FwdBwd(PhaseIbs phaseIbs) {
        PhaseData phaseData = phaseIbs.phaseData();
        this.phaseIbs = phaseIbs;
        int nMarkers = phaseData.targGT().nMarkers();
        this.maxStates = phaseData.par().phase_states()/2;
        this.states = new PhaseStates(phaseIbs, maxStates);
        this.pRecomb = phaseData.pRecomb();
        this.nMismatches = new byte[nMarkers][maxStates];
        this.bwd = new float[maxStates];
        float pErr = phaseData.err();
        float pNoErr = 1.0f - pErr;
        this.emProb = new float[] {pNoErr, pErr};
    }

    /**
     * Stores the HMM reference haplotypes and states probabilities for the
     * specified target haplotype, and returns the number of HMM states
     * per marker.  The contract for this method is undefined
     * if the number of elements in each row of the specified arrays is not
     * greater than or equal to {@code this.maxStates()}.
     *
     * @param hap a target haplotype index
     * @param refHaps the array in which the reference haplotypes for each
     * hidden state will be stored
     * @param stateProbs the array in which estimated probabilities for each
     * hidden state will be stored
     * @return the number of hidden states at each marker
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.phaseIbs().phaseData().targGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code refHaps.length < this.phaseIbs().phaseData.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code stateProbs.length < this.phaseIbs().phaseData.targGT().nMarkers()}
     * @throws NullPointerException if
     * {@code refHaps == null || stateProbs == null}
     */
    public int run(int hap, int[][] refHaps, float[][] stateProbs) {
        int nStates = states.ibsStates(hap, refHaps, nMismatches);
        runFwd(stateProbs, nStates);
        runBwd(stateProbs, nStates);
        return nStates;
    }

    private void runFwd(float[][] stateProbs, int nStates) {
        float lastSum = 0.0f;
        for (int j=0; j<nStates; ++j) {
            stateProbs[0][j] = emProb[nMismatches[0][j]];
            lastSum += stateProbs[0][j];
        }
        for (int m=1; m<stateProbs.length; ++m) {
            int mM1 = m - 1;
            float pRec = pRecomb.get(m);
            float shift = pRec/nStates;
            float scale = (1.0f - pRec)/lastSum;
            lastSum = 0.0f;
            for (int j=0; j<nStates; ++j) {
                float em = emProb[nMismatches[m][j]];
                stateProbs[m][j] = em*(scale*stateProbs[mM1][j] + shift);
                lastSum += stateProbs[m][j];
            }
        }
    }

    private void runBwd(float[][] stateProbs, int nStates) {
        int inclEnd = stateProbs.length - 1;
        Arrays.fill(bwd, 0, nStates, 1.0f/nStates);
        for (int m=inclEnd-1; m>=0; --m) {
            int mP1 = m + 1;
            float sum = 0.0f;
            for (int j=0; j<nStates; ++j) {
                bwd[j] *= emProb[nMismatches[mP1][j]];
                sum += bwd[j];
            }
            float pRec = pRecomb.get(mP1);
            float scale = (1.0f - pRec)/sum;
            float shift = pRec/nStates;
            sum = 0.0f;
            for (int j=0; j<nStates; ++j) {
                bwd[j] = scale*bwd[j] + shift;
                stateProbs[m][j] *= bwd[j];
                sum += stateProbs[m][j];
            }
            for (int j=0; j<nStates; ++j) {
                stateProbs[m][j] /= sum;
            }
        }
    }

    /**
     * Return the haplotype IBS
     * @return the haplotype IBS
     */
    public PhaseIbs phaseIbs() {
        return phaseIbs;
    }

    /**
     * Returns the maximum number of HMM states.
     * @return the maximum number of HMM states
     */
    public int maxStates() {
        return maxStates;
    }
}
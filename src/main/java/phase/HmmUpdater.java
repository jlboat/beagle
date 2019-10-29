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

/**
 * <p>Instances of class {@code HmmUpdater} have methods for one-step
 * updates of forward or backward HMM values.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HmmUpdater {

    private final float[] pMismatch;

    /**
     * Constructs a new {@code HmmUpdater} instance from the specified
     * data.
     * @param pErr the probability that the observed allele does not
     * match the state allele
     */
    public HmmUpdater(float pErr) {
        this.pMismatch = new float[] {1.0f - pErr, pErr};
    }

    /**
     * Updates the forward values and returns the sum of the updated forward
     * values.
     * @param fwd the array of forward values that will be updated
     * @param pSwitch the probability of jumping to a random state
     * @param sum the sum of forward values in the specified array
     * @param mismatch the number of allele mismatches (0 or 1) for
     * each HMM state
     * @param nStates the number of states
     * @return the sum of the updated forward values
     * @throws IndexOutOfBoundsException if
     * {@code fwd.length < nStates || mismatch.length < nStates}
     * @throws NullPointerException if
     * {@code fwd == null || mismatch == null}
     */
    public float fwdUpdate(float[] fwd, float sum, float pSwitch,
            byte[] mismatch, int nStates) {
        float shift = pSwitch/nStates;
        float scale = (1.0f - pSwitch)/sum;
        sum = 0.0f;
        for (int k=0; k<nStates; ++k) {
            fwd[k] = pMismatch[mismatch[k]]*(scale*fwd[k] + shift);
            sum += fwd[k];
        }
        return sum;
    }

    /**
     * Updates the backward values.
     * @param bwd the array of backward values that will be updated
     * @param pSwitch the probability of jumping to a random state
     * @param mismatch the number of allele mismatches (0 or 1) for
     * each HMM state
     * @param nStates the number of states
     * @throws IndexOutOfBoundsException if
     * {@code bwd.length < nStates || mismatch.length < nStates}
     * @throws NullPointerException if
     * {@code bwd == null || mismatch == null}
     */
    public void bwdUpdate(float[] bwd, float pSwitch, byte[] mismatch,
            int nStates) {
        float sum = 0.0f;
        for (int k=0; k<nStates; ++k) {
            bwd[k] *= pMismatch[mismatch[k]];
            sum += bwd[k];
        }
        float shift = pSwitch/nStates;
        float scale = (1.0f - pSwitch)/sum;
        for (int k=0; k<nStates; ++k) {
            bwd[k] = scale*bwd[k] + shift;
        }
    }
}

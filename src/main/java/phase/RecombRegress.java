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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import math.Regress;

/**
 * <p>Class {@code RecombRegress} uses linear regression to estimate the
 * recombination factor for a haploid Li and Stephens hidden Markov model.</p>
 *
 * <p>Instances of class {@code RecombRegress} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RecombRegress {

    private final PhaseData phaseData;
    private final EstPhase estPhase;
    private final int nMarkers;
    private final FloatArray genDist;
    private final FloatArray pRecomb;
    private final PhaseStates states;
    private final List<int[]> refAl;
    private final byte[][] alMatch1;
    private final byte[][] alMatch2;

    private final int[] hap1;
    private final int[] hap2;
    private final float[] fwd;
    private final float[] bwd;
    private final float[] fwdM1;
    private final float[][] savedBwd;
    private final float[] pMismatch;
    private final HmmUpdater hmmUpdater;
    private final Regress regress;

    private double sumY = 0.0f;
    private int nStates = -1;

    /**
     * Creates a {@code RecombRegress} instance for the specified data.
     *
     * @param phaseIbs the IBS haplotype segments
     * @param regress the object for storing data for regressing the
     * probability of switching reference haplotypes between consecutive
     * markers on genetic distance
     * @throws NullPointerException if
     * {@code (phaseIBS == null || regress == null)}
     */
    public RecombRegress(PhaseIbs phaseIbs, Regress regress) {
        if (regress==null) {
            throw new NullPointerException(Regress.class.toString());
        }
        this.phaseData = phaseIbs.phaseData();
        this.estPhase = phaseData.estPhase();
        this.nMarkers = phaseData.targGT().nMarkers();
        this.genDist = phaseData.map().genDist();
        this.pRecomb = phaseData.pRecomb();
        int maxStates = phaseData.par().phase_states();
        this.states = new PhaseStates(phaseIbs, maxStates);
        this.refAl = new ArrayList<>(0);
        this.alMatch1 = new byte[nMarkers][maxStates];
        this.alMatch2 = new byte[nMarkers][maxStates];

        this.hap1 = new int[nMarkers];
        this.hap2 = new int[nMarkers];
        this.savedBwd = new float[nMarkers][maxStates];
        this.fwd = new float[maxStates];
        this.bwd = new float[maxStates];
        this.fwdM1 = new float[maxStates];
        float pErr = phaseData.err();
        this.pMismatch = new float[] {1.0f - pErr, pErr};
        this.hmmUpdater = new HmmUpdater(pErr);
        this.regress = regress;
    }

    /**
     * Returns the current input data for genotype phasing.
     * @return the current input data for genotype phasing
     */
    public PhaseData phaseData() {
        return phaseData;
    }

    /**
     * Stores regression data for regressing the probability of switching
     * reference haplotypes between consecutive markers on genetic distance
     * for the specified sample.
     * @param sample a sample index
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= phaseData().targGT().nSamples()}
     */
    public void update(int sample) {
        estPhase.getHaps(sample, hap1, hap2);
        this.nStates = states.ibsStates(sample, hap1, hap2, refAl,
                alMatch1, alMatch2);
        updateRegress(alMatch1);
        updateRegress(alMatch2);
    }

    private void updateRegress(byte[][] alMatch) {
        Arrays.fill(bwd, 0, nStates, 1.0f);
        Arrays.fill(savedBwd[nMarkers-1], 0, nStates, 1.0f);
        for (int m=nMarkers-2; m>=0; --m) {
            int mP1 = m + 1;
            hmmUpdater.bwdUpdate(bwd, pRecomb.get(mP1), alMatch[mP1], nStates);
            System.arraycopy(bwd, 0, savedBwd[m], 0, nStates);
        }
        float hFactor = (float) (nStates / (nStates - 1.0));
        Arrays.fill(fwd, 0, nStates, 1.0f/nStates);
        float sum = 1.0f;
        sum = hmmUpdater.fwdUpdate(fwd, sum, pRecomb.get(0), alMatch[0], nStates);
        for (int m=1; m<nMarkers; ++m) {
            sum = fwdUpdate(m, alMatch[m], sum, hFactor);
        }
    }

    private float fwdUpdate(int m, byte[] alMatch, float lastSum, float hFactor) {
        float pSwitch = pRecomb.get(m);
        float f = ((1.0f - pSwitch) + (pSwitch/nStates))/lastSum;
        float partNumer = 0.0f;
        float denom = 0.0f;

        float[] storedBwd = savedBwd[m];
        System.arraycopy(fwd, 0, fwdM1, 0, nStates);
        float sum = hmmUpdater.fwdUpdate(fwd, lastSum, pRecomb.get(m), alMatch, nStates);
        for (int k=0; k<nStates; ++k) {
            partNumer += pMismatch[alMatch[k]]*f*fwdM1[k]*storedBwd[k];
            denom += fwd[k]*storedBwd[k];
        }
        float num = denom - partNumer;
        float xVal = genDist.get(m);
        float yVal = (hFactor*num/denom);
        regress.add(xVal, yVal);
        sumY += yVal;
        return sum;
    }

    /**
     * Returns the sum of the probabilities of switching reference
     * haplotypes between consecutive markers stored by {@code this.update()}.
     * @return the sum of the probabilities of switching reference haplotypes
     * between consecutive markers
     */
    public double sumY() {
        return sumY;
    }
}

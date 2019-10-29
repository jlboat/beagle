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
package main;

import phase.EstPhase;
import phase.PhaseLS;
import java.util.Random;
import phase.FixedPhaseData;
import phase.PhaseData;
import vcf.GT;
import vcf. RefGT;

/**
 * <p>Class {@code MainHelper} is an auxiliary class with methods called by
 * the {@code main.Main} class.
 * </p>
 * <p>Instances of class {@code MainHelper} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MainHelper {

    private final Par par;
    private final RunStats runStats;
    private final Random rand;

    private float recombFactor;

    /**
     * Constructs a new {@code MainHelper} instance.
     * @param par the command line parameters
     * @param runStats the class for collecting and printing run-time statistics
     * @param seed random number seed
     * @throws NullPointerException
     * if {@code runStarts == null}
     */
    public MainHelper(Par par, RunStats runStats, long seed) {
        this.par = par;
        this.runStats = runStats;
        this.rand = new Random(seed);
    }

    /**
     * Phases the current window of genotype data.
     * @param fpd the input data for the current window
     * @return the phased genotype data
     * @throws IllegalArgumentException if {@code fpd.par()!=this.par()}
     * @throws NullPointerException if {@code fpd == null}
     */
    public GT phase(FixedPhaseData fpd) {
        if (fpd.par()!=par) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (fpd.targGT().isPhased()) {
            return fpd.targGT();
        }
        EstPhase hiFreqEstPhase = phaseHiFreqVariants(fpd);
        return phaseLowFreqVariants(fpd, hiFreqEstPhase);
    }

    private EstPhase phaseHiFreqVariants(FixedPhaseData fpd) {
        recombFactor = (float) 0.04f*par.ne()/fpd.nHaps();
        EstPhase hiFreqPhase = createEstPhase(fpd);
        int nIts = par.burnin() + par.iterations();
        for (int it=0; it<nIts; ++it) {
            long t0 = System.nanoTime();
            PhaseData phaseData = new PhaseData(fpd, hiFreqPhase,
                    recombFactor, it, rand.nextLong());
            PhaseLS.runStage1(phaseData, updateRecombFactor(it));
            recombFactor = phaseData.recombFactor();
            long nanos = System.nanoTime() - t0;
            runStats.printPhasingItUpdate(it+1, it<par.burnin(), nanos);
        }
        return hiFreqPhase;
    }

    private EstPhase createEstPhase(FixedPhaseData fpd) {
        GT targGT = fpd.hiFreqTargGT();
        RefGT refGT = fpd.hiFreqRefGT();
        int overlap = fpd.hiFreqOverlap();
        return new EstPhase(targGT, refGT, overlap, rand.nextLong());
    }

    private boolean updateRecombFactor(int it) {
        int nBurninIts = par.burnin();
        return (it==(nBurninIts-1) || it==nBurninIts);
    }

    private GT phaseLowFreqVariants(FixedPhaseData fpd, EstPhase estPhase) {
        if (fpd.hiFreqIndices().size() < fpd.targGT().nMarkers()) {
            int nIts = par.burnin() + par.iterations();
            PhaseData phaseData = new PhaseData(fpd, estPhase, recombFactor,
                    nIts, rand.nextLong());
            return PhaseLS.runStage2(phaseData);
        }
        else {
            return estPhase.hapsGT();
        }
    }
}

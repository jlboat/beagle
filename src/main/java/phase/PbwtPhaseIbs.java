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

import beagleutil.PbwtUpdater;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import main.Par;

/**
 * <p>Class {@code PbwtPhaseIBS} uses the Positional Burrows-Wheeler
 * Transform (PBWT) to find long IBS haplotypes for each sample that
 * contain a specified small genomic interval.</p>
 *
 * <p>Instances of class {@code PBWT} are thread-safe.</p>
 *
 * <p>Reference: Durbin, R. 2014. Bioinformatics 30(9):1266–1272.
 * doi:10.1093/bioinformatics/btu014</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PbwtPhaseIbs implements PhaseIbs {

    private static final int BURNIN_CANDIDATES = 100;
    private static final int MAX_PHASE_CANDIDATES = 90;
    private static final int MIN_PHASE_CANDIDATES = 5;

    private final PhaseData phaseData;
    private final CodedSteps codedSteps;
    private final int[][] ibsHaps;  //[step][targ hap]

    /**
     * Constructs a new {@code PbwtPhaseIBS} instance from the
     * specified data.
     * @param phaseData the current input data for genotype phasing
     * @param useBwd {@code true} if reverse order PBWT should be used
     * @throws NullPointerException if {@code phaseData == null}
     */
    public PbwtPhaseIbs(PhaseData phaseData, boolean useBwd) {
        Par par = phaseData.par();
        this.phaseData = phaseData;
        this.codedSteps = phaseData.codedSteps();
        int nSteps = codedSteps.nSteps();
        int initNBatches = Math.min(nSteps, par.nthreads());
        int stepsPerBatch = (codedSteps.nSteps()/initNBatches) + 1;
        int nBatches = nSteps / stepsPerBatch;
        if (nBatches * stepsPerBatch < nSteps) {
            ++nBatches;
        }
        int nBufferSteps = (int) Math.rint(par.buffer()/par.phase_step());

        if (useBwd) {
            this.ibsHaps = IntStream.range(0, nBatches)
                    .parallel()
                    .mapToObj(j -> bwdIbsHaps(phaseData, j, nBufferSteps, stepsPerBatch))
                    .flatMap(a -> Arrays.stream(a))
                    .toArray(int[][]::new);
        }
        else {
            this.ibsHaps = IntStream.range(0, nBatches)
                    .parallel()
                    .mapToObj(j -> fwdIbsHaps(phaseData, j, nBufferSteps, stepsPerBatch))
                    .flatMap(a -> Arrays.stream(a))
                    .toArray(int[][]::new);
        }
    }

    private static int nCandidates(PhaseData phaseData) {
        int nCandidates = BURNIN_CANDIDATES;
        int it = phaseData.it();
        Par par = phaseData.par();
        if (it>=par.burnin()) {
            double nItsRemaining = par.burnin() + par.iterations() - it;
            double p = (double) nItsRemaining / par.iterations();
            nCandidates = (int) Math.round(p*MAX_PHASE_CANDIDATES);
            nCandidates = Math.max(nCandidates, MIN_PHASE_CANDIDATES);
        }
        return Math.min(nCandidates, phaseData.nHaps());
    }

    private static int[][] bwdIbsHaps(PhaseData phaseData, int batch,
            int nBufferSteps, int stepsPerBatch) {
        CodedSteps codedSteps = phaseData.codedSteps();
        int nCandidates = nCandidates(phaseData);
        int nSteps = codedSteps.nSteps();
        int startStep = batch*stepsPerBatch;
        int endStep = Math.min(startStep + stepsPerBatch, nSteps);
        int bufferEnd = Math.min(endStep + nBufferSteps, nSteps);
        assert startStep < nSteps;

        int[][] ibsHaps = new int[endStep - startStep][];
        int nHaps = codedSteps.nHaps();
        PbwtUpdater pbwt = new PbwtUpdater(nHaps);
        int[] a = IntStream.range(0, nHaps).toArray();
        int[] d = IntStream.range(0, nHaps+1).map(j -> (bufferEnd-1)).toArray(); // last entry is sentinal

        for (int step=(bufferEnd-1); step>=endStep; --step) {
            int nAlleles = codedSteps.get(step).valueSize();
            pbwt.bwdUpdate(codedSteps.get(step), nAlleles, step, a, d);
        }
        for (int step=(endStep-1); step>=startStep; --step) {
            int nAlleles = codedSteps.get(step).valueSize();
            pbwt.bwdUpdate(codedSteps.get(step), nAlleles, step, a, d);
            ibsHaps[step-startStep] = getBwdIbsHaps(phaseData, step,
                    a, d, nCandidates);
        }
        return ibsHaps;
    }

    private static int[][] fwdIbsHaps(PhaseData phaseData, int batch,
            int nBufferSteps, int stepsPerBatch) {
        CodedSteps codedSteps = phaseData.codedSteps();
        int nCandidates = nCandidates(phaseData);
        int nSteps = codedSteps.nSteps();
        int startStep = batch*stepsPerBatch;
        int endStep = Math.min(startStep + stepsPerBatch, nSteps);
        int bufferStart = Math.max(0, startStep - nBufferSteps);
        assert startStep < nSteps;

        int[][] ibsHaps = new int[endStep - startStep][];
        int nHaps = codedSteps.nHaps();
        PbwtUpdater pbwt = new PbwtUpdater(nHaps);
        int[] a = IntStream.range(0, nHaps).toArray();
        int[] d = IntStream.range(0, nHaps+1).map(j -> bufferStart).toArray(); // last entry is sentinal

        for (int step=bufferStart; step<startStep; ++step) {
            int nAlleles = codedSteps.get(step).valueSize();
            pbwt.fwdUpdate(codedSteps.get(step), nAlleles, step, a, d);
        }
        for (int step=startStep; step<endStep; ++step) {
            int nAlleles = codedSteps.get(step).valueSize();
            pbwt.fwdUpdate(codedSteps.get(step), nAlleles, step, a, d);
            ibsHaps[step-startStep] = getfwdIbsHaps(phaseData, step,
                    a, d, nCandidates);
        }
        return ibsHaps;
    }

    private static int[] getBwdIbsHaps(PhaseData phaseData,
            int step, int[] a, int[] d, int maxCandidates) {
        CodedSteps codedSteps = phaseData.codedSteps();
        Random rand = new Random(phaseData.seed() + step);
        int nTargHaps = codedSteps.nTargHaps();
        int mStart = codedSteps.stepStart(step);
        int mInclEnd = codedSteps.stepEnd(step) - 1;
        int[] selectedHaps = new int[nTargHaps];
        d[0] = d[a.length] = step - 2;  // set sentinals
        // no need to save and restore old d[0], d[a.length] values
        for (int i=0; i<a.length; ++i) {
            if (a[i]<nTargHaps) {
                int u = i;          // inclusive start
                int v = i + 1;      // exclusive end
                int uMatchEnd = d[u];
                int vMatchEnd = d[v];
                while ((v - u)<maxCandidates && (step<=uMatchEnd || step<=vMatchEnd)) {
                    if (uMatchEnd<=vMatchEnd) {
                        vMatchEnd = Math.min(d[++v], vMatchEnd);
                    }
                    else {
                        uMatchEnd = Math.min(d[--u], uMatchEnd);
                    }
                }
                selectedHaps[a[i]] = getMatch(phaseData, mStart,
                        mInclEnd, i, u, v, a, rand);
            }
        }
        return selectedHaps;
    }

    private static int[] getfwdIbsHaps(PhaseData phaseData,
            int step, int[] a, int[] d, int maxCandidates) {
        CodedSteps codedSteps = phaseData.codedSteps();
        Random rand = new Random(phaseData.seed() + step);
        int nTargHaps = codedSteps.nTargHaps();
        int mStart = codedSteps.stepStart(step);
        int mInclEnd = codedSteps.stepEnd(step) - 1;
        int[] selectedHaps = new int[nTargHaps];
        d[0] = d[a.length] = step + 2;  // set sentinals
        // no need to save and restore old d[0], d[a.length] values
        for (int i=0; i<a.length; ++i) {
            if (a[i]<nTargHaps) {
                int u = i;          // inclusive start
                int v = i + 1;      // exclusive end
                int uMatchStart = d[u];
                int vMatchStart = d[v];
                while ((v - u)<maxCandidates && (uMatchStart<=step || vMatchStart<=step)) {
                    if (vMatchStart<=uMatchStart) {
                        vMatchStart = Math.max(d[++v], vMatchStart);
                    }
                    else {
                        uMatchStart = Math.max(d[--u], uMatchStart);
                    }
                }
                selectedHaps[a[i]] = getMatch(phaseData, mStart,
                        mInclEnd, i, u, v, a, rand);
            }
        }
        return selectedHaps;
    }

    private static int getMatch(PhaseData phaseData, int mStart, int mInclEnd,
            int i, int iStart, int iEnd, int[] a, Random rand) {
        int iLength = iEnd - iStart;
        if (iLength==1) {
            return -1;
        }
        Ibs2 ibs2 = phaseData.ibs2();
        int sample = a[i]>>1;
        int match = -1;
        int index = iStart + rand.nextInt(iLength);
        for (int j=0; j<iLength && match==-1; ++j) {
            int sample2 = a[index]>>1;
            if (ibs2.areIbs2(sample, sample2, mStart)==false
                    && ibs2.areIbs2(sample, sample2, mInclEnd)==false) {
                match = a[index];
            }
            if (++index==iEnd) {
                index = iStart;
            }
        }
        return match;
    }

    @Override
    public PhaseData phaseData() {
        return phaseData;
    }

    @Override
    public void addIbsHaps(int hap, int step, PhaseStates phaseStates) {
        if (ibsHaps[step][hap]>=0) {
            phaseStates.updateFields(ibsHaps[step][hap], step);
        }
    }
}

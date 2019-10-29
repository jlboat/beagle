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

import blbutil.Const;
import blbutil.Utilities;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import math.Regress;
import vcf.GT;

/**
 * <p>Class {@code PhaseLS} estimated genotypes phase using
 * a haploid Li and Stephens hidden Markov model.</p>
 *
 * <p>Instances of class {@code PhaseLS} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseLS {

    private PhaseLS() {
        // private constructor to prevent instantiation
    }

    /**
     * Estimates phased haplotypes at high-frequency markers in the
     * target samples.
     * @param phaseData the current input data for updating genotype phase
     * estimates
     * @param updateRecombFactor {@code true} if
     * {@code phaseData.recombFactor()} should be updated
     * @throws NullPointerException if {@code phaseData == null}
     */
    public static void runStage1(PhaseData phaseData, boolean updateRecombFactor) {
        long t0 = System.nanoTime();
        boolean useBwd = (phaseData.it() & 1)==0;
        PhaseIbs phaseIbs = new PbwtPhaseIbs(phaseData, useBwd);

        ExecutorService es = createExecService(phaseData.par().nthreads());
        if (updateRecombFactor) {
            updateRecombFactor(es, phaseData, phaseIbs);
        }
        phase(es, phaseData, phaseIbs);
        shutdownExecService(es);
    }

    private static void updateRecombFactor(ExecutorService es,
            PhaseData phaseData, PhaseIbs phaseIbs) {
        int nThreads = phaseData.par().nthreads();
        double maxSumY = Math.max(5000.0/nThreads, 200.0);
        Regress regress = new Regress();
        int nSamples = phaseData.targGT().nSamples();
        CountDownLatch latch = new CountDownLatch(nThreads);
        for (int j=0; j<nThreads; ++j) {
            Random rand = new Random(phaseData.seed() +  j);
            es.submit(() -> {
                try {
                    RecombRegress rr = new RecombRegress(phaseIbs, regress);
                    while (rr.sumY()<maxSumY) { // thread-confined for reproducability
                        rr.update(rand.nextInt(nSamples));
                    }
                    latch.countDown();
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        await(latch);
        float recombFactor = (float) regress.beta();
        if (recombFactor>0.0 && Float.isFinite(recombFactor)) {
            phaseData.setRecombFactor(recombFactor);
        }
        else {
            System.out.println(Const.nl + "WARNING: no recombFactor update: "
                    + recombFactor + Const.nl);
        }
    }

    private static void phase(ExecutorService es, PhaseData phaseData,
            PhaseIbs phaseIbs) {
        int nThreads = phaseData.par().nthreads();
        AtomicInteger sampleIndices = new AtomicInteger(0);
        int nSamples = phaseData.targGT().nSamples();
        CountDownLatch latch = new CountDownLatch(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    PhaseBaum1 baum = new PhaseBaum1(phaseIbs);
                    int sample = sampleIndices.getAndIncrement();
                    while (sample>=0 && sample<nSamples) {
                        baum.phase(sample);
                        sample = sampleIndices.getAndIncrement();
                    }
                    latch.countDown();
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        await(latch);
    }

    /**
     * Returns phased genotypes at all markers.
     * @param phaseData the current input data for updating genotype phase
     * estimates
     * @return phased genotypes at all markers
     * @throws NullPointerException if {@code fpd == null || phaseData == null}
     */
    public static GT runStage2(PhaseData phaseData) {
        GT gt = phaseData.fpd().targGT();
        int nThreads = phaseData.par().nthreads();
        int nSamples = gt.nSamples();
        HapImputer imputableHaps = new HapImputer(gt.markers(), gt.samples());
        PhaseIbs phaseIbs = new LowFreqPhaseIbs(phaseData);

        AtomicInteger sampleIndex = new AtomicInteger(0);
        ExecutorService es = createExecService(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(() -> {
                try {
                    ImputeBaum baum = new ImputeBaum(phaseIbs, imputableHaps);
                    int sample = sampleIndex.getAndIncrement();
                    while (sample>=0 && sample<nSamples) {
                        baum.phase(sample);
                        sample = sampleIndex.getAndIncrement();
                    }
                }
                catch (Throwable t) {
                    Utilities.exit(t);
                }
            } );
        }
        shutdownExecService(es);
        GT imputedGT = imputableHaps.imputedHaps();
        return imputedGT;
    }

    private static void await(CountDownLatch latch) {
        try {
            latch.await();
        }
        catch (InterruptedException e) {
            Utilities.exit("ERROR", e);
        }
    }

    private static ExecutorService createExecService(int nThreads) {
        return Executors.newFixedThreadPool(nThreads);
    }

    private static void shutdownExecService(ExecutorService es) {
        try {
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (InterruptedException e) {
            Utilities.exit("ERROR", e);
        }
    }
}

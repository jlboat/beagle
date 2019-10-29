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

import vcf.HapsGT;
import blbutil.FloatArray;
import main.Par;
import vcf.GT;
import vcf.MarkerMap;
import vcf. RefGT;

/**
 * <p>Class {@code PhaseData} contains the current input data for updating
 * genotype phase.
 * </p>
 * <p>Instances of class {@code PhaseData} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PhaseData {

    private final FixedPhaseData fpd;
    private final Par par;
    private final EstPhase estPhase;
    private final GT targGT;
    private final HapsGT phasedTargHaps;
    private final RefGT refGT;
    private final MarkerMap map;
    private final Ibs2 ibs2;
    private final CodedSteps codedSteps;
    private final int nTargHaps;
    private final int nHaps;
    private final float err;
    private volatile float recombFactor;
    private volatile FloatArray pRecomb;

    private final int it;
    private final int nItsRemaining;
    private final long seed;

    /**
     * Constructs a new {@code PhaseData} instance from the specified
     * data.
     * @param fpd the input data for phasing that is the same in each iteration
     * @param estPhase the current estimate of phased target genotypes
     * @param recombFactor the factor multiplied by genetic distance to
     * obtain the probability of transitioning to a random HMM state
     * @param it the 0-based iteration
     * @param seed seed for random numbers
     *
     * @throws IllegalArgumentException if
     * {@code recombFactor < 0 || Double.isFinite(recombFactor) == false}
     * @throws IllegalArgumentException if
     * {@code fpd.hiFreqTargGT().equals(estPhase.targGT()) == false}
     * @throws IllegalArgumentException if
     * {@code it < 0 || it > (fpd.par().burnin() + fpd.par().iterations())}
     * @throws NullPointerException if {@code fpd == null || estPhase == null}
     */
    public PhaseData(FixedPhaseData fpd, EstPhase estPhase,
            float recombFactor, int it, long seed) {
        if (recombFactor < 0 || Double.isFinite(recombFactor)==false) {
            throw new IllegalArgumentException(String.valueOf(recombFactor));
        }
        if (fpd.hiFreqTargGT().equals(estPhase.targGT())==false) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (it<0 || it > (fpd.par().burnin() + fpd.par().iterations())) {
            throw new IllegalArgumentException(String.valueOf(it));
        }
        this.fpd = fpd;
        this.par = fpd.par();
        this.estPhase = estPhase;
        this.targGT = fpd.hiFreqTargGT();
        this.phasedTargHaps = estPhase.hapsGT();
        this.refGT = fpd.hiFreqRefGT();
        this.map = fpd.hiFreqMap();
        this.ibs2 = fpd.hiFreqIbs2();
        this.codedSteps = new CodedSteps(phasedTargHaps, refGT, map,
                par.phase_step(), par.scaleFactor(), seed);
        this.nTargHaps = targGT.nHaps();
        this.nHaps = fpd.nHaps();
        this.err = fpd.err();
        this.recombFactor = recombFactor;
        this.pRecomb = map.pRecomb(recombFactor);
        this.it = it;
        this.nItsRemaining = par.burnin() + par.iterations() - it;
        this.seed = seed;
    }

    /**
     * Returns the input data for phasing that is the same in each iteration.
     * @return the input data for phasing that is the same in each iteration
     */
    public FixedPhaseData fpd() {
        return fpd;
    }
    /**
     * Returns the command line parameters
     * @return the command line parameters
     */
    public Par par() {
        return par;
    }

    /**
     * Returns the estimated phase.
     * @return the estimated phase
     */
    public EstPhase estPhase() {
        return estPhase;
    }

    /**
     * Returns the allele carried by the specified high-frequency  marker
     * and haplotype. The first reference haplotype index is
     * this.targGT().nHaps().
     * @param marker the high-frequency marker index
     * @param hap the haplotype index
     * @return the allele carried by the specified marker and haplotype
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.targGT().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.nHaps()}
     */
    public int allele(int marker, int hap) {
        if (hap>=nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        if (hap < nTargHaps) {
            return phasedTargHaps.allele(marker, hap);
        }
        else {
            return refGT.allele(marker, hap - nTargHaps);
        }
    }

    /**
     * Returns the estimated phased target haplotypes at the high-frequency
     * markers.
     * @return the estimated phased target haplotypes at the high-frequency
     * markers
     */
    public GT phasedTarg() {
        return phasedTargHaps;
    }

    /**
     * Return the total number of target and reference haplotypes.
     * @return the total number of target and reference haplotypes
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Return the allele mismatch probability.
     * @return the allele mismatch probability
     */
    public float err() {
        return err;
    }

    /**
     * Returns the input genotypes for the target samples at the
     * high-frequency target data markers.
     * @return the input genotypes for the target samples at the
     * high-frequency target data markers
     */
    public GT targGT() {
        return targGT;
    }

    /**
     * Returns the phased, nonmissing reference genotypes at the
     * high-frequency markeres, or {@code null} if there are
     * no reference samples.
     * @return the phased, nonmissing reference genotypes at the
     * high-frequency markeres, or {@code null} if there are
     * no reference samples
     */
    public GT refGT() {
        return refGT;
    }

    /**
     * Returns the ibs2 status at each high-frequency marker.
     * @return the ibs2 status at each high-frequency marker
     */
    public Ibs2 ibs2() {
        return ibs2;
    }

    /**
     * Returns the coded steps.
     * @return the coded steps
     */
    public CodedSteps codedSteps() {
        return codedSteps;
    }

    /**
     * Return the marker map that stores genetic positions of each
     * high-frequency marker.
     * @return the marker map that stores genetic positions of each
     * high-frequency marker
     */
    public MarkerMap map() {
        return map;
    }

    /**
     * Returns the factor multiplied by genetic distance to
     * obtain the probability of transitioning to a random HMM state.
     * @return the recombination factor
     */
    public float recombFactor() {
        return recombFactor;
    }

    /**
     * Return a {@code FloatArray} of size {@code this.targGT().nMarkers()}
     * whose {@code k}-th element is the probability of transitioning to a
     * random HMM state between the {@code k}-th high-frequency target marker
     * and the previous high-frequency marker.
     * @return a {@code FloatArray} of size {@code this.targGT().nMarkers()}
     * whose {@code k}-th element is the probability of transitioning to a
     * random HMM state between the {@code k}-th high-frequency target marker
     * and the previous high-frequency marker
     */
    public FloatArray pRecomb() {
        return pRecomb;
    }

    /**
     * Returns the proportion of unphased heterozygotes that should be
     * left unphased at the end of iteration {@code this.it()} for the
     * specified sample.
     * @param sample the sample index
     * @return the proportion of unphased heterozygotes that should be
     * left unphased at the end of iteration {@code this.it()} for the
     * specified sample
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.targGT().nSamples()}
     */
    public double leaveUnphasedProp(int sample) {
        if (sample<0 || sample>=this.targGT.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        if (it<par.burnin()) {
            return 1.0;
        }
        else if (nItsRemaining==1) {
            return 0.0;
        }
        else {
            return Math.pow(estPhase.unphased(sample).size(), -1.0/nItsRemaining);
        }
    }

    /**
     * Returns the current iteration (0 based).
     * @return the current iteration
     */
    public int it() {
        return it;
    }

    /**
     * Returns the seed for generating random numbers.
     * @return the seed for generating random numbers
     */
    public long seed() {
        return seed;
    }

    /**
     * Sets the recombination factor to the specified value.
     * @param recombFactor the new value for the recombination factor
     * @throws IllegalArgumentException if
     * {@code recombFactor < 0 || Double.isFinite(recombFactor) == false}
     */
    public void setRecombFactor(float recombFactor) {
        if (recombFactor < 0 || Double.isFinite(recombFactor)==false) {
            throw new IllegalArgumentException(String.valueOf(recombFactor));
        }
        this.recombFactor = recombFactor;
        this.pRecomb =  map.pRecomb(recombFactor);
    }
}

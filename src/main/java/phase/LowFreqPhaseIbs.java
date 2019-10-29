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

import java.util.stream.IntStream;

/**
 * <p>Class {@code LowFreqPhaseIbs} identifies haplotypes that share a long
 * IBS segment or a low frequency variant with a specified haplotype
 * in a specified genomic interval.</p>
 *
 * <p>Instances of {@code LowFreqPhaseIbs} are immutable.</p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowFreqPhaseIbs implements PhaseIbs {

    private final PhaseData phaseData;
    private final CodedSteps codedSteps;

    private final PhaseIbs fwdPhaseIbs;
    private final PhaseIbs bwdPhaseIbs;

    private final int[][] fwdMatch;    //[step][targ hap]
    private final int[][] bwdMatch;    //[step][targ hap]

    /**
     * Constructs a new {@code LowFreqPhaseIbs} object from the specified data.
     * @param phaseData the current input data for updating genotype phase
     * estimates at high-frequency markers
     *
     * @throws IllegalArgumentException if
     * {@code fpd.hiFreqTargGT().markers().equals(phaseData.targGT().markers()) == false}
     * @throws NullPointerException if {@code fpd == null || phaseData == null}
     */
    public LowFreqPhaseIbs(PhaseData phaseData) {
        this.phaseData = phaseData;
        this.codedSteps = phaseData.codedSteps();
        this.fwdPhaseIbs = new PbwtPhaseIbs(phaseData, false);
        this.bwdPhaseIbs = new PbwtPhaseIbs(phaseData, true);

        BestMatch bestMatch = new BestMatch(phaseData);
        this.fwdMatch = IntStream.range(0, codedSteps.nSteps())
                .parallel()
                .mapToObj(step -> bestMatch.fwdMatch(step))
                .toArray(int[][]::new);
        this.bwdMatch = IntStream.range(0, codedSteps.nSteps())
                .parallel()
                .mapToObj(step -> bestMatch.bwdMatch(step))
                .toArray(int[][]::new);
    }

    @Override
    public PhaseData phaseData() {
        return phaseData;
    }

    @Override
    public void addIbsHaps(int hap, int step, PhaseStates phaseStates) {
        fwdPhaseIbs.addIbsHaps(hap, step, phaseStates);
        bwdPhaseIbs.addIbsHaps(hap, step, phaseStates);
        if (fwdMatch[step][hap] != -1) {
            phaseStates.updateFields(fwdMatch[step][hap], step);
        }
        if (bwdMatch[step][hap] != -1) {
            phaseStates.updateFields(bwdMatch[step][hap], step);
        }
    }
}

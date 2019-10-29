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
 * <p>Instances of interface {@code PhaseIbs} identify haplotypes that are
 * identical by state with a specified haplotype in a specified genomic
 * interval.</p>
 *
 * <p>All instances of {@code PhaseIbs} are required to be immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface PhaseIbs {

    /**
     * Returns the input data for the next phase update.
     * @return the input data for the next phase update
     */
    public PhaseData phaseData();

    /**
     * Adds haplotypes that are IBS with the specified haplotype to the HMM
     * state space.
     * @param hap a haplotype index
     * @param step an index of a genomic interval
     * @param phaseStates the object for constructing the HMM state space
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.phaseData().targGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.codedSteps().nSteps()}
     * @throws NullPointerException if {@code phaseStates == null}
     */
    void addIbsHaps(int hap, int step, PhaseStates phaseStates);
}

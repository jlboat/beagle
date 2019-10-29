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

/**
 * <p>Class {@code CompHapSegment} represents a copied haplotype segment
 * in a composite reference haplotype.</p>
 *
 * <p>Instances of class {@code CompHapSegment} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CompHapSegment implements Comparable<CompHapSegment> {

    private int hap;
    private int step;
    private final int compHapIndex;

    /**
     * Constructs a new {@code CompHapSegment} from the specified data.
     * @param hap the haplotype
     * @param step the step
     * @param compHapIndex the composite haplotype index
     */
    public CompHapSegment(int hap, int step, int compHapIndex) {
        this.hap = hap;
        this.step = step;
        this.compHapIndex = compHapIndex;
    }

    /**
     * Update the copied haplotype to the specified haplotype
     * @param hap the new haplotype
     */
    public void updateHap(int hap) {
        this.hap = hap;
    }

    /**
     * Updates the step to the specified value
     * @param step the new step value
     */
    public void updateStep(int step) {
        this.step = step;
    }

    /**
     * Returns the copied haplotype.
     * @return the copied haplotype
     */
    public int hap() {
        return hap;
    }

    /**
     * Returns the last recorded step for {@code this.hap()}.
     * @return the last recorded step for {@code this.hap()}
     */
    public int step() {
        return step;
    }

    /**
     * Returns the composite haplotype index.
     * @return the composite haplotype index
     */
    public int compHapIndex() {
        return compHapIndex;
    }

    /**
     * Compares the specified segment to {@code this} for order.  Returns
     * -1, 0, or 1 according to whether {@code this.end()} is less than,
     * equal, or greater than {@code seg.end()}.
     * @param seg the object to be compared
     * @return -1, 0, or 1 according to whether {@code this.end()} is less
     * than, equal, or greater than {@code seg.end()}
     */
    @Override
    public int compareTo(CompHapSegment seg) {
        if (this.step!=seg.step) {
            return this.step<seg.step ? -1 : 1;
        } else {
            return 0;
        }
    }

}

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
package vcf;

import ints.IntArray;
import ints.WrappedIntArray;
import java.io.Closeable;
import main.Pedigree;

/**
 * Interface {@code Data} represents a sliding window of target VCF records
 * or a sliding window of reference and target VCF records.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Data extends Closeable {

    public static final IntArray ZERO_FREQ_ARRAY = new WrappedIntArray(new int[0]);
    public static final IntArray HIGH_FREQ_ARRAY = new WrappedIntArray(new int[0]);

    /**
     * Returns the pedigree.
     * @return the pedigree
     */
    Pedigree ped();

    /**
     * Returns the genetic map.
     * @return the genetic map
     */
    GeneticMap genMap();

    /**
     * Returns {@code true} if the current window of VCF records is the last
     * window for the chromosome and returns {@code false} otherwise.
     * @return {@code true} if the current window of VCF records is the last
     * window for the chromosome
     */
    boolean lastWindowOnChrom();

    /**
     * Returns {@code true} if the sliding window of VCF records can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding window of VCF records can advance
     */
    boolean canAdvanceWindow();

    /**
     * Advances the sliding window of VCF records.
     *
     * @throws IllegalArgumentException if a format error in the input data
     * is detected
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow() == false}
     */
    void advanceWindowCm();

    /**
     * Returns the current window index.  The first window has index 1.
     * @return the current window index
     */
    public int windowIndex();

    /**
     * Returns the number of target VCF records in the union of the
     * current window and all previous windows.
     * @return the number of target VCF records in the union of the
     * current window and all previous windows
     */
    int nTargMarkersSoFar();

    /**
     * Returns the number of markers in the current window.
     * @return the number of markers in the current window
     */
     int nMarkers();

    /**
     * Returns the number of markers in the union of the current window
     * and all previous windows.
     * @return the number of markers in the union of the current window
     * and all previous windows
     */
    int nMarkersSoFar();

    /**
     * Returns the genotype likelihoods for the target samples
     * restricted to the target data markers in the current window.
     * The returned {@code GL} instance will contain no markers if
     * {@code this.advanceWindow()} has not yet been invoked.
     * @return the genotype likelihoods for the target samples
     * restricted to the target data markers in the current window
     */
    GT targGT();

    /**
     * Returns the  phased, nonmissing reference genotype data
     * for the current window, or {@code null} if there are no reference data
     * @return the reference genotype data for the current window or
     * {@code null} if there are no reference data
     */
    RefGT refGT();

    /**
     * Returns the  phased, nonmissing reference genotype data
     * for the target data markers in the current window.
     * Returns {@code null} if there are no reference data
     * @return the reference genotype data for the target data markers
     * in the current window or {@code null} if there are no reference data
     */
    RefGT restrictRefGT();

    /**
     * <p>Returns the indices of the reference and target carriers for each
     * low-frequency allele at the target data markers.  The reference sample
     * indices will be shifted by the number of target samples so that the
     * first reference sample will have an index equal to the number of target
     * samples. An element of the returned array will be empty and equal to
     * {@code Data.ZERO_FREQ_ARRAY} if the allele has no carriers, and the
     * the element will be empty and equal to {@code Data.HIGH_FREQ_ARRAY}
     * if the number of carriers of the allele exceeds the specified
     * maximum number of carriers.</p>
     *
     * <p>The list of carriers for the {@code k}-th allele of the {@code j}-th
     * target marker are stored in entry {@code (j, k)} of the returned array.
     * if the number of carriers is less than or equal to the specified
     * maximum number of carriers.</p>
     * @param maxCarriers the maximum number of carriers in any list
     * of the returned array.
     * @return the indices of the reference and target carriers for each
     * low-frequency allele
     */
    IntArray[][] carriers(int maxCarriers);

    /**
     * Return a {@code MarkerIndices} instance which stores the overlap
     * with the current marker window and adjacent marker windows
     * and the mappings between marker indices and target marker indices.
     * @return a {@code MarkerIndices} instance
     */
    MarkerIndices markerIndices();

    /**
     * Releases any I/O resources controlled by this object.
     */
    @Override
    void close();
}

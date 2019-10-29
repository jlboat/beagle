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

import ints.IntList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code MarkerIndices} stores the overlap with adjacent marker
 * windows and the mappings between marker indices and the target marker
 * indices.
 * </p>
 * <p>Instances of class {@code MarkerIndices} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MarkerIndices {

    private final int prevSplice;
    private final int nextOverlap;
    private final int nextSplice;

    private final int[] targMarkerToMarker;
    private final int[] markerToTargMarker;

    private final int prevTargSplice;
    private final int nextTargOverlap;
    private final int nextTargSplice;

    /**
     * Constructs a {@code MarkerIndices} instance from the specified data.
     * @param inTarg an array whose {@code k}-th element is {@code true}
     * if the {@code k}-th marker is present in the target genotype data
     * @param prevOverlap the ending marker index (exclusive) for the overlap
     * with the previous window
     * @param nextOverlap the starting marker index (inclusive) for the overlap
     * with the next window
     *
     * @throws IndexOutOfBoundsException if
     * {@code prevOverlap < 0 || prevOverlap > inTarg.length}
     * @throws IndexOutOfBoundsException if
     * {@code nextOverlap < 0 || nextOverlap > inTarg.length}
     * @throws NullPointerException if {@code inTarg==null}
     */
    public MarkerIndices(boolean[] inTarg, int prevOverlap, int nextOverlap) {
        if (prevOverlap < 0 || prevOverlap > inTarg.length) {
            throw new IndexOutOfBoundsException(String.valueOf(prevOverlap));
        }
        if (nextOverlap < 0 || nextOverlap > inTarg.length) {
            throw new IndexOutOfBoundsException(String.valueOf(nextOverlap));
        }
        this.prevSplice = prevOverlap/2;
        this.nextOverlap = nextOverlap;
        this.nextSplice = (inTarg.length + nextOverlap) >>> 1;

        this.targMarkerToMarker = targMarkerToMarker(inTarg);
        this.markerToTargMarker = markerToTargMarker(targMarkerToMarker,
                inTarg.length);

        this.prevTargSplice = targIndex(targMarkerToMarker, prevSplice);
        this.nextTargOverlap = targIndex(targMarkerToMarker, nextOverlap);
        this.nextTargSplice = targIndex(targMarkerToMarker, nextSplice);
    }

    /**
     * Constructs a {@code MarkerIndices} instance from the specified data
     * @param prevOverlap the ending marker index (exclusive) for the overlap
     * with the previous window
     * @param nextOverlap the starting marker index (inclusive) for the overlap
     * with the next window
     * @param nMarkers the number of markers
     *
     * @throws IllegalArgumentException if {@code nMarkers < 0}
     * @throws IndexOutOfBoundsException if
     * {@code prevOverlap < 0 || prevOverlap > nMarkers}
     * @throws IndexOutOfBoundsException if
     * {@code nextOverlap < 0 || nextOverlap > nMarkers}
     */
    public MarkerIndices(int prevOverlap, int nextOverlap, int nMarkers) {
        if (nMarkers < 0) {
            throw new IllegalArgumentException(String.valueOf(nMarkers));
        }
        if (prevOverlap < 0 || prevOverlap > nMarkers) {
            throw new IndexOutOfBoundsException(String.valueOf(prevOverlap));
        }
        if (nextOverlap < 0 || nextOverlap > nMarkers) {
            throw new IndexOutOfBoundsException(String.valueOf(nextOverlap));
        }
        this.prevSplice = prevOverlap/2;
        this.nextOverlap = nextOverlap;
        this.nextSplice = (nMarkers + nextOverlap) >>> 1;

        int[] identityMap = IntStream.range(0, nMarkers)
                .parallel()
                .toArray();
        this.targMarkerToMarker = identityMap;
        this.markerToTargMarker = identityMap;

        this.prevTargSplice = prevSplice;
        this.nextTargOverlap = nextOverlap;
        this.nextTargSplice = nextSplice;
    }

    /* return first target marker on or after specified marker */
    private static int targIndex(int[] targMarkerToMarker, int marker) {
        int insPt = Arrays.binarySearch(targMarkerToMarker, marker);
        if (insPt<0) {
            insPt = -insPt - 1;
        }
        return insPt;
    }

    private static int[] targMarkerToMarker(boolean[] inTarg) {
        IntList il = new IntList(inTarg.length/64 + 1);
        for (int j=0; j<inTarg.length; ++j) {
            if (inTarg[j]) {
                il.add(j);
            }
        }
        return il.toArray();
    }

    private int[] markerToTargMarker(int[] targMarkerToMarker, int nMarkers) {
        int[] ia = IntStream.range(0, nMarkers)
                .parallel()
                .map(i -> -1)
                .toArray();
        for (int j=0; j<targMarkerToMarker.length; ++j) {
            ia[targMarkerToMarker[j]] = j;
        }
        return ia;
    }

    /**
     * Returns the first marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} the next marker window is from
     * a different chromosome.
     * @return the first marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextOverlap() {
        return nextOverlap;
    }

    /**
     * Returns the first target marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} if there is no overlap or if there are
     * no target markers in the overlap.
     * @return the first target marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextTargOverlap() {
        return nextTargOverlap;
    }

    /**
     * Returns the first marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first marker index after the splice point with
     * the previous marker window
     */
    public int prevSplice() {
        return prevSplice;
    }

    /**
     * Returns the first target marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first target marker index after the splice point with
     * the previous marker window
     */
    public int prevTargSplice() {
        return prevTargSplice;
    }

    /**
     * Returns the first marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nMarkers()} if there is no overlap or if there are
     * no markers after the splice point.
     * @return the first marker index after the next splice point
     */
    public int nextSplice() {
        return nextSplice;
    }

    /**
     * Returns the first target marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nTargMarkers()} if there is no overlap or if there are
     * no target markers after the splice point
     * @return the first target marker index after the next splice point
     */
    public int nextTargSplice() {
        return nextTargSplice;
    }

    /**
     * Returns the index of the specified marker in the reference data markers.
     * @param targetMarker index of a marker in the list of target data markers
     * @return the index of the specified marker in the reference data markers
     * @throws IndexOutOfBoundsException if
     * {@code targetMarker < 0 || targetMarker >= this.nTargMarkers()}
     */
    public int targMarkerToMarker(int targetMarker) {
        return targMarkerToMarker[targetMarker];
    }

    /**
     * Returns an array of length {@code this.nTargMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers.
     * @return an array of length {@code this.nTargMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers
     */
    public int[] targMarkerToMarker() {
        return targMarkerToMarker.clone();
    }

    /**
     * Returns the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data.
     * @param marker index of a marker in the reference data
     * @return the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}.
     */
     public int markerToTargMarker(int marker) {
         return markerToTargMarker[marker];
     }

    /**
     * Returns an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data.
     * @return an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data
     */
     public int[] markerToTargMarker() {
         return markerToTargMarker.clone();
     }
}

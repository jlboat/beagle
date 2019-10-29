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

import beagleutil.Samples;
import blbutil.SampleFileIt;
import ints.IntArray;
import ints.IntList;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import main.Par;
import main.Pedigree;

/**
 * <p>Class {@code AllData} represents a sliding window of reference and
 * target VCF records.
 * </p>
 * <p>Instances of class {@code AllData} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AllData implements Data {

    private final Pedigree ped;
    private Window<RefGTRec> refWindow;
    private RefGT refGT;
    private RefGT restrictRefGT;
    private GTRec[] targRecs;  // missing markers as null entries
    private GT targGT;
    private MarkerIndices markerIndices;
    private int window = 0;

    private final WindowIt<RefGTRec> refWindowIt;
    private final RestrictedVcfWindow targWindow;

    private int cumMarkerCnt = 0;

    /**
     * Constructs and returns a new {@code AllData} instance from VCF
     * recList returned by the specified {@code SampleFileIt} objects.
     *
     * @param supplier an object to supply the reference file iterator
     * @param targIt an iterator that returns target VCF recList
     * @param par the command line parameters
     * @return a new {@code AllData} instance
     *
     * @throws IllegalArgumentException if either the reference data or
     * target data contain no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if
     * {@code refIt == null || targetIt == null || genMap == null}
     */
    public static AllData allData(Supplier<SampleFileIt<RefGTRec>> supplier,
            SampleFileIt<GTRec> targIt, Par par) {
        GeneticMap genMap = GeneticMap.geneticMap(par.map(), par.chromInt());
        WindowIt<RefGTRec> refWindow = WindowIt.newInstance(supplier, genMap,
                par.window(), par.overlap());
        if (refWindow.samples().nSamples()==0 || targIt.samples().nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        RestrictedVcfWindow targetWindow = new RestrictedVcfWindow(targIt);
        AllData allData = new AllData(par, refWindow, targetWindow);
        assert allData.canAdvanceWindow();
        allData.advanceWindowCm();
        return allData;
    }

    private AllData(Par par, WindowIt<RefGTRec> refWind,
            RestrictedVcfWindow targWind) {
        this.ped = new Pedigree(targWind.samples(), par.ped());
        this.refWindowIt = refWind;
        this.targWindow = targWind;

        this.refWindow = null;
        this.refGT = null;
        this.restrictRefGT = null;
        this.targRecs = new GTRec[0];
        this.targGT = targGT(targWind.samples(), targRecs, ped);
        this.markerIndices = null;
    }

    private static GT targGT(Samples samples, GTRec[] targData,
            Pedigree ped) {
        return new BasicGT(samples, targData);
    }

    @Override
    public Pedigree ped() {
        return ped;
    }

    @Override
    public GeneticMap genMap() {
        return refWindowIt.genMap();
    }

    @Override
    public boolean lastWindowOnChrom() {
        return refWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return refWindowIt.hasNext();
    }

    @Override
    public void advanceWindowCm() {
        refWindow = refWindowIt.next();
        cumMarkerCnt += (refWindow.nMarkers() - refWindow.prevOverlap());
        RefGTRec[] refRecs = refWindow.recList().toArray(new RefGTRec[0]);
        refGT = new RefGT(refRecs);
        targRecs = targWindow.advanceWindow(refGT.markers());
        this.markerIndices = new MarkerIndices(inTarg(targRecs),
                refWindow.prevOverlap(), refWindow.nextOverlap());
        int[] refIndices = markerIndices.targMarkerToMarker();
        targGT = targGTWindow(targWindow.samples(), targRecs, refIndices, ped);
        restrictRefGT = restrictRecs(targGT.markers(), refRecs, refIndices);
        ++window;
    }

    @Override
    public IntArray[][] carriers(int maxCarriers) {
        return IntStream.range(0, targGT.nMarkers())
                .parallel()
                .mapToObj(i -> carriers(i, maxCarriers))
                .toArray(IntArray[][]::new);
    }

    private IntArray[] carriers(int marker, int maxCarriers) {
        int nAlleles = targGT.marker(marker).nAlleles();
        IntList[] carriers = IntStream.range(0, nAlleles)
                .mapToObj(i -> new IntList(16))
                .toArray(IntList[]::new);
        int nTargSamples = targGT.nSamples();
        int nRefSamples = (restrictRefGT!=null) ? restrictRefGT.nSamples() : 0;
        for (int s=0; s<nTargSamples; ++s) {
            int a1 = targGT.allele1(marker, s);
            int a2 = targGT.allele2(marker, s);
            if (a1>=0 && carriers[a1].size()<=maxCarriers) {
                carriers[a1].add(s);
            }
            if (a2>=0 && a2!=a1 && carriers[a2].size()<=maxCarriers) {
                carriers[a2].add(s);
            }
        }
        for (int s=0; s<nRefSamples; ++s) {
            int a1 = restrictRefGT.allele1(marker, s);
            int a2 = restrictRefGT.allele2(marker, s);
            if (a1>=0 && carriers[a1].size()<=maxCarriers) {
                carriers[a1].add(nTargSamples + s);
            }
            if (a2>=0 && a2!=a1 && carriers[a2].size()<=maxCarriers) {
                carriers[a2].add(nTargSamples + s);
            }
        }
        return Arrays.stream(carriers)
                .map(list -> {
                    if (list.isEmpty()) {
                        return Data.ZERO_FREQ_ARRAY;
                    }
                    else if (list.size() <= maxCarriers) {
                        return new WrappedIntArray(list);
                    }
                    else {
                        return Data.HIGH_FREQ_ARRAY;
                    }
                })
                .toArray(IntArray[]::new);
    }

    @Override
    public int windowIndex() {
        return window;
    }

    public static boolean[] inTarg(GTRec[] recs) {
        boolean[] inTarg = new boolean[recs.length];
        for (int j=0; j<recs.length; ++j) {
            if (recs[j]!=null) {
                inTarg[j] = true;
            }
        }
        return inTarg;
    }

    private static GT targGTWindow(Samples samples, GTRec[] recs,
            int[] refMarkerIndex, Pedigree ped) {
        GTRec[] restricted = new GTRec[refMarkerIndex.length];
        for (int j=0; j<refMarkerIndex.length; ++j) {
            restricted[j] = recs[refMarkerIndex[j]];
        }
        return new BasicGT(samples, restricted);
    }

    private RefGT restrictRecs(Markers targMarkers,
            RefGTRec[] refData, int[] targToRef) {
        assert targMarkers.nMarkers()==targToRef.length;
        RefGTRec[] restricted = new RefGTRec[targToRef.length];
        for (int j=0; j<targToRef.length; ++j) {
            restricted[j] = refData[targToRef[j]];
        }
        return new RefGT(targMarkers, refWindowIt.samples(), restricted);
    }

    @Override
    public int nTargMarkersSoFar() {
        return targWindow.cumMarkerCnt();
    }

    @Override
    public int nMarkers() {
        return refGT.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public GT targGT() {
       return targGT;
    }

    @Override
    public RefGT refGT() {
        return refGT;
    }

    @Override
    public RefGT restrictRefGT() {
        return restrictRefGT;
    }

    @Override
    public void close() {
        refWindowIt.close();
        targWindow.close();
    }

    @Override
    public MarkerIndices markerIndices() {
        return markerIndices;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(20);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}

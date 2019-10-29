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
 * <p>Class {@code TargetData} represents a sliding window of
 * target VCF records.
 * </p>
 * <p>Instances of class {@code TargetData} are not thread-safe.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TargetData implements Data {

    private final Pedigree ped;
    private int window = 0;
    private Window<GTRec> currentWindow;
    private GTRec[] recs;
    private GT gt;
    private MarkerIndices markerIndices;

    private final WindowIt<GTRec> targWindowIt;
    private int cumMarkerCnt = 0;

    /**
     * Constructs and returns a new {@code TargetData} instance from
     * VcfRecords returned by the specified {@code SampleFileIt} objects.
     *
     * @param supplier a supplier for the sample file iterator
     * @param par the command line parameters
     * @return a new {@code TargetData} instance
     *
     * @throws IllegalArgumentException if the data returned by
     * the specified iterator contains no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws IllegalArgumentException if
     * {@code overlap < 0 || Float.isFinite(overlap) == false}
     * @throws IllegalArgumentException if
     * {@code window <= overlap || Float.isFinite(window) == false}
     * @throws NullPointerException if
     * {@code it == null || ped == null || genMap == null}
     */
    public static TargetData targetData(Par par,
            Supplier<SampleFileIt<GTRec>> supplier) {
        GeneticMap genMap = GeneticMap.geneticMap(par.map(), par.chromInt());
        WindowIt<GTRec> targWindow = WindowIt.newInstance(supplier, genMap,
                par.window(), par.overlap());
        TargetData targetData = new TargetData(par, targWindow);
        assert targetData.canAdvanceWindow();
        targetData.advanceWindowCm();
        return targetData;
    }

    private TargetData(Par par, WindowIt<GTRec> targWindow) {
        this.ped = new Pedigree(targWindow.samples(), par.ped());
        this.targWindowIt = targWindow;
        this.currentWindow = null;
        this.recs = new GTRec[0];
        this.gt = targGT(targWindow.samples(), recs, ped);
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
        return targWindowIt.genMap();
    }

    @Override
    public boolean lastWindowOnChrom() {
        return currentWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return targWindowIt.hasNext();
    }

    @Override
    public void advanceWindowCm() {
        currentWindow = targWindowIt.next();
        cumMarkerCnt += (currentWindow.nMarkers() - currentWindow.prevOverlap());
        recs = currentWindow.recList().toArray(new GTRec[0]);
        gt = targGT(targWindowIt.samples(), recs, ped);
        markerIndices = new MarkerIndices(currentWindow.prevOverlap(),
                currentWindow.nextOverlap(), gt.nMarkers());
        ++window;
    }

    @Override
    public IntArray[][] carriers(int maxCarriers) {
        return Arrays.stream(recs)
                .parallel()
                .map(rec -> carriers(rec, maxCarriers))
                .toArray(IntArray[][]::new);
    }

    private static IntArray[] carriers(GTRec rec, int maxCarriers) {
        int nAlleles = rec.marker().nAlleles();
        IntList[] carriers = IntStream.range(0, nAlleles)
                .mapToObj(i -> new IntList(16))
                .toArray(IntList[]::new);
        for (int s=0, n=rec.nSamples(); s<n; ++s) {
            int a1 = rec.allele1(s);
            int a2 = rec.allele2(s);
            if (a1>=0 && carriers[a1].size()<=maxCarriers) {
                carriers[a1].add(s);
            }
            if (a2>=0 && a2!=a1 && carriers[a2].size()<=maxCarriers) {
                carriers[a2].add(s);
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

    @Override
    public int nTargMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public int nMarkers() {
        return gt.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return cumMarkerCnt;
    }

    @Override
    public GT targGT() {
        return gt;
    }

    @Override
    public RefGT refGT() {
        // no reference genotype data to return
        return null;
    }

    @Override
    public RefGT restrictRefGT() {
        // no reference genotype data to return
        return null;
    }

    @Override
    public void close() {
       targWindowIt.close();
    }

    @Override
    public MarkerIndices markerIndices() {
        return markerIndices;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("vcf.NonRefData");
        return sb.toString();
    }
}

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

import blbutil.Utilities;
import ints.IndexArray;
import ints.IntArray;
import ints.IntList;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.Marker;
import vcf.MarkerMap;

/**
 * <p>Class {@code CodedSteps} divides phased genotype data
 * into non-overlapping intervals (the steps), indexes the unique
 * allele sequences in each interval, and stores a map of haplotype
 * index to allele sequence index for each interval.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CodedSteps {

    private final int nMarkers;
    private final int nHaps;
    private final int nTargHaps;
    private final int[] stepStarts;
    private final IndexArray[] codedSteps;

    /**
     * Constructs a new {@code CodedSteps} instance from the specified data.
     * @param targGT the phased target genotype data
     * @param refGT the phased phased reference genotype data or {@code null}
     * if there is no reference data
     * @param map the genetic map
     * @param step the step length in cM
     * @param scaleFactor factor by which to scale the number of steps
     * @param seed the random seed
     * @throws IllegalArgumentException if
     * {@code map.genDist().size()!=targGT.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code  refGT != null && targGT.markers().equals(refGT.markers()) == false}
     * @throws IllegalArgumentException if
     * {@code step <= 0.0 || Double.isFinite(step) == false}
     * @throws NullPointerException if {@code targGT == null || map == null}
     */
    public CodedSteps(GT targGT, GT refGT, MarkerMap map, double step,
            float scaleFactor, long seed) {
        if (map.genDist().size()!=targGT.nMarkers()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (refGT!=null && targGT.markers().equals(refGT.markers())==false) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (step<=0.0 || Double.isFinite(step)==false) {
            throw new IllegalArgumentException(String.valueOf(step));
        }
        Random rand = new Random(seed);
        HapData hapData = new HapData(targGT, refGT);
        this.nMarkers = hapData.nMarkers();
        this.nHaps = hapData.nHaps();
        this.nTargHaps = targGT.nHaps();
        int[] stepStarts0 = stepStarts(map, step, rand);
        IndexArray[] codedSteps0 = IntStream.range(0, stepStarts0.length)
                .parallel()
                .mapToObj(j -> codeStep(hapData, stepStarts0, j))
                .toArray(IndexArray[]::new);
        if (scaleFactor==1.0f) {
            this.stepStarts = stepStarts0;
            this.codedSteps = codedSteps0;
        }
        else {
            int[] resizedIndices = resizedIndices(stepStarts0.length,
                    scaleFactor, rand);
            this.stepStarts = Arrays.stream(resizedIndices)
                    .parallel()
                    .map(j -> stepStarts0[j])
                    .toArray();
            this.codedSteps = Arrays.stream(resizedIndices)
                    .parallel()
                    .mapToObj(j -> codedSteps0[j])
                    .toArray(IndexArray[]::new);
        }
    }

    private static int[] resizedIndices(int size, float scaleFactor, Random rand) {
        int newLength = (int) Math.ceil(size * scaleFactor);
        if (newLength < 40) {
            newLength = 40;
        }
        int quotient = newLength / size;
        int remainder = newLength - quotient*size;

        int[] indices = IntStream.range(0, size).toArray();
        int[] resizedIndices = new int[newLength];
        int start = 0;
        for (int j=0; j<quotient; ++j) {
            System.arraycopy(indices, 0, resizedIndices, start, indices.length);
            start += indices.length;
        }
        if (remainder>0) {
            Utilities.shuffle(indices, remainder, rand);
            System.arraycopy(indices, 0, resizedIndices, start, remainder);
        }
        Arrays.sort(resizedIndices);
        return resizedIndices;
    }

    private static int[] stepStarts(MarkerMap map, double step, Random rand) {
        double[] pos = map.genPos().toArray();
        IntList indices = new IntList(pos.length/10);

        indices.add(0);
        double nextPos =  pos[0] + rand.nextDouble()*step;
        int index = nextIndex(pos, 0, nextPos);
        while (index < pos.length) {
            indices.add(index);
            nextPos = pos[index] + step;
            index = nextIndex(pos, index, nextPos);
        }
        return indices.toArray();
    }

    private static int nextIndex(double[] pos, int start, double targetPos) {
        int nextIndex = Arrays.binarySearch(pos, start, pos.length, targetPos);
        return (nextIndex<0) ? -nextIndex-1 : nextIndex;
    }

    private static IndexArray codeStep(HapData hapData, int[] starts,
            int startsIndex) {
        int start = starts[startsIndex];
        int end = (startsIndex+1) < starts.length ?  starts[startsIndex+1]
                : hapData.nMarkers();
        if ((end-start)==1) {
            return indexArray(hapData, start);
        }
        else {
            return indexArray(hapData, start, end);
        }
    }

    private static IndexArray indexArray(HapData hapData, int start,
            int end) {
        int nTargHaps = hapData.nTargHaps();
        int nHaps = hapData.nHaps();
        int[] hapToSeq = IntStream.range(0, nHaps).map(i -> 1).toArray();
        int seqCnt = 2; // seq 0 is reserved for sequences not found in target
        for (int m=start; m<end; ++m) {
            int nAlleles = hapData.marker(m).nAlleles();
            int[] seqMap = new int[seqCnt*nAlleles];
            seqCnt = 1;
            for (int h=0; h<nTargHaps; ++h) {
                int allele = hapData.allele(m, h);
                int index = nAlleles*hapToSeq[h] + allele;
                if (seqMap[index] == 0) {
                    seqMap[index] = seqCnt++;
                }
                hapToSeq[h] = seqMap[index];
            }
            for (int h=nTargHaps; h<nHaps; ++h) {
                if (hapToSeq[h] != 0) {
                    int allele = hapData.allele(m, h);
                    int index = hapToSeq[h]*nAlleles + allele;
                    hapToSeq[h] = seqMap[index];
                }
            }
        }
        IntArray intArray = IntArray.create(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }

    private static IndexArray indexArray(HapData hapData, int marker) {
        int nAlleles = hapData.marker(marker).nAlleles();
        IntArray intArray =  new IntArray() {
            @Override
            public int size() {
                return hapData.nHaps();
            }

            @Override
            public int get(int index) {
                return hapData.allele(marker, index);
            }
        };
        return new IndexArray(intArray, nAlleles);
    }

    /**
     * Returns the number of markers
     * @return the number of markeres
     */
    public int nMarkers() {
        return nMarkers;
    }

    /**
     * Returns the number of target and reference haplotypes.
     * @return the number of target and reference haplotypes
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Returns the number of target haplotypes.
     * @return the number of target haplotypes
     */
    public int nTargHaps() {
        return nTargHaps;
    }

    /**
     * Returns the number of steps.
     * @return the number of steps
     */
    public int nSteps() {
        return stepStarts.length;
    }

    /**
     * Returns the first marker index in the specified step.
     * @param step a step index
     * @return the first marker index in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int stepStart(int step) {
        return stepStarts[step];
    }

    /**
     * Returns the last marker index (exclusive) in the specified step.
     * @param step a step index
     * @return the lastt marker index (exclusive) in the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int stepEnd(int step) {
        return (step+1 < stepStarts.length) ? stepStarts[step+1] : nMarkers;
    }

    /**
     * Returns a map from haplotype index to allele sequence index
     * for the specified step
     * @param step a step index
     * @return a map from haplotype index to allele sequence index
     * for the specified step
     * @throws IllegalArgumentException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public IndexArray get(int step) {
        return codedSteps[step];
    }

    private static class HapData {
        private final int nTargHaps;
        private final int nHaps;
        private final GT targGT;
        private final GT refGT;

        public HapData(GT targGT, GT refGT) {
            this.nTargHaps = targGT.nHaps();
            this.nHaps = nTargHaps + (refGT==null ? 0 : refGT.nHaps());
            this.targGT = targGT;
            this.refGT = refGT;
        }

        public int nTargHaps() {
            return targGT.nHaps();
        }

        public int nHaps() {
            return nHaps;
        }

        public int nMarkers() {
            return targGT.nMarkers();
        }

        public Marker marker(int m) {
            return targGT.marker(m);
        }

        public int allele (int marker, int hap) {
            if (hap>=nHaps) {
                throw new IndexOutOfBoundsException(String.valueOf(hap));
            }
            if (hap<nTargHaps) {
                return targGT.allele(marker, hap);
            }
            else {
                return refGT.allele(marker, hap - nTargHaps);
            }
        }
    }
}

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

package imp;

import ints.IndexArray;
import ints.IntArray;
import ints.IntList;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code CodedSteps} divides phased genotype data
 * into non-overlapping intervals (the steps), indexes the unique
 * allele sequences in each interval, and stores a map of haplotype
 * index to allele sequence index for each interval.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CodedSteps {

    private final ImpData impData;
    private final int[] stepStarts;
    private final IndexArray[] codedSteps;

    /**
     * Constructs a new {@code CodedSteps} instance from the specified data.
     * @param impData input data for genotype imputation
     * @throws NullPointerException if {@code impData == null}
     */
    public CodedSteps(ImpData impData) {
        this.impData = impData;
        this.stepStarts = stepStarts(impData);
        this.codedSteps = IntStream.range(0, stepStarts.length)
                .parallel()
                .mapToObj(j -> codeStep(impData, stepStarts, j))
                .toArray(IndexArray[]::new);
    }

    private static int[] stepStarts(ImpData impData) {
        double[] pos = impData.pos();
        double step = impData.par().imp_step();
        IntList indices = new IntList(pos.length/10);
        indices.add(0);
        double nextPos =  pos[0] + step/2;  // make first step be half-length
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

    private static IndexArray codeStep(ImpData impData, int[] starts,
            int startIndex) {
        int nRefHaps = impData.nRefHaps();
        int nHaps = impData.nHaps();
        int[] hapToSeq = IntStream.range(0, nHaps).map(i -> 1).toArray();
        int start = starts[startIndex];
        int end = (startIndex+1) < starts.length ?  starts[startIndex+1]
                : impData.nClusters();

        int seqCnt = 2; // seq 0 is reserved for sequences not found in target
        for (int m=start; m<end; ++m) {
            IndexArray h2s = impData.hapToSeq(m);
            IntArray codedHaps = h2s.intArray();
            int nAlleles = h2s.valueSize();
            int[] seqMap = new int[seqCnt*nAlleles];

            seqCnt = 1;
            for (int h=nRefHaps; h<nHaps; ++h) {
                int index = nAlleles*hapToSeq[h] + codedHaps.get(h);
                if (seqMap[index] == 0) {
                    seqMap[index] = seqCnt++;
                }
                hapToSeq[h] = seqMap[index];
            }
            for (int h=0; h<nRefHaps; ++h) {
                if (hapToSeq[h] != 0) {
                    int index = hapToSeq[h]*nAlleles + codedHaps.get(h);
                    hapToSeq[h] = seqMap[index];
                }
            }
        }
        IntArray intArray = IntArray.create(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }

    /**
     * Return the input data for genotype imputation used to construct
     * {@code this}.
     * @return the input data for genotype imputation used to construct
     * {@code this}
     */
    public ImpData impData() {
        return impData;
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
}

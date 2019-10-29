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

import ints.IntArray;
import ints.IntList;
import blbutil.Utilities;
import ints.IndexArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;
import main.Par;

/**
 * <p>Class {@code ImpIbs} identifies haplotypes that share a long
 * IBS segment with a specified haplotype.
 * </p>
 * <p>Instances of {@code ImpIbs} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImpIbs {

    private final ImpData impData;
    private final int nRefHaps;
    private final long seed;
    private final int nStates;
    private final int nSteps;
    private final int nHapsPerStep;

    private final CodedSteps codedSteps;
    private final int[][][] ibsHaps; //[window][targ hap][ibs_set]

    /**
     * Constructs a new {@code ImpIbs} object from the specified data.
     * @param impData the input data for genotype imputation
     *
     * @throws NullPointerException if {@code impData == null}
     */
    public ImpIbs(ImpData impData) {
        Par par = impData.par();
        this.impData = impData;
        this.seed = par.seed();
        this.nRefHaps = impData.nRefHaps();
        this.nStates = par.imp_states();
        this.nSteps = par.imp_nsteps();
        int nStepsPerSegment = Math.round(par.imp_segment()/par.imp_step());
        this.nHapsPerStep = (par.imp_states() / nStepsPerSegment);

        this.codedSteps = new CodedSteps(impData);
        this.ibsHaps = IntStream.range(0, codedSteps.nSteps())
                .parallel()
                .mapToObj(j -> getIbsHaps(codedSteps, j))
                .toArray(int[][][]::new);
    }

    private int[][] getIbsHaps(CodedSteps codedSteps, int index) {
        int nTargHaps = impData.nTargHaps();
        int[][] results = new int[nTargHaps][];
        int nStepsToMerge = Math.min(nSteps, codedSteps.nSteps() - index);
        List<IntList> children = initPartition(codedSteps.get(index));
        List<IntList> nextParents = new ArrayList<>(children.size());
        initUpdateResults(children, nextParents, results);
        for (int i=1; i<nStepsToMerge; ++i) {
            int initCapacity = Math.min(nTargHaps, 2*nextParents.size());
            List<IntList> parents = nextParents;
            nextParents = new ArrayList<>(initCapacity);
            IndexArray codedStep = codedSteps.get(index+i);
            for (int j=0, nj=parents.size(); j<nj; ++j) {
                IntList parent = parents.get(j);
                children = partition(parent, codedStep);
                updateResults(parent, children, nextParents, results);
            }
        }
        finalUpdateResults(nextParents, results);
        return results;
    }

    private List<IntList> initPartition(IndexArray codedStep) {
        IntList[] list = new IntList[codedStep.valueSize()];
        IntArray hap2Seq = codedStep.intArray();
        int nHaps = hap2Seq.size();
        List<IntList> children = new ArrayList<>();
        for (int h=nRefHaps; h<nHaps; ++h) {
            int seq = hap2Seq.get(h);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
        }
        for (int h=0; h<nHaps; ++h) {
            int seq = hap2Seq.get(h);
            if (list[seq]!=null) {
                list[seq].add(h);
            }
        }
        return children;
    }

    private List<IntList> partition(IntList parent, IndexArray codedStep) {
        IntList[] list = new IntList[codedStep.valueSize()];
        IntArray hap2Seq = codedStep.intArray();
        int nParentHaps = parent.size();
        List<IntList> children = new ArrayList<>();
        int targStart = insPt(parent, nRefHaps);
        for (int k=targStart; k<nParentHaps; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]==null) {
                list[seq] = new IntList();
                children.add(list[seq]);
            }
        }
        for (int k=0; k<nParentHaps; ++k) {
            int hap = parent.get(k);
            int seq = hap2Seq.get(hap);
            if (list[seq]!=null) {
                list[seq].add(hap);
            }
        }
        return children;
    }

    private void initUpdateResults(List<IntList> children,
            List<IntList> nextParents, int[][] result) {
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList hapList = children.get(j);
            int nRef = insPt(hapList, nRefHaps);
            if (nRef <= nHapsPerStep) {
                setResult(hapList, nRef, hapList.copyOf(nRef), result);
            }
            else {
                nextParents.add(hapList);
            }
        }
    }

    private void updateResults(IntList parent, List<IntList> children, List<IntList> nextIbs, int[][] results) {
        for (int k=0, nk=children.size(); k<nk; ++k) {
            IntList child = children.get(k);
            int nChildRef = insPt(child, nRefHaps);
            if (nChildRef <= nHapsPerStep) {
                int[] ibsList = ibsHaps(parent, child, nChildRef);
                setResult(child, nChildRef, ibsList, results);
                child.clear();
            }
            else {
                nextIbs.add(child);
            }
        }
    }

    private int[] ibsHaps(IntList parent, IntList child, int nChildRef) {
        IntList combined = new IntList(nHapsPerStep);
        for (int j=0; j<nChildRef; ++j) {
            combined.add(child.get(j));
        }
        int size = nHapsPerStep - nChildRef;
        Random rand = new Random(seed + parent.get(0));
        IntList uniqToParent = uniqToParent(parent, child, nChildRef);
        int[] randSubset = randomSubset(uniqToParent, size, rand);
        for (int i : randSubset) {
            combined.add(i);
        }
        int[] ia = combined.toArray();
        Arrays.sort(ia);    // is this needed?
        return ia;
    }

    private IntList uniqToParent(IntList parent, IntList child, int nChildRef) {
        int nChildRefM1 = nChildRef - 1;
        int nParentRef = insPt(parent, nRefHaps);
        IntList uniqToParent = new IntList(parent.size());
        int c = 0;
        int cVal = child.get(c);
        for (int p=0; p<nParentRef; ++p) {
            int pVal = parent.get(p);
            while (cVal < pVal && c < nChildRefM1) {
                cVal = child.get(++c);
            }
            if (pVal != cVal) {
                uniqToParent.add(pVal);
            }
        }
        return uniqToParent;
    }

    private static int[] randomSubset(IntList list, int size, Random rand) {
        if (list.size() < size) {
            size = list.size();
        }
        int[] ia = list.toArray();
        for (int j=0; j<size; ++j) {
            int x = rand.nextInt(ia.length-j);
            int tmp = ia[j];
            ia[j] = ia[j+x];
            ia[j+x] = tmp;
        }
        return Arrays.copyOf(ia, size);
    }

    private void finalUpdateResults(List<IntList> children, int[][] results) {
        for (int j=0, n=children.size(); j<n; ++j) {
            IntList child = children.get(j);
            int nRef = insPt(child, nRefHaps);
            int[] ibsList = child.copyOf(nRef);
            if (nHapsPerStep < ibsList.length) {
                Random rand = new Random(seed + child.get(0));
                Utilities.shuffle(ibsList, rand);
                ibsList = Arrays.copyOf(ibsList, nHapsPerStep);
                Arrays.sort(ibsList);
            }
            setResult(child, nRef, ibsList, results);
        }
    }

    private void setResult(IntList child, int firstTargIndex,
            int[] ibsHaps, int[][] result) {
        for (int j=firstTargIndex, nj=child.size(); j<nj; ++j) {
            result[child.get(j) - nRefHaps] = ibsHaps;
        }
    }

    private static int insPt(IntList il, int nRefHaps) {
        int index = il.binarySearch(nRefHaps);
        return index >= 0 ? index : -index - 1;
    }

    /**
     * Returns an array containing reference haplotype indices that
     * that are IBS with the specified target haplotype in an interval
     * beginning with the specified step. The returned array will contain fewer
     * than {@code this.nHapsPerStep()} haplotypes if the number of reference
     * haplotypes that are IBS with specified target haplotype in the specified
     * step is less than {@code this.nHapsPerStep()}.
     * @param hap a haplotype index
     * @param step a step index
     * @return an array containing reference haplotype indices that
     * that are IBS with the specified target haplotype
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.hapPairs().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code step < 0 || step >= this.nSteps()}
     */
    public int[] ibsHaps(int hap, int step) {
        return ibsHaps[step][hap].clone();
    }

    /**
     * Return the data for genotype imputation in the marker window.
     * @return the data for genotype imputation in the marker window
     */
    public ImpData impData() {
        return impData;
    }

    /**
     * Returns the coded steps.  The coded steps stores the starting marker
     * for each step along with the index of the allele sequence carried by
     * each haplotype in each step.
     * @return the coded steps
     */
    public CodedSteps codedSteps() {
        return codedSteps;
    }
}

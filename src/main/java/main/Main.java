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

import vcf.GeneticMap;
import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import imp.ImpData;
import imp.ImpLS;
import java.io.File;
import java.util.HashSet;
import java.util.Locale;
import java.util.Set;
import vcf.AllData;
import bref.Bref3It;
import imp.StateProbs;
import ints.LongArray;
import java.util.concurrent.atomic.AtomicReferenceArray;
import java.util.function.Supplier;
import phase.FixedPhaseData;
import vcf.VcfIt;
import vcf.Data;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.FilterUtil;
import vcf.GT;
import vcf.Markers;
import vcf.TargetData;
import vcf.RefIt;
import vcf.GTRec;
import vcf.HapsGT;
import vcf.RefGTRec;

/**
 * Class {@code Main} is the entry class for the Beagle program.
 * See {@code Par.usage()} and online program documentation for usage
 * instructions.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Main {

    /**
     * The program name and version.
     */
    public static final String VERSION = "(version 5.1)";
    public static final String PROGRAM = "beagle.21Sep19.ec3.jar";
    public static final String COMMAND = "java -jar beagle.21Sep19.ec3.jar";

    /**
     * The copyright string.
     */
    public static final String COPYRIGHT = "Copyright (C) 2014-2019 Brian L. Browning";

    /**
     * The program name and a brief help message.
     */
    public static final String SHORT_HELP = Main.PROGRAM + " " + VERSION
            + Const.nl + Main.COPYRIGHT
            + Const.nl + "Enter \"java -jar beagle.21Sep19.ec3.jar\" to "
            + "list command line argument";

    private final Par par;
    private final GeneticMap genMap;
    private final Data data;
    private final RunStats runStats;
    private final WindowWriter windowWriter;

    /**
     * Entry point to Beagle program.  See {@code Parameters.usage()} and
     * online program documentation for usage instructions.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        if (args.length==0) {
            System.out.println(PROGRAM + " " + VERSION);
            System.out.println(COPYRIGHT);
            System.out.println(Par.usage());
            System.exit(0);
        }
        Par par = parameters(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nthreads()));
        RunStats runStats = new RunStats(par);
        runStats.printStartInfo();

        try (Data data = data(par, runStats);
                WindowWriter winOut = new WindowWriter(data.targGT().samples(),
                        par.out())) {
            Main main = new Main(par, data, winOut, runStats);
            main.phaseData();
            runStats.printSummaryAndClose(data.nTargMarkersSoFar(),
                    data.nMarkersSoFar());
        }
    }

    private Main(Par par, Data data, WindowWriter windWriter, RunStats runStats) {
        assert par!=null;
        assert data!=null;
        assert windWriter!=null;
        assert runStats!=null;
        this.par = par;
        this.genMap = data.genMap();
        this.data = data;
        this.runStats = runStats;
        this.windowWriter = windWriter;
    }

    private static Data data(Par par, RunStats runStats) {
        Filter<String> sFilter = FilterUtil.sampleFilter(par.excludesamples());
        Filter<Marker> mFilter = FilterUtil.markerFilter(par.excludemarkers());
        ChromInterval chromInterval = par.chromInt();
        if (par.ref()==null) {
            Supplier<SampleFileIt<GTRec>> targSupplier =
                    () -> targIt(par, mFilter, sFilter, chromInterval);
            return TargetData.targetData(par, targSupplier);
        }
        else {
            SampleFileIt<GTRec> targIt = targIt(par, mFilter, sFilter, chromInterval);
            Supplier<SampleFileIt<RefGTRec>> refSupplier = refSupplier(par,
                    mFilter, sFilter, chromInterval, runStats);
            return AllData.allData(refSupplier, targIt, par);
        }
    }

    /*
     * Phases the data, imputes ungenotyped markers, and performed IBD segment
     * detection.
     */
    private void phaseData() {
        runStats.printSampleSummary(data);
        GT overlap = null;
        int window = 0;
        do {
            if (++window > 1) {
                data.advanceWindowCm();
            }
            runStats.printWindowUpdate(data);
            GT phasedTarg = phasedTarg(overlap, (par.seed() + window));
            printOutput(data, phasedTarg);
            overlap = phasedOverlap(data, phasedTarg);
        } while (data.canAdvanceWindow());
    }

    private GT phasedTarg(GT overlap, long seed) {
        if (data.targGT().isPhased()) {
            return data.targGT();
        }
        else {
            FixedPhaseData fpd = new FixedPhaseData(par, genMap, data, overlap);
            MainHelper mh = new MainHelper(par, runStats, seed);
            return mh.phase(fpd);
        }
    }

    private void printOutput(Data data, GT phasedTarg) {
        assert par.gt()!=null;
        int nThreads = par.nthreads();
        if (data.nMarkers()==data.targGT().nMarkers() || par.impute()==false) {
            int targStart = data.markerIndices().prevTargSplice();
            int targEnd = data.markerIndices().nextTargSplice();
            windowWriter.printPhased(phasedTarg, targStart, targEnd, nThreads);
        }
        else {
            long t0 = System.nanoTime();
            ImpData impData = new ImpData(par, data, phasedTarg, genMap);
            AtomicReferenceArray<StateProbs> stateProbs = ImpLS.stateProbs(impData);
            int refStart = data.markerIndices().prevSplice();
            int refEnd = data.markerIndices().nextSplice();
            windowWriter.printImputed(impData, refStart, refEnd, stateProbs);
            runStats.imputationNanos(System.nanoTime() - t0);
            runStats.printImputationUpdate();
        }
    }

    private GT phasedOverlap(Data data, GT phasedTarg) {
        assert phasedTarg.isPhased();
        int nextOverlap = data.markerIndices().nextTargOverlap();
        int nextSplice = data.markerIndices().nextTargSplice();
        int nMarkers = nextSplice - nextOverlap;
        if (nMarkers==0) {
            return null;
        }
        else {
            int nHaps = phasedTarg.nHaps();
            Markers markers = phasedTarg.markers().restrict(nextOverlap, nextSplice);
            LongArray[] haps = new LongArray[nHaps];
            int[] hap = new int[nMarkers];
            for (int h=0; h<nHaps; ++h) {
                for (int m=0; m<nMarkers; ++m) {
                    hap[m] = phasedTarg.allele(nextOverlap + m, h);
                }
                haps[h] = markers.allelesToBits(hap);
            }
            return new HapsGT(markers, phasedTarg.samples(), haps);
        }
    }

    private static SampleFileIt<GTRec> targIt(Par par,
            Filter<Marker> markerFilter, Filter<String> sampleFilter,
            ChromInterval chromInterval) {
        FileIt<String> it = InputIt.fromGzipFile(par.gt());
        SampleFileIt<GTRec> targIt = VcfIt.create(it, sampleFilter,
                markerFilter,  VcfIt.toBitSetGT);
        if (chromInterval!=null) {
            targIt = new IntervalVcfIt<>(targIt, chromInterval);
        }
        return targIt;
    }

    private static Supplier<SampleFileIt<RefGTRec>> refSupplier(Par par,
            Filter<Marker> mFilter, Filter<String> sFilter,
            ChromInterval chromInt, RunStats runStats) {
        Filter<Marker> mFilter2 = updateFilter(par, mFilter, sFilter, chromInt);
        return () -> {
            SampleFileIt<RefGTRec> refIt;
            String filename = par.ref().toString();
            if (filename.endsWith(".bref")) {
                String s = Const.nl + "ERROR: bref format (.bref) is not supported"
                         + Const.nl + "       Reference files should be in bref3 format (.brer3)" ;
                Utilities.exit(s);
            }
            if (filename.endsWith(".bref3")) {
                refIt = new Bref3It(par.ref(), mFilter2);
            }
            else {
                if (filename.endsWith(".vcf")==false
                        && filename.endsWith(".vcf.gz")==false) {
                    runStats.println(Const.nl
                            + "WARNING: unrecognized reference file type "
                            + "(expected \".bref3\", \".vcf\", or \".vcf.gz\")"
                            + Const.nl);
                }
                FileIt<String> it = InputIt.fromGzipFile(par.ref());
                refIt = RefIt.create(it, sFilter, mFilter2,
                        RefIt.MAX_EM_BUFFER_SIZE);
            }
            if (chromInt!=null) {
                refIt = new IntervalVcfIt<>(refIt, chromInt);
            }
            return refIt;
        } ;
    }

    private static Filter<Marker> updateFilter(Par par, Filter<Marker> mFilter,
            Filter<String> sFilter, ChromInterval chromInt) {
        if (par.impute() && par.gt()!=null) {
            return mFilter;
        }
        else {
            Set<Marker> includedMarkers = new HashSet<>(50000);
            try (SampleFileIt<GTRec> vcfIt = targIt(par, mFilter, sFilter,
                    chromInt)) {
                while (vcfIt.hasNext()) {
                    includedMarkers.add(vcfIt.next().marker());
                }
            }
            return Filter.includeFilter(includedMarkers);
        }
    }

    /*
     * Checks that certain parameters are consistent, and prints error
     * message and exits if parameters are inconsistent.
     *
     * @param args the command line arguments.
     */
    private static Par parameters(String[] args) {
        // warnings are printed in RunStats.startInfo() method
        Par par = new Par(args);
        checkOutputPrefix(par);
        if (1.1*par.overlap() >= par.window()) {
            String s = SHORT_HELP + Const.nl
                    + Const.nl + "ERROR: The \"window\" parameter must be at least "
                    + "1.1 times the \"overlap\" parameter"
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
        return par;
    }

    private static void checkOutputPrefix(Par par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(Par.usage() + s);
        }

        File vcfOut = new File(par.out() + ".vcf.gz");
        if (vcfOut.equals(par.ref())) {
            String s = "ERROR: VCF output file equals input file: " + par.ref();
            Utilities.exit(Par.usage() + s);
        }
        if (vcfOut.equals(par.gt())) {
            String s = "ERROR: VCF output file equals input file: " + par.gt();
            Utilities.exit(Par.usage() + s);
        }
    }
}

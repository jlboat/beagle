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

import beagleutil.ChromInterval;
import blbutil.Const;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code Parameters} represents the parameters for a Beagle analysis.
 * </p>
 * <p>Instances of class {@code Parameters} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Par {

    private final String[] args;

    // data parameters
    private final File gt;
    private final File ref;
    private final String out;
    private final File ped;
    private final File map;
    private final ChromInterval chromInt;
    private final File excludesamples;
    private final File excludemarkers;

    // phasing parameters
    private final int burnin;
    private final int iterations;
    private final int phase_states;
    private final float phase_step;
    private final float rare;

    // imputation parameters
    private final boolean impute;
    private final int imp_states;
    private final float imp_segment;
    private final float imp_step;
    private final int imp_nsteps;
    private final float cluster;
    private final boolean ap;
    private final boolean gp;

    // general parameters
    private final float ne;
    private final float err;
    private final float window;
    private final float overlap;
    private final long seed;
    private final int nthreads;
    private final float buffer;

    // undocumented parameters
    private final File truth;

    // default phasing parameters
    private static final int D_BURNIN = 6;
    private static final int D_ITERATIONS = 12;
    private static final int D_PHASE_STATES = 280;
    private static final float D_PHASE_STEP = 0.006f;
    private static final float D_RARE = 0.0015f;

    // default imputation parameters
    private static final boolean D_IMPUTE = true;
    private static final int D_IMP_STATES = 1600;
    private static final float D_IMP_SEGMENT = 6.0f;
    private static final float D_IMP_STEP = 0.1f;
    private static final int D_IMP_NSTEPS = 7;
    private static final float D_CLUSTER = 0.005f;
    private static final boolean D_AP = false;
    private static final boolean D_GP = false;

    // default general parameters
    private static final int D_NE = 1_000_000;
    private static final float D_ERR = -Float.MIN_VALUE;
    private static final float D_WINDOW = 40.0f;
    private static final float D_OVERLAP = 4.0f;
    private static final int D_SEED = -99999;
    private static final int D_NTHREADS = Integer.MAX_VALUE;
    private static final float D_BUFFER = 0.6f;


    /**
     * Constructs a new {@code Parameters} instance from the specified
     * command line arguments.
     * @param args the command line arguments
     * @throws IllegalArgumentException if a command line argument
     * is incorrectly specified
     * @throws NumberFormatException if a numeric value for a parameter
     * is incorrectly specified
     * @throws NullPointerException if {@code args ==  null}
     */
    public Par(String[] args) {
        int IMAX = Integer.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;

        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        // data input/output parameters
        gt = Validate.getFile(
                Validate.stringArg("gt", argsMap, false, null, null));
        ref = Validate.getFile(
                Validate.stringArg("ref", argsMap, false, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        ped = Validate.getFile(
                Validate.stringArg("ped", argsMap, false, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, false, null, null));
        chromInt = parseChromInt(Validate.stringArg("chrom", argsMap, false, null, null));
        excludesamples = Validate.getFile(
                Validate.stringArg("excludesamples", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));

        // phasing parameters
        burnin = Validate.intArg("burnin", argsMap, false, D_BURNIN, 1, IMAX);
        iterations = Validate.intArg("iterations", argsMap, false, D_ITERATIONS, 1, IMAX);
        phase_states = Validate.intArg("phase-states", argsMap, false, D_PHASE_STATES, 1, IMAX);
        phase_step = Validate.floatArg("phase-step", argsMap, false, D_PHASE_STEP, FMIN, FMAX);
        rare = Validate.floatArg("rare", argsMap, false, D_RARE, FMIN, 0.5f);

        // imputation parameters
        impute = Validate.booleanArg("impute", argsMap, false, D_IMPUTE);
        imp_states = Validate.intArg("imp-states", argsMap, false, D_IMP_STATES, 1, IMAX);
        imp_segment = Validate.floatArg("imp-segment", argsMap, false, D_IMP_SEGMENT, FMIN, FMAX);
        imp_step = Validate.floatArg("imp-step", argsMap, false, D_IMP_STEP, FMIN, FMAX);
        imp_nsteps = Validate.intArg("imp-nsteps", argsMap, false, D_IMP_NSTEPS, 1, IMAX);
        cluster = Validate.floatArg("cluster", argsMap, false, D_CLUSTER, 0.0f, FMAX);
        ap = Validate.booleanArg("ap", argsMap, D_AP, false);
        gp = Validate.booleanArg("gp", argsMap, D_GP, false);

        // general parameters
        ne = Validate.floatArg("ne", argsMap, false, D_NE, FMIN, FMAX);
        err = Validate.floatArg("err", argsMap, false, D_ERR, -FMIN, FMAX);
        window = Validate.floatArg("window", argsMap, false, D_WINDOW, FMIN, IMAX);
        overlap = Validate.floatArg("overlap", argsMap, false, D_OVERLAP, FMIN, IMAX);
        buffer = Validate.floatArg("overlap", argsMap, false, D_BUFFER, FMIN, IMAX);
        seed = Validate.longArg("seed", argsMap, false, D_SEED, LMIN, LMAX);
        nthreads = modNthreads(Validate.intArg("nthreads", argsMap, false, D_NTHREADS, 1, IMAX));

        // undocumented parameters
        truth =  Validate.getFile(Validate.stringArg("truth", argsMap, false, null, null));

        Validate.confirmEmptyMap(argsMap);
    }

    /**
     * Returns the command line arguments.
     * @return the command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a description of the Beagle command line arguments.
     * @return a description of the Beagle command line arguments
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Usage: " + Main.COMMAND + " [arguments]" + nl
                + nl
                + "data parameters ..." + nl
                + "  gt=<VCF file: use GT field>                        (optional)" + nl
                + "  ref=<bref3 or VCF file with phased genotypes>      (optional)" + nl
                + "  out=<output file prefix>                           (required)" + nl
//                + "  ped=<linkage format pedigree file>                 (optional)" + nl
                + "  map=<PLINK map file with cM units>                 (optional)" + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)" + nl
                + "  excludesamples=<file with 1 sample ID per line>    (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>    (optional)" + nl + nl

                + "phasing parameters ..." + nl
                + "  burnin=<number of burnin iterations>               (default=" + D_BURNIN + ")" + nl
                + "  iterations=<number of phasing iterations>          (default=" + D_ITERATIONS + ")" + nl
                + "  phase-states=<model states for phasing>            (default=" + D_PHASE_STATES + ")" + nl + nl

                + "imputation parameters ..." + nl
                + "  impute=<impute ungenotyped markers (true/false)>   (default=" + D_IMPUTE + ")" + nl
                + "  imp-states=<model states for imputation>           (default=" + D_IMP_STATES + ")" + nl
                + "  imp-segment=<min haplotype segment length (cM)>    (default=" + D_IMP_SEGMENT + ")" + nl
                + "  imp-step=<IBS step length (cM)>                    (default=" + D_IMP_STEP + ")" + nl
                + "  imp-nsteps=<number of IBS steps>                   (default=" + D_IMP_NSTEPS + ")" + nl
                + "  cluster=<max cM in a marker cluster>               (default=" + D_CLUSTER + ")" + nl
                + "  ap=<print posterior allele probabilities>          (default=" + D_AP + ")" + nl
                + "  gp=<print posterior genotype probabilities>        (default=" + D_GP + ")" + nl + nl

                + "general parameters ..." + nl
                + "  ne=<effective population size>                     (default=" + D_NE + ")" + nl
                + "  err=<allele mismatch rate>                         (default: data dependent)" + nl
                + "  window=<window length in cM>                       (default=" + D_WINDOW + ")" + nl
                + "  overlap=<window overlap in cM>                     (default=" + D_OVERLAP + ")" + nl
                + "  seed=<random seed>                                 (default=" + D_SEED + ")" + nl
                + "  nthreads=<number of threads>                       (default: machine dependent)" + nl + nl;

    }

    private static ChromInterval parseChromInt(String str) {
        ChromInterval ci = ChromInterval.parse(str);
        if (str!=null && str.length()>0 && ci==null) {
            throw new IllegalArgumentException("Invalid chrom parameter: " + str);
        }
        return ci;
    }

    /**
     * Returns the nthreads parameter, which is equal to
     * {@code Runtime.getRuntime().availableProcessors()} if
     * {@code nthreads == Integer.MAX_VALUE}.
     * @return the nthreads parameter
     */
    private static int modNthreads(int nthreads) {
        if (nthreads==Integer.MAX_VALUE) {
            return Runtime.getRuntime().availableProcessors();
        }
        else {
            return nthreads;
        }
    }

    // data parameters

    /**
     * Returns the gt parameter or {@code null} if no gt parameter was
     * specified.
     * @return the gt parameter or {@code null} if no gt parameter was
     * specified
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the ref parameter or {@code null} if no ref parameter was
     * specified.
     * @return the ref parameter or {@code null} if no ref parameter was
     * specified
     */
    public File ref() {
        return ref;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the ped parameter or {@code null}
     * if no ped parameter was specified.
     *
     * @return the ped parameter or {@code null}
     * if no ped parameter was specified
     */
    public File ped() {
        return null;
//        return ped;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
     */
    public File map() {
        return map;
    }

    /**
     * Returns the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     *
     * @return the chromosome interval or {@code null} if no chrom
     * parameter was specified.
     */
    public ChromInterval chromInt() {
        return chromInt;
    }
    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    // phasing parameters

    /**
     * Returns the burnin parameter.
     * @return the burnin parameter
     */
    public int burnin() {
        return burnin;
    }

    /**
     * Returns the iterations parameter.
     * @return the iterations parameter
     */
    public int iterations() {
        return iterations;
    }

    /**
     * Returns the phase-states parameter.
     * @return the phase-states parameter
     */
    public int phase_states() {
        return phase_states;
    }

    /**
     * Return the ratio of the phase-states parameter and the default
     * phase-states parameter
     * @return the ratio of the phase-states parameter and the default
     * phase-states parameter
     */
    public float scaleFactor() {
        return (phase_states==D_PHASE_STATES) ? 1.0f
                : (float) phase_states / D_PHASE_STATES;
    }

    /**
     * Returns the phase-step parameter.
     * @return the phase-step parameter
     */
    public float phase_step() {
        return phase_step;
    }

    /**
     * Returns the rare parameter
     * @return the rare parameter
     */
    public float rare() {
        return rare;
    }

    // imputation parameters

    /**
     * Returns the impute parameter.
     * @return the impute parameter
     */
    public boolean impute() {
        return impute;
    }

    /**
     * Returns the imp-states parameter.
     * @return the imp-states parameter
     */
    public int imp_states() {
        return imp_states;
    }

    /**
     * Returns the imp-segment parameter.
     * @return the imp-segment parameter
     */
    public float imp_segment() {
        return imp_segment;
    }

    /**
     * Returns the imp-step parameter.
     * @return the imp-step parameter
     */
    public float imp_step() {
        return imp_step;
    }

    /**
     * Returns the imp-nsteps parameter.
     * @return the imp-nsteps parameter
     */
    public int imp_nsteps() {
        return imp_nsteps;
    }

    /**
     * Returns the cluster parameter.
     * @return the cluster parameter
     */
    public float cluster() {
        return cluster;
    }

    /**
     * Returns the ap parameter.
     * @return the ap parameter
     */
    public boolean ap() {
        return ap;
    }

    /**
     * Returns the gp parameter.
     * @return the gp parameter
     */
    public boolean gp() {
        return gp;
    }

    // general parameters

    /**
     * Returns the ne parameter
     * @return the ne parameter
     */
    public float ne() {
        return ne;
    }

    /**
     * Returns the allele mismatch parameter for the specified number of samples.
     * @param nHaps the number of reference and target haplotypes
     * @return the allele mismatch parameter for the specified number of samples
     */
    public float err(int nHaps) {
        if (nHaps <= 0) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        return (err>=0) ? err : liStephensErr(nHaps);
    }

    /**
     * Return an approximation to the allele mismatch probability suggested by
     * Li and Stephens.  The approximation uses a Riemann sum approximation
     * of the natural log function.
     * Refs:
     * Li N, Stephens M. Genetics 2003 Dec;165(4):2213-33
     * Marchini J, Howie B. Myers S, McVean G, Donnelly P. 2007:39(7)906-13.
     *
     * @param nHaps the number of haplotypes
     * @return the allele mismatch probability suggested by Li and Stepehens
     */
    private static float liStephensErr(int nHaps) {
        double theta = 1/((Math.log(nHaps) + 0.5));
        return (float) (theta/(2*(theta + nHaps)));
    }

    /**
     * Returns the window parameter.
     * @return the window parameter
     */
    public float window() {
        return window;
    }

    /**
     * Return the overlap parameter.
     * @return the overlap parameter.
     */
    public float overlap() {
        return overlap;
    }

    /**
     * Return the buffer parameter.
     * @return the buffer parameter.
     */
    public float buffer() {
        return buffer;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    // undocumented parameters

    /**
     * Returns the truth parameter
     * @return the truth
     */
    public File truth() {
        return truth;
    }
}

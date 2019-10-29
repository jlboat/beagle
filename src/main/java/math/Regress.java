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
package math;

import java.util.concurrent.atomic.DoubleAdder;
import java.util.concurrent.atomic.LongAdder;

/**
 * <p>Class {@code Regress} estimates a regression coefficient.
 * </p>
 * <p>Instances of class {@code Regress} are not thread-safe, but
 * concurrent updates are permitted, and invocation of the {@code beta()}
 * method in the absence of concurrent updates returns an accurate result.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Regress {

    private final LongAdder cnt;
    private final DoubleAdder sumX;
    private final DoubleAdder sumY;
    private final DoubleAdder sumXX;
    private final DoubleAdder sumXY;

    /**
     * Constructs a new {@code Regress} instance.
     */
    public Regress() {
        this.cnt = new LongAdder();
        this.sumX = new DoubleAdder();
        this.sumY = new DoubleAdder();
        this.sumXX = new DoubleAdder();
        this.sumXY = new DoubleAdder();
    }

    /**
     * Records the specified values of the dependent and independent variables.
     * @param x the value of the independent variable
     * @param y the value of the dependent variable
     */
    public void add(double x, double y) {
        cnt.increment();
        sumX.add(x);
        sumY.add(y);
        sumXX.add(x*x);
        sumXY.add(x*y);
    }

    /**
     * Records the specified values of the dependent and independent variables.
     * @param regress recorded values of the independent and dependent variables
     * @throws NullPointerException if {@code regress == null}t
     */
    public void add(Regress regress) {
        this.sumX.add(regress.sumX());
        this.sumY.add(regress.sumY());
        this.sumXX.add(regress.sumXX());
        this.sumXY.add(regress.sumXY());
        this.cnt.add(regress.cnt());
    }

    /**
     * Returns the number of recorded values of the independent variable.
     * The returned value is NOT an atomic snapshot. An accurate result is
     * guaranteed only if no concurrent updates occur during method invocation.
     *
     * @return the number of recorded values of the independent variable
     */
    public long cnt() {
        return cnt.sum();
    }

    /**
     * Returns the sum of recorded values of the independent variable.
     * The returned value is NOT an atomic snapshot. An accurate result is
     * guaranteed only if no concurrent updates occur during method invocation.
     * @return the sum of recorded values of the independent variable
     */
    public double sumX() {
        return sumX.sum();
    }

    /**
     * Returns the sum of recorded values of the dependent variable.
     * The returned value is NOT an atomic snapshot. An accurate result is
     * guaranteed only if no concurrent updates occur during method invocation.
     * @return the sum of recorded values of the dependent variable
     */
    public double sumY() {
        return sumY.sum();
    }

    /**
     * Returns the sum of the squared recorded values of the independent
     * variable.
     * The returned value is NOT an atomic snapshot. An accurate result is
     * guaranteed only if no concurrent updates occur during method invocation.
     * @return the sum of the squared recorded values of the independent variable
     */
    public double sumXX() {
        return sumXX.sum();
    }

    /**
     * Returns the sum of the products of the recorded values of the
     * independent and dependent variables.
     * The returned value is NOT an atomic snapshot. An accurate result is
     * guaranteed only if no concurrent updates occur during method invocation.
     * @return the sum of the products of the recorded values of the
     * independent and dependent variables
     */
    public double sumXY() {
        return sumXY.sum();
    }

    /**
     * Deletes all recorded values from {@code this}.
     */
    public void reset() {
        cnt.reset();
        sumX.reset();
        sumY.reset();
        sumXX.reset();
        sumXY.reset();
    }

    /**
     * Returns the regression coefficient for the recorded values of the
     * independent and dependent variables. The returned value is NOT an
     * atomic snapshot. An accurate result is guaranteed only if no concurrent
     * updates occur during method invocation.
     * @return the regression coefficient for the recorded values of the
     * independent and dependent variables
     */
    public double beta() {
        long n = cnt.sum();
        double sx = sumX.sum();
        double sy = sumY.sum();
        double sxx = sumXX.sum();
        double sxy = sumXY.sum();
        return (n*sxy - (sx*sy))/(n*sxx - (sx*sx));
    }

//    public static void main(String[] args) {
//        double[] x = new double[] {1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83};
//        double[] y = new double[] {52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46};
//
//        Regress regress = new Regress();
//        for (int j=0; j<x.length; ++j) {
//            regress.add(x[j], y[j]);
//        }
//        System.out.println("Beta from Wikipedia Simple Linear Regression article: 61.272");
//        System.out.println("Calculated Beta:                                      " + regress.beta());
//    }
}

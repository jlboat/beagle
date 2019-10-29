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

import beagleutil.IntInterval;
import blbutil.Const;
import java.util.Comparator;

/**
 * <p>Class {@code SampleSeg} represents a haplotype shared with a sample.
 * </p>
 * Instances of class {@code SampleSeg} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SampleSeg implements Comparable<SampleSeg>, IntInterval {

    private final int sample;
    private final int start;
    private final int end;

    /**
     * Constructs a new {@code SampleSeg} instance from the specified data.
     * @param sample the sample index
     * @param start the start of the segment (inclusive)
     * @param end the end of the segment (inclusive)
     * @throws IllegalArgumentException if {@code start > end}
     */
    public SampleSeg(int sample, int start, int end) {
        if (start > end) {
            throw new IllegalArgumentException(String.valueOf(start));
        }
        this.sample = sample;
        this.start = start;
        this.end = end;
    }

    /**
     * Returns the sample index.
     * @return the sample index
     */
    public int sample() {
        return sample;
    }

    /**
     * Returns the start of the segment (inclusive).
     * @return the start of the segment (inclusive)
     */
    @Override
    public int start() {
        return start;
    }

    /**
     * Returns the end of the segment (inclusive).
     * @return the end of the segment (inclusive)
     */
    @Override
    public int inclEnd() {
        return end;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(70);
        sb.append("[");
        sb.append(sample);
        sb.append(": ");
        sb.append(start);
        sb.append(Const.hyphen);
        sb.append(end);
        sb.append("]");
        return sb.toString();
    }

    /**
     * <p>Returns the hash code value for this object. The hash code is defined
     * by the following calculation:
     * </p>
     * <pre>
     *  int hash = 5;
     *  hash = 89 * hash + this.hap();
     *  hash = 89 * hash + this.start();
     *  hash = 89 * hash + this.end();
     </pre>
     * @return the hash code value for this object
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 89*hash + this.sample;
        hash = 89*hash + this.start;
        hash = 89*hash + this.end;
        return hash;
    }

    /**
     * Compares the specified object with this {@code SampleSeg} for
     * equality. Returns {@code true} if the specified object is a
     * {@code SampleSeg} instance and if this {@code SampleSeg} is
     * equal to the specified {@code SampleSeg}, and returns
     * {@code false}  otherwise.  Two {@code SampleSeg}  instances
     * are equal if they have equal haplotype indices,
     * equal starting marker indices, and equal ending marker indices.
     * @param o the reference object with which to compare.
     * @return {@code true} if the specified object is an
     * {@code SampleSeg} instance and if this {@code SampleSeg} is
     * equal to the specified {@code SampleSeg}
     */
    @Override
    public boolean equals(Object o) {
        if (o==null) {
            return false;
        }
        if (getClass()!=o.getClass()) {
            return false;
        }
        final SampleSeg other=(SampleSeg) o;
        return (this.sample==other.sample && this.start==other.start
                && this.end==other.end);
    }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     * {@code SampleSeg} instances are ordered first by
     * {@code this.start()}, then by {@code this.end()},
     * and finally by {@code this.sample()}.
     * @param hs the {@code SampleSeg} to be compared
     * @return a negative integer, zero, or a positive integer as this
     * {@code SampleSeg} is less than, equal to, or greater than the
     * specified {@code SampleSeg}
     * @throws NullPointerException if {@code o == null}
     */
    @Override
    public int compareTo(SampleSeg hs) {
        if (this.start != hs.start) {
            return (this.start < hs.start) ? -1 : 1;
        }
        else if (this.end != hs.end) {
            return (this.end < hs.end) ? -1 : 1;
        }
        if (this.sample != hs.sample) {
            return (this.sample < hs.sample) ? -1 : 1;
        }
        return 0;
    }

    /**
     * Returns a comparator that orders first by {@code this.sample()}, then
     * by {@code this.start()}, and finally by {@code this.end()}.
     * @return a comparator that orders first by {@code this.sample()}, then
     * by {@code this.start()}, and finally by {@code this.end()}
     */
    public static Comparator<SampleSeg> sampleComp() {
        return (SampleSeg t1, SampleSeg t2) -> {
            if (t1.sample != t2.sample) {
                return (t1.sample < t2.sample) ? -1 : 1;
            }
            if (t1.start != t2.start) {
                return (t1.start < t2.start) ? -1 : 1;
            }
            else if (t1.end != t2.end) {
                return (t1.end < t2.end) ? -1 : 1;
            }
            return 0;
        } ;
    }
}

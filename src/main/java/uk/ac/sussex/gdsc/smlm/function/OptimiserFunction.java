/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.function;

import java.util.Arrays;

import gnu.trove.list.array.TDoubleArrayList;

/**
 * Allow optimisation using Apache Commons Math 3 Optimiser
 */
public abstract class OptimiserFunction
{
    /** The x. */
    protected TDoubleArrayList x = null;
    /** The y. */
    protected TDoubleArrayList y = null;

    /**
     * Adds the point.
     *
     * @param x
     *            the x
     * @param y
     *            the y
     */
    public void addPoint(double x, double y)
    {
        if (this.x == null)
        {
            this.x = new TDoubleArrayList();
            this.y = new TDoubleArrayList();
        }
        this.x.add(x);
        this.y.add(y);
    }

    /**
     * Adds the data.
     *
     * @param x
     *            the x
     * @param y
     *            the y
     */
    public void addData(float[] x, float[] y)
    {
        this.x = new TDoubleArrayList();
        this.y = new TDoubleArrayList();
        for (int i = 0; i < x.length; i++)
        {
            this.x.add(x[i]);
            this.y.add(y[i]);
        }
    }

    /**
     * Adds the data.
     *
     * @param x
     *            the x
     * @param y
     *            the y
     */
    public void addData(double[] x, double[] y)
    {
        this.x = new TDoubleArrayList();
        this.y = new TDoubleArrayList();
        for (int i = 0; i < x.length; i++)
        {
            this.x.add(x[i]);
            this.y.add(y[i]);
        }
    }

    /**
     * Gets the x data
     *
     * @return the x
     */
    public double[] getX()
    {
        return x.toArray();
    }

    /**
     * Gets the y data
     *
     * @return the y
     */
    public double[] getY()
    {
        return y.toArray();
    }

    /**
     * Gets the weights array. This is an array filled with ones.
     *
     * @return the weights
     */
    public double[] getWeights()
    {
        final double[] w = new double[y.size()];
        Arrays.fill(w, 1);
        return w;
    }

    /**
     * Get the size.
     *
     * @return the size
     */
    public int size()
    {
        return x.size();
    }
}

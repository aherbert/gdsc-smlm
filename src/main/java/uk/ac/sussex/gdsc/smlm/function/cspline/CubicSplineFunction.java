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
package uk.ac.sussex.gdsc.smlm.function.cspline;

import uk.ac.sussex.gdsc.core.math.interpolation.CustomTricubicFunction;
import uk.ac.sussex.gdsc.smlm.function.Gradient2Function;

/**
 * Represent a cubic spline function. N splines are drawn into a target region.
 * <p>
 * The parameters are [Background + n *{Intensity, X, Y, Z}]. The spline can be scaled-down before sampling (i.e.
 * drawing on the target region). Only one sample is taken per index in the target region.
 */
public abstract class CubicSplineFunction implements Gradient2Function
{
    /** Index of the background in the parameters array. */
    public static final int BACKGROUND = 0;
    /** Index of the signal intensity in the parameters array. */
    public static final int SIGNAL = 1;
    /** Index of the x-position in the parameters array */
    public static final int X_POSITION = 2;
    /** Index of the y-position in the parameters array */
    public static final int Y_POSITION = 3;
    /** Index of the z-position in the parameters array */
    public static final int Z_POSITION = 4;
    /** Index of the x-standard deviation in the parameters array */

    /** The number of parameters per spline. */
    public static final int PARAMETERS_PER_PEAK = 4;

    /**
     * Gets the name of the parameter assuming a 2D Gaussian function.
     *
     * @param index
     *            the index (zero or above)
     * @return the name
     */
    public static String getName(int index)
    {
        final int i = 1 + (index - 1) % PARAMETERS_PER_PEAK;
        switch (i)
        {
            //@formatter:off
			case BACKGROUND: return "Background";
			case SIGNAL: return "Signal";
			case X_POSITION: return "X";
			case Y_POSITION: return "Y";
			case Z_POSITION: return "Z";
			default: return "Unknown: "+index;
			//@formatter:on
        }
    }

    /**
     * Gets the peak number (zero-based index) of the parameter assuming a cubic spline function.
     *
     * @param index
     *            the index (zero or above)
     * @return the peak number
     */
    public static int getPeak(int index)
    {
        if (index < 1)
            return 0;
        return (index - 1) / PARAMETERS_PER_PEAK;
    }

    /**
     * Gets the index of the parameter in a multi-peak parameter array assuming a cubic spline function.
     *
     * @param peak
     *            the peak number (zero-based index)
     * @param parameterIndex
     *            the parameter index for a single peak (this can use the class constants, e.g.
     *            {@link CubicSplineFunction#SIGNAL})
     * @return the index
     */
    public static int getIndex(int peak, int parameterIndex)
    {
        if (parameterIndex < 1)
            return 0;
        return peak * PARAMETERS_PER_PEAK + parameterIndex;
    }

    /**
     * Gets the name of the gradient parameter.
     *
     * @param index
     *            the index (must be within the array returned from {@link #gradientIndices()})
     * @return the name
     */
    public String getGradientParameterName(int index)
    {
        return getName(gradientIndices()[index]);
    }

    /**
     * Locate the index within the gradient indices for the specified parameter.
     *
     * @param parameterIndex
     *            the parameter index
     * @return the gradient index (or -1 if not present)
     */
    public int findGradientIndex(int parameterIndex)
    {
        final int[] gradientIndices = gradientIndices();
        for (int i = 0; i < gradientIndices.length; i++)
            if (gradientIndices[i] == parameterIndex)
                return i;
        return -1;
    }

    /**
     * Internal class to control visiting the correct cubic spline node for each [x][y] index in the
     * target range {@code [0 <= x < maxx]}, {@code [0 <= y < maxy]}.
     */
    protected abstract class TargetSpline
    {
        /** The id. */
        int id;

        /** The offset used for derivatives. */
        int offset;

        /** Working space for first order gradients. */
        double[] dfda = new double[3];

        /** Working space for second order gradients. */
        double[] d2fda2 = new double[3];

        /**
         * The x-index within the xy splines for x=0 in the target region.
         * It is offset by the scale for faster iteration with pre-increment loops.
         * This may be negative indicating that the spline does not overlap the target at x=0.
         */
        int ix0;
        /**
         * The y-index within the xy splines for y=0 in the target region.
         * It is offset by the scale for faster iteration with pre-increment loops.
         * This may be negative indicating that the spline does not overlap the target at y=0.
         */
        int iy0;

        /** The current y-index during iteration. */
        int yindex;

        /** The index of the current (X,Y) index within the xy splines. Used during forEach iteration. */
        int yxindex;

        /** The xy splines for the target z-position. */
        CustomTricubicFunction[] xySplines;

        /** Flag for each x-index to indicate if the spline overlaps the target region. */
        boolean[] activeX = new boolean[maxx];

        /** The target intensity multiplied by the scale^2 to normalise the integral. */
        double tI_by_s2;
        /** The target intensity multiplied by the scale^3 and negated. Used to scale the first order gradients. */
        double neg_tI_by_s3;
        /** The target intensity multiplied by the scale^4. Used to scale the second order gradients. */
        double tI_by_s4;

        /**
         * Initialise the target. This checks if the spline, shifted to centre at the given XYZ coordinates, will
         * overlap the target region. If true then initialisation is performed for function evaluation.
         *
         * @param id
         *            the id (used to write the correct derivatives)
         * @param tI
         *            the target Intensity
         * @param tX
         *            the target X
         * @param tY
         *            the target Y
         * @param tZ
         *            the target Z
         * @param order
         *            the derivative order
         * @return true, if the spline partially overlaps with the target region
         */
        public boolean initialise(int id, double tI, double tX, double tY, double tZ, int order)
        {
            // Map z to a position in the spline
            // We want 0 to be in the centre.
            // Note: Scale up the input parameter to the spline scale.
            final double z = cz + scale * tZ;

            if (z < 0 || z > maxSz)
                return false;

            // Shift the scaled XY spline bounds by the target centre
            final double x1 = lx + tX;
            final double x2 = ux + tX;
            final double y1 = ly + tY;
            final double y2 = uy + tY;

            // Check if it is within the region
            if (!(x2 > 0 && y2 > 0 && x1 < maxx && y1 < maxy))
                return false;

            // Convert the lower bounds to integer grid in the target region,
            // i.e. we sample the region at x=0,1,2,...
            // We want the first integer that the function overlaps,
            // i.e. can interpolate a value for so it must be above the lower bounds
            final int ix1 = (int) Math.ceil(x1);
            final int iy1 = (int) Math.ceil(y1);

            // How far into the unscaled function is the first point.
            // i.e. x=1 may be 0.6 above the scaled lower bound (0.4) but that would require
            // the first sample to be taken at spline[1] @ 0.2 if the scale is 2.
            final double x = scale * (ix1 - x1);
            final double y = scale * (iy1 - y1);

            // This is the first index for the spline sample
            final int ix = (int) x;
            final int iy = (int) y;
            int iz = (int) z;
            if (iz == maxSz)
                // Special edge case. Interpolation uses the node below with a (z-iz) value of 1
                iz--;

            // Get the spline index position for 0,0 offset by the scale (for pre-increment loops)
            ix0 = ix - scale * ix1 - scale;
            iy0 = iy - scale * iy1 - scale;
            // Store the xy splines for the z position
            xySplines = splines[iz];

            // Set the working flag for all x
            for (int i = 0, xindex = ix0; i < maxx; i++)
            {
                xindex += scale;
                // Note that in theory we could interpolate if xindex==maxSx
                // but this requires a new power table with (x-ix)==1 and previous xindex.
                // For speed this situation is ignored to avoid computing additional
                // power tables.
                activeX[i] = xindex >= 0 && xindex < maxSx;
            }

            computePowerTable(x - ix, y - iy, z - iz, order);

            // The scale is the increment we sample the PSF.
            // In order to have the same integral we adjust the intensity.
            this.tI_by_s2 = tI * scale2;
            this.id = id;
            if (order > 0)
            {
                this.offset = 1 + id * 4;
                this.neg_tI_by_s3 = -tI_by_s2 * scale;
                if (order == 2)
                    this.tI_by_s4 = tI_by_s2 * scale2;
            }

            return true;
        }

        /**
         * Reset for iteration through YX-order.
         */
        public void reset()
        {
            yindex = iy0;
        }

        /**
         * Checks if is the next Y-index is active. If true then it initialises the worker for iteration through x.
         *
         * @return true, if is next Y active
         */
        public boolean isNextYActive()
        {
            // pre-increment yindex
            yindex += scale;
            // Note that in theory we could interpolate if yindex==maxSy
            // but this requires a new power table with (y-iy)==1 and previous yindex.
            // For speed this situation is ignored to avoid computing additional
            // power tables.
            if (yindex >= 0 && yindex < maxSy)
            {
                // The y-index is inside the XY spline data
                // Reset the yx-index for iteration
                yxindex = yindex * maxSx + ix0;
                return true;
            }
            return false;
        }

        /**
         * Checks if is the next Y-index is active. If true then it initialises the worker for iteration through x.
         * Otherwise it resets the gradients.
         *
         * @param df_da
         *            the df da
         * @return true, if is next Y active
         */
        public boolean isNextYActive(double[] df_da)
        {
            // pre-increment yindex
            yindex += scale;
            if (yindex >= 0 && yindex < maxSy)
            {
                // The y-index is inside the XY spline data
                // Reset the yx-index for iteration
                yxindex = yindex * maxSx + ix0;
                return true;
            }
            // Zero gradients
            df_da[offset] = 0;
            df_da[offset + 1] = 0;
            df_da[offset + 2] = 0;
            df_da[offset + 3] = 0;
            return false;
        }

        /**
         * Checks if is the next Y-index is active. If true then it initialises the worker for iteration through x.
         * Otherwise it resets the gradients.
         *
         * @param df_da
         *            the df da
         * @param d2f_da2
         *            the d 2 f da 2
         * @return true, if is next Y active
         */
        public boolean isNextYActive(double[] df_da, double[] d2f_da2)
        {
            // pre-increment yindex
            yindex += scale;
            if (yindex >= 0 && yindex < maxSy)
            {
                // The y-index is inside the XY spline data
                // Reset the yx-index for iteration
                yxindex = yindex * maxSx + ix0;
                return true;
            }
            // Zero gradients
            df_da[offset] = 0;
            df_da[offset + 1] = 0;
            df_da[offset + 2] = 0;
            df_da[offset + 3] = 0;
            d2f_da2[offset + 1] = 0;
            d2f_da2[offset + 2] = 0;
            d2f_da2[offset + 3] = 0;
            return false;
        }

        /**
         * Compute the power tables for the given spline position and derivative order.
         *
         * @param x
         *            the x (range 0-1)
         * @param y
         *            the y (range 0-1)
         * @param z
         *            the z (range 0-1)
         * @param order
         *            the order
         */
        abstract public void computePowerTable(double x, double y, double z, int order);

        /**
         * Compute the value at the given x-index.
         * Assumes that the current y-index has been set with a call to #{@link TargetSpline#isNextYActive()}.
         *
         * @param x
         *            the x
         * @return the value
         */
        public double value(int x)
        {
            yxindex += scale; // pre-increment
            return (activeX[x]) ? tI_by_s2 * computeValue(xySplines[yxindex]) : 0;
        }

        /**
         * Compute the value.
         *
         * @param customTricubicFunction
         *            the custom tricubic function
         * @return the value
         */
        abstract public double computeValue(CustomTricubicFunction customTricubicFunction);

        /**
         * Compute the value and derivatives at the given x-index.
         * Assumes that the current y-index has been set with a call to {@link TargetSpline#isNextYActive(double[])}.
         *
         * @param x
         *            the x
         * @param df_da
         *            the df da
         * @return the value
         */
        public double value(int x, double[] df_da)
        {
            yxindex += scale; // pre-increment
            if (activeX[x])
            {
                final double v = computeValue1(xySplines[yxindex]);
                // Copy the gradients into the correct position and account for the intensity.
                // Negate the gradients as a shift of the position moves the spline the
                // other direction. Also scale the gradients appropriately.
                df_da[offset] = scale2 * v;
                df_da[offset + 1] = neg_tI_by_s3 * dfda[0];
                df_da[offset + 2] = neg_tI_by_s3 * dfda[1];
                df_da[offset + 3] = neg_tI_by_s3 * -dfda[2];
                return tI_by_s2 * v;
            }
            // Zero gradients
            df_da[offset] = 0;
            df_da[offset + 1] = 0;
            df_da[offset + 2] = 0;
            df_da[offset + 3] = 0;
            return 0;
        }

        /**
         * Compute the value and first-order derivatives. The derivatives are stored in {@link TargetSpline#dfda}.
         *
         * @param customTricubicFunction
         *            the custom tricubic function
         * @return the value
         */
        abstract public double computeValue1(CustomTricubicFunction customTricubicFunction);

        /**
         * Compute the value and derivatives at the given x-index.
         * Assumes that the current y-index has been set with a call to
         * #{@link TargetSpline#isNextYActive(double[],double[])}.
         *
         * @param x
         *            the x
         * @param df_da
         *            the df da
         * @param d2f_da2
         *            the d 2 f da 2
         * @return the value
         */
        public double value(int x, double[] df_da, double[] d2f_da2)
        {
            yxindex += scale; // pre-increment
            if (activeX[x])
            {
                final double v = computeValue2(xySplines[yxindex]);
                // Copy the gradients into the correct position and account for the intensity.
                // Negate the gradients as a shift of the position moves the spline the
                // other direction. Also scale the gradients appropriately.
                df_da[offset] = scale2 * v;
                df_da[offset + 1] = neg_tI_by_s3 * dfda[0];
                df_da[offset + 2] = neg_tI_by_s3 * dfda[1];
                df_da[offset + 3] = neg_tI_by_s3 * -dfda[2];
                d2f_da2[offset + 1] = tI_by_s4 * d2fda2[0];
                d2f_da2[offset + 2] = tI_by_s4 * d2fda2[1];
                d2f_da2[offset + 3] = tI_by_s4 * d2fda2[2];
                return tI_by_s2 * v;
            }
            // Zero gradients
            df_da[offset] = 0;
            df_da[offset + 1] = 0;
            df_da[offset + 2] = 0;
            df_da[offset + 3] = 0;
            d2f_da2[offset + 1] = 0;
            d2f_da2[offset + 2] = 0;
            d2f_da2[offset + 3] = 0;
            return 0;
        }

        /**
         * Compute the value, first- and second-order derivatives. The derivatives are stored in
         * {@link TargetSpline#dfda} and {@link TargetSpline#d2fda2}.
         *
         * @param customTricubicFunction
         *            the custom tricubic function
         * @return the value
         */
        abstract public double computeValue2(CustomTricubicFunction customTricubicFunction);

        /**
         * Checks if the power table is at the boundary of the cubic polynomial.
         *
         * @param dimension
         *            the dimension
         * @return true, if is node boundary
         */
        abstract public boolean isNodeBoundary(int dimension);
    }

    /**
     * Double precision computation of the target spline.
     */
    protected class DoubleTargetSpline extends TargetSpline
    {

        /** The table 1. */
        double[] table1 = new double[64];

        /** The table 2. */
        double[] table2;

        /** The table 3. */
        double[] table3;

        /** The table 6. */
        double[] table6;

        /** {@inheritDoc} */
        @Override
        public void computePowerTable(double x, double y, double z, int order)
        {
            table1 = CustomTricubicFunction.computePowerTable(x, y, z);

            // Create tables for derivatives
            if (order == 1)
            {
                table2 = multiply(2, table2);
                table3 = multiply(3, table3);
            }
            else if (order == 2)
            {
                table2 = multiply(2, table2);
                table3 = multiply(3, table3);
                table6 = multiply(6, table6);
            }
        }

        /**
         * Multiply.
         *
         * @param n
         *            the n
         * @param table
         *            the table
         * @return the double[]
         */
        private double[] multiply(int n, double[] table)
        {
            if (table == null)
                table = new double[64];
            for (int i = 0; i < 64; i++)
                table[i] = table1[i] * n;
            return table;
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1);
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue1(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1, table2, table3, dfda);
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue2(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1, table2, table3, table6, dfda, d2fda2);
        }

        /** {@inheritDoc} */
        @Override
        public boolean isNodeBoundary(int dimension)
        {
            return CustomTricubicFunction.isBoundary(dimension, table1);
        }
    }

    /**
     * Single precision computation of the target spline.
     */
    protected class FloatTargetSpline extends TargetSpline
    {

        /** The table 1. */
        float[] table1 = new float[64];

        /** The table 2. */
        float[] table2;

        /** The table 3. */
        float[] table3;

        /** The table 6. */
        float[] table6;

        /** {@inheritDoc} */
        @Override
        public void computePowerTable(double x, double y, double z, int order)
        {
            table1 = CustomTricubicFunction.computeFloatPowerTable(x, y, z);

            // Create tables for derivatives
            if (order == 1)
            {
                table2 = multiply(2, table2);
                table3 = multiply(3, table3);
            }
            else if (order == 2)
            {
                table2 = multiply(2, table2);
                table3 = multiply(3, table3);
                table6 = multiply(6, table6);
            }
        }

        /**
         * Multiply.
         *
         * @param n
         *            the n
         * @param table
         *            the table
         * @return the float[]
         */
        private float[] multiply(int n, float[] table)
        {
            if (table == null)
                table = new float[64];
            for (int i = 0; i < 64; i++)
                table[i] = table1[i] * n;
            return table;
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1);
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue1(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1, table2, table3, dfda);
        }

        /** {@inheritDoc} */
        @Override
        public double computeValue2(CustomTricubicFunction customTricubicFunction)
        {
            return customTricubicFunction.value(table1, table2, table3, table6, dfda, d2fda2);
        }

        /** {@inheritDoc} */
        @Override
        public boolean isNodeBoundary(int dimension)
        {
            return CustomTricubicFunction.isBoundary(dimension, table1);
        }
    }

    /** The scale to reduce the size of the spline before mapping to the target range (maxx * maxy). */
    protected int scale = 1;

    /** The scale squared (stored for convenience). */
    private int scale2;

    /** The x centre of the spline (unscaled). */
    protected double cx;

    /** The y centre of the spline (unscaled). */
    protected double cy;

    /** The z centre of the spline (unscaled). */
    protected double cz;

    /** The scaled lower x bound of the spline function with the centre at x=0,y=0 in the target region */
    protected double lx;

    /** The scaled lower y bound of the spline function with the centre at x=0,y=0 in the target region */
    protected double ly;

    /** The scaled upper x bound of the spline function with the centre at x=0,y=0 in the target region */
    protected double ux;

    /** The scaled upper y bound of the spline function with the centre at x=0,y=0 in the target region */
    protected double uy;

    /** Max size of spline data in the x-dimension. */
    protected final int maxSx;

    /** Max size of spline data in the y-dimension. */
    protected final int maxSy;

    /** Max size of spline data in the z-dimension. */
    protected final int maxSz;

    /** The tricubic spline packed as Z * YX arrays. */
    protected final CustomTricubicFunction[][] splines;

    /** The target range in the x-dimension. */
    protected final int maxx;

    /** The target range in the y-dimension. */
    protected final int maxy;

    /** The target background. */
    protected double tB;

    /**
     * Instantiates a new cubic spline function.
     *
     * @param splineData
     *            the spline data
     * @param maxx
     *            The maximum x value of the 2-dimensional data
     * @param maxy
     *            The maximum y value of the 2-dimensional data
     * @throws IllegalArgumentException
     *             If the function does not have an integer grid spacing from the origin
     */
    public CubicSplineFunction(CubicSplineData splineData, int maxx, int maxy) throws IllegalArgumentException
    {
        this.splines = splineData.splines;
        this.maxx = (maxx < 1) ? 1 : maxx;
        this.maxy = (maxy < 1) ? 1 : maxy;
        maxSx = splineData.maxx;
        maxSy = splineData.maxy;
        maxSz = splines.length;
        // Centre in the middle, assuming the min is zero
        cx = (maxSx / 2.0);
        cy = (maxSy / 2.0);
        cz = (maxSz / 2.0);
        setScale(1);
    }

    /**
     * Instantiates a new cubic spline function.
     *
     * @param splineData
     *            the spline data
     * @param maxx
     *            The maximum x value of the 2-dimensional data
     * @param maxy
     *            The maximum y value of the 2-dimensional data
     * @param cx
     *            the x centre of the spline data
     * @param cy
     *            the y centre of the spline data
     * @param cz
     *            the z centre of the spline data
     * @param scale
     *            the scale of the spline data
     * @throws IllegalArgumentException
     *             If the function does not have an integer grid spacing from the origin
     */
    public CubicSplineFunction(CubicSplineData splineData, int maxx, int maxy, double cx, double cy, double cz,
            int scale) throws IllegalArgumentException
    {
        this.splines = splineData.splines;
        this.maxx = (maxx < 1) ? 1 : maxx;
        this.maxy = (maxy < 1) ? 1 : maxy;
        maxSx = splineData.maxx;
        maxSy = splineData.maxy;
        maxSz = splines.length;
        this.cx = cx;
        this.cy = cy;
        this.cz = cz;
        setScale(scale);
    }

    /**
     * Update function bounds.
     */
    private void updateFunctionBounds()
    {
        // Store the bounds of the cubic spline if it were positioned at 0,0
        lx = -cx / scale;
        ux = (maxSx - cx) / scale;
        ly = -cy / scale;
        uy = (maxSy - cy) / scale;
        // Check if the centre was within the function
        if (lx > 0 || ly > 0 || ux < 0 || uy < 0 || cz < 0 || cz > maxSz)
            throw new IllegalArgumentException("Require the centre within the cubic spline");
    }

    /**
     * Gets the maximum x value of the 2-dimensional data.
     *
     * @return the maximum x value of the 2-dimensional data.
     */
    public int getMaxX()
    {
        return maxx;
    }

    /**
     * Gets the maximum y value of the 2-dimensional data.
     *
     * @return the maximum y value of the 2-dimensional data.
     */
    public int getMaxY()
    {
        return maxy;
    }

    /**
     * Gets the scale to map the cubic spline function to the integer grid. E.g. set a scale of 2 to render the spline
     * at half its size.
     *
     * @return the scale
     */
    public int getScale()
    {
        return scale;
    }

    /**
     * Sets the scale to map the cubic spline function to the integer grid. E.g. set a scale of 2 to render the spline
     * at half its size.
     *
     * @param scale
     *            the new scale
     * @throws IllegalArgumentException
     *             If the scale is not strictly positive
     */
    public void setScale(int scale) throws IllegalArgumentException
    {
        if (scale < 1)
            throw new IllegalArgumentException();
        this.scale = scale;
        updateFunctionBounds();
        scale2 = scale * scale;
    }

    /**
     * Gets the number of splines to draw.
     *
     * @return the number of splines to draw
     */
    abstract public int getN();

    /**
     * Gets the centre X.
     *
     * @return the centre X
     */
    public double getCentreX()
    {
        return cx;
    }

    /**
     * Sets the centre X.
     *
     * @param cx
     *            the new centre X
     */
    public void setCentreX(double cx)
    {
        this.cx = cx;
        updateFunctionBounds();
    }

    /**
     * Gets the centre Y.
     *
     * @return the centre Y
     */
    public double getCentreY()
    {
        return cy;
    }

    /**
     * Sets the centre Y.
     *
     * @param cy
     *            the new centre Y
     */
    public void setCentreY(double cy)
    {
        this.cy = cy;
        updateFunctionBounds();
    }

    /**
     * Gets the centre Z.
     *
     * @return the centre Z
     */
    public double getCentreZ()
    {
        return cz;
    }

    /**
     * Sets the centre Z.
     *
     * @param cz
     *            the new centre Z
     */
    public void setCentreZ(double cz)
    {
        this.cz = cz;
        updateFunctionBounds();
    }

    // The following properties may be overridden by optimised versions (e.g. no background computation)

    /**
     * @return True if the function can evaluate the background gradient.
     */
    public boolean evaluatesBackground()
    {
        return true;
    }

    /**
     * @return True if the function can evaluate the signal gradient.
     */
    public boolean evaluatesSignal()
    {
        return true;
    }

    /**
     * @return True if the function can evaluate the XY-position gradient.
     */
    public boolean evaluatesPosition()
    {
        return evaluatesX() && evaluatesY();
    }

    /**
     * @return True if the function can evaluate the X-position gradient.
     */
    public boolean evaluatesX()
    {
        return true;
    }

    /**
     * @return True if the function can evaluate the Y-position gradient.
     */
    public boolean evaluatesY()
    {
        return true;
    }

    /**
     * @return True if the function can evaluate the Z-position gradient.
     */
    public boolean evaluatesZ()
    {
        return true;
    }

    /** {@inheritDoc} */
    @Override
    public int size()
    {
        return maxx * maxy;
    }

    /** {@inheritDoc} */
    @Override
    public void initialise(double[] a)
    {
        initialise(a, 0);
    }

    /** {@inheritDoc} */
    @Override
    public void initialise0(double[] a)
    {
        initialise(a, 0);
    }

    /** {@inheritDoc} */
    @Override
    public void initialise1(double[] a)
    {
        initialise(a, 1);
    }

    /** {@inheritDoc} */
    @Override
    public void initialise2(double[] a)
    {
        initialise(a, 2);
    }

    /**
     * Initialise.
     *
     * @param a
     *            the a
     * @param order
     *            the order
     */
    protected abstract void initialise(double[] a, int order);

    /**
     * Checks if the gradient parameter is on a cubic spline node boundary. If true then the second order derivative
     * will not be smooth as they are not constant across spline points.
     * <p>
     * This only applies to XYZ gradients.
     *
     * @param gradientIndex
     *            the gradient index
     * @return true, if the peak is on a node boundary
     */
    abstract public boolean isNodeBoundary(int gradientIndex);
}

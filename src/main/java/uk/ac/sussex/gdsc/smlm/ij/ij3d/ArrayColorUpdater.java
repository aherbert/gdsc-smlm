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
package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;

/**
 * Class to create the colour array for a set of vertex coordinates with the same colour.
 *
 * @author Alex Herbert
 */
public abstract class ArrayColorUpdater
{
    /** The colour array. Used for GeometryArrays to colour coordinates with the same colour. */
    final float[] pointColor;

    /**
     * Instantiates a new array color updater.
     *
     * @param n
     *            the number of coordinates
     * @param withAlpha
     *            Set to true to use colours with alpha
     */
    ArrayColorUpdater(int n, boolean withAlpha)
    {
        pointColor = new float[((withAlpha) ? 4 : 3) * n];
    }

    /**
     * The size of the point colour array.
     *
     * @return the size
     */
    public int size()
    {
        return pointColor.length;
    }

    /**
     * Gets the number of coordinates.
     *
     * @return the n
     */
    public abstract int getN();

    /**
     * Set to true to use colours with alpha.
     *
     * @return true, if using colours with alpha
     */
    public abstract boolean withAlpha();

    /**
     * Gets the size of a single colour (either 3 or 4).
     *
     * @return the color size
     */
    public abstract int getColorSize();

    /**
     * Set the colour and alpha value and gets the current color array.
     *
     * @param color4f
     *            the color 4 f
     * @return the color
     */
    public abstract float[] getColors(Color4f color4f);

    /**
     * Set the colour and alpha value and gets the current color array.
     *
     * @param color
     *            the color
     * @param alpha
     *            the alpha
     * @return the color
     */
    public abstract float[] getColors(Color3f color, float alpha);

    /**
     * Set the colour value and gets the current color array.
     *
     * @param color
     *            the color
     * @return the color
     */
    public abstract float[] getColors(Color3f color);

    /**
     * Set the alpha value and gets the current color array.
     *
     * @param alpha
     *            the alpha
     * @return the color
     */
    public abstract float[] getColors(float alpha);

    /**
     * Class to update single colour coordinates with alpha.
     */
    private static class SingleArrayColorUpdater4 extends ArrayColorUpdater
    {
        SingleArrayColorUpdater4()
        {
            super(1, true);
        }

        @Override
        public int getN()
        {
            return 1;
        }

        @Override
        public boolean withAlpha()
        {
            return true;
        }

        @Override
        public int getColorSize()
        {
            return 4;
        }

        @Override
        public float[] getColors(Color4f color)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            pointColor[3] = color.w;
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color, float alpha)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            pointColor[3] = alpha;
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            return pointColor;
        }

        @Override
        public float[] getColors(float alpha)
        {
            pointColor[3] = alpha;
            return pointColor;
        }
    }

    /**
     * Class to update multiple colour coordinates with alpha.
     */
    private static class MultiArrayColorUpdater4 extends ArrayColorUpdater
    {
        MultiArrayColorUpdater4(int n)
        {
            super(n, true);
        }

        @Override
        public int getN()
        {
            return size() / 4;
        }

        @Override
        public boolean withAlpha()
        {
            return true;
        }

        @Override
        public int getColorSize()
        {
            return 4;
        }

        @Override
        public float[] getColors(Color4f color)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                pointColor[i++] = color.w;
            }
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color, float alpha)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                pointColor[i++] = alpha;
            }
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                i++;
            }
            return pointColor;
        }

        @Override
        public float[] getColors(float alpha)
        {
            for (int i = 3; i < pointColor.length; i += 4)
                pointColor[i] = alpha;
            return pointColor;
        }
    }

    /**
     * Class to update single colour coordinates with no alpha.
     */
    private static class SingleArrayColorUpdater3 extends ArrayColorUpdater
    {
        SingleArrayColorUpdater3()
        {
            super(1, false);
        }

        @Override
        public int getN()
        {
            return 1;
        }

        @Override
        public boolean withAlpha()
        {
            return false;
        }

        @Override
        public int getColorSize()
        {
            return 3;
        }

        @Override
        public float[] getColors(Color4f color)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            // Ignore alpha
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color, float alpha)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            // Ignore alpha
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color)
        {
            pointColor[0] = color.x;
            pointColor[1] = color.y;
            pointColor[2] = color.z;
            return pointColor;
        }

        @Override
        public float[] getColors(float alpha)
        {
            // Ignore alpha
            return pointColor;
        }
    }

    /**
     * Class to update multiple colour coordinates with no alpha.
     */
    private static class MultiArrayColorUpdater3 extends ArrayColorUpdater
    {
        MultiArrayColorUpdater3(int n)
        {
            super(n, false);
        }

        @Override
        public int getN()
        {
            return size() / 3;
        }

        @Override
        public boolean withAlpha()
        {
            return false;
        }

        @Override
        public int getColorSize()
        {
            return 3;
        }

        @Override
        public float[] getColors(Color4f color)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                // Ignore alpha
            }
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color, float alpha)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                // Ignore alpha
            }
            return pointColor;
        }

        @Override
        public float[] getColors(Color3f color)
        {
            for (int i = 0; i < pointColor.length;)
            {
                pointColor[i++] = color.x;
                pointColor[i++] = color.y;
                pointColor[i++] = color.z;
                // Ignore alpha
            }
            return pointColor;
        }

        @Override
        public float[] getColors(float alpha)
        {
            // Ignore alpha
            return pointColor;
        }
    }

    /**
     * Creates the array color updater.
     *
     * @param n
     *            the number of coordinates
     * @param withAlpha
     *            Set to true to use colours with alpha
     * @return the array color updater
     */
    public static ArrayColorUpdater create(int n, boolean withAlpha)
    {
        if (withAlpha)
            return (n == 1) ? new SingleArrayColorUpdater4() : new MultiArrayColorUpdater4(n);
        return (n == 1) ? new SingleArrayColorUpdater3() : new MultiArrayColorUpdater3(n);
    }
}

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
package uk.ac.sussex.gdsc.smlm.ij.utils;

import java.util.Arrays;

import org.jtransforms.fft.FloatFFT_3D;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import ij.ImageStack;
import ij.process.FloatProcessor;
import uk.ac.sussex.gdsc.core.utils.FloatEquality;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.smlm.function.StandardFloatValueProcedure;
import uk.ac.sussex.gdsc.smlm.function.gaussian.Gaussian2DFunction;
import uk.ac.sussex.gdsc.smlm.function.gaussian.GaussianFunctionFactory;
import uk.ac.sussex.gdsc.smlm.function.gaussian.QuadraticAstigmatismZModel;

@SuppressWarnings({ "javadoc" })
public class FloatDHT3DTest
{
    static int size = 16;
    static double centre = (size - 1) / 2.0;

    final static double gamma = 2;
    final static int zDepth = 5;
    static QuadraticAstigmatismZModel zModel = new QuadraticAstigmatismZModel(gamma, zDepth);

    private static FloatDHT3D createData(double cx, double cy, double cz)
    {
        final Gaussian2DFunction f = GaussianFunctionFactory.create2D(1, size, size,
                GaussianFunctionFactory.FIT_ASTIGMATISM, zModel);
        final int length = size * size;
        final float[] data = new float[size * length];
        final double[] a = new double[1 + Gaussian2DFunction.PARAMETERS_PER_PEAK];
        a[Gaussian2DFunction.SIGNAL] = 1;
        a[Gaussian2DFunction.X_POSITION] = cx;
        a[Gaussian2DFunction.Y_POSITION] = cy;
        a[Gaussian2DFunction.X_SD] = 1;
        a[Gaussian2DFunction.Y_SD] = 1;
        final StandardFloatValueProcedure p = new StandardFloatValueProcedure();
        for (int z = 0; z < size; z++)
        {
            a[Gaussian2DFunction.Z_POSITION] = z - cz;
            p.getValues(f, a, data, z * length);
        }
        return new FloatDHT3D(size, size, size, data, false);
    }

    private static FloatDHT3D createData()
    {
        return createData(centre, centre, centre);
    }

    private static FloatDHT3D createOctants(int w, int h, int d)
    {
        return new FloatDHT3D(createOctantsStack(w, h, d));
    }

    static ImageStack createOctantsStack(int w, int h, int d)
    {
        final int w_2 = w / 2;
        final int h_2 = h / 2;
        final int d_2 = d / 2;
        final ImageStack stack = new ImageStack(w, h, d);
        final FloatProcessor fp = new FloatProcessor(w, h);
        final float[] pixels = (float[]) fp.getPixels();
        fill(fp, w_2, 0, w_2, h_2, 1);
        fill(fp, 0, 0, w_2, h_2, 2);
        fill(fp, 0, h_2, w_2, h_2, 3);
        fill(fp, w_2, h_2, w_2, h_2, 4);
        for (int z = 0; z < d_2; z++)
            stack.setPixels(pixels.clone(), 1 + z);
        fill(fp, w_2, 0, w_2, h_2, 5);
        fill(fp, 0, 0, w_2, h_2, 6);
        fill(fp, 0, h_2, w_2, h_2, 7);
        fill(fp, w_2, h_2, w_2, h_2, 8);
        for (int z = d_2; z < d; z++)
            stack.setPixels(pixels.clone(), 1 + z);
        return stack;
    }

    static void fill(FloatProcessor fp, int x, int y, int w, int h, double value)
    {
        fp.setRoi(x, y, w, h);
        fp.setValue(value);
        fp.fill();
    }

    @Test
    public void canSwapOctants()
    {
        FloatDHT3D dht;

        // Simple test
        final float[] data = new float[] { 2, 1, 3, 4, 6, 5, 7, 8 };
        dht = new FloatDHT3D(2, 2, 2, data.clone(), false);
        dht.swapOctants();
        checkOctants(data, dht.getData());

        final int[] test = new int[] { 2, 4, 6 };
        for (final int w : test)
            for (final int h : test)
                for (final int d : test)
                {
                    dht = createOctants(w, h, d);

                    final float[] in = dht.getData().clone();

                    // This just tests that the swap of the DHT and the stack matches
                    final ImageStack stack = dht.getImageStack();
                    //uk.ac.sussex.gdsc.core.ij.Utils.display("Test", stack);
                    dht.swapOctants();
                    FloatDHT3D.swapOctants(stack);

                    final float[] e = new FloatDHT3D(stack).getData();
                    final float[] o = dht.getData();

                    checkOctants(in, o);

                    Assertions.assertArrayEquals(e, o);
                }
    }

    private static void checkOctants(float[] in, float[] out)
    {
        final int[] swap = new int[9];
        swap[1] = 7;
        swap[2] = 8;
        swap[3] = 5;
        swap[4] = 6;
        swap[5] = 3;
        swap[6] = 4;
        swap[7] = 1;
        swap[8] = 2;
        for (int i = 0; i < in.length; i++)
            Assertions.assertEquals(in[i], swap[(int) out[i]]);
    }

    @Test
    public void canConvolveAndDeconvolve()
    {
        final FloatDHT3D dht = createData();
        final float[] pixels = dht.getData().clone();
        dht.transform();

        final FloatDHT3D copy = dht.copy();
        copy.initialiseFastMultiply();

        final FloatDHT3D convolved = dht.multiply(dht);
        final FloatDHT3D deconvolved = convolved.divide(dht);

        final FloatDHT3D convolved2 = dht.multiply(copy);
        final FloatDHT3D deconvolved2 = convolved.divide(copy);

        Assertions.assertArrayEquals(convolved.getData(), convolved2.getData());
        Assertions.assertArrayEquals(deconvolved.getData(), deconvolved2.getData());

        float[] e = dht.getData();
        float[] o = deconvolved.getData();
        for (int i = 0; i < e.length; i++)
            Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-6f, 1e-6f));

        deconvolved.inverseTransform();

        // Test after reverse transform
        e = pixels;
        o = deconvolved.getData();

        for (int i = 0; i < e.length; i++)
            Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], 1e-7f, 1e-7f));
    }

    @Test
    public void canCorrelate()
    {
        final FloatDHT3D dht = createData();
        dht.transform();

        final FloatDHT3D copy = dht.copy();
        copy.initialiseFastMultiply();

        // Centre of power spectrum
        final int icentre = size / 2;

        for (int z = -1; z <= 1; z++)
            for (int y = -1; y <= 1; y++)
                for (int x = -1; x <= 1; x++)
                {
                    final FloatDHT3D dht2 = createData(centre + x, centre + y, centre + z);
                    dht2.transform();

                    final FloatDHT3D correlation = dht2.conjugateMultiply(dht);
                    final FloatDHT3D correlation2 = dht2.conjugateMultiply(copy);
                    Assertions.assertArrayEquals(correlation.getData(), correlation2.getData());

                    correlation.inverseTransform();
                    correlation.swapOctants();

                    final float[] pixels = correlation.getData();

                    final int i = SimpleArrayUtils.findMaxIndex(pixels);
                    final int[] xyz = correlation.getXYZ(i);

                    // This is how far dht has to move to align with dht2.
                    // To align dht2 with dht would be the opposite sign.
                    final int ox = xyz[0] - icentre;
                    final int oy = xyz[1] - icentre;
                    final int oz = xyz[2] - icentre;
                    //logger.fine(TestLog.getSupplier("Shift [%d,%d,%d], centre [%d,%d,%d]", x, y, z, xyz[0], xyz[1], xyz[2]);
                    Assertions.assertEquals(x, ox);
                    Assertions.assertEquals(y, oy);
                    Assertions.assertEquals(z, oz);
                }
    }

    @Test
    public void canConvertToDFT()
    {
        final FloatDHT3D dht = createData();
        final float[] input = dht.getData().clone();
        dht.transform();

        final FloatImage3D[] result = dht.toDFT(null, null);

        final float rel = 1e-6f;
        final float abs = 1e-6f;

        // Test reverse transform
        final FloatDHT3D dht2 = FloatDHT3D.fromDFT(result[0], result[1], null);

        final float[] e = dht.getData();
        final float[] o = dht2.getData();
        for (int i = 0; i < e.length; i++)
            Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(e[i], o[i], rel, abs));

        // Test verses full forward transform
        final FloatFFT_3D fft = new FloatFFT_3D(dht.ns, dht.nr, dht.nc);
        final float[] dft = Arrays.copyOf(input, 2 * e.length);
        fft.realForwardFull(dft);

        final float[] or = result[0].getData();
        final float[] oi = result[1].getData();
        for (int i = 0, j = 0; i < dft.length; i += 2, j++)
        {
            Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i], or[j], rel, abs));
            Assertions.assertTrue(FloatEquality.almostEqualRelativeOrAbsolute(dft[i + 1], oi[j], rel, abs));
        }
    }
}

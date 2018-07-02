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
package gdsc.smlm.ij.results;

import java.awt.Rectangle;
import java.util.Arrays;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.TurboList;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.CalibrationWriter;
import gdsc.smlm.data.config.PSFHelper;
import gdsc.smlm.data.config.PSFProtos.PSF;
import gdsc.smlm.data.config.PSFProtos.PSFType;
import gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import gdsc.smlm.results.Gaussian2DPeakResultHelper;
import gdsc.smlm.results.PeakResult;
import gdsc.test.TestSettings;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Test the IJImagePeakResults functionality.
 */
public class IJImagePeakResultsTest
{
	//private RandomGenerator rand = TestSettings.getRandomGenerator(System.currentTimeMillis() + System.identityHashCode(this));
	static Calibration calibration;
	static PSF psf;
	static
	{
		psf = PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D);
		CalibrationWriter cw = new CalibrationWriter();
		cw.setDistanceUnit(DistanceUnit.PIXEL);
		cw.setIntensityUnit(IntensityUnit.COUNT);
		calibration = cw.getCalibration();
	}

	private static final String title = "Test";
	Rectangle bounds = new Rectangle(0, 0, 3, 5);

	@Test
	public void canAddToSinglePixels()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		add(fp, r, 1, 1, 1);
		add(fp, r, 1, 2, 4);
		add(fp, r, 0, 1, 2);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canAddToSinglePixelsWithInvalidPositions()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		add(fp, r, 1, 1, 1);
		add(fp, r, 1, 2, 4);
		add(fp, r, 0, 1, 2);
		for (int x : new int[] { -1, 0, 1, bounds.width, bounds.width + 1 })
			for (int y : new int[] { -1, 0, 1, bounds.height, bounds.height + 1 })
				add(fp, r, x, y, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canAddToZeroSizeImage()
	{
		Rectangle zeroBounds = new Rectangle();
		canAddToZeroSizeImage(zeroBounds);
		zeroBounds.width = 1;
		canAddToZeroSizeImage(zeroBounds);
		zeroBounds.width = 0;
		zeroBounds.height = 1;
		canAddToZeroSizeImage(zeroBounds);
		zeroBounds.width = -1;
		zeroBounds.height = -1;
		canAddToZeroSizeImage(zeroBounds);
	}

	private static void canAddToZeroSizeImage(Rectangle bounds)
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		begin(r);
		for (int x : new int[] { -1, 0, 1, bounds.width, bounds.width + 1 })
			for (int y : new int[] { -1, 0, 1, bounds.height, bounds.height + 1 })
				r.add(x, y, 1);
		r.end();
		float[] expecteds = new float[1];
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	private void add(FloatProcessor fp, IJImagePeakResults r, int x, int y, float value)
	{
		addValue(r, x, y, value);
		fp.putPixelValue(x, y, fp.getPixelValue(x, y) + value);
	}

	@Test
	public void canInterpolateInMiddleOfPixel()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 1.5f, 1);
		fp.putPixelValue(1, 1, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInXAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1f, 1.5f, 2);
		fp.putPixelValue(0, 1, 1);
		fp.putPixelValue(1, 1, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInXAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 2f, 1.5f, 2);
		fp.putPixelValue(1, 1, 1);
		fp.putPixelValue(2, 1, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInYAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 1f, 2);
		fp.putPixelValue(1, 0, 1);
		fp.putPixelValue(1, 1, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInYAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 2f, 2);
		fp.putPixelValue(1, 1, 1);
		fp.putPixelValue(1, 2, 1);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInX()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.25f, 1.5f, 2);
		fp.putPixelValue(0, 1, 0.5f);
		fp.putPixelValue(1, 1, 1.5f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInX()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.75f, 1.5f, 2);
		fp.putPixelValue(1, 1, 1.5f);
		fp.putPixelValue(2, 1, 0.5f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInY()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 1.25f, 2);
		fp.putPixelValue(1, 0, 0.5f);
		fp.putPixelValue(1, 1, 1.5f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInY()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 1.75f, 2);
		fp.putPixelValue(1, 1, 1.5f);
		fp.putPixelValue(1, 2, 0.5f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInXYAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1f, 1f, 4);
		fp.putPixelValue(0, 0, 1f);
		fp.putPixelValue(0, 1, 1f);
		fp.putPixelValue(1, 0, 1f);
		fp.putPixelValue(1, 1, 1f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInXYAtPixelEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 2f, 2f, 4);
		fp.putPixelValue(1, 1, 1f);
		fp.putPixelValue(2, 1, 1f);
		fp.putPixelValue(1, 2, 1f);
		fp.putPixelValue(2, 2, 1f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateDownInXY()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.25f, 1.25f, 1);
		fp.putPixelValue(0, 0, 0.25f * 0.25f);
		fp.putPixelValue(0, 1, 0.75f * 0.25f);
		fp.putPixelValue(1, 0, 0.75f * 0.25f);
		fp.putPixelValue(1, 1, 0.75f * 0.75f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canInterpolateUpInXY()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.75f, 1.75f, 1);
		fp.putPixelValue(1, 1, 0.75f * 0.75f);
		fp.putPixelValue(2, 1, 0.75f * 0.25f);
		fp.putPixelValue(1, 2, 0.75f * 0.25f);
		fp.putPixelValue(2, 2, 0.25f * 0.25f);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void noInterpolateDownInXAtImageEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 0.5f, 1.5f, 2);
		fp.putPixelValue(0, 1, 2);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void noInterpolateUpInXAtImageEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 2.5f, 1.5f, 2);
		fp.putPixelValue(2, 1, 2);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void noInterpolateDownInYAtImageEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 0.5f, 2);
		fp.putPixelValue(1, 0, 2);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void noInterpolateUpInYAtImageEdge()
	{
		IJImagePeakResults r = new IJImagePeakResults(title, bounds, 1);
		r.setDisplayFlags(IJImagePeakResults.DISPLAY_WEIGHTED);
		FloatProcessor fp = new FloatProcessor(bounds.width, bounds.height);
		begin(r);
		addValue(r, 1.5f, 2.5f, 2);
		fp.putPixelValue(1, 2, 2);
		r.end();
		float[] expecteds = getImage(fp);
		float[] actuals = getImage(r);
		Assert.assertArrayEquals(expecteds, actuals, 0);
	}

	@Test
	public void canAddUsingDifferentMethods()
	{
		canAddUsingDifferentMethods(0);
	}

	@Test
	public void canAddUsingDifferentMethodsEqualized()
	{
		canAddUsingDifferentMethods(IJImagePeakResults.DISPLAY_EQUALIZED);
	}

	@Test
	public void canAddUsingDifferentMethodsEqualizedWeighted()
	{
		canAddUsingDifferentMethods(IJImagePeakResults.DISPLAY_EQUALIZED | IJImagePeakResults.DISPLAY_WEIGHTED);
	}

	@Test
	public void canAddUsingDifferentMethodsMax()
	{
		canAddUsingDifferentMethods(IJImagePeakResults.DISPLAY_MAX);
	}

	@Test
	public void canAddUsingDifferentMethodsReplace()
	{
		canAddUsingDifferentMethods(IJImagePeakResults.DISPLAY_REPLACE);
	}

	@Test
	public void canAddUsingDifferentMethodsWeighted()
	{
		canAddUsingDifferentMethods(IJImagePeakResults.DISPLAY_WEIGHTED);
	}

	private void canAddUsingDifferentMethods(int displayFlags)
	{
		RandomGenerator rand = TestSettings.getRandomGenerator();
		displayFlags |= IJImagePeakResults.DISPLAY_SIGNAL;

		IJImagePeakResults[] r = new IJImagePeakResults[8];
		PSF psf = PSFHelper.create(PSFType.ONE_AXIS_GAUSSIAN_2D);
		for (int i = 0; i < r.length; i++)
		{
			r[i] = new IJImagePeakResults(title + i, bounds, 1);
			r[i].setDisplayFlags(displayFlags);
			r[i].setPSF(psf);
			begin(r[i]);
		}

		int size = 20;
		int[] t = new int[size];
		float[] x = new float[size];
		float[] y = new float[size];
		float[] v = new float[size];
		for (int i = 0; i < size; i++)
		{
			t[i] = i;
			x[i] = (rand.nextFloat() * bounds.width);
			y[i] = (rand.nextFloat() * bounds.height);
			v[i] = (rand.nextFloat());

			addPeakResult(r[0], x[i], y[i], v[i]);
			addPeakResult(r[1], t[i], x[i], y[i], v[i]);
			addValue(r[2], x[i], y[i], v[i]);
			addValue(r[3], t[i], x[i], y[i], v[i]);
		}

		addPeakResults(r[4], x, y, v);
		addPeakResults(r[5], t, x, y, v);
		addValues(r[6], x, y, v);
		addValues(r[7], t, x, y, v);

		float[][] image = new float[r.length][];
		for (int i = 0; i < r.length; i++)
		{
			r[i].end();
			image[i] = getImage(r[i]);
			TestSettings.debug("[%d] = %s\n", i, Arrays.toString(image[i]));
		}

		// Test single value adds
		float[] expecteds = image[0];
		for (int i = 1; i < 4; i++)
		{
			float[] actuals = image[i];
			Assert.assertArrayEquals("Image" + i, expecteds, actuals, 0);
		}

		// Test multi value adds
		expecteds = image[4];
		for (int i = 5; i < image.length; i++)
		{
			float[] actuals = image[i];
			Assert.assertArrayEquals("Image" + i, expecteds, actuals, 0);
		}

		// Test they are roughly the same (differences occur due to floating point summation
		Assert.assertArrayEquals("Single != Multi", expecteds, image[0], 1e-5f);
	}

	private static void begin(IJImagePeakResults r)
	{
		r.setPSF(psf);
		r.setDisplayImage(false);
		r.setCalibration(calibration);
		r.begin();
	}

	private static void addPeakResult(IJImagePeakResults r, float x, float y, float v)
	{
		r.add(new PeakResult(x, y, v));
	}

	private static void addPeakResult(IJImagePeakResults r, int t, float x, float y, float v)
	{
		r.add(new PeakResult(t, 0, 0, 0, 0, 0, 0, createParams(x, y, v), null));
	}

	private static float[] createParams(float x, float y, float v)
	{
		return Gaussian2DPeakResultHelper.createOneAxisParams(0, v, x, y, 0, 1);
	}

	private static void addPeakResults(IJImagePeakResults r, float[] x, float[] y, float[] v)
	{
		TurboList<PeakResult> results = new TurboList<PeakResult>(x.length);
		for (int i = 0; i < x.length; i++)
			results.add(new PeakResult(x[i], y[i], v[i]));
		r.addAll(results);
	}

	private static void addPeakResults(IJImagePeakResults r, int[] t, float[] x, float[] y, float[] v)
	{
		TurboList<PeakResult> results = new TurboList<PeakResult>(x.length);
		for (int i = 0; i < x.length; i++)
			results.add(new PeakResult(t[i], 0, 0, 0, 0, 0, 0, createParams(x[i], y[i], v[i]), null));
		r.addAll(results);
	}

	private static void addValue(IJImagePeakResults r, float x, float y, float v)
	{
		r.add(x, y, v);
	}

	private static void addValue(IJImagePeakResults r, int t, float x, float y, float v)
	{
		r.add(t, x, y, v);
	}

	private static void addValues(IJImagePeakResults r, float[] x, float[] y, float[] v)
	{
		r.add(x, y, v);
	}

	private static void addValues(IJImagePeakResults r, int[] t, float[] x, float[] y, float[] v)
	{
		r.add(t, x, y, v);
	}

	private static float[] getImage(IJImagePeakResults r)
	{
		return getImage(r.getImagePlus().getProcessor());
	}

	private static float[] getImage(ImageProcessor ip)
	{
		return (float[]) ip.convertToFloat().getPixels();
	}
}

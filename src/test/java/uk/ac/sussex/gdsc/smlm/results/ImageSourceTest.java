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
package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

import uk.ac.sussex.gdsc.core.utils.Random;
import uk.ac.sussex.gdsc.test.TestSettings;

@SuppressWarnings({ "javadoc", "unused" })
public class ImageSourceTest
{
	@Test
	public void canConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final MemoryImageSource source = new MemoryImageSource(w, h, createData(w, h, n));
		Assert.assertEquals(w, source.getWidth());
		Assert.assertEquals(h, source.getHeight());
		Assert.assertEquals(n, source.getFrames());
	}

	@Test(expected = RuntimeException.class)
	public void nullDataThrowsConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = null;
		new MemoryImageSource(w, h, data);
	}

	@Test(expected = RuntimeException.class)
	public void nullArrayDataThrowsConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		data[2] = null;
		new MemoryImageSource(w, h, data);
	}

	@Test(expected = RuntimeException.class)
	public void invalidLengthArrayDataThrowsConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		data[2] = new float[w * h + 1];
		new MemoryImageSource(w, h, data);
	}

	@Test
	public void invalidWidthThrowsConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		int ok = 0;
		new MemoryImageSource(w, h, data);
		if (canConstruct(w + 1, h, data))
			ok++;
		if (canConstruct(w - 1, h, data))
			ok++;
		if (canConstruct(-1, h, data))
			ok++;
		if (canConstruct(0, h, data))
			ok++;
		Assert.assertEquals(0, ok);
	}

	@Test
	public void invalidHeightThrowsConstructMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		int ok = 0;
		new MemoryImageSource(w, h, data);
		if (canConstruct(w, h + 1, data))
			ok++;
		if (canConstruct(w, h - 1, data))
			ok++;
		if (canConstruct(w, -1, data))
			ok++;
		if (canConstruct(w, 0, data))
			ok++;
		Assert.assertEquals(0, ok);
	}

	private static boolean canConstruct(int w, int h, float[][] data)
	{
		try
		{
			new MemoryImageSource(w, h, data);
			return true;
		}
		catch (final RuntimeException e)
		{
			return false;
		}
	}

	@Test
	public void memoryImageSourceCanReturnDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final MemoryImageSource source = new MemoryImageSource(w, h, data);
		int i = 0;
		float[] d = null;
		Assert.assertTrue(source.open());
		while ((d = source.next()) != null)
			Assert.assertArrayEquals(data[i++], d, 0);
		Assert.assertEquals(i, data.length);
	}

	@Test
	public void memoryImageSourceCanReturnDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final MemoryImageSource source = new MemoryImageSource(w, h, data);

		final int[] frames = new int[data.length];
		for (int i = 0; i < data.length; i++)
			frames[i] = i + 1;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		Assert.assertTrue(source.open());
		for (int i = 0; i < data.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue(source.isValid(frame));
			Assert.assertArrayEquals(data[frame - 1], source.get(frame), 0);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void memoryImageSourceCanReturnCroppedDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final float[][] data = createData(w, h, 15);
		final MemoryImageSource source = new MemoryImageSource(w, h, data);
		int i = 0;
		float[] d = null;
		Assert.assertTrue(source.open());
		while ((d = source.next(bounds)) != null)
		{
			Assert.assertEquals(bounds.width * bounds.height, d.length);
			Assert.assertArrayEquals(crop(data[i], w, bounds), d, 0);
			i++;
		}
		Assert.assertEquals(i, data.length);
	}

	@Test
	public void memoryImageSourceCanReturnCroppedDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final float[][] data = createData(w, h, 15);
		final MemoryImageSource source = new MemoryImageSource(w, h, data);

		final int[] frames = new int[data.length];
		for (int i = 0; i < data.length; i++)
			frames[i] = i + 1;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		Assert.assertTrue(source.open());
		for (int i = 0; i < data.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue(source.isValid(frame));
			final float[] d = source.get(frame, bounds);
			Assert.assertEquals(bounds.width * bounds.height, d.length);
			Assert.assertArrayEquals(crop(data[frame - 1], w, bounds), d, 0);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test(expected = RuntimeException.class)
	public void memoryImageSourceThrowsWithInvalidBounds()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		data[2] = new float[w * h + 1];
		final MemoryImageSource source = new MemoryImageSource(w, h, data);
		source.next(new Rectangle(-1, 0, 1, 1));
	}

	@Test
	public void canConstructInterlacedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final MemoryImageSource source = new MemoryImageSource(w, h, createData(w, h, n));
		final InterlacedImageSource iSource = new InterlacedImageSource(source, start, size, skip);
		Assert.assertEquals(w, iSource.getWidth());
		Assert.assertEquals(h, iSource.getHeight());
		Assert.assertEquals(12 * 2 / 3, iSource.getFrames());
		Assert.assertEquals(start, iSource.getStart());
		Assert.assertEquals(size, iSource.getSize());
		Assert.assertEquals(skip, iSource.getSkip());
	}

	@Test(expected = RuntimeException.class)
	public void nullImageSourceThrowsConstructInterlacedImageSource()
	{
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		new InterlacedImageSource(null, start, size, skip);
	}

	@Test(expected = RuntimeException.class)
	public void invalidStartThrowsConstructInterlacedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int start = 0;
		final int size = 2;
		final int skip = 1;
		new InterlacedImageSource(new MemoryImageSource(w, h, createData(w, h, n)), start, size, skip);
	}

	@Test(expected = RuntimeException.class)
	public void invalidSizeThrowsConstructInterlacedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int start = 4;
		final int size = 0;
		final int skip = 1;
		new InterlacedImageSource(new MemoryImageSource(w, h, createData(w, h, n)), start, size, skip);
	}

	@Test(expected = RuntimeException.class)
	public void invalidSkipThrowsConstructInterlacedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int start = 4;
		final int size = 2;
		final int skip = -1;
		new InterlacedImageSource(new MemoryImageSource(w, h, createData(w, h, n)), start, size, skip);
	}

	@Test
	public void interlacedImageSourceCanReturnDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip);

		Assert.assertTrue(source.open());

		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };
		int i = 0;
		float[] d = null;
		while ((d = source.next()) != null)
		{
			final int startFrame = source.getStartFrameNumber();
			final int endFrame = source.getEndFrameNumber();
			//System.out.printf("Read %d - %d\n", startFrame, endFrame);
			Assert.assertEquals("Start and end frames do not match", startFrame, endFrame);
			Assert.assertEquals(expected[i], startFrame);
			Assert.assertArrayEquals(data[startFrame - 1], d, 0);
			i++;
		}
		Assert.assertEquals(i, source.getFrames());
	}

	@Test
	public void interlacedImageSourceCanReturnDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip);

		final int[] frames = new int[data.length];
		for (int i = 0; i < data.length; i++)
			frames[i] = i + 1;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };
		Assert.assertTrue(source.open());
		for (int i = 0; i < data.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue(source.isValid(frame));
			final float[] d = source.get(frame);
			if (isExpected(frame, expected))
				Assert.assertArrayEquals(data[frame - 1], d, 0);
			else
				Assert.assertNull(d);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void interlacedImageSourceCanReturnCroppedDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final ImageSource source = new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip);

		Assert.assertTrue(source.open());

		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };
		int i = 0;
		float[] d = null;
		while ((d = source.next(bounds)) != null)
		{
			final int startFrame = source.getStartFrameNumber();
			final int endFrame = source.getEndFrameNumber();
			Assert.assertEquals("Start and end frames do not match", startFrame, endFrame);
			Assert.assertEquals(expected[i], startFrame);
			Assert.assertEquals(bounds.width * bounds.height, d.length);
			Assert.assertArrayEquals(crop(data[startFrame - 1], w, bounds), d, 0);
			i++;
		}
		Assert.assertEquals(i, source.getFrames());
	}

	@Test
	public void interlacedImageSourceCanReturnCroppedDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final float[][] data = createData(w, h, 15);
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final ImageSource source = new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip);

		final int[] frames = new int[data.length];
		for (int i = 0; i < data.length; i++)
			frames[i] = i + 1;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };
		Assert.assertTrue(source.open());
		for (int i = 0; i < data.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue(source.isValid(frame));
			final float[] d = source.get(frame, bounds);
			if (isExpected(frame, expected))
			{
				Assert.assertEquals(bounds.width * bounds.height, d.length);
				Assert.assertArrayEquals(crop(data[frame - 1], w, bounds), d, 0);
			}
			else
				Assert.assertNull(d);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void canConstructAggregatedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int aggregate = 3;
		final MemoryImageSource source = new MemoryImageSource(w, h, createData(w, h, n));
		final AggregatedImageSource iSource = new AggregatedImageSource(source, aggregate);
		Assert.assertEquals(w, iSource.getWidth());
		Assert.assertEquals(h, iSource.getHeight());
		Assert.assertEquals(n / aggregate, iSource.getFrames());
		Assert.assertEquals(aggregate, iSource.getAggregate());
	}

	@Test(expected = RuntimeException.class)
	public void nullImageSourceThrowsConstructAggregatedImageSource()
	{
		final int aggregate = 3;
		new AggregatedImageSource(null, aggregate);
	}

	@Test(expected = RuntimeException.class)
	public void invalidAggregateThrowsConstructAggregatedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final int aggregate = 1;
		new AggregatedImageSource(new MemoryImageSource(w, h, createData(w, h, n)), aggregate);
	}

	@Test
	public void aggregatedImageSourceCanReturnDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final int aggregate = 3;
		final float[][] data = createData(w, h, 15);
		final ImageSource source = new AggregatedImageSource(new MemoryImageSource(w, h, data), aggregate);

		Assert.assertTrue(source.open());

		int i = 1;
		int ii = 0;
		float[] d = null;
		while ((d = source.next()) != null)
		{
			ii++;
			Assert.assertEquals(i, source.getStartFrameNumber());
			Assert.assertEquals(i + 2, source.getEndFrameNumber());
			final float[] all = combine(data[i - 1], data[i], data[i + 1]);
			Assert.assertArrayEquals(all, d, 0);
			i += 3;
		}
		Assert.assertEquals(ii, source.getFrames());
	}

	@Test
	public void aggregatedImageSourceCanReturnDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final int aggregate = 3;
		final float[][] data = createData(w, h, 15);
		final ImageSource source = new AggregatedImageSource(new MemoryImageSource(w, h, data), aggregate);

		final int[] frames = new int[data.length / 3];
		for (int i = 0, frame = 1; i < frames.length; i++, frame += 3)
			frames[i] = frame;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		Assert.assertTrue(source.open());
		for (int i = 0; i < frames.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue("Invalid frame " + frame, source.isValid(frame));
			final float[] d = source.get(frame);
			Assert.assertEquals(frame, source.getStartFrameNumber());
			Assert.assertEquals(frame + 2, source.getEndFrameNumber());
			final float[] all = combine(data[frame - 1], data[frame], data[frame + 1]);
			Assert.assertArrayEquals("Invalid frame data " + frame, all, d, 0);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void aggregatedImageSourceCanReturnCroppedDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final int aggregate = 3;
		final float[][] data = createData(w, h, 15);
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final ImageSource source = new AggregatedImageSource(new MemoryImageSource(w, h, data), aggregate);

		Assert.assertTrue(source.open());

		int i = 1;
		int ii = 0;
		float[] d = null;
		while ((d = source.next(bounds)) != null)
		{
			ii++;
			Assert.assertEquals(i, source.getStartFrameNumber());
			Assert.assertEquals(i + 2, source.getEndFrameNumber());
			final float[] all = combine(crop(data[i - 1], w, bounds), crop(data[i], w, bounds), crop(data[i + 1], w, bounds));
			Assert.assertArrayEquals(all, d, 0);
			i += 3;
		}
		Assert.assertEquals(ii, source.getFrames());
	}

	@Test
	public void aggregatedImageSourceCanReturnCroppedDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final int aggregate = 3;
		final float[][] data = createData(w, h, 15);
		final Rectangle bounds = new Rectangle(2, 1, 3, 1);
		final ImageSource source = new AggregatedImageSource(new MemoryImageSource(w, h, data), aggregate);

		final int[] frames = new int[data.length / 3];
		for (int i = 0, frame = 1; i < frames.length; i++, frame += 3)
			frames[i] = frame;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		Assert.assertTrue(source.open());
		for (int i = 0; i < frames.length; i++)
		{
			final int frame = frames[i];
			Assert.assertTrue(source.isValid(frame));
			final float[] d = source.get(frame, bounds);
			Assert.assertEquals(frame, source.getStartFrameNumber());
			Assert.assertEquals(frame + 2, source.getEndFrameNumber());
			final float[] all = combine(crop(data[frame - 1], w, bounds), crop(data[frame], w, bounds),
					crop(data[frame + 1], w, bounds));
			Assert.assertArrayEquals("Invalid frame data " + frame, all, d, 0);
		}
		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void canConstructAggregatedInterlacedImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final float[][] data = createData(w, h, n);
		final int aggregate = 3;
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip);
		final AggregatedImageSource iSource = new AggregatedImageSource(source, aggregate);
		Assert.assertEquals(w, iSource.getWidth());
		Assert.assertEquals(h, iSource.getHeight());
		Assert.assertEquals(12 * 2 / 3, source.getFrames());
		Assert.assertEquals((int) Math.ceil((double) source.getFrames() / aggregate), iSource.getFrames());
		Assert.assertEquals(aggregate, iSource.getAggregate());
	}

	@Test(expected = RuntimeException.class)
	public void constructInterlacedAggregatedImageSourceThrows()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final float[][] data = createData(w, h, n);
		final int aggregate = 3;
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new AggregatedImageSource(new MemoryImageSource(w, h, data), aggregate);
		new InterlacedImageSource(source, start, size, skip);
	}

	@Test
	public void aggregatedInterlacedImageSourceCanReturnDataWithNext()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final float[][] data = createData(w, h, n);
		final int aggregate = 3;
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new AggregatedImageSource(
				new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip), aggregate);

		Assert.assertTrue(source.open());

		// Set the expected frames returned by the interlacing
		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };

		int i = 0, ii = 0;
		float[] d = null;
		while ((d = source.next()) != null)
		{
			// Get the range for the data
			final int endE = FastMath.min(i + 2, expected.length - 1);
			final int startFrame = expected[i];
			final int endFrame = expected[endE];

			// Check the correct range is returned
			Assert.assertEquals(startFrame, source.getStartFrameNumber());
			Assert.assertEquals(endFrame, source.getEndFrameNumber());

			// Check the data is collated correctly
			final float[] all = new float[data[0].length];
			for (int e = i; e <= endE; e++)
			{
				final int frame = expected[e] - 1;
				for (int j = 0; j < all.length; j++)
					all[j] += data[frame][j];
			}
			Assert.assertArrayEquals(all, d, 0);
			i += 3;
			ii++;
		}
		Assert.assertEquals(ii, source.getFrames());
	}

	@Test
	public void aggregatedInterlacedImageSourceCanReturnDataWithGet()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final float[][] data = createData(w, h, n);
		final int aggregate = 3;
		final int start = 4;
		final int size = 2;
		final int skip = 1;
		final ImageSource source = new AggregatedImageSource(
				new InterlacedImageSource(new MemoryImageSource(w, h, data), start, size, skip), aggregate);

		// Set the expected frames returned by the interlacing
		final int[] expected = new int[] { 4, 5, 7, 8, 10, 11, 13, 14 };

		// Randomly pick points from the positions used by next()
		final int[] frames = new int[source.getFrames()];
		for (int i = 0, ii = 0; ii < expected.length; i++, ii += 3)
			frames[i] = ii;
		final RandomGenerator rg = TestSettings.getRandomGenerator();
		Random.shuffle(frames, rg);

		Assert.assertTrue(source.open());
		for (int i = 0; i < frames.length; i++)
		{
			// Get the range for the data
			int e = frames[i];
			final int endE = FastMath.min(e + 2, expected.length - 1);
			final int startFrame = expected[e];
			final int endFrame = expected[endE];

			// Get the data
			final float[] d = source.get(startFrame);

			// Check the correct range is returned
			Assert.assertEquals(startFrame, source.getStartFrameNumber());
			Assert.assertEquals(endFrame, source.getEndFrameNumber());

			// Check the data is collated correctly
			final float[] all = new float[data[0].length];
			for (; e <= endE; e++)
			{
				final int frame = expected[e] - 1;
				for (int j = 0; j < all.length; j++)
					all[j] += data[frame][j];
			}
			Assert.assertArrayEquals(all, d, 0);
		}

		// Check all the data is valid but skipped interlaced points return null
		for (int i = 0; i < data.length; i++)
		{
			final int frame = i + 1;
			Assert.assertTrue(source.isValid(frame));
			if (!isExpected(frame, expected))
				Assert.assertNull(source.get(frame));
		}

		Assert.assertFalse(source.isValid(0));
		Assert.assertFalse(source.isValid(data.length + 1));
	}

	@Test
	public void canSerialiseMemoryImageSource()
	{
		final int w = 5;
		final int h = 3;
		final int n = 15;
		final String name = "canSerialiseMemoryImageSource";

		final MemoryImageSource source = new MemoryImageSource(w, h, createData(w, h, n));
		source.setName(name);
		source.setFreeMemoryOnClose(true);

		final String xml = source.toXML();
		//System.out.println(xml);

		final MemoryImageSource source2 = (MemoryImageSource) ImageSource.fromXML(xml);

		Assert.assertEquals(w, source2.getWidth());
		Assert.assertEquals(h, source2.getHeight());
		Assert.assertEquals(n, source2.getFrames());
		Assert.assertEquals(name, source2.getName());
		Assert.assertEquals(true, source2.isFreeMemoryOnClose());

		float[] data;
		while ((data = source.next()) != null)
		{
			final float[] data2 = source2.next();
			Assert.assertArrayEquals(data, data2, 1e-6f);
		}
	}

	/**
	 * Create data using the specified dimensions and the number of frames. Each frame will have a different base
	 * number and each index will be unique in the frame.
	 *
	 * @param w
	 *            width
	 * @param h
	 *            height
	 * @param n
	 *            The number of frames
	 * @return The data
	 */
	private static float[][] createData(int w, int h, int n)
	{
		final List<float[]> data = new ArrayList<>(n);
		float f = 1;
		for (int i = 0; i < n; i++)
		{
			data.add(createData(w, h, f));
			f *= 2;
		}
		return data.toArray(new float[0][0]);
	}

	/**
	 * Create data using the specified dimensions and the base number. Each index will be unique in the frame.
	 *
	 * @param w
	 *            width
	 * @param h
	 *            height
	 * @param f
	 *            the base number
	 * @return The data
	 */
	private static float[] createData(int w, int h, float f)
	{
		final float[] data = new float[w * h];
		Arrays.fill(data, f);
		for (int i = 0; i < data.length; i++)
			data[i] += (double) i / data.length;
		return data;
	}

	/**
	 * Check if the frame is contained in the expected frames array.
	 *
	 * @param frame
	 *            the frame
	 * @param expected
	 *            the expected frames
	 * @return true if within the array
	 */
	private static boolean isExpected(int frame, int[] expected)
	{
		for (final int e : expected)
			if (e == frame)
				return true;
		return false;
	}

	/**
	 * Crop the data to the specified bounds.
	 *
	 * @param data
	 *            the data
	 * @param width
	 *            the width
	 * @param bounds
	 *            the bounds
	 * @return the cropped data
	 */
	float[] crop(float[] data, int width, Rectangle bounds)
	{
		final float[] newData = new float[bounds.width * bounds.height];
		for (int y = bounds.y, ii = 0; y < bounds.y + bounds.height; y++)
			for (int x = bounds.x; x < bounds.x + bounds.width; x++, ii++)
			{
				final int index = y * width + x;
				newData[ii] = data[index];
			}
		return newData;
	}

	/**
	 * Sum all the input arrays.
	 *
	 * @param f
	 *            the arrays
	 * @return The summed array
	 */
	private static float[] combine(float[]... f)
	{
		final float[] all = Arrays.copyOf(f[0], f[0].length);
		for (int i = 1; i < f.length; i++)
			for (int j = 0; j < all.length; j++)
				all[j] += f[i][j];
		return all;
	}
}
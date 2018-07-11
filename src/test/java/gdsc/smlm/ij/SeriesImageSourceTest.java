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
package gdsc.smlm.ij;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Random;
import gdsc.core.utils.SimpleArrayUtils;
import gdsc.smlm.results.ImageSource.ReadHint;
import gdsc.test.TestSettings;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.io.TiffEncoder;
import ij.measure.Calibration;

@SuppressWarnings({ "javadoc" })
public class SeriesImageSourceTest
{
	int w = 10, h = 5, d = 7;

	@Test
	public void canReadBigTIFFSequentially() throws IOException
	{
		canReadBigTIFFSequentially(false, true);
	}

	@Test
	public void canReadBigTIFFSequentiallyInMemory() throws IOException
	{
		canReadBigTIFFSequentially(true, true);
	}

	@Test
	public void canReadBigTIFFSequentiallyBE() throws IOException
	{
		canReadBigTIFFSequentially(false, false);
	}

	@Test
	public void canReadBigTIFFSequentiallyInMemoryBE() throws IOException
	{
		canReadBigTIFFSequentially(true, false);
	}

	private void canReadBigTIFFSequentially(boolean inMemory, boolean intelByteOrder) throws IOException
	{
		int n = 2;
		String[] filenames = createFilenames(n);
		ImageStack[] stacks = createSeries(filenames, intelByteOrder);
		SeriesImageSource source = new SeriesImageSource("Test", filenames);
		if (!inMemory)
			source.setBufferLimit(0); // To force standard reading functionality
		source.setReadHint(ReadHint.SEQUENTIAL);
		source.open();
		Assert.assertEquals(w, source.getWidth());
		Assert.assertEquals(h, source.getHeight());
		Assert.assertEquals(d * n, source.getFrames());
		for (int i = 0; i < stacks.length; i++)
		{
			for (int j = 0; j < d; j++)
			{
				float[] e = (float[]) stacks[i].getPixels(j + 1);
				float[] o = source.next();
				Assert.assertArrayEquals(e, o, 0);
			}
		}
		Assert.assertNull(source.next());
		source.close();
	}

	@Test
	public void canReadBigTIFFNonSequentially() throws IOException
	{
		canReadBigTIFFNonSequentially(false, true);
	}

	@Test
	public void canReadBigTIFFNonSequentiallyInMemory() throws IOException
	{
		canReadBigTIFFNonSequentially(true, true);
	}

	@Test
	public void canReadBigTIFFNonSequentiallyBE() throws IOException
	{
		canReadBigTIFFNonSequentially(false, false);
	}

	@Test
	public void canReadBigTIFFNonSequentiallyInMemoryBE() throws IOException
	{
		canReadBigTIFFNonSequentially(true, false);
	}

	private void canReadBigTIFFNonSequentially(boolean inMemory, boolean intelByteOrder) throws IOException
	{
		int n = 2;
		String[] filenames = createFilenames(n);
		ImageStack[] stacks = createSeries(filenames, intelByteOrder);
		SeriesImageSource source = new SeriesImageSource("Test", filenames);
		if (!inMemory)
			source.setBufferLimit(0); // To force standard reading functionality
		source.setReadHint(ReadHint.NONSEQUENTIAL);
		source.open();
		Assert.assertEquals(w, source.getWidth());
		Assert.assertEquals(h, source.getHeight());
		Assert.assertEquals(d * n, source.getFrames());
		float[][] pixels = new float[n * d][];
		for (int i = 0, k = 0; i < stacks.length; i++)
		{
			for (int j = 0; j < d; j++)
			{
				pixels[k++] = (float[]) stacks[i].getPixels(j + 1);
			}
		}

		RandomGenerator r = TestSettings.getRandomGenerator();
		for (int i = 0; i < 3; i++)
		{
			int[] random = Random.sample(pixels.length / 2, pixels.length, r);
			for (int frame : random)
			{
				//System.out.printf("[%d] frame = %d\n", i, frame);
				float[] e = pixels[frame];
				float[] o = source.get(frame + 1); // 1-base index on the frame
				Assert.assertArrayEquals(e, o, 0);
			}
		}
	}

	private String[] createFilenames(int n) throws IOException
	{
		String[] filenames = new String[n];
		for (int i = 0; i < n; i++)
		{
			File path = File.createTempFile(this.getClass().getSimpleName(), ".tif");
			path.deleteOnExit();
			filenames[i] = path.getCanonicalPath();
		}
		return filenames;
	}

	private ImageStack[] createSeries(String[] filenames, boolean intelByteOrder) throws IOException
	{
		int n = filenames.length;
		ImageStack[] stacks = new ImageStack[n];
		int index = 0;
		int length = w * h;
		for (int i = 0; i < n; i++)
		{
			ImageStack stack = new ImageStack(w, h);
			for (int j = 0; j < d; j++)
			{
				stack.addSlice(null, SimpleArrayUtils.newArray(length, index, 1f));
				index += length;
			}
			ImagePlus imp = new ImagePlus(null, stack);
			// Add a calibration with origin
			Calibration c = imp.getCalibration();
			c.xOrigin = 4;
			c.yOrigin = 5;
			saveAsTiff(imp, filenames[i], intelByteOrder);
			stacks[i] = stack;
		}
		return stacks;
	}

	private void saveAsTiff(ImagePlus imp, String path, boolean intelByteOrder) throws IOException
	{
		// IJ.saveAsTiff(imp, path);

		FileInfo fi = imp.getFileInfo();
		fi.nImages = imp.getStackSize();
		ij.Prefs.intelByteOrder = intelByteOrder;
		DataOutputStream out = null;
		try
		{
			TiffEncoder file = new TiffEncoder(fi);
			out = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(path)));
			file.write(out);
			out.close();
		}
		finally
		{
			if (out != null)
				try
				{
					out.close();
				}
				catch (IOException e)
				{
				}
		}
	}
}

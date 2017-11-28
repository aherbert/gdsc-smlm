package gdsc.smlm.ij;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.Random;
import gdsc.core.utils.SimpleArrayUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;

public class SeriesImageSourceTest
{
	int w = 10, h = 5, d = 7;

	@Test
	public void canReadBigTIFFSequentially() throws IOException
	{
		canReadBigTIFFSequentially(false);
	}

	@Test
	public void canReadBigTIFFSequentiallyInMemory() throws IOException
	{
		canReadBigTIFFSequentially(true);
	}

	private void canReadBigTIFFSequentially(boolean inMemory) throws IOException
	{
		int n = 2;
		String[] filenames = createFilenames(n);
		ImageStack[] stacks = createSeries(filenames);
		SeriesImageSource source = new SeriesImageSource("Test", filenames);
		if (!inMemory)
			source.setBufferLimit(0); // To force standard reading functionality
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
	}

	@Test
	public void canReadBigTIFFNonSequentially() throws IOException
	{
		canReadBigTIFFNonSequentially(false);
	}

	@Test
	public void canReadBigTIFFNonSequentiallyInMemory() throws IOException
	{
		canReadBigTIFFNonSequentially(true);
	}
	
	private void canReadBigTIFFNonSequentially(boolean inMemory) throws IOException
	{
		int n = 2;
		String[] filenames = createFilenames(n);
		ImageStack[] stacks = createSeries(filenames);
		SeriesImageSource source = new SeriesImageSource("Test", filenames);
		if (!inMemory)
			source.setBufferLimit(0); // To force standard reading functionality
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

		RandomGenerator r = new Well19937c(30051977);
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

	private ImageStack[] createSeries(String[] filenames)
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
			IJ.saveAsTiff(new ImagePlus(null, stack), filenames[i]);
			stacks[i] = stack;
		}
		return stacks;
	}
}

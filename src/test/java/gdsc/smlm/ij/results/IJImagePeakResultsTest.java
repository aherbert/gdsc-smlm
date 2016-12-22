package gdsc.smlm.ij.results;

import java.awt.Rectangle;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.TurboList;
import gdsc.smlm.function.gaussian.Gaussian2DFunction;
import gdsc.smlm.results.PeakResult;

/**
 * Test the IJImagePeakResults functionality.
 */
public class IJImagePeakResultsTest
{
	private RandomGenerator rand = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

	String title = "Test";
	Rectangle bounds = new Rectangle(0, 0, 3, 3);

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
		IJImagePeakResults[] r = new IJImagePeakResults[8];
		for (int i = 0; i < r.length; i++)
		{
			r[i] = new IJImagePeakResults(title, bounds, 1);
			r[i].setDisplayFlags(displayFlags);
			r[i].begin();
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
			image[i] = (float[]) r[i].getImagePlus().getProcessor().convertToFloat().getPixels();
		}

		float[] expecteds = image[0];
		for (int i = 1; i < r.length; i++)
		{
			float[] actuals = image[i];
			Assert.assertArrayEquals(expecteds, actuals, 0);
		}
	}

	private void addPeakResult(IJImagePeakResults r, float x, float y, float v)
	{
		r.add(new PeakResult(x, y, 1, v));
	}

	private void addPeakResult(IJImagePeakResults r, int t, float x, float y, float v)
	{
		r.add(new PeakResult(t, 0, 0, 0, 0, 0, createParams(x, y, v), null));
	}

	private float[] createParams(float x, float y, float v)
	{
		float[] params = new float[7];
		params[Gaussian2DFunction.X_POSITION] = x;
		params[Gaussian2DFunction.Y_POSITION] = y;
		params[Gaussian2DFunction.SIGNAL] = v;
		return params;
	}

	private void addPeakResults(IJImagePeakResults r, float[] x, float[] y, float[] v)
	{
		TurboList<PeakResult> results = new TurboList<PeakResult>(x.length);
		for (int i = 0; i < x.length; i++)
			results.add(new PeakResult(x[i], y[i], 1, v[i]));
		r.addAll(results);
	}

	private void addPeakResults(IJImagePeakResults r, int[] t, float[] x, float[] y, float[] v)
	{
		TurboList<PeakResult> results = new TurboList<PeakResult>(x.length);
		for (int i = 0; i < x.length; i++)
			results.add(new PeakResult(t[i], 0, 0, 0, 0, 0, createParams(x[i], y[i], v[i]), null));
		r.addAll(results);
	}

	private void addValue(IJImagePeakResults r, float x, float y, float v)
	{
		r.add(x, y, v);
	}

	private void addValue(IJImagePeakResults r, int t, float x, float y, float v)
	{
		r.add(t, x, y, v);
	}

	private void addValues(IJImagePeakResults r, float[] x, float[] y, float[] v)
	{
		r.add(x, y, v);
	}

	private void addValues(IJImagePeakResults r, int[] t, float[] x, float[] y, float[] v)
	{
		r.add(t, x, y, v);
	}
}

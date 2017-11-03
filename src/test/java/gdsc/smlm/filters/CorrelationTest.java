package gdsc.smlm.filters;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.jtransforms.fft.FloatFFT_2D;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.utils.SimpleArrayUtils;
import ij.plugin.filter.EDM;
import ij.process.ByteProcessor;
import ij.process.FHT2;
import ij.process.FloatProcessor;

public class CorrelationTest
{
	@Test
	public void canCrossCorrelate()
	{
		int size = 16;
		int ex = 5, ey = 7;
		int ox = 1, oy = 2;
		RandomGenerator r = new Well19937c(30051977);
		FloatProcessor fp1 = createProcessor(size, ex, ey, 4, 4, r);
		// This must be offset from the centre
		FloatProcessor fp2 = createProcessor(size, size / 2 + ox, size / 2 + oy, 4, 4, r);

		float[] i1 = ((float[]) fp1.getPixels()).clone();
		float[] i2 = ((float[]) fp2.getPixels()).clone();

		FHT2 fht1 = new FHT2(fp1);
		fht1.transform();
		FHT2 fht2 = new FHT2(fp2);
		fht2.transform();

		FHT2 fhtE = fht1.conjugateMultiply(fht2);
		fhtE.inverseTransform();
		fhtE.swapQuadrants();

		float[] e = (float[]) fhtE.getPixels();
		int max = SimpleArrayUtils.findMaxIndex(e);
		int x = max % 16;
		int y = max / 16;

		Assert.assertEquals(ex, x + ox);
		Assert.assertEquals(ey, y + oy);

		// Do this with JTransform
		// TODO - Does this need a full transform realForwardFull() ?
		
		FloatFFT_2D fft = new FloatFFT_2D(size, size);
		fft.realForward(i1);
		fft.realForward(i2);

		// Conjugate multiply
		float[] o = new float[i1.length];

		for (int i = 0; i < o.length; i += 2)
		{
			float a = i1[i];
			float b = i2[i];
			float c = i1[i + 1];
			float d = -i2[i + 1]; // Complex conjugate has opposite imaginary part
			// Complex multiplication:
			// https://en.wikipedia.org/wiki/Complex_number#Multiplication_and_division
			o[i] = (a * c - b * d);
			o[i + 1] = (b * c + a * d);
		}

		// Transform back
		fft.realInverse(o, true);
		FloatProcessor op = new FloatProcessor(size, size, o);
		fht1.swapQuadrants(op);

		int max2 = SimpleArrayUtils.findMaxIndex(o);

		Assert.assertEquals(max, max2);
		for (int i = 0; i < o.length; i++)
			System.out.printf("[%d] %f vs %f = %f\n", i, e[i], o[i], o[i] / e[i]);
		//Assert.assertArrayEquals(e, o, 0);
	}

	private FloatProcessor createProcessor(int size, int x, int y, int w, int h, RandomGenerator r)
	{
		ByteProcessor bp = new ByteProcessor(size, size);
		bp.setColor(255);
		bp.fillOval(x, y, w, h);
		EDM e = new EDM();
		FloatProcessor fp = e.makeFloatEDM(bp, 0, true);
		if (r != null)
		{
			float[] d = (float[]) fp.getPixels();
			for (int i = 0; i < d.length; i++)
				d[i] += r.nextFloat() * 0.01;
		}
		return fp;
	}
}

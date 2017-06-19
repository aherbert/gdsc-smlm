package gdsc.smlm.results;

import org.junit.Assert;
import org.junit.Test;

public class PeakResultHelperTest
{
	@Test
	public void canConvertLocalBackgroundToNoise()
	{
		double gain = 6;

		double[] photons = { 0, 1, 2, 4, 10, 50, 100 };

		// CCD
		for (double p : photons)
		{
			// Assuming a Poisson distribution N photons should have a noise of sqrt(N).
			// However the input and output are in ADU counts so we apply the gain.
			double n = PeakResultHelper.localBackgroundToNoise(p * gain, gain, false);
			Assert.assertEquals("CCD " + p, Math.sqrt(p) * gain, n, 0);
		}

		// EM-CCD
		for (double p : photons)
		{
			// Assuming a Poisson distribution N photons should have a noise of sqrt(N * 2)
			// (due to the EM-CCD noise factor of 2).
			// However the input and output are in ADU counts so we apply the gain.
			double n = PeakResultHelper.localBackgroundToNoise(p * gain, gain, true);
			Assert.assertEquals("EM-CCD " + p, Math.sqrt(2 * p) * gain, n, 0);
		}
	}

	@Test
	public void canConvertLocalBackgroundToNoiseAndBack()
	{
		double gain = 6;

		double[] photons = { 0, 1, 2, 4, 10, 50, 100 };

		for (boolean emCCD : new boolean[] { false, true })
		{
			for (double p : photons)
			{
				double b = p * gain;
				double n = PeakResultHelper.localBackgroundToNoise(b, gain, emCCD);
				double b2 = PeakResultHelper.noiseToLocalBackground(n, gain, emCCD);
				Assert.assertEquals(emCCD + " " + p, b, b2, 1e-6);
			}
		}
	}

}

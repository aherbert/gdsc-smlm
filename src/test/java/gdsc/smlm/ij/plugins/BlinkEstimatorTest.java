package gdsc.smlm.ij.plugins;

import gdsc.smlm.fitting.function.Gaussian2DFunction;
import gdsc.smlm.fitting.utils.DoubleEquality;
import gdsc.smlm.model.ActivationEnergyImageModel;
import gdsc.smlm.model.CompoundMoleculeModel;
import gdsc.smlm.model.FluorophoreSequenceModel;
import gdsc.smlm.model.ImageModel;
import gdsc.smlm.model.LocalisationModel;
import gdsc.smlm.model.MoleculeModel;
import gdsc.smlm.model.SpatialDistribution;
import gdsc.smlm.model.SpatialIllumination;
import gdsc.smlm.model.UniformDistribution;
import gdsc.smlm.model.UniformIllumination;
import gdsc.smlm.results.Calibration;
import gdsc.smlm.results.MemoryPeakResults;
import gdsc.smlm.utils.Statistics;
import gdsc.smlm.utils.StoredDataStatistics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

public class BlinkEstimatorTest
{
	private RandomGenerator rand = new Well19937c(System.currentTimeMillis() + System.identityHashCode(this));

	// Set to sensible simulation parameters
	double diffusionRate = 0.25; // pixels^2/sec
	double pixelPitch = 107;
	int msPerFrame = 1000;
	double photons = 1000;
	float psfWidth = 1.2f;

	final double relativeError = 0.2;
	final double minPhotons = 20;
	final double pDelete = 0;
	final double pAdd = 0;

	int LOW = 0;
	int MEDIUM = 1;
	int HIGH = 2;

	int MIN_FITTED_POINTS = 5;
	int MAX_FITTED_POINTS = 15;
	double[] nBlinks = { 0.5, 1.5, 4 };
	double[] tOn = { 1.5, 3, 8 };
	double[] tOff = { 2.5, 5, 10 };

	// If true then test against the real population statistics. 
	// If false then test against the sampled statistics (i.e. using integer frames).
	// Note: When false the success rate is very low so the estimation method does actually 
	// account for the integer frame sampling and get the population statistics.
	boolean usePopulationStatistics = true;

	@Test
	public void canEstimateBlinkingFromSimulationWithLowNBlinksAndMediumOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[LOW], tOn[MEDIUM], tOff[MEDIUM], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndMediumOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[MEDIUM], tOn[MEDIUM], tOff[MEDIUM], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithHighNBlinksAndMediumOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[HIGH], tOn[MEDIUM], tOff[MEDIUM], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithLowNBlinksAndHighOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[LOW], tOn[HIGH], tOff[HIGH], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndHighOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[MEDIUM], tOn[HIGH], tOff[HIGH], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithHighNBlinksAndHighOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[HIGH], tOn[HIGH], tOff[HIGH], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithLowNBlinksAndLowOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[LOW], tOn[LOW], tOff[LOW], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithMediumNBlinksAndLowOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[MEDIUM], tOn[LOW], tOff[LOW], particles, fixedFraction, false, true);
	}

	@Test
	public void canEstimateBlinkingFromSimulationWithHighNBlinksAndLowOnOffTimesWithFixedMolecules()
	{
		int particles = 1000;
		double fixedFraction = 1;
		estimateBlinking(nBlinks[HIGH], tOn[LOW], tOff[LOW], particles, fixedFraction, false, true);
	}

	@Test
	public void findOptimalFittedPoints()
	{
		int particles = 1000;
		double fixedFraction = 1;

		for (boolean timeAtLowerBound : new boolean[] { false })
		{
			int[] count = new int[MAX_FITTED_POINTS + 1];
			int tests = 0;
			for (int run = 0; run < 3; run++)
			{
				for (double n : nBlinks)
				{
					for (int i = 0; i < tOn.length; i++)
					{
						tests++;
						TreeSet<Integer> ok = estimateBlinking(n, tOn[i], tOff[i], particles, fixedFraction,
								timeAtLowerBound, false);
						for (int fit : ok)
							count[fit]++;
					}
				}
			}
			System.out.printf("Time@LowerBound = %b\n", timeAtLowerBound);
			for (int nFittedPoints = MIN_FITTED_POINTS; nFittedPoints <= MAX_FITTED_POINTS; nFittedPoints++)
			{
				System.out.printf("%2d = %2d/%2d |", nFittedPoints, count[nFittedPoints], tests);
				for (int i = 0; i < count[nFittedPoints]; i++)
					System.out.printf("-");
				System.out.printf("\n");
			}
		}
	}

	private TreeSet<Integer> estimateBlinking(double nBlinks, double tOn, double tOff, int particles,
			double fixedFraction, boolean timeAtLowerBound, boolean doAssert)
	{
		SpatialIllumination activationIllumination = new UniformIllumination(100);
		int totalSteps = 100;
		double eAct = totalSteps * 0.3 * activationIllumination.getAveragePhotons();

		ImageModel imageModel = new ActivationEnergyImageModel(eAct, activationIllumination, tOn, tOff, nBlinks, 0, 0);
		imageModel.setRandomGenerator(rand);

		double[] max = new double[] { 256, 256, 32 };
		double[] min = new double[3];
		SpatialDistribution distribution = new UniformDistribution(min, max, rand.nextInt());
		List<CompoundMoleculeModel> compounds = new ArrayList<CompoundMoleculeModel>(1);
		CompoundMoleculeModel c = new CompoundMoleculeModel(1, 0, 0, 0, Arrays.asList(new MoleculeModel(0, 0, 0, 0)));
		c.setDiffusionRate(diffusionRate);
		compounds.add(c);

		List<CompoundMoleculeModel> molecules = imageModel.createMolecules(compounds, particles, distribution, false);

		// Activate fluorophores
		List<? extends FluorophoreSequenceModel> fluorophores = imageModel.createFluorophores(molecules, totalSteps);

		totalSteps = checkTotalSteps(totalSteps, fluorophores);

		imageModel.setUseGridWalk(false);
		List<LocalisationModel> localisations = imageModel.createImage(compounds, fixedFraction, totalSteps, photons,
				0.5, false);

		//		// Remove localisations to simulate missed counts. 
		//		List<LocalisationModel> newLocalisations = new ArrayList<LocalisationModel>(localisations.size());
		//		boolean[] id = new boolean[fluorophores.size() + 1];
		//		Statistics photonStats = new Statistics();
		//		for (LocalisationModel l : localisations)
		//		{
		//			photonStats.add(l.getIntensity());
		//			// Remove by intensity threshold and optionally at random.
		//			if (l.getIntensity() < minPhotons || rand.nextDouble() < pDelete)
		//				continue;
		//			newLocalisations.add(l);
		//			id[l.getId()] = true;
		//		}
		//		localisations = newLocalisations;
		//		System.out.printf("Photons = %f\n", photonStats.getMean());
		//
		//		List<FluorophoreSequenceModel> newFluorophores = new ArrayList<FluorophoreSequenceModel>(fluorophores.size());
		//		for (FluorophoreSequenceModel f : fluorophores)
		//		{
		//			if (id[f.getId()])
		//				newFluorophores.add(f);
		//		}
		//		fluorophores = newFluorophores;

		MemoryPeakResults results = new MemoryPeakResults();
		results.setCalibration(new Calibration(pixelPitch, 1, msPerFrame));
		for (LocalisationModel l : localisations)
		{
			// Remove by intensity threshold and optionally at random.
			if (l.getIntensity() < minPhotons || rand.nextDouble() < pDelete)
				continue;
			float[] params = new float[7];
			params[Gaussian2DFunction.X_POSITION] = (float) l.getX();
			params[Gaussian2DFunction.Y_POSITION] = (float) l.getY();
			params[Gaussian2DFunction.X_WIDTH] = params[Gaussian2DFunction.Y_WIDTH] = psfWidth;
			params[Gaussian2DFunction.AMPLITUDE] = (float) (l.getIntensity() / (2 * Math.PI * psfWidth * psfWidth));
			results.add(l.getTime(), 0, 0, 0, 0, 0, params, null);
		}

		// Add random localisations
		for (int i = (int) (localisations.size() * pAdd); i-- > 0;)
		{
			float[] params = new float[7];
			params[Gaussian2DFunction.X_POSITION] = (float) (rand.nextDouble() * max[0]);
			params[Gaussian2DFunction.Y_POSITION] = (float) (rand.nextDouble() * max[1]);
			params[Gaussian2DFunction.X_WIDTH] = params[Gaussian2DFunction.Y_WIDTH] = psfWidth;
			// Intensity doesn't matter at the moment for tracing
			params[Gaussian2DFunction.AMPLITUDE] = (float) (photons / (2 * Math.PI * psfWidth * psfWidth));
			results.add(1 + rand.nextInt(totalSteps), 0, 0, 0, 0, 0, params, null);
		}

		// Get actual simulated stats ...
		Statistics statsNBlinks = new Statistics();
		Statistics statsTOn = new Statistics();
		Statistics statsTOff = new Statistics();
		Statistics statsSampledNBlinks = new Statistics();
		Statistics statsSampledTOn = new Statistics();
		StoredDataStatistics statsSampledTOff = new StoredDataStatistics();
		for (FluorophoreSequenceModel f : fluorophores)
		{
			statsNBlinks.add(f.getNumberOfBlinks());
			statsTOn.add(f.getOnTimes());
			statsTOff.add(f.getOffTimes());
			int[] on = f.getSampledOnTimes();
			statsSampledNBlinks.add(on.length);
			statsSampledTOn.add(on);
			statsSampledTOff.add(f.getSampledOffTimes());
		}

		System.out.printf("N = %d (%d), N-blinks = %f, tOn = %f, tOff = %f, Fixed = %f\n", fluorophores.size(),
				localisations.size(), nBlinks, tOn, tOff, fixedFraction);
		System.out.printf("Actual N-blinks = %f (%f), tOn = %f (%f), tOff = %f (%f), 95%% = %f, max = %f\n",
				statsNBlinks.getMean(), statsSampledNBlinks.getMean(), statsTOn.getMean(), statsSampledTOn.getMean(),
				statsTOff.getMean(), statsSampledTOff.getMean(), statsSampledTOff.getStatistics().getPercentile(95),
				statsSampledTOff.getStatistics().getMax());
		System.out.printf("-=-=--=-\n");

		BlinkEstimator be = new BlinkEstimator();
		be.maxDarkTime = (int) (tOff * 10);
		be.msPerFrame = msPerFrame;
		be.relativeDistance = false;
		double d = ImageModel.getRandomMoveDistance(diffusionRate);
		be.searchDistance = (fixedFraction < 1) ? Math.sqrt(2 * d * d) * 3 : 0;
		be.timeAtLowerBound = timeAtLowerBound;
		be.showPlots = false;

		//Assert.assertTrue("Max dark time must exceed the dark time of the data (otherwise no plateau)",
		//		be.maxDarkTime > statsSampledTOff.getStatistics().getMax());

		int nMolecules = fluorophores.size();
		if (usePopulationStatistics)
		{
			nBlinks = statsNBlinks.getMean();
			tOff = statsTOff.getMean();
		}
		else
		{
			nBlinks = statsSampledNBlinks.getMean();
			tOff = statsSampledTOff.getMean();
		}

		// See if any fitting regime gets a correct answer
		TreeSet<Integer> ok = new TreeSet<Integer>();
		for (int nFittedPoints = MIN_FITTED_POINTS; nFittedPoints <= MAX_FITTED_POINTS; nFittedPoints++)
		{
			be.nFittedPoints = nFittedPoints;
			be.computeBlinkingRate(results, true);

			double moleculesError = DoubleEquality.relativeError(nMolecules, be.getNMolecules());
			double blinksError = DoubleEquality.relativeError(nBlinks, be.getNBlinks());
			double offError = DoubleEquality.relativeError(tOff * msPerFrame, be.getTOff());
			System.out.printf("Error %d: N = %f, blinks = %f, tOff = %f : %f\n", nFittedPoints, moleculesError,
					blinksError, offError, (moleculesError + blinksError + offError) / 3);

			if (moleculesError < relativeError && blinksError < relativeError && offError < relativeError)
			{
				ok.add(nFittedPoints);
				System.out.printf("-=-=--=-\n");
				System.out.printf("*** Correct at %d fitted points ***\n", nFittedPoints);
				if (doAssert)
					break;
			}

			//if (!be.isIncreaseNFittedPoints())
			//	break;
		}

		System.out.printf("-=-=--=-\n");

		if (doAssert)
			Assert.assertFalse(ok.isEmpty());

		//Assert.assertEquals("Invalid N-blinks", nBlinks, be.getNBlinks(), nBlinks * relativeError);
		//Assert.assertEquals("Invalid N-molecules", fluorophores.size(), be.getNMolecules(), fluorophores.size() *	relativeError);
		//Assert.assertEquals("Invalid t-off", tOff * msPerFrame, be.getTOff(), tOff * msPerFrame * relativeError);
		return ok;
	}

	private int checkTotalSteps(int totalSteps, List<? extends FluorophoreSequenceModel> fluorophores)
	{
		for (FluorophoreSequenceModel f : fluorophores)
		{
			if (totalSteps < f.getEndTime())
				totalSteps = (int) (f.getEndTime() + 1);
		}
		return totalSteps;
	}
}

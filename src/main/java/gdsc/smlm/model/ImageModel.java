package gdsc.smlm.model;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.utils.Maths;
import gdsc.smlm.utils.StoredDataStatistics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

/**
 * Contains a model for an image of blinking fluorophores.
 * <p>
 * Based on the work of Coltharp et al (2012) Accurate Construction of photoactivated localization microscopy images for
 * quantitative measurements. PLOS One 7, Issue 12, pp 1-15
 */
public abstract class ImageModel
{
	protected double tOn;
	protected double tOff, tOff2;
	protected double nBlinks, nBlinks2;

	/**
	 * Specifies the maximum number of frames that will be simulated in {@link #createFluorophores(List, int)}. Any
	 * activation time above this limit returned by {@link #createActivationTime(double[])} will be ignored and there
	 * will be no call to {@link #createFluorophore(int, double[], double)}.
	 */
	protected int frameLimit;

	private RandomGenerator random;
	private RandomDataGenerator randomGenerator;
	private RealDistribution photonDistribution;
	private SpatialDistribution confinementDistribution = null;
	private int confinementAttempts = 5;
	private boolean useGridWalk = true;
	private boolean useGeometricDistribution = false;
	private boolean photonBudgetPerFrame = true;
	private boolean rotation2D = false;

	/**
	 * Construct a new image model
	 * 
	 * @param tOn
	 *            Average on-state time
	 * @param tOff
	 *            Average off-state time for the first dark state
	 * @param tOff2
	 *            Average off-state time for the second dark state
	 * @param nBlinks
	 *            Average number of blinks int the first dark state (used for each burst between second dark states)
	 * @param nBlinks2
	 *            Average number of blinks into the second dark state
	 */
	public ImageModel(double tOn, double tOff, double tOff2, double nBlinks, double nBlinks2)
	{
		init(tOn, tOff, tOff2, nBlinks, nBlinks2,
				new Well19937c(System.currentTimeMillis() + System.identityHashCode(this)));
	}

	private void init(double tOn, double tOff, double tOff2, double nBlinks, double nBlinks2, RandomGenerator rand)
	{
		checkParameter("tOn", tOn);
		checkParameter("tOff", tOff);
		checkParameter("tOff2", tOff2);
		checkParameter("nBlinks", nBlinks);
		checkParameter("nBlinks2", nBlinks2);
		this.tOn = tOn;
		this.tOff = tOff;
		this.tOff2 = tOff2;
		this.nBlinks = nBlinks;
		this.nBlinks2 = nBlinks2;
		setRandomGenerator(rand);
	}

	protected void checkParameter(String param, double value)
	{
		if (value < 0)
			throw new IllegalArgumentException(param + " cannot be less than 0");
	}

	/**
	 * @return the tOn
	 */
	public double gettOn()
	{
		return tOn;
	}

	/**
	 * @return the tOff
	 */
	public double gettOff()
	{
		return tOff;
	}

	/**
	 * @return the tOff2
	 */
	public double gettOff2()
	{
		return tOff2;
	}

	/**
	 * @return the nBlinks
	 */
	public double getnBlinks()
	{
		return nBlinks;
	}

	/**
	 * @return the nBlinks2
	 */
	public double getnBlinks2()
	{
		return nBlinks2;
	}

	/**
	 * @return The random data generator
	 */
	protected RandomDataGenerator getRandom()
	{
		return randomGenerator;
	}

	/**
	 * Creates N particles by sampling from the distribution. If the distribution returns null (no coordinates) then no
	 * more attempts to generate particles is made.
	 * <p>
	 * The compound's molecule coordinates will be updated by calling {@link CompoundMoleculeModel#centre()}.
	 * 
	 * @param compounds
	 * @param particles
	 * @param distribution
	 * @param rotate
	 * @return
	 */
	public List<CompoundMoleculeModel> createMolecules(List<CompoundMoleculeModel> compounds, int particles,
			SpatialDistribution distribution, boolean rotate)
	{
		if (compounds == null || compounds.isEmpty())
			throw new IllegalArgumentException("Compounds cannot be null or empty");

		// Sum the fractions
		double total = 0;
		for (CompoundMoleculeModel c : compounds)
		{
			total += c.getFraction();
			c.centre();
		}

		ArrayList<CompoundMoleculeModel> molecules = new ArrayList<CompoundMoleculeModel>(particles);
		//int[] count = new int[compounds.size()];
		for (int i = 1; i <= particles; i++)
		{
			// Get the next random compound
			double d = random.nextDouble() * total;
			//int j=0;
			for (CompoundMoleculeModel c : compounds)
			{
				d -= c.getFraction();
				if (d <= 0)
				{
					final double[] xyz = distribution.next();
					if (xyz == null || xyz.length < 3)
						return molecules;

					//count[j]++;
					CompoundMoleculeModel m = new CompoundMoleculeModel(i, xyz, copyMolecules(c), false);
					m.setDiffusionRate(c.getDiffusionRate());
					if (rotate)
						rotate(m);
					molecules.add(m);
					//System.out.printf("XYZ = %f,%f,%f\n", m.getX(), m.getY(), m.getZ());
					break;
				}
				//j++;
			}
		}
		//for (int j=0; j<count.length; j++)
		//	System.out.printf("Compound %d = %d\n", j+1, count[j]);
		return molecules;
	}

	private List<MoleculeModel> copyMolecules(CompoundMoleculeModel c)
	{
		int n = c.getSize();
		List<MoleculeModel> list = new ArrayList<MoleculeModel>(n);
		while (n-- > 0)
		{
			MoleculeModel m = c.getMolecule(n);
			list.add(new MoleculeModel(n + 1, m.getMass(), Arrays.copyOf(m.getCoordinates(), 3)));
		}
		return list;
	}

	private void rotate(CompoundMoleculeModel m)
	{
		if (rotation2D)
			m.rotateRandomAngle(CompoundMoleculeModel.Z_AXIS, 180, random);
		else
			m.rotateRandom(180, random);
	}

	/**
	 * Generates fluorophores for the molecules spatial positions. Only simulate up to the given maximum number of
	 * frames.
	 * <p>
	 * Replace the CompoundMoleculeModel objects in the list with new compounds containing the generated
	 * FluorophoreSequenceModel objects. If no fluorophores can be generated for a compound then it is removed. Since
	 * the fluorophores are part of a compound their coordinates are relative to the compound centre of mass.
	 * <p
	 * Note that the activation energy is sampled at the spatial position without movement. This is therefore an
	 * approximation of the energy the molecule would receive if it were moving.
	 * 
	 * @param molecules
	 * @param frames
	 * @return
	 */
	public List<? extends FluorophoreSequenceModel> createFluorophores(List<CompoundMoleculeModel> molecules, int frames)
	{
		frameLimit = frames;
		ArrayList<FluorophoreSequenceModel> list = new ArrayList<FluorophoreSequenceModel>(molecules.size());
		for (int i = 0; i < molecules.size();)
		{
			CompoundMoleculeModel c = molecules.get(i);
			List<FluorophoreSequenceModel> fluorophores = new ArrayList<FluorophoreSequenceModel>(c.getSize());
			List<MoleculeModel> removed = new ArrayList<MoleculeModel>(c.getSize());
			for (int n = c.getSize(); n-- > 0;)
			{
				// Molecule Id is ignored since we renumber the sorted collection at the end
				final double[] xyz = c.getRelativeCoordinates(n);
				final double tAct = createActivationTime(xyz);
				FluorophoreSequenceModel f = null;
				if (frames == 0 || tAct < frames)
					f = createFluorophore(0, xyz, tAct);
				if (f != null && f.getCoordinates() != null && f.getEndTime() > f.getStartTime())
				{
					list.add(f);
					fluorophores.add(f);
				}
				else
				{
					// If we remove the molecule from the compound then the rotation will not be around
					// the true COM but the new COM of the reduced size compound. Keep these molecules
					// so the rotation is valid.
					// Note: Fluorophores have mass zero. If this changes then this will have to be updated. 
					removed.add(new MoleculeModel(0, c.getRelativeCoordinates(n)));
				}
			}
			if (fluorophores.isEmpty())
			{
				molecules.remove(i);
				// No increment as the list has been reduced in size
			}
			else
			{
				removed.addAll(fluorophores);
				CompoundMoleculeModel newC = new CompoundMoleculeModel(c.getId(), c.getCoordinates(), removed, false);
				newC.setDiffusionRate(c.getDiffusionRate());
				molecules.set(i, newC);

				// Debug the fluorophore position
				//c = molecules.get(i);
				//System.out.printf("fluorophores XYZ = %f,%f,%f (%d)\n", c.getX(), c.getY(), c.getZ(), fluorophores.size());
				//for (int j = 0; j < c.getSize(); j++)
				//{
				//	double[] xyz = c.getCoordinates(j);
				//	System.out.printf("fluorophore %d XYZ = %f,%f,%f\n", j + 1, xyz[0], xyz[1], xyz[2]);
				//}

				i++;
			}
		}
		// Sort by time
		Collections.sort(list);
		// Renumber
		int id = 1;
		for (FluorophoreSequenceModel f : list)
			f.setId(id++);
		return list;
	}

	/**
	 * Create a fluorophore activation time for the given position. Derived classes can check the {@link #frameLimit}
	 * variable to determine the limit of the simulation. If non-zero only activation times below this will be used to
	 * generate fluorophores.
	 * 
	 * @param xyz
	 * @return the activation time
	 */
	protected abstract double createActivationTime(double[] xyz);

	/**
	 * Create a fluorophore with the given id, position and activation time
	 * 
	 * @param id
	 * @param xyz
	 * @param tAct
	 *            the activation time (generated using {@link #createActivationTime(double[])})
	 * @return a fluorophore
	 */
	protected abstract FluorophoreSequenceModel createFluorophore(int id, double[] xyz, double tAct);

	/**
	 * Add to the list but link up the continuous pulses with previous/next pointers
	 * 
	 * @param localisations
	 * @param models
	 */
	private void linkLocalisations(List<LocalisationModel> localisations, LocalisationModel[] models)
	{
		LocalisationModel previous = null;
		for (int i = 0; i < models.length; i++)
		{
			if (models[i] != null)
			{
				models[i].setPrevious(previous);
				localisations.add(models[i]);
			}
			previous = models[i];
		}
	}

	private void generateOnTimes(int maxFrames, final double frameInterval, List<double[]> bursts, int sequenceStartT,
			double[] onTime, int[] state)
	{
		for (int i = 0; i < bursts.size(); i++)
		{
			double[] burst = bursts.get(i);
			int on = (int) (burst[0] / frameInterval) - sequenceStartT;
			int off = (int) (burst[1] / frameInterval) - sequenceStartT;

			// If the off-time is truncated the modulus function used below is
			// not valid.
			boolean offTimeTrucated = false;

			if (on >= onTime.length)
				break;
			if (off >= onTime.length)
			{
				burst[1] = maxFrames * frameInterval;
				off = onTime.length - 1;
				offTimeTrucated = true;
			}

			if (on == off)
			{
				onTime[on] += (burst[1] - burst[0]);
				state[on] |= LocalisationModel.SINGLE;
			}
			else
			{
				for (int t = on; t <= off; t++)
				{
					if (t > on)
						state[t] |= LocalisationModel.PREVIOUS;
					if (t < off)
						state[t] |= LocalisationModel.NEXT;

					// Get signal on-time
					if (t == on)
					{
						// Signal turned on this frame
						onTime[t] += frameInterval - (burst[0] % frameInterval);
					}
					else if (t == off)
					{
						// Signal turned off this frame
						if (offTimeTrucated)
						{
							// The off-time was truncated so the burst is on for
							// the entire duration
							onTime[t] += 1;
						}
						else
						{
							onTime[t] += (burst[1] % frameInterval);
						}
					}
					else
					{
						// Continuous signal
						onTime[t] = frameInterval;
						state[t] |= LocalisationModel.CONTINUOUS;
					}
				}
			}
		}
	}

	private void sortByTime(List<LocalisationModel> localisations)
	{
		// Sort by time-only (not signal as well which is the default).
		// This means the localisations are in the same order as the input
		// fluorophores.
		Collections.sort(localisations, new Comparator<LocalisationModel>()
		{
			public int compare(LocalisationModel o1, LocalisationModel o2)
			{
				if (o1.getTime() < o2.getTime())
					return -1;
				if (o1.getTime() > o2.getTime())
					return 1;
				return 0;
			}
		});
	}

	/**
	 * Simulate an image of fluorophores. The total amount of time a fluorophore is on (i.e. sum of tOn) is used to
	 * create the signal strength using the specified correlation.
	 * <p>
	 * A second set of fluorophores, independent of the first, are generated using
	 * {@link #createFluorophores(List, int)}. The correlated on-times will be created by combining the times using the
	 * provided correlation (r):
	 * 
	 * <pre>
	 * // X1 : Fluorophore total on-times
	 * // X2 : Independently generated fluorophore total on-times
	 * // r  : correlation
	 * a = sqrt(1 - r * r)
	 * newX = (r * X1 + a * X2) / (r + a)
	 * </pre>
	 * 
	 * The signal is proportional to newly generated on-times (newX) with an average of the specified photon budget.
	 * <p>
	 * The photon budget can either be distributed evenly over the fluorophore lifetime or per frame (see
	 * {@link #isPhotonBudgetPerFrame()}). Each frame signal output will be subject to Poisson sampling.
	 * <p>
	 * If the input correlation is zero then the number of photons will be sampled from the configured photon
	 * distribution or, if this is null, will be uniform.
	 * <p>
	 * A random fraction of the fluorophores will move using a random walk with the diffusion coefficient defined in the
	 * compound.
	 * 
	 * @param compoundFluorophores
	 *            The compounds containing the fluorophores
	 * @param fixedFraction
	 *            The fraction of molecules that are fixed
	 * @param maxFrames
	 *            The maximum frame for the simulation
	 * @param photonBudget
	 *            the average number of photons per frame/lifetime; see {@link #isPhotonBudgetPerFrame()}
	 * @param correlation
	 *            The correlation between the number of photons and the total on-time of the fluorophore
	 * @param rotate
	 *            Rotate the molecule per frame
	 * @return the localisations
	 */
	public List<LocalisationModel> createImage(List<CompoundMoleculeModel> compoundFluorophores, double fixedFraction,
			int maxFrames, double photonBudget, double correlation, boolean rotate)
	{
		List<LocalisationModel> localisations = new ArrayList<LocalisationModel>();
		// Extract the fluorophores in all the compounds
		ArrayList<FluorophoreSequenceModel> fluorophores = new ArrayList<FluorophoreSequenceModel>(
				compoundFluorophores.size());
		for (CompoundMoleculeModel c : compoundFluorophores)
		{
			for (int i = c.getSize(); i-- > 0;)
				if (c.getMolecule(i) instanceof FluorophoreSequenceModel)
				{
					fluorophores.add((FluorophoreSequenceModel) c.getMolecule(i));
				}
		}
		final int nMolecules = fluorophores.size();

		// Check the correlation bounds. 
		// Correlation for photons per frame verses total on time should be negative.
		// Correlation for total photons verses total on time should be positive.
		double r;
		if (photonBudgetPerFrame)
			r = (correlation < -1) ? -1 : (correlation > 0) ? 0 : correlation;
		else
			r = (correlation < 0) ? 0 : (correlation > 1) ? 1 : correlation;

		// Create a photon budget for each fluorophore
		double[] photons = new double[nMolecules];

		// Generate a second set of on times using the desired correlation
		if (r != 0)
		{
			// Observations on real data show:
			// - there is a weak positive correlation between total photons and time
			// - There is a weak negative correlation between photons per frame and total on-time

			// How to generate a random correlation:
			// http://www.uvm.edu/~dhowell/StatPages/More_Stuff/CorrGen.html
			// http://stats.stackexchange.com/questions/13382/how-to-define-a-distribution-such-that-draws-from-it-correlate-with-a-draw-from

			// Used for debugging the correlation
			//StoredDataStatistics onTime = new StoredDataStatistics();
			//StoredDataStatistics allPhotons = new StoredDataStatistics();
			//PearsonsCorrelation c = new PearsonsCorrelation();

			// Create a second set of fluorophores. This is used to generate the correlated photon data
			List<FluorophoreSequenceModel> fluorophores2 = new ArrayList<FluorophoreSequenceModel>();
			while (fluorophores2.size() < fluorophores.size())
			{
				FluorophoreSequenceModel f = createFluorophore(0, new double[] { 0, 0, 0 }, maxFrames);
				if (f != null)
					fluorophores2.add(f);
			}

			final double a = Math.sqrt(1 - r * r);

			// Q - How to generate a negative correlation?
			// Generate a positive correlation then invert the data and shift to the same range
			boolean invert = (r < 0);
			r = Math.abs(r);

			StoredDataStatistics correlatedOnTime = new StoredDataStatistics();
			for (int i = 0; i < nMolecules; i++)
			{
				final double X = getTotalOnTime(fluorophores.get(i));
				final double Z = getTotalOnTime(fluorophores2.get(i));
				//onTime.add(X);
				correlatedOnTime.add((r * X + a * Z) / (r + a));
			}

			if (invert)
			{
				final double min = correlatedOnTime.getStatistics().getMin();
				final double range = correlatedOnTime.getStatistics().getMax() - min;
				StoredDataStatistics newCorrelatedOnTime = new StoredDataStatistics();
				for (double d : correlatedOnTime.getValues())
				{
					newCorrelatedOnTime.add(range - d + min);
				}
				correlatedOnTime = newCorrelatedOnTime;
			}

			// Get the average on time from the correlated sample. 
			// Using the population value (tOn * (1+nBlinks)) over-estimates the on time.
			final double averageTotalTOn = correlatedOnTime.getMean();

			// Now create the localisations
			double[] correlatedOnTimes = correlatedOnTime.getValues();
			for (int i = 0; i < nMolecules; i++)
			{
				// Generate photons using the correlated on time data
				double p = photonBudget * correlatedOnTimes[i] / averageTotalTOn;

				//final double X = getTotalOnTime(fluorophores.get(i));
				//double p = (X > 0) ? photonBudget * correlatedOnTimes[i] / X : 0;
				//double p = (X > 0) ? randomGenerator.nextGamma(photonBudget, correlatedOnTimes[i] / X) : 0;
				//double p = randomGenerator.nextGamma(photonBudget, correlatedOnTimes[i] / X);
				//allPhotons.add(p);

				photons[i] = p;
			}

			//System.out.printf("t = %f, p = %f, R = %f\n", averageTotalTOn, allPhotons.getMean(),
			//		c.correlation(onTime.getValues(), allPhotons.getValues()));
		}
		else
		{
			// Sample from the provided distribution. Do not over-write the class level distribution to allow 
			// running again with a different shape parameter / photon budget.
			if (photonDistribution != null)
			{
				// Ensure the custom distribution is scaled to the correct photon budget
				final double photonScale = photonBudget / photonDistribution.getNumericalMean();

				// Generate photons sampling from the photon budget
				for (int i = 0; i < nMolecules; i++)
				{
					photons[i] = photonDistribution.sample() * photonScale;
				}
			}
			else
			{
				// No distribution so use the same number for all
				Arrays.fill(photons, photonBudget);
			}
		}

		int photonIndex = 0;
		for (CompoundMoleculeModel c : compoundFluorophores)
		{
			final boolean fixed = (random.nextDouble() < fixedFraction);
			photonIndex += createLocalisation(c, localisations, fixed, maxFrames, photons, photonIndex, !fixed &&
					rotate);
		}

		sortByTime(localisations);

		return localisations;
	}

	private double getTotalOnTime(FluorophoreSequenceModel f)
	{
		return Maths.sum(f.getOnTimes());
	}

	/**
	 * Create localisations for each time frame using the fluorophores in the compound. The compound is moved using the
	 * specified diffusion rate.
	 * 
	 * @param compound
	 *            Compound containing one or more fluorophores
	 * @param localisations
	 *            Output list to put all localisations
	 * @param fixed
	 * @param maxFrames
	 *            The maximum number of frames. Fluorophores still active beyond this frame will not be represented
	 * @param photons
	 *            Array of the photon budget for all the fluorophores
	 * @param photonIndex
	 *            The current index within the photon budget array
	 * @param rotate
	 *            Rotate at each step
	 * @return The number of fluorophores that were drawn
	 */
	private int createLocalisation(CompoundMoleculeModel compound, List<LocalisationModel> localisations,
			boolean fixed, int maxFrames, double[] photons, int photonIndex, boolean rotate)
	{
		// -=-=-
		// TODO: 
		// Account for the variable nature of fluorophore intensities in real data?
		// Perhaps add quantum efficiency simulation using a Binomial distribution?
		// -=-=-
		final double frameInterval = 1;

		// Extract the fluorophores
		ArrayList<FluorophoreSequenceModel> fluorophores = new ArrayList<FluorophoreSequenceModel>();
		for (int i = compound.getSize(); i-- > 0;)
			if (compound.getMolecule(i) instanceof FluorophoreSequenceModel)
			{
				fluorophores.add((FluorophoreSequenceModel) compound.getMolecule(i));
			}
		final int nFluorophores = fluorophores.size();

		// Get the duration of each fluorophore sequence
		List<List<double[]>> bursts = new ArrayList<List<double[]>>(nFluorophores);
		int[] sequenceStartT = new int[nFluorophores];
		int[] sequenceEndT = new int[nFluorophores];

		for (int i = 0; i < nFluorophores; i++)
		{
			FluorophoreSequenceModel fluorophore = fluorophores.get(i);
			bursts.add(fluorophore.getBurstSequence());
			sequenceStartT[i] = (int) (fluorophore.getStartTime() / frameInterval);
			sequenceEndT[i] = (int) (fluorophore.getEndTime() / frameInterval);
		}

		// Get the maximum and minimum start and end times
		int sequenceStart = Maths.min(sequenceStartT);
		int sequenceEnd = Maths.max(sequenceEndT);

		// Time-range check. Note that the final frames are 1-based
		if (sequenceStart > maxFrames - 1)
			return nFluorophores;
		if (sequenceEnd > maxFrames - 1)
			sequenceEnd = maxFrames - 1;

		// For each fluorophore build a graph of on-time in frame intervals
		final int nSteps = sequenceEnd - sequenceStart + 1;
		double[][] onTime = new double[nFluorophores][nSteps];
		int[][] state = new int[nFluorophores][nSteps];
		double[] totalOnTime = new double[nFluorophores];
		double[] photonBudget = new double[nFluorophores];
		for (int i = 0; i < nFluorophores; i++)
		{
			generateOnTimes(maxFrames, frameInterval, bursts.get(i), sequenceStart, onTime[i], state[i]);
			totalOnTime[i] = Maths.sum(onTime[i]);
			photonBudget[i] = photons[i + photonIndex];
		}

		// Convert diffusion co-efficient into the standard deviation for the random walk
		final double diffusionRate = (fixed) ? 0 : getRandomMoveDistance(compound.getDiffusionRate());
		//System.out.printf("N=%d, rate=%f\n", compound.getSize(), diffusionRate);

		// For the duration of the sequence, move the molecule and output
		// position and intensity
		LocalisationModel[][] models = new LocalisationModel[nFluorophores][nSteps];

		if (confinementDistribution != null)
			confinementDistribution.initialise(compound.getCoordinates());
		//int pass=0,fail=0;

		for (int t = sequenceStart, step = 0; t <= sequenceEnd; t++, step++)
		{
			for (int i = 0; i < nFluorophores; i++)
			{
				if (onTime[i][step] > 0 && photonBudget[i] > 0)
				{
					final double intensity;
					if (photonBudgetPerFrame)
					{
						// -=-=-
						// Photon budget is the rate per frame 
						// -=-=-
						// This suffers from the exponential distribution of the tOn which destroys the 
						// shape of the gamma distribution when tOn is short, i.e. may onTime values are less than 1
						//intensity = randomGenerator.nextPoisson(photonBudget * onTime[i]);

						// Ignore the time that the fluorophore is on and simply use the random sample. 
						// This allows the simulation to match observed experimental data.
						intensity = randomGenerator.nextPoisson(photonBudget[i]);
					}
					else
					{
						// -=-=-
						// When photon budget is for the entire lifetime then allocate proportionately 
						// using the time fraction
						// -=-=-
						intensity = randomGenerator.nextPoisson(photonBudget[i] * onTime[i][step] / totalOnTime[i]);
					}

					// TODO - Is this useful? Add orientation brightness
					//intensity *= randomGenerator.nextUniform(0.7, 1.3);

					//rawPhotons.add(intensity);
					models[i][step] = new LocalisationModel(fluorophores.get(i).getId(), t + 1,
							compound.getCoordinates(i), intensity, state[i][step]);
				}
			}

			// Move the compound
			if (diffusionRate > 0)
			{
				if (confinementDistribution != null)
				{
					// TODO
					// Confined diffusion only asks if the molecule moved into a region that is allowed.
					// A better model would be to ask if the move went out of bounds. If yes then the molecules 
					// should be reflected back into the allowed region using the vector of the original move.  

					// This only works because the coordinates are a reference
					double[] xyz = compound.getCoordinates();
					double[] originalXyz = Arrays.copyOf(xyz, 3);
					for (int n = confinementAttempts; n-- > 0;)
					{
						if (useGridWalk)
							compound.walk(diffusionRate, random);
						else
							compound.move(diffusionRate, random);
						if (!confinementDistribution.isWithin(compound.getCoordinates()))
						{
							//fail++;
							// Reset position
							for (int j = 0; j < 3; j++)
								xyz[j] = originalXyz[j];
						}
						else
						{
							//pass++;
							// The move was allowed
							break;
						}
					}
				}
				else
				{
					if (useGridWalk)
						compound.walk(diffusionRate, random);
					else
						compound.move(diffusionRate, random);
				}
			}

			if (rotate)
				rotate(compound);
		}
		//System.out.printf("Pass = %d, fail = %d\n", pass, fail);

		for (int i = 0; i < nFluorophores; i++)
		{
			linkLocalisations(localisations, models[i]);
		}

		return nFluorophores;
	}

	/**
	 * Set the random generator for creating the image
	 * 
	 * @param random
	 */
	public void setRandomGenerator(RandomGenerator random)
	{
		if (random == null)
			throw new NullPointerException("Random generator must not be null");
		this.random = random;
		this.randomGenerator = new RandomDataGenerator(random);
	}

	/**
	 * Convert diffusion co-efficient (D) into the average step size required for a random diffusion. The step size is
	 * per dimension.
	 * 
	 * @see {@link gdsc.smlm.model.MoleculeModel#move(double, RandomDataGenerator) }
	 * @see {@link gdsc.smlm.model.MoleculeModel#walk(double, RandomGenerator) }
	 * @param diffusionRateInPixelsPerStep
	 * @return
	 */
	public static double getRandomMoveDistance(double diffusionRateInPixelsPerStep)
	{
		// Convert diffusion co-efficient into the standard deviation for the random move in each dimension
		// For 1D diffusion: sigma^2 = 2D
		//                   sigma = sqrt(2D)
		return Math.sqrt(2 * diffusionRateInPixelsPerStep);
	}

	/**
	 * @return true if the random diffusion will use a grid walk
	 */
	public boolean isUseGridWalk()
	{
		return useGridWalk;
	}

	/**
	 * Set to true if the random diffusion will use a grid walk. The alternative is an off grid random move along a unit
	 * vector.
	 * 
	 * @param useGridWalk
	 *            true if the random diffusion will use a grid walk
	 */
	public void setUseGridWalk(boolean useGridWalk)
	{
		this.useGridWalk = useGridWalk;
	}

	/**
	 * @return the photon distribution used for the fluorophore
	 */
	public RealDistribution getPhotonDistribution()
	{
		return photonDistribution;
	}

	/**
	 * Set the distribution used to generate the photon budget of a fluorophore
	 * 
	 * @param photonDistribution
	 *            the photon distribution to set
	 */
	public void setPhotonDistribution(RealDistribution photonDistribution)
	{
		this.photonDistribution = photonDistribution;
	}

	/**
	 * @return the useGeometricDistribution
	 */
	public boolean isUseGeometricDistribution()
	{
		return useGeometricDistribution;
	}

	/**
	 * @param useGeometricDistribution
	 *            the useGeometricDistribution to set
	 */
	public void setUseGeometricDistribution(boolean useGeometricDistribution)
	{
		this.useGeometricDistribution = useGeometricDistribution;
	}

	/**
	 * @return if true the image will be created using the photon budget per frame
	 */
	public boolean isPhotonBudgetPerFrame()
	{
		return photonBudgetPerFrame;
	}

	/**
	 * Set to true if the photon budget is per frame. The default is for the lifetime of the fluorophore.
	 * 
	 * @param photonBudgetPerFrame
	 *            if true the image will be created using the photon budget per frame
	 */
	public void setPhotonBudgetPerFrame(boolean photonBudgetPerFrame)
	{
		this.photonBudgetPerFrame = photonBudgetPerFrame;
	}

	/**
	 * Set the distribution used to confine any diffusing molecules
	 * 
	 * @param confinementDistribution
	 */
	public void setConfinementDistribution(SpatialDistribution confinementDistribution)
	{
		this.confinementDistribution = confinementDistribution;
	}

	/**
	 * @return the confinementAttempts
	 */
	public int getConfinementAttempts()
	{
		return confinementAttempts;
	}

	/**
	 * Set the number of attempts to move a diffusing molecule. If the molecule cannot be successfully moved within the
	 * confinement distribution then it remains fixed.
	 * 
	 * @param confinementAttempts
	 *            the confinementAttempts to set
	 */
	public void setConfinementAttempts(int confinementAttempts)
	{
		this.confinementAttempts = confinementAttempts;
	}

	/**
	 * Set to true if rotation is around the z-axis (i.e. rotation in the 2D plane)
	 * 
	 * @return True if using 2D rotation
	 */
	public boolean isRotation2D()
	{
		return rotation2D;
	}

	/**
	 * Set to true to rotate around the z-axis (i.e. rotation in the 2D plane)
	 * 
	 * @param rotation2D
	 *            True if using 2D rotation
	 */
	public void setRotation2D(boolean rotation2D)
	{
		this.rotation2D = rotation2D;
	}
}

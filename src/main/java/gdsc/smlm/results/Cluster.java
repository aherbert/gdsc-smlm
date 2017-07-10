package gdsc.smlm.results;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.util.FastMath;

import gdsc.core.data.utils.TypeConverter;
import gdsc.core.utils.Maths;
import gdsc.smlm.data.config.UnitConfig.DistanceUnit;

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

/**
 * Define a cluster of localisations
 */
public class Cluster implements Comparable<Cluster>
{
	public enum CentroidMethod
	{
		STANDARD, SIGNAL_WEIGHTED
	}

	protected ArrayList<PeakResult> results = new ArrayList<PeakResult>(2);
	private float[] centroid = null;
	private int id;

	public Cluster()
	{
	}

	public Cluster(PeakResult result)
	{
		add(result);
	}

	public int size()
	{
		return results.size();
	}

	public ArrayList<PeakResult> getPoints()
	{
		return results;
	}

	public void add(PeakResult result)
	{
		results.add(result);
		centroid = null;
	}

	public float[] getCentroid(CentroidMethod method)
	{
		if (centroid == null && !results.isEmpty())
		{
			switch (method)
			{
				case SIGNAL_WEIGHTED:
					float[] weights = new float[results.size()];
					int i = 0;
					for (PeakResult result : results)
					{
						weights[i] = Math.abs(result.getSignal());
					}
					// Normalise weights?
					return getCentroid(results, weights);

				case STANDARD:
				default:
					return getCentroid();
			}
		}
		return centroid;
	}

	private float[] getCentroid(ArrayList<PeakResult> results, float[] weights)
	{
		centroid = new float[2];
		double sum = 0;
		int i = 0;
		for (PeakResult result : results)
		{
			final float w = weights[i++];
			sum += w;
			centroid[0] += result.getXPosition() * w;
			centroid[1] += result.getYPosition() * w;
		}
		centroid[0] /= sum;
		centroid[1] /= sum;
		return centroid;
	}

	public float[] getCentroid()
	{
		if (centroid == null && !results.isEmpty())
		{
			centroid = new float[2];
			for (PeakResult result : results)
			{
				centroid[0] += result.getXPosition();
				centroid[1] += result.getYPosition();
			}
			centroid[0] /= results.size();
			centroid[1] /= results.size();
		}
		return centroid;
	}

	/**
	 * Remove the calculated centroid from memory (forces a refresh of the centroid)
	 */
	public void resetCentroid()
	{
		centroid = null;
	}

	/**
	 * @return The standard deviation of the distances from the centroid
	 */
	public float getStandardDeviation()
	{
		if (results.size() < 2)
			return 0;
		getCentroid();
		double ssx = 0;
		for (PeakResult result : results)
		{
			final double dx = result.getXPosition() - centroid[0];
			final double dy = result.getYPosition() - centroid[1];
			final double d2 = dx * dx + dy * dy;
			ssx += d2;
		}
		return (float) Math.sqrt(ssx / (results.size() - 1));
	}

	/**
	 * Calculate the weighted localisation precision using the PC-PALM formula of Sengupta, et al (2013) Nature
	 * Protocols 8, 345.
	 * <p>
	 * Also sets the centroid if it has not been calculated using the signal weighted centre-of-mass
	 * <p>
	 * Note that the PeakResult must have valid values in the precision field, otherwise a value of 1 is used.
	 *
	 * @param converter the converter to convert the distances to nm
	 * @return The weighted localisation precision of the group peak (in nm)
	 */
	public double getLocalisationPrecision(TypeConverter<DistanceUnit> converter)
	{
		if (converter == null || converter.to() != DistanceUnit.NM)
			return 0;
		
		final int n = size();
		if (n == 0)
		{
			centroid = null;
			return 0;
		}

		if (n == 1)
		{
			PeakResult result = results.get(0);
			if (centroid == null)
			{
				centroid = new float[] { result.getXPosition(), result.getYPosition() };
			}
			return checkPrecision(result.getPrecision());
		}

		float[] photons = new float[results.size()];
		int i = 0;
		for (PeakResult result : results)
		{
			photons[i++] = Math.abs(result.getSignal());
		}

		double sumNi = 0;
		i = 0;
		double xm = 0, ym = 0;
		for (PeakResult result : results)
		{
			final float Ni = photons[i++];
			sumNi += Ni;
			xm += result.getXPosition() * Ni;
			ym += result.getYPosition() * Ni;
		}
		xm /= sumNi;
		ym /= sumNi;

		if (centroid == null)
		{
			centroid = new float[] { (float) xm, (float) ym };
		}

		i = 0;
		double sumXi2Ni = 0, sumYi2Ni = 0, sumS2 = 0;
		for (PeakResult result : results)
		{
			final float Ni = photons[i++];

			double dx = converter.convert(result.getXPosition() - xm);
			double dy = converter.convert(result.getYPosition() - ym);

			sumXi2Ni += dx * dx * Ni;
			sumYi2Ni += dy * dy * Ni;
			sumS2 += Maths.pow2(checkPrecision(result.getPrecision())) * Ni;
		}

		double sumNin = sumNi * n;
		double sumS2_sumNin = sumS2 / sumNin;
		double sxm = Math.sqrt(sumXi2Ni / sumNin + sumS2_sumNin) / 1.414213562;
		double sym = Math.sqrt(sumYi2Ni / sumNin + sumS2_sumNin) / 1.414213562;

		double sPeak = FastMath.max(sxm, sym);

		return sPeak;
	}

	private static double checkPrecision(double p)
	{
		return (p > 0 && p < Double.MAX_VALUE) ? p : 1;

	}

	/**
	 * @return The first PeakResult in the cluster (or null)
	 */
	public PeakResult getHead()
	{
		if (results.isEmpty())
			return null;
		return results.get(0);
	}

	/**
	 * @return The last PeakResult in the cluster (or null)
	 */
	public PeakResult getTail()
	{
		if (results.isEmpty())
			return null;
		return results.get(results.size() - 1);
	}

	/**
	 * @return The total signal
	 */
	public double getSignal()
	{
		double sum = 0;
		for (PeakResult result : results)
		{
			sum += result.getSignal();
		}
		return sum;
	}

	/**
	 * Sort in time order
	 */
	public void sort()
	{
		Collections.sort(results);
	}

	public int getId()
	{
		return id;
	}

	public void setId(int id)
	{
		this.id = id;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Cluster that)
	{
		// Sort by ID ascending
		return this.id - that.id;
	}

	/**
	 * Expand any localisations that have a different start and end frame into a series.
	 * Note that this will increase the size of the cluster.
	 * <p>
	 * The results are copies save for the end frame. This makes analysis of the signal invalid as it will have been
	 * increased n-fold for each localisation that spans n frames. The original multi-frame result is removed.
	 */
	public void expandToSingles()
	{
		ArrayList<PeakResult> extra = null;
		ArrayList<PeakResult> remove = null;
		for (PeakResult result : results)
		{
			if (result.getFrame() != result.getEndFrame())
			{
				if (extra == null)
				{
					extra = new ArrayList<PeakResult>();
					remove = new ArrayList<PeakResult>();
				}
				remove.add(result);
				for (int peak = result.getFrame(); peak <= result.getEndFrame(); peak++)
					extra.add(new ExtendedPeakResult(peak, result.origX, result.origY, result.origValue, result.error,
							result.noise, result.params, result.paramStdDevs, peak, result.getId()));
			}
		}
		if (extra == null)
			return;
		for (PeakResult result : remove)
			results.remove(result);
		for (PeakResult result : extra)
			add(result);
	}

	/**
	 * Remove the first and last result. If the size is 2 or less then the new size will be zero.
	 */
	public void removeEnds()
	{
		if (size() <= 2)
		{
			results.clear();
		}
		else
		{
			results = new ArrayList<PeakResult>(results.subList(1, size() - 1));
		}
		resetCentroid();
	}

	/**
	 * @return The mean-squared displacement between adjacent localisations
	 */
	public double getMSD()
	{
		if (size() < 2)
			return 0;
		double msd = 0;
		PeakResult last = null;
		for (PeakResult result : results)
		{
			if (last != null)
			{
				msd += last.distance2(result);
			}
			last = result;
		}
		return msd / (size() - 1);
	}

	/**
	 * @return The mean displacement between adjacent localisations
	 */
	public double getMeanPerFrame()
	{
		if (size() < 2)
			return 0;
		double msd = 0;
		PeakResult last = null;
		for (PeakResult result : results)
		{
			if (last != null)
			{
				msd += last.distance(result);
			}
			last = result;
		}
		return msd / (size() - 1);
	}
}

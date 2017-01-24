package gdsc.smlm.ij.settings;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Contain the settings for the clustering algorithm
 * 
 * @author Alex Herbert
 */
public class OPTICSSettings implements Cloneable
{
	/**
	 * Options for displaying the clustering image
	 */
	public enum ImageMode
	{
		//@formatter:off
		CLUSTER_ID {
			@Override
			public String getName() { return "Cluster Id"; };
			@Override
			public float getValue(float value, int clusterId) 
			{ 
				return clusterId; 
			}
			@Override
			public boolean isMapped()
			{
				return true;
			}
		},
		CLUSTER_DEPTH {
			@Override
			public String getName() { return "Cluster Depth"; };
			@Override
			public float getValue(float value, int clusterId) 
			{ 
				return clusterId; 
			}
			@Override
			public boolean isMapped()
			{
				return true;
			}
		},
		VALUE {
			@Override
			public String getName() { return "Value"; };
			@Override
			public boolean canBeWeighted() { return true; }
			@Override
			public float getValue(float value, int clusterId) { return value; }
		},
		COUNT {
			@Override
			public String getName() { return "Count"; };
			@Override
			public boolean canBeWeighted() { return true; }
			@Override
			public float getValue(float value, int clusterId) { return 1f; }
		},
		NONE {
			@Override
			public String getName() { return "None"; };
			@Override
			public float getValue(float value, int clusterId) { return 0; }
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		/**
		 * Return the value to draw.
		 *
		 * @param value
		 *            The value of the cluster point
		 * @param clusterId
		 *            The cluster Id of the cluster point
		 * @return The value
		 */
		abstract public float getValue(float value, int clusterId);

		/**
		 * Can be weighted.
		 *
		 * @return true, if successful
		 */
		public boolean canBeWeighted()
		{
			return false;
		}
		
		/**
		 * Checks if is a mapped value.
		 *
		 * @return true, if is mapped
		 */
		public boolean isMapped()
		{
			return false;
		}

		@Override
		public String toString()
		{
			return getName();
		}
	}

	/**
	 * Options for plotting the OPTICS algorithm
	 */
	public enum OPTICSMode
	{
		//@formatter:off
		FAST_OPTICS {
			@Override
			public String getName() { return "FastOPTICS"; };
		},
		OPTICS {
			@Override
			public String getName() { return "OPTICS"; };
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		@Override
		public String toString()
		{
			return getName();
		}
	}

	/**
	 * Options for plotting the OPTICS results
	 */
	public enum ClusteringMode
	{
		//@formatter:off
		XI {
			@Override
			public String getName() { return "Xi"; };
		},
		DBSCAN {
			@Override
			public String getName() { return "pseudo-DBSCAN"; };
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		@Override
		public String toString()
		{
			return getName();
		}
	}

	/**
	 * Options for plotting the OPTICS results
	 */
	public enum PlotMode
	{
		//@formatter:off
		ON {
			@Override
			public String getName() { return "On"; };
		},
		WITH_CLUSTERS {
			@Override
			public String getName() { return "With clusters"; };
			@Override
			public boolean isDrawClusters() { return true; }
		},
		HIGHLIGHTED {
			@Override
			public String getName() { return "Highlighted"; };
			@Override
			public boolean isHighlightProfile() { return true; }
		},
		HIGHLIGHTED_WITH_CLUSTERS {
			@Override
			public String getName() { return "Highlighted with clusters"; };
			@Override
			public boolean isHighlightProfile() { return true; }
			@Override
			public boolean isDrawClusters() { return true; }
		},
		COLOURED_WITH_CLUSTERS {
			@Override
			public String getName() { return "Coloured with clusters"; };
			@Override
			public boolean isHighlightProfile() { return true; }
			@Override
			public boolean isColourProfile() { return true; }
			@Override
			public boolean isDrawClusters() { return true; }
		},
		COLOURED_BY_DEPTH_WITH_CLUSTERS {
			@Override
			public String getName() { return "Coloured by depth with clusters"; };
			@Override
			public boolean isHighlightProfile() { return true; }
			@Override
			public boolean isColourProfile() { return true; }
			@Override
			public boolean isColourProfileByDepth() { return true; }
			@Override
			public boolean isDrawClusters() { return true; }
		},
		OFF {
			@Override
			public String getName() { return "Off"; };
		};
		//@formatter:on

		/**
		 * Gets the name.
		 *
		 * @return the name
		 */
		abstract public String getName();

		/**
		 * @return True if the profile should be highlighted for top-cluster regions
		 */
		public boolean isHighlightProfile()
		{
			return false;
		}

		/**
		 * @return True if the profile should be coloured using the cluster
		 */
		public boolean isColourProfile()
		{
			return false;
		}

		/**
		 * @return True if the profile should be coloured using the cluster depth
		 */
		public boolean isColourProfileByDepth()
		{
			return false;
		}

		/**
		 * @return If clusters should be drawn on the plot
		 */
		public boolean isDrawClusters()
		{
			return false;
		}

		/**
		 * @return True if the clusters are needed
		 */
		public boolean requiresClusters()
		{
			return isDrawClusters() || isHighlightProfile() || isColourProfile();
		}

		@Override
		public String toString()
		{
			return getName();
		}
	}

	// Affect creating the OPTICS manager

	/**
	 * The input results dataset to use
	 */
	public String inputOption = "";

	// Affect running OPTICS

	/**
	 * The OPTICS algorithm to use.
	 */
	private OPTICSMode opticsMode = OPTICSMode.FAST_OPTICS;

	/** The number of splits to compute (if below 1 it will be auto-computed using the size of the data) */
	public int numberOfSplitSets = 0;

	/**
	 * The generating distance, i.e. the distance to search for neighbours of a point. Set to zero to auto-calibrate
	 * using the expected density of uniformly spread random points.
	 */
	public double generatingDistance = 0;

	/**
	 * The minimum number of neighbours to define a core point.
	 * <p>
	 * Note that the minimum cardinality (i.e. count of the number of neighbours) in the paper discussing Generalised
	 * DBSCAN is recommended to be 2 x dimensions, so 4 for a 2D dataset.
	 */
	public int minPoints = 4;

	// OPTICS clustering

	/**
	 * The clustering mode to use on the OPTICS results.
	 */
	private ClusteringMode clusteringMode = ClusteringMode.XI;

	// Affect running OPTICS Xi

	/**
	 * The steepness parameter for the OPTICS hierarchical clustering algorithm using the reachability profile.
	 */
	public double xi = 0.03;
	/**
	 * Set to true to only show the top-level clusters, i.e. child clusters will be merged into their parents.
	 */
	public boolean topLevel = false;

	// Affect DBSCAN clustering

	/**
	 * The number of samples to take for the k-distance plot. This should be 1-10% of the data.
	 */
	public int samples = 100;
	/**
	 * The fraction of the data to sample for the k-distance plot. Recommended to be 1-10%.
	 */
	public double sampleFraction = 0.05;
	/**
	 * The fraction of noise in the k-distance plot. The clustering distance is set as the next distance after noise has
	 * been ignored.
	 */
	public double fractionNoise = 0.05;
	/**
	 * The clustering distance for DBSCAN.
	 */
	public double clusteringDistance = 0;
	/**
	 * Set to true to only include core point in clusters. Note: Non-core points can be assigned arbitrarily to clusters
	 * if they are on the border of two clusters due to the arbitrary processing order of input points.
	 */
	public boolean core = false;

	// Affect display of results

	/**
	 * The magnification scale of the output image
	 */
	public double imageScale = 2;
	/**
	 * The output image mode
	 */
	private ImageMode imageMode = ImageMode.CLUSTER_ID;
	/**
	 * Set to true to weight the image data over nearest neighbour pixels
	 */
	public boolean weighted = false;
	/**
	 * Set to true to equalise the image histogram (allowing viewing high dynamic range data)
	 */
	public boolean equalised = false;

	/**
	 * The plot mode for the reachability distance profile
	 */
	private PlotMode plotMode = PlotMode.COLOURED_WITH_CLUSTERS;

	/**
	 * Set to true to draw the convex hull of each cluster as an outline
	 */
	public boolean outline = true;

	/**
	 * Set to true to draw the spanning tree (connections between each point and its parent)
	 */
	public boolean spanningTree = true;

	/**
	 * Set to true to draw the overlay coloured using the depth of the OPTICS hierarchy (default is the cluster ID)
	 */
	public boolean overlayColorByDepth = true;

	public OPTICSMode getOPTICSMode()
	{
		return opticsMode;
	}

	public int getOPTICSModeOridinal()
	{
		if (opticsMode == null)
			return 0;
		return opticsMode.ordinal();
	}

	public void setOPTICSMode(OPTICSMode mode)
	{
		opticsMode = mode;
	}

	public void setOPTICSMode(int mode)
	{
		OPTICSMode[] values = OPTICSMode.values();
		if (mode < 0 || mode >= values.length)
			mode = 0;
		this.opticsMode = values[mode];
	}

	public ImageMode getImageMode()
	{
		return imageMode;
	}

	public int getImageModeOridinal()
	{
		if (imageMode == null)
			return 0;
		return imageMode.ordinal();
	}

	public void setImageMode(ImageMode mode)
	{
		imageMode = mode;
	}

	public void setImageMode(int mode)
	{
		ImageMode[] values = ImageMode.values();
		if (mode < 0 || mode >= values.length)
			mode = 0;
		this.imageMode = values[mode];
	}

	public ClusteringMode getClusteringMode()
	{
		return clusteringMode;
	}

	public int getClusteringModeOridinal()
	{
		if (clusteringMode == null)
			return 0;
		return clusteringMode.ordinal();
	}

	public void setClusteringMode(ClusteringMode mode)
	{
		clusteringMode = mode;
	}

	public void setClusteringMode(int mode)
	{
		ClusteringMode[] values = ClusteringMode.values();
		if (mode < 0 || mode >= values.length)
			mode = 0;
		this.clusteringMode = values[mode];
	}

	public PlotMode getPlotMode()
	{
		return plotMode;
	}

	public int getPlotModeOridinal()
	{
		if (plotMode == null)
			return 0;
		return plotMode.ordinal();
	}

	public void setPlotMode(PlotMode mode)
	{
		plotMode = mode;
	}

	public void setPlotMode(int mode)
	{
		PlotMode[] values = PlotMode.values();
		if (mode < 0 || mode >= values.length)
			mode = 0;
		this.plotMode = values[mode];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public OPTICSSettings clone()
	{
		try
		{
			return (OPTICSSettings) super.clone();
		}
		catch (CloneNotSupportedException ex)
		{
			return null;
		}
	}
}

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
			public float getValue(float value, int clusterId) { return clusterId; }
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
		 * Return the value to draw
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
		 * @return True if the profile should be coloured using the cluster colour for top-cluster regions
		 */
		public boolean isColourProfile()
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
	public String inputOption = "";

	// Affect running OPTICS
	public double generatingDistance = 0;
	public int minPoints = 5;

	// OPTICS clustering
	public ClusteringMode clusteringMode = ClusteringMode.XI;

	// Affect running OPTICS Xi
	public double xi = 0.03;
	public boolean topLevel = false;

	// Affect DBSCAN clustering
	public double clusteringDistance = 0;
	public boolean core = false;

	// Affect display of results
	public double imageScale = 2;
	private ImageMode imageMode = ImageMode.CLUSTER_ID;
	public boolean weighted = false;
	public boolean equalised = false;

	private PlotMode plotMode = PlotMode.COLOURED_WITH_CLUSTERS;

	public boolean outline = true;

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

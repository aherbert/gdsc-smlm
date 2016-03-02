package gdsc.common.analytics;

import java.util.ArrayList;
import java.util.List;

/*
 * <ul>
 * <li>Copyright (c) 2016 Alex Herbert
 * </ul>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Base class for parameters to allow storing custom dimensions and metrics.
	 * Note that custom dimensions have to be created for your site before they can be used in analytics reports.
	 * 
	 * @see https://support.google.com/analytics/answer/2709829
 * 
 * @see https://developers.google.com/analytics/devguides/collection/protocol/v1/parameters
 */
public class Parameters
{
	/**
	 * Stores a custom dimension
	 */
	public class CustomDimension
	{
		final int index;
		final String value;

		public CustomDimension(int index, String value)
		{
			this.index = index;
			this.value = value;
		}
	}
	/**
	 * Stores a custom metric
	 */
	public class CustomMetric
	{
		final int index;
		final int value;

		public CustomMetric(int index, int value)
		{
			this.index = index;
			this.value = value;
		}
	}

	private List<CustomDimension> customDimensions = null;
	private List<CustomMetric> customMetrics = null;

	/**
	 * Default constructor
	 */
	public Parameters()
	{
	}

	/**
	 * Add a custom dimension
	 * <p>
	 * Note that custom dimensions have to be created for your site before they can be used in analytics reports.
	 * 
	 * @see https://support.google.com/analytics/answer/2709829
	 * 
	 * @param index The dimension index (1-20 or 1-200 for premium accounts)
	 * @param value The dimension value (must not be null)
	 */
	public void addCustomDimension(int index, String value)
	{
		if (index < 1 || index > 200)
			return;
		if (value == null)
			return;
		if (customDimensions == null)
			customDimensions = new ArrayList<CustomDimension>(1);
		customDimensions.add(new CustomDimension(index, value));
	}

	/**
	 * @return The custom dimensions
	 */
	public List<CustomDimension> getCustomDimensions()
	{
		return customDimensions;
	}

	/**
	 * @return The number of customer dimensions
	 */
	public int getNumberOfCustomDimensions()
	{
		return (customDimensions == null) ? 0 : customDimensions.size();
	}
	
	/**
	 * @return True if there are custom dimensions
	 */
	public boolean hasCustomeDimensions()
	{
		return customDimensions != null;
	}

	/**
	 * Add a custom metric
	 * <p>
	 * Note that custom metrics have to be created for your site before they can be used in analytics reports.
	 * 
	 * @see https://support.google.com/analytics/answer/2709829
	 * 
	 * @param index The dimension index (1-20 or 1-200 for premium accounts)
	 * @param value The metric value
	 */
	public void addCustomMetric(int index, int value)
	{
		if (index < 1 || index > 200)
			return;
		if (customMetrics == null)
			customMetrics = new ArrayList<CustomMetric>(1);
		customMetrics.add(new CustomMetric(index, value));
	}

	/**
	 * @return The custom metrics
	 */
	public List<CustomMetric> getCustomMetrics()
	{
		return customMetrics;
	}

	/**
	 * @return The number of customer metrics
	 */
	public int getNumberOfCustomMetrics()
	{
		return (customMetrics == null) ? 0 : customMetrics.size();
	}
	
	/**
	 * @return True if there are custom metrics
	 */
	public boolean hasCustomeMetrics()
	{
		return customMetrics != null;
	}
}
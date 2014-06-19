package gdsc.smlm.ij.utils;

import gdsc.smlm.fitting.function.OptimiserFunction;
import ij.IJ;

/**
 * Allow progress tracking of the Apache Commons Math 3 Optimiser in ImageJ
 */
public abstract class LoggingOptimiserFunction extends OptimiserFunction
{
	private boolean logging = false;
	private int evalCount = 0;
	protected String name = "Optimiser";

	public LoggingOptimiserFunction(String name)
	{
		this.name = name;
	}

	/**
	 * Log the count of evaluations to the ImageJ status bar
	 * 
	 * @param b
	 */
	public void setLogging(boolean b)
	{
		logging = b;
		if (b)
			evalCount = 0;
	}

	public void increment()
	{
		evalCount++;
		if (logging)
			IJ.showStatus(name + " Evaluation " + evalCount);
	}

	/**
	 * @return The function name
	 */
	public String getName()
	{
		return name;
	}
}
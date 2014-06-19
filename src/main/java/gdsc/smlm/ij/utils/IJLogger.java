package gdsc.smlm.ij.utils;

import ij.IJ;

/**
 * Log to the ImageJ log window
 */
public class IJLogger implements gdsc.smlm.fitting.logging.Logger
{
	/* (non-Javadoc)
	 * @see gdsc.smlm.fitting.logging.Logger#info(java.lang.String)
	 */
	public void info(String message)
	{
		IJ.log(message);
	}

	/* (non-Javadoc)
	 * @see gdsc.smlm.fitting.logging.Logger#info(java.lang.String, java.lang.Object[])
	 */
	public void info(String format, Object... args)
	{
		IJ.log(String.format(format, args));
	}
}
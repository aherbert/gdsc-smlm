/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 * 
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package gdsc.smlm;

public class TestSettings
{
	/**
	 * Set this to true to run all the speed tests. Tests will run and output the speed
	 * difference between two methods. Set to false to ignore tests.
	 */
	public static final boolean RUN_SPEED_TESTS = false;

	/**
	 * Set this to true to assert the speed tests, i.e. check if something is faster than another.
	 * Note that in some cases the speed difference between two methods is very close and so the test
	 * can fail a fraction of the time it is run. Leaving at false will ensure the package can be built.
	 */
	public static final boolean ASSERT_SPEED_TESTS = false;

	/**
	 * Set this to true to run fitting tests. Certain tests may fail the strict pass criteria due to the use of random
	 * seeds creating poor data. These test can be disabled.
	 */
	public static final boolean RUN_FITTING_TESTS = false;

	// This functionality would be better served with log4j (or another logging framework)
	// TODO - Add a logging framework for the tests. 

	// Low log levels are prioritised

	public enum LogLevel
	{
		SILENT
		{
			@Override
			int getLevel()
			{
				return 0;
			}
		},
		ERROR
		{
			@Override
			int getLevel()
			{
				return 1;
			}
		},
		INFO
		{
			@Override
			int getLevel()
			{
				return 2;
			}
		},
		DEBUG
		{
			@Override
			int getLevel()
			{
				return 3;
			}
		};

		abstract int getLevel();
	}

	/** The verbosity for logging. */
	private static LogLevel logLevel;
	private static int verbosity;
	static
	{
		setLogLevel(LogLevel.INFO);
	}

	public static void setLogLevel(LogLevel logLevel)
	{
		if (logLevel == null)
			logLevel = LogLevel.INFO;
		TestSettings.logLevel = logLevel;
		// Do not log if silent
		verbosity = (logLevel == LogLevel.SILENT) ? Integer.MIN_VALUE : logLevel.getLevel();
	}

	public static LogLevel getLogLevel()
	{
		return logLevel;
	}

	/**
	 * Provide messages for logging.
	 */
	public interface MessageProvider
	{
		String getMessage();
	}

	/**
	 * Provide messages for logging. To be used when providing the arguments for the message is computationally intense.
	 */
	public static abstract class Message
	{
		final String format;

		public Message(String format)
		{
			this.format = format;
		}

		/**
		 * Wrap the variable length arguments list.
		 *
		 * @param args
		 *            the arguments for the message
		 * @return the arguments for the message as an array
		 */
		public Object[] wrap(Object... args)
		{
			return args;
		}

		/**
		 * Gets the arguments for the message.
		 *
		 * @return the arguments
		 */
		public abstract Object[] getArgs();
	}

	public static void info(String format, Object... args)
	{
		log(LogLevel.INFO, format, args);
	}

	public static void info(MessageProvider msg)
	{
		log(LogLevel.INFO, msg);
	}

	public static void info(Message msg)
	{
		log(LogLevel.INFO, msg);
	}

	public static void log(LogLevel level, String format, Object... args)
	{
		if (allow(level))
			System.out.printf(format, args);
	}

	public static void log(LogLevel level, MessageProvider msg)
	{
		if (allow(level))
			System.out.print(msg.getMessage());
	}

	public static void log(LogLevel level, Message msg)
	{
		if (allow(level))
			System.out.printf(msg.format, msg.getArgs());
	}

	private static boolean allow(LogLevel level)
	{
		return level.getLevel() <= verbosity;
	}
}

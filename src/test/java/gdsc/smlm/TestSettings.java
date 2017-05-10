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
	private static LogLevel logLevel = LogLevel.INFO;

	private static int verbosity;

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

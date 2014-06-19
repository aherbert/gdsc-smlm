package gdsc.smlm;

public class TestSettings
{
	/**
	 * Set this to true to run all the speed tests. Tests will run and output the speed
	 * difference between two methods. Set to false to ignore tests.
	 */
	public static final boolean RUN_SPEED_TESTS = false;
	
	/**
	 * Set this to true assert the speed tests, i.e. check if something is faster than another.
	 * Note that in some cases the speed difference between two methods is very close and so the test 
	 * can fail a fraction of the time it is run. Leaving at false will ensure the package can be built. 
	 */
	public static final boolean ASSERT_SPEED_TESTS = false;
}

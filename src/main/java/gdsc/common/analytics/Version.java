package gdsc.common.analytics;

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
 * Store the version of the code
 */
public class Version
{
	/**
	 * The major version
	 */
	public static final int MAJOR = 1;
	/**
	 * The minor version
	 */
	public static final int MINOR = 0;
	/**
	 * The release version
	 */
	public static final int RELEASE = 0;

	/**
	 * The major version string
	 */
	public static final String VERSION_X;
	/**
	 * The major.minor version string
	 */
	public static final String VERSION_X_X;
	/**
	 * The major.minor.release version string
	 */
	public static final String VERSION_X_X_X;

	static
	{
		VERSION_X = getVersion(1);
		VERSION_X_X = getVersion(2);
		VERSION_X_X_X = getVersion(3);
	}

	/**
	 * Get the version as a string. The string is built as major.minor.release using the specified number of levels.
	 * 
	 * @param levels
	 *            The number of levels (1-3).
	 * @return The version
	 */
	public static String getVersion(int levels)
	{
		String version = Integer.toString(MAJOR);
		if (levels > 1)
			version += '.' + Integer.toString(MINOR);
		if (levels > 2)
			version += '.' + Integer.toString(RELEASE);
		return version;
	}
}
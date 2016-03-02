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
 * Represent a session
 */
public class Session
{
	/**
	 * Google sessions timeout after 30 minutes of inactivity
	 */
	private long timeout = 30 * 60;
	/**
	 * Timestamp of the session.
	 */
	private long now;

	/**
	 * Create a new session
	 */
	public Session()
	{
		now = 0;
	}

	/**
	 * Get the number of seconds since the epoch (midnight, January 1, 1970 UTC)
	 * 
	 * @return The timestamp in seconds
	 */
	public static long timestamp()
	{
		return System.currentTimeMillis() / 1000L;
	}

	/**
	 * Check if the session is new (i.e. has not been initialised, has timed out, or been reset).
	 * Calling this refreshes the current session to prevent timeout.
	 */
	public boolean isNew()
	{
		// Get the current session expire time
		final long expires = now + timeout;
		// Get the current time.
		now = timestamp();
		// Check if the session has expired
		return (now > expires);
	}

	/**
	 * Reset and start a new session
	 */
	public void reset()
	{
		now = 0;
	}

	/**
	 * @return the timeout
	 */
	public long getTimeout()
	{
		return timeout;
	}

	/**
	 * @param timeout the timeout to set
	 */
	public void setTimeout(long timeout)
	{
		this.timeout = timeout;
	}
}
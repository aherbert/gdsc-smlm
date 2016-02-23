package gdsc.analytics;

/*
 * <ul>
 * <li>Copyright (c) 2010 Daniel Murphy
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
 *
 * @see https://code.google.com/archive/p/jgoogleanalyticstracker/
 */

import java.security.SecureRandom;

/**
 * Represent the Google Analytics session data for a visitor.
 * 
 * @see http://www.cardinalpath.com/ga-basics-the-structure-of-cookie-values/
 * @author Alex Herbert
 */
public class SessionData
{
	private int visitorId;
	private long initial;
	private long previous;
	private long current;
	private int sessionNumber;
	private int timeout = 30 * 60;

	private SessionData(int visitorId, int sessionNumber)
	{
		this.visitorId = visitorId;
		this.initial = 0;
		this.previous = 0;
		this.current = 0;
		this.sessionNumber = sessionNumber;
	}

	/**
	 * Initializes a new session data, with new random visitor id
	 */
	public static SessionData newSessionData()
	{
		final int visitorId = (new SecureRandom().nextInt() & 0x7FFFFFFF);
		return new SessionData(visitorId, 1);
	}

	public int getVisitorId()
	{
		return visitorId;
	}

	public long getInitial()
	{
		return initial;
	}

	public long getPrevious()
	{
		return previous;
	}

	/**
	 * Get the current time.
	 * <p>
	 * This should be called before {@link #getInitial()} and {@link #getPrevious()} as they are updated.
	 * 
	 * @return The current time
	 */
	public long getCurrent()
	{
		final long now = System.currentTimeMillis() / 1000L;
		if (initial == 0)
		{
			// If this is the first check of the time then initialise the session
			initial = previous = now;
		}
		else
		{
			// Else update the previous time
			previous = current;
			// Check the timeout and start a new session if necessary
			if (now > current + timeout)
				newSession();
		}
		return current = now;
	}

	public int getSessionNumber()
	{
		return sessionNumber;
	}

	/**
	 * Increment the session number to start a new session
	 */
	public void newSession()
	{
		// Do not update the timestamps here. 
		// The next call to getCurrent() will do that anyway.
		this.sessionNumber++;
	}

	/**
	 * Set the session timeout. After this amount of time the session number will increment.
	 * 
	 * @param timeout
	 *            The timeout in seconds
	 */
	public void setTimeout(int timeout)
	{
		this.timeout = timeout;
	}

	/**
	 * @return The timeout in seconds
	 */
	public int getTimeout()
	{
		return timeout;
	}
}
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

/**
 * Common client data. Allows caching of the client component of the Google Analytics URL.
 */
public class ClientData
{
	private final String trackingCode;
	private String encoding = "UTF-8";
	private String screenResolution = null;
	private String userLanguage = null;
	private String hostName = null;
	private SessionData sessionData;
	private String url = null;

	/**
	 * constructs with the tracking code and a new session data.
	 * 
	 * @param trackingCode
	 */
	public ClientData(String trackingCode)
	{
		this(trackingCode, SessionData.newSessionData());
	}

	/**
	 * constructs with the tracking code using the provided session data.
	 * 
	 * @param trackingCode
	 */
	public ClientData(String trackingCode, SessionData sessionData)
	{
		if (trackingCode == null)
		{
			throw new RuntimeException("Tracking code cannot be null");
		}
		this.trackingCode = trackingCode;
		this.sessionData = sessionData;
	}

	/**
	 * @return the trackingCode
	 */
	public String getTrackingCode()
	{
		return trackingCode;
	}

	/**
	 * @return the encoding
	 */
	public String getEncoding()
	{
		return encoding;
	}

	/**
	 * @return the screenResolution
	 */
	public String getScreenResolution()
	{
		return screenResolution;
	}

	/**
	 * @return the userLanguage
	 */
	public String getUserLanguage()
	{
		return userLanguage;
	}

	/**
	 * @return The hostname
	 */
	public String getHostName()
	{
		return hostName;
	}

	/**
	 * @return the session data, used to track unique sessions
	 */
	public SessionData getSessionData()
	{
		return sessionData;
	}

	/**
	 * Sets the character encoding of the client. like UTF-8
	 * 
	 * @param argEncoding
	 *            the encoding to set
	 */
	public void setEncoding(String argEncoding)
	{
		encoding = argEncoding;
	}

	/**
	 * Sets the screen resolution, like "1280x800".
	 * 
	 * @param argScreenResolution
	 *            the screenResolution to set
	 */
	public void setScreenResolution(String argScreenResolution)
	{
		screenResolution = argScreenResolution;
	}

	/**
	 * Sets the user language, like "EN-us"
	 * 
	 * @param argUserLanguage
	 *            the userLanguage to set
	 */
	public void setUserLanguage(String argUserLanguage)
	{
		userLanguage = argUserLanguage;
	}

	/**
	 * Set the hostname
	 * 
	 * @param hostName
	 *            the hostname
	 */
	public void setHostName(String hostName)
	{
		this.hostName = hostName;
	}

	/**
	 * @return The client component of the URL
	 */
	public String getUrl()
	{
		return url;
	}

	/**
	 * Set the client component of the Google Analytics URL. This can be used to cache part of the URL.
	 * 
	 * @param url
	 *            The client component of the URL
	 */
	public void setUrl(String url)
	{
		this.url = url;
	}
}
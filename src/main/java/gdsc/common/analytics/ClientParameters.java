package gdsc.common.analytics;

import java.util.regex.Pattern;

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
 * Common client data. Allows caching of the client component of the Google Analytics URL.
 * 
 * @see https://developers.google.com/analytics/devguides/collection/protocol/v1/parameters
 */
public class ClientParameters extends Parameters
{
	private final String trackingId;
	private final String clientId;
	private final String applicationName;
	private final Session session;

	private String screenResolution = null;
	private String userLanguage = null;
	private String hostName = null;
	private String userAgent = null;
	private String applicationId = null;
	private String applicationVersion = null;
	private boolean anonymized = false;

	private String url = null;

	/**
	 * Constructs with the tracking Id. If the client Id is null or empty then a new UUID will be created.
	 * 
	 * @param trackingId
	 *            Tracking Id (must not be null)
	 * @param clientId
	 *            Client Id (optional)
	 * @param applicationName
	 *            Application name (must not be null)
	 */
	public ClientParameters(String trackingId, String clientId, String applicationName)
	{
		if (trackingId == null || trackingId.length() == 0)
			throw new IllegalArgumentException("Tracking code cannot be null");
		if (!Pattern.matches("[A-Z]+-[0-9]+-[0-9]+", trackingId))
			System.out.printf("Tracking code appears invalid: %s\n", trackingId);
		if (clientId == null || clientId.length() == 0)
			clientId = java.util.UUID.randomUUID().toString();
		if (applicationName == null || applicationName.length() == 0)
			throw new IllegalArgumentException("Application name cannot be null");
		this.trackingId = trackingId;
		this.clientId = clientId;
		this.applicationName = applicationName;
		this.session = new Session();
	}

	/**
	 * @return the tracking Id
	 */
	public String getTrackingId()
	{
		return trackingId;
	}

	/**
	 * @return The client Id
	 */
	public String getClientId()
	{
		return clientId;
	}

	/**
	 * @return the application name
	 */
	public String getApplicationName()
	{
		return applicationName;
	}

	/**
	 * Check if this is a new session
	 */
	public boolean isNewSession()
	{
		return session.isNew();
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
	 * @return the user agent
	 */
	public String getUserAgent()
	{
		return userAgent;
	}

	/**
	 * @return the application Id
	 */
	public String getApplicationId()
	{
		return applicationId;
	}

	/**
	 * @return the application version
	 */
	public String getApplicationVersion()
	{
		return applicationVersion;
	}

	/**
	 * @return True if the IP address of the sender will be anonymized
	 */
	public boolean isAnonymized()
	{
		return anonymized;
	}

	/**
	 * @return The client component of the URL
	 */
	public String getUrl()
	{
		return url;
	}

	/**
	 * Sets the screen resolution, like "1280x800".
	 * 
	 * @param argScreenResolution
	 *            the screenResolution to set
	 */
	public void setScreenResolution(String argScreenResolution)
	{
		this.url = null;
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
		this.url = null;
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
		this.url = null;
		this.hostName = hostName;
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

	/**
	 * @param userAgent
	 *            the user agent to set
	 */
	public void setUserAgent(String userAgent)
	{
		this.url = null;
		this.userAgent = userAgent;
	}

	/**
	 * @param applicationId
	 *            the application Id to set
	 */
	public void setApplicationId(String applicationId)
	{
		this.url = null;
		this.applicationId = applicationId;
	}

	/**
	 * @param applicationVersion
	 *            the application version to set
	 */
	public void setApplicationVersion(String applicationVersion)
	{
		this.url = null;
		this.applicationVersion = applicationVersion;
	}

	/**
	 * Set the state of IP anonymization
	 * 
	 * @param anonymize
	 *            True if the IP address of the sender will be anonymized
	 */
	public void setAnonymized(boolean anonymized)
	{
		this.url = null;
		this.anonymized = anonymized;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.analytics.Parameters#addCustomDimension(int, java.lang.String)
	 */
	@Override
	public void addCustomDimension(int index, String value)
	{
		this.url = null;
		super.addCustomDimension(index, value);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.analytics.Parameters#addCustomMetric(int, int)
	 */
	@Override
	public void addCustomMetric(int index, int value)
	{
		this.url = null;
		super.addCustomMetric(index, value);
	}

	/**
	 * Reset the session (i.e. start a new session)
	 */
	public void resetSession()
	{
		this.session.reset();
	}
}
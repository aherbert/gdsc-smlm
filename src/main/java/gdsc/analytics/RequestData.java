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
 * Google Analytics request data
 * 
 * @author Alex Herbert
 */
public class RequestData
{
	/**
	 * The type of request
	 * 
	 * @author Alex Herbert
	 */
	public enum Type
	{
		EVENT("event"), TRANSACTION("transaction"), ITEM("item"), CUSTOM_VARIABLE("custom_variable"), PAGE("page");

		private String name;

		private Type(String name)
		{
			this.name = name;
		}

		@Override
		public String toString()
		{
			return name;
		}
	}

	private Type type = Type.PAGE;
	private String value = null;
	private String pageTitle = null;
	private String pageURL = null;

	/**
	 * @param type
	 *            The type of request (defaults to page)
	 */
	public void setType(Type type)
	{
		this.type = type;
	}

	/**
	 * @return The type of request
	 */
	public String getType()
	{
		return type.toString();
	}

	/**
	 * @param value
	 *            the value (used for events or custom variables)
	 */
	public void setValue(String value)
	{
		this.value = value;
	}

	/**
	 * @return the value (used for events or custom variables)
	 */
	public String getValue()
	{
		return value;
	}

	/**
	 * Sets the page title, which will be the Content Title in Google Analytics
	 * 
	 * @param argContentTitle
	 *            the contentTitle to set
	 */
	public void setPageTitle(String argContentTitle)
	{
		pageTitle = argContentTitle;
	}

	/**
	 * @return the page title
	 */
	public String getPageTitle()
	{
		return pageTitle;
	}

	/**
	 * The page url, which is required. Traditionally this is of the form "/content/page.html", but you can
	 * put anything here (like "/com/dmurph/test.java").
	 * 
	 * @param argPageURL
	 *            the pageURL to set
	 */
	public void setPageURL(String argPageURL)
	{
		pageURL = argPageURL;
	}

	/**
	 * @return the page URL
	 */
	public String getPageURL()
	{
		return pageURL;
	}
}

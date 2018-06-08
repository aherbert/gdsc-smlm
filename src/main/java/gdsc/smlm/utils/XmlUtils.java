package gdsc.smlm.utils;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

/**
 * Extend the XML Utilities using XStream
 */
public class XmlUtils extends gdsc.core.utils.XmlUtils
{
	private static XStream xs = null;

	/**
	 * Convert an object to XML
	 * 
	 * @param obj
	 * @return XML string representation
	 */
	public static String toXML(Object obj)
	{
		init();
		if (xs != null)
		{
			try
			{
				return xs.toXML(obj);
			}
			catch (XStreamException ex)
			{
				//ex.printStackTrace();
			}
		}
		return "";
	}

	/**
	 * Load an object from the XML string representation
	 * 
	 * @param xml
	 * @return the object
	 */
	public static Object fromXML(String xml)
	{
		init();
		if (xs != null)
		{
			try
			{
				return xs.fromXML(xml);
			}
			catch (XStreamException ex)
			{
				//ex.printStackTrace();
			}
		}
		return null;
	}

	private static void init()
	{
		if (xs == null)
		{
			xs = new XStream(new DomDriver());
			XStream.setupDefaultSecurity(xs); // to be removed after 1.5
			xs.allowTypesByWildcard(new String[] { "gdsc.smlm.**" });
		}
	}
}

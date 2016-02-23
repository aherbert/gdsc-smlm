package gdsc.smlm.ij.plugins;

import gdsc.smlm.ij.ImageJTracker;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.settings.FilterSettings;
import gdsc.smlm.ij.settings.GlobalSettings;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.ij.utils.Utils;
import gdsc.smlm.results.filter.AndFilter;
import gdsc.smlm.results.filter.Filter;
import gdsc.smlm.results.filter.OrFilter;
import gdsc.smlm.results.filter.PrecisionFilter;
import gdsc.smlm.results.filter.SNRFilter;
import gdsc.smlm.results.filter.WidthFilter;
import gdsc.smlm.utils.TextUtils;
import gdsc.smlm.utils.XmlUtils;
import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.awt.Checkbox;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactoryConfigurationError;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Creates an XML file of configured filters from a template.
 */
public class CreateFilters implements PlugIn, ItemListener
{
	private static final String TITLE = "Create Filters";

	private static boolean enumerateEarly = true;
	private GlobalSettings gs;
	private FilterSettings filterSettings;

	private Pattern pattern = Pattern.compile("(\\S+)=\"(\\S+):(\\S+):(\\S+)\"(\\S*)");

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		ImageJTracker.recordPlugin(TITLE, arg);
		
		if (!showDialog())
			return;

		// This method assumes valid XML elements have been input with no root node.
		// The output will be an expanded set of the same XML elements with specific 
		// attributes updated.

		// Each top level element is processed. All the attributes are scanned and if 
		// they contain a 'min:max:increment' token the attribute is enumerated to 
		// the output XML, replicating the element.

		// Add a dummy root element to allow the XML to be loaded as a document
		String xml = "<root>" + filterSettings.filterTemplate + "</root>";

		// Load the XML as a document
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

		IJ.showStatus("Creating filters");

		try
		{
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document dom = db.parse(new InputSource(new StringReader(xml)));
			Element docElement = dom.getDocumentElement();

			StringWriter sw = new StringWriter();
			sw.write("<linked-list>");

			// For each element
			int childCount = docElement.getChildNodes().getLength();
			int total = 0;
			for (int c = 0; c < childCount; c++)
			{
				Node node = docElement.getChildNodes().item(c);
				if (node.getNodeType() != Node.ELEMENT_NODE)
					continue;

				total += processElement(sw, node);
			}

			sw.write("</linked-list>");

			if (total > 0)
				saveFilters(sw, total);
			else
				IJ.error(TITLE, "No filters created");
		}
		catch (Exception e)
		{
			IJ.error(TITLE, "Unable to load the input XML:\n" + e.getMessage());
			IJ.showStatus("");
		}
	}

	private int processElement(StringWriter sw, Node node) throws TransformerFactoryConfigurationError,
			TransformerException
	{
		// Get entire element as a string
		String xmlString = XmlUtils.getString(node, false);

		ArrayList<StringBuilder> out = new ArrayList<StringBuilder>();

		// Process through the XML appending to the current output.
		String[] tokens = xmlString.split("\\s+");
		for (String token : tokens)
		{
			// Any attributes with enumerations should be expanded.
			Matcher match = pattern.matcher(token);
			if (match.find())
			{
				final String prefix = " " + match.group(1) + "=\"";
				// Use Big Decimal for the enumeration to preserve the precision of the input text
				// (i.e. using doubles for an enumeration can lose precision and fail to correctly enumerate)
				BigDecimal min, max, inc;
				try
				{
					min = new BigDecimal(match.group(2));
					max = new BigDecimal(match.group(3));
					inc = new BigDecimal(match.group(4));

					if (min.compareTo(max) > 0 || inc.compareTo(BigDecimal.ZERO) <= 0)
						throw new RuntimeException("Invalid 'min:max:increment' attribute: " + token);
				}
				catch (NumberFormatException e)
				{
					throw new RuntimeException("Invalid 'min:max:increment' attribute: " + token, e);
				}
				final String suffix = "\"" + match.group(5);

				// Enumerate the attribute
				ArrayList<String> attributeText = new ArrayList<String>();
				for (BigDecimal bd = min; bd.compareTo(max) <= 0; bd = bd.add(inc))
					attributeText.add(prefix + bd.toString() + suffix);

				// Add to the current output
				ArrayList<StringBuilder> out2 = new ArrayList<StringBuilder>(out.size() * attributeText.size());

				if (enumerateEarly)
				{
					// Enumerate earlier attributes first
					for (String text : attributeText)
					{
						for (StringBuilder sb : out)
						{
							out2.add(new StringBuilder(sb.toString()).append(text));
						}
					}
				}
				else
				{
					// Enumerate later attributes first
					for (StringBuilder sb : out)
					{
						final String current = sb.toString();
						for (String text : attributeText)
						{
							out2.add(new StringBuilder(current).append(text));
						}
					}
				}

				out = out2;
			}
			else
			{
				if (out.isEmpty())
					out.add(new StringBuilder(token));
				else
				{
					for (StringBuilder sb : out)
						sb.append(" ").append(token);
				}
			}
		}

		if (out.size() > 0)
		{
			sw.write("<FilterSet name=\"");
			sw.write(getName(out.get(0)));
			sw.write("\"><filters class=\"linked-list\">");
			for (StringBuilder sb : out)
				sw.write(sb.toString());
			sw.write("</filters></FilterSet>");
		}
		return out.size();
	}

	private String getName(StringBuilder sb)
	{
		Filter f = Filter.fromXML(sb.toString());
		if (f != null)
			return f.getType().replaceAll("&", "&amp;");
		return "";
	}

	private void saveFilters(StringWriter sw, int total)
	{
		// Save the output to file
		IJ.showStatus("Saving filters");
		String filename = Utils.getFilename("Filter_File", filterSettings.filterSetFilename);
		if (filename != null)
		{
			OutputStreamWriter out = null;
			try
			{
				filterSettings.filterSetFilename = filename;
				// Append .xml if no suffix
				if (filename.lastIndexOf('.') < 0)
					filterSettings.filterSetFilename += ".xml";

				FileOutputStream fos = new FileOutputStream(filterSettings.filterSetFilename);
				out = new OutputStreamWriter(fos, "UTF-8");
				out.write(XmlUtils.prettyPrintXml(sw.toString()));
				SettingsManager.saveSettings(gs);
				IJ.showStatus(total + " filters: " + filterSettings.filterSetFilename);
			}
			catch (Exception e)
			{
				IJ.log("Unable to save the filter sets to file: " + e.getMessage());
			}
			finally
			{
				if (out != null)
				{
					try
					{
						out.close();
					}
					catch (IOException e)
					{
						// Ignore
					}
				}
			}

		}
	}

	private boolean showDialog()
	{
		GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage("Create a set of filters for use in the Filter Analysis plugin.\nAttributes will be enumerated if they are of the form 'min:max:increment'");

		gs = SettingsManager.loadSettings();
		filterSettings = gs.getFilterSettings();

		gd.addTextAreas(filterSettings.filterTemplate, null, 20, 80);
		gd.addCheckbox("Enumerate_early attributes first", enumerateEarly);
		gd.addCheckbox("Show_demo_filters", false);

		if (!(java.awt.GraphicsEnvironment.isHeadless() || IJ.isMacro()))
		{
			Checkbox cb = (Checkbox) gd.getCheckboxes().get(1);
			cb.addItemListener(this);
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		filterSettings.filterTemplate = gd.getNextText();
		enumerateEarly = gd.getNextBoolean();
		boolean demoFilters = gd.getNextBoolean();

		if (demoFilters)
		{
			logDemoFilters();
			return false;
		}

		return SettingsManager.saveSettings(gs);
	}

	public void itemStateChanged(ItemEvent e)
	{
		// When the checkbox is clicked, output the list of available filters to the ImageJ log

		Checkbox cb = (Checkbox) e.getSource();
		if (cb.getState())
		{
			cb.setState(false);

			logDemoFilters();
		}
	}

	private void logDemoFilters()
	{
		comment(TITLE + " example template filters");
		IJ.log("");
		comment("Filters are described using XML");
		comment("Filter attibutes that take the form 'min:max:increment' will be enumerated");
		IJ.log("");

		comment("Single filters");
		IJ.log("");
		demo(new SNRFilter(10), new String[] { "10:20:1" });
		demo(new PrecisionFilter(30), new String[] { "30:50:2" });
		IJ.log("");

		comment("Combined filters");
		IJ.log("");
		demo(new AndFilter(new SNRFilter(10), new WidthFilter(2)), new String[] { "10:20:1", "1.5:2.5:0.2" });
		demo(new OrFilter(new PrecisionFilter(30), new AndFilter(new SNRFilter(10), new WidthFilter(2))), new String[] {
				"30:40:2", "10:20:1", "1.5:2.5:0.2" });
		IJ.log("");
	}

	private void demo(Filter filter, String... attributeSubstitutions)
	{
		// Create the filter XML
		String xml = filter.toXML();
		if (attributeSubstitutions != null)
		{
			// Process the XML substituting attributes in the order they occur using a SAX parser. 
			// Write the new XMl to a buffer.
			StringBuilder sb = new StringBuilder();
			try
			{
				SAXParserFactory factory = SAXParserFactory.newInstance();
				SAXParser saxParser = factory.newSAXParser();
				saxParser.parse(new InputSource(new StringReader(xml)), new AttributeSubstitutionHandler(sb,
						attributeSubstitutions));
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}

			IJ.log(XmlUtils.prettyPrintXml(sb.toString()));
		}
		else
		{
			IJ.log(xml);
		}
	}

	/*
	 * Inner class for the Callback Handlers. This class replaces the attributes in the XML with the given
	 * substitutions. Attributes are processed in order. Substitutions are ignored (skipped) if they are null.
	 */
	class AttributeSubstitutionHandler extends DefaultHandler
	{
		StringBuilder sb;
		String[] attributeSubstitutions;
		int substitutionCount = 0;

		public AttributeSubstitutionHandler(StringBuilder sb, String[] attributeSubstitutions)
		{
			this.sb = sb;
			this.attributeSubstitutions = attributeSubstitutions;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.xml.sax.helpers.DefaultHandler#startElement(java.lang.String, java.lang.String, java.lang.String,
		 * org.xml.sax.Attributes)
		 * 
		 * Only start elements have attributes so this is where the substitutions are made
		 */
		@Override
		public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException
		{
			sb.append("<").append(qName);
			for (int attribute = 0; attribute < attributes.getLength(); attribute++)
			{
				sb.append(" ");
				if (substitutionCount < attributeSubstitutions.length &&
						!attributes.getLocalName(attribute).equals("class"))
				{
					String nextSubstitution = attributeSubstitutions[substitutionCount++];
					if (nextSubstitution != null)
					{
						sb.append(attributes.getLocalName(attribute)).append("=\"").append(nextSubstitution)
								.append("\"");
						continue;
					}
				}

				sb.append(attributes.getLocalName(attribute)).append("=\"").append(attributes.getValue(attribute))
						.append("\"");
			}
			sb.append(">");
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see org.xml.sax.helpers.DefaultHandler#endElement(java.lang.String, java.lang.String, java.lang.String)
		 * 
		 * We must respect the end elements since combined filters require them. They can be later stripped using a
		 * pretty print XML method.
		 */
		@Override
		public void endElement(String uri, String localName, String qName) throws SAXException
		{
			sb.append("</").append(qName).append(">");
		}
	}

	private void comment(String text)
	{
		IJ.log(TextUtils.wrap("<!-- " + text + " -->", 80));
	}
}

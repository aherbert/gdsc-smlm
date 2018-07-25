/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.ij.plugins;

import java.awt.Checkbox;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.FileOutputStream;
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
import javax.xml.transform.TransformerFactoryConfigurationError;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.helpers.DefaultHandler;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.GUIFilterSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.filter.AndFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.Filter;
import uk.ac.sussex.gdsc.smlm.results.filter.OrFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.PrecisionFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.SNRFilter;
import uk.ac.sussex.gdsc.smlm.results.filter.WidthFilter;

/**
 * Creates an XML file of configured filters from a template.
 */
public class CreateFilters implements PlugIn, ItemListener
{
	private static final String TITLE = "Create Filters";

	private static boolean enumerateEarly = true;
	private GUIFilterSettings.Builder filterSettings;

	private final Pattern pattern = Pattern.compile("(\\S+)=\"(\\S+):(\\S+):(\\S+)\"(\\S*)");

	/*
	 * (non-Javadoc)
	 *
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	@Override
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (!showDialog())
			return;

		// This method assumes valid XML elements have been input with no root node.
		// The output will be an expanded set of the same XML elements with specific
		// attributes updated.

		// Each top level element is processed. All the attributes are scanned and if
		// they contain a 'min:max:increment' token the attribute is enumerated to
		// the output XML, replicating the element.

		// Add a dummy root element to allow the XML to be loaded as a document
		final String xml = "<root>" + filterSettings.getFilterTemplate() + "</root>";

		// Load the XML as a document
		final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

		IJ.showStatus("Creating filters");

		try
		{
			final DocumentBuilder db = dbf.newDocumentBuilder();
			final Document dom = db.parse(new InputSource(new StringReader(xml)));
			final Element docElement = dom.getDocumentElement();

			final StringWriter sw = new StringWriter();
			sw.write("<linked-list>");

			// For each element
			final int childCount = docElement.getChildNodes().getLength();
			int total = 0;
			for (int c = 0; c < childCount; c++)
			{
				final Node node = docElement.getChildNodes().item(c);
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
		catch (final Exception e)
		{
			IJ.error(TITLE, "Unable to load the input XML:\n" + e.getMessage());
			IJ.showStatus("");
		}
	}

	private int processElement(StringWriter sw, Node node)
			throws TransformerFactoryConfigurationError
	{
		// Get entire element as a string
		final String xmlString = uk.ac.sussex.gdsc.core.utils.XmlUtils.getString(node, false);

		ArrayList<StringBuilder> out = new ArrayList<>();

		// Process through the XML appending to the current output.
		final String[] tokens = xmlString.split("\\s+");
		for (final String token : tokens)
		{
			// Any attributes with enumerations should be expanded.
			final Matcher match = pattern.matcher(token);
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
				catch (final NumberFormatException e)
				{
					throw new RuntimeException("Invalid 'min:max:increment' attribute: " + token, e);
				}
				final String suffix = "\"" + match.group(5);

				// Enumerate the attribute
				final ArrayList<String> attributeText = new ArrayList<>();
				for (BigDecimal bd = min; bd.compareTo(max) <= 0; bd = bd.add(inc))
					attributeText.add(prefix + bd.toString() + suffix);

				// Add to the current output
				final ArrayList<StringBuilder> out2 = new ArrayList<>(out.size() * attributeText.size());

				if (enumerateEarly)
					// Enumerate earlier attributes first
					for (final String text : attributeText)
						for (final StringBuilder sb : out)
							out2.add(new StringBuilder(sb.toString()).append(text));
				else
					// Enumerate later attributes first
					for (final StringBuilder sb : out)
					{
						final String current = sb.toString();
						for (final String text : attributeText)
							out2.add(new StringBuilder(current).append(text));
					}

				out = out2;
			}
			else if (out.isEmpty())
				out.add(new StringBuilder(token));
			else
				for (final StringBuilder sb : out)
					sb.append(" ").append(token);
		}

		if (out.size() > 0)
		{
			sw.write("<FilterSet name=\"");
			sw.write(getName(out.get(0)));
			sw.write("\"><filters class=\"linked-list\">");
			for (final StringBuilder sb : out)
				sw.write(sb.toString());
			sw.write("</filters></FilterSet>");
		}
		return out.size();
	}

	private static String getName(StringBuilder sb)
	{
		final Filter f = Filter.fromXML(sb.toString());
		if (f != null)
			return f.getType().replaceAll("&", "&amp;");
		return "";
	}

	private void saveFilters(StringWriter sw, int total)
	{
		// Save the output to file
		IJ.showStatus("Saving filters");
		final String filename = Utils.getFilename("Filter_File", filterSettings.getFilterSetFilename());
		if (filename != null)
		{
			filterSettings.setFilterSetFilename(filename);
			// Append .xml if no suffix
			if (filename.lastIndexOf('.') < 0)
				filterSettings.setFilterSetFilename(filterSettings.getFilterSetFilename() + ".xml");

			try (OutputStreamWriter out = new OutputStreamWriter(
					new FileOutputStream(filterSettings.getFilterSetFilename()), "UTF-8"))
			{
				out.write(uk.ac.sussex.gdsc.core.utils.XmlUtils.prettyPrintXml(sw.toString()));
				SettingsManager.writeSettings(filterSettings.build());
				IJ.showStatus(total + " filters: " + filterSettings.getFilterSetFilename());
			}
			catch (final Exception e)
			{
				IJ.log("Unable to save the filter sets to file: " + e.getMessage());
			}
		}
	}

	private boolean showDialog()
	{
		final GenericDialog gd = new GenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addMessage(
				"Create a set of filters for use in the Filter Analysis plugin.\nAttributes will be enumerated if they are of the form 'min:max:increment'");

		filterSettings = SettingsManager.readGUIFilterSettings(0).toBuilder();

		gd.addTextAreas(filterSettings.getFilterTemplate(), null, 20, 80);
		gd.addCheckbox("Enumerate_early attributes first", enumerateEarly);
		gd.addCheckbox("Show_demo_filters", false);

		if (Utils.isShowGenericDialog())
		{
			final Checkbox cb = (Checkbox) gd.getCheckboxes().get(1);
			cb.addItemListener(this);
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		filterSettings.setFilterTemplate(gd.getNextText());
		enumerateEarly = gd.getNextBoolean();
		final boolean demoFilters = gd.getNextBoolean();

		if (demoFilters)
		{
			logDemoFilters();
			return false;
		}

		return SettingsManager.writeSettings(filterSettings.build());
	}

	@Override
	public void itemStateChanged(ItemEvent e)
	{
		// When the checkbox is clicked, output the list of available filters to the ImageJ log

		final Checkbox cb = (Checkbox) e.getSource();
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
		comment("Note: This example is a subset. All filters are described in the user manual");
		IJ.log("");

		comment("Single filters");
		IJ.log("");
		demo(new SNRFilter(10), new String[] { "10:20:1" });
		demo(new PrecisionFilter(30), new String[] { "30:50:2" });
		IJ.log("");

		comment("Combined filters");
		IJ.log("");
		demo(new AndFilter(new SNRFilter(10), new WidthFilter(2)), new String[] { "10:20:1", "1.5:2.5:0.2" });
		demo(new OrFilter(new PrecisionFilter(30), new AndFilter(new SNRFilter(10), new WidthFilter(2))),
				new String[] { "30:40:2", "10:20:1", "1.5:2.5:0.2" });
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
			final StringBuilder sb = new StringBuilder();
			try
			{
				final SAXParserFactory factory = SAXParserFactory.newInstance();
				final SAXParser saxParser = factory.newSAXParser();
				saxParser.parse(new InputSource(new StringReader(xml)),
						new AttributeSubstitutionHandler(sb, attributeSubstitutions));
			}
			catch (final Exception e)
			{
				e.printStackTrace();
			}

			//IJ.log(xml);
			//IJ.log(sb.toString());
			xml = sb.toString();
		}
		IJ.log(uk.ac.sussex.gdsc.core.utils.XmlUtils.prettyPrintXml(xml));
	}

	/*
	 * Inner class for the Callback Handlers. This class replaces the attributes in the XML with the given
	 * substitutions. Attributes are processed in order. Substitutions are ignored (skipped) if they are null.
	 */
	private class AttributeSubstitutionHandler extends DefaultHandler
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
		public void startElement(String uri, String localName, String qName, Attributes attributes)
		{
			sb.append("<").append(qName);
			for (int attribute = 0; attribute < attributes.getLength(); attribute++)
			{
				sb.append(" ");
				final String name = attributes.getQName(attribute);
				if (substitutionCount < attributeSubstitutions.length && !name.equals("class"))
				{
					final String nextSubstitution = attributeSubstitutions[substitutionCount++];
					if (nextSubstitution != null)
					{
						sb.append(name).append("=\"").append(nextSubstitution).append("\"");
						continue;
					}
				}

				sb.append(name).append("=\"").append(attributes.getValue(attribute)).append("\"");
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
		public void endElement(String uri, String localName, String qName)
		{
			sb.append("</").append(qName).append(">");
		}
	}

	private static void comment(String text)
	{
		IJ.log(TextUtils.wrap("<!-- " + text + " -->", 80));
	}
}

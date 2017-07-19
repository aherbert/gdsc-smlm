package gdsc.smlm.ij.plugins;

import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpression;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import com.thoughtworks.xstream.XStream;
import com.thoughtworks.xstream.XStreamException;
import com.thoughtworks.xstream.io.xml.DomDriver;

import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.CalibrationProtos.Calibration;
import gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import gdsc.smlm.engine.FitConfiguration;

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

import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.ij.settings.BatchRun;
import gdsc.smlm.ij.settings.BatchSettings;
import gdsc.smlm.ij.settings.ParameterSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.ExtendedGenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/**
 * Runs the Peak Fit plugin in a batch.
 * <p>
 * The batch specifies the set of images to process. For each image the batch can specify a set of values for each of
 * the fitting parameters. The Peak Fit plugin is then run for each combination of parameters and the results of each
 * run saved to file.
 * 
 * @deprecated This should be updated to use JSON and methods from the Google Proto Buffers library
 */
@Deprecated
public class BatchPeakFit implements PlugIn
{
	private static final String TITLE = "Batch Peak Fit";

	private static String configFilename = "";

	private XStream xs;
	private TextField configFilenameText;

	/**
	 * Default constructor
	 */
	public BatchPeakFit()
	{
		xs = createXStream();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		SMLMUsageTracker.recordPlugin(this.getClass(), arg);

		if (!showDialog())
			return;

		runBatch(configFilename);
	}

	/**
	 * Reads the batch configuration file. For each parameter variation, create a fitting configuration. Then run the
	 * fit engine on each input image using each configuration.
	 * 
	 * @param configurationFilename
	 */
	private void runBatch(String configurationFilename)
	{
		BatchSettings settings = loadSettings(configurationFilename);
		if (settings == null || settings.parameters.isEmpty())
		{
			IJ.log("No settings for the fitting engine");
			return;
		}

		if (!new File(settings.resultsDirectory).exists())
		{
			if (!new File(settings.resultsDirectory).mkdirs())
			{
				IJ.log("Unable to create the results directory: " + settings.resultsDirectory);
				return;
			}
		}

		Document doc = getDefaultSettingsXmlDocument();
		if (doc == null)
			return;

		// Create XML for each variation
		ArrayList<String> xmlSettings = new ArrayList<String>();
		setParameters(settings.parameters, 0, doc, xmlSettings);

		// Run all the variants on the input images
		for (String imageFilename : settings.images)
		{
			ImagePlus imp = IJ.openImage(imageFilename);
			if (imp == null)
			{
				IJ.log("Unable to load image: " + imageFilename);
				continue;
			}

			processImage(settings, imp, xmlSettings);
		}
	}

	private Document getDefaultSettingsXmlDocument()
	{
		// Create an XML document of the default settings
		Document doc = null;
		try
		{
			String configXml = xs.toXML(new FitEngineConfiguration());
			doc = loadDocument(configXml);
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		return doc;
	}

	/**
	 * Modify the XML document using the specified values for the given parameter. For each value
	 * call the method recursively for the next parameter. If there are no more parameters
	 * then add the XML document to the xmlSettings.
	 * 
	 * @param parameters
	 *            The list of parameters
	 * @param i
	 *            Parameter to process
	 * @param doc
	 *            The XML document containing all the current parameters
	 * @param xmlSettings
	 *            A list of XML parameter settings
	 */
	private void setParameters(ArrayList<ParameterSettings> parameters, int i, Document doc,
			ArrayList<String> xmlSettings)
	{
		// For all the fields within Calibration, PSF, FitEngineSettings
		// print out the field name and default value. The FitEngineSettings has to
		// extract the FitSettings separately.
		// The user can then set any field to a list of comma-separated values and the 
		// plugin will produce a functional fit configuration for each combination. 
		//FitSettings  s = FitSettings.getDefaultInstance();
		//Descriptor d = s.getDescriptorForType();
		//FieldDescriptor fd = d.findFieldByNumber(0);
		//fd.

		if (i < parameters.size())
		{
			ParameterSettings param = parameters.get(i);
			NodeList nodes = doc.getElementsByTagName(param.name);
			if (nodes.getLength() == 1)
			{
				// For each value, set the parameter and move to the next
				String[] values = param.value.split(",");
				for (String value : values)
				{
					Node node = nodes.item(0);
					node.setTextContent(value);
					setParameters(parameters, i + 1, doc, xmlSettings);
				}
			}
			else
			{
				// Just move to the next parameter
				setParameters(parameters, i + 1, doc, xmlSettings);
			}
		}
		else
		{
			// Add the final XML to the parameters to run
			TransformerFactory tf = TransformerFactory.newInstance();
			Transformer transformer;
			try
			{
				transformer = tf.newTransformer();
				transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
				StringWriter writer = new StringWriter();
				transformer.transform(new DOMSource(doc), new StreamResult(writer));
				xmlSettings.add(writer.getBuffer().toString());
			}
			catch (TransformerConfigurationException e)
			{
				e.printStackTrace();
			}
			catch (TransformerException e)
			{
				e.printStackTrace();
			}
		}
	}

	private BatchSettings loadSettings(String configurationFilename)
	{
		BatchSettings settings = null;
		FileInputStream fs = null;
		try
		{
			fs = new FileInputStream(configurationFilename);
			settings = (BatchSettings) xs.fromXML(fs);
		}
		catch (ClassCastException ex)
		{
			//ex.printStackTrace();
		}
		catch (FileNotFoundException ex)
		{
			ex.printStackTrace();
		}
		catch (XStreamException ex)
		{
			ex.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
		return settings;
	}

	private Document loadDocument(String xml)
	{
		DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
		DocumentBuilder docBuilder;
		try
		{
			docBuilder = docFactory.newDocumentBuilder();
			return docBuilder.parse(new ByteArrayInputStream(xml.getBytes()));
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		return null;
	}

	private void processImage(BatchSettings settings, ImagePlus imp, ArrayList<String> xmlSettings)
	{
		String imageFilename = imp.getOriginalFileInfo().directory + imp.getOriginalFileInfo().fileName;
		String basename = getBaseName(imp.getOriginalFileInfo().fileName);
		String format = settings.resultsDirectory + File.separatorChar + basename + ".%05d";

		String statusSuffix = String.format(" / %d: %s", xmlSettings.size(), imp.getOriginalFileInfo().fileName);

		Calibration calibration = settings.getCalibration();
		for (int i = 0; i < xmlSettings.size(); i++)
		{
			IJ.showStatus((i + 1) + statusSuffix);

			// Create the configuration
			FitEngineConfiguration fitConfig = null;
			try
			{
				fitConfig = (FitEngineConfiguration) xs.fromXML(xmlSettings.get(i));
			}
			catch (XStreamException e)
			{
				// Ignore
			}
			if (fitConfig == null)
				continue;

			// No need to skip settings that do not make sense as we will catch exceptions.
			// This relies on the fit engine to throw exceptions for invalid settings.

			// Update the configuration
			fitConfig.getFitConfiguration().setCalibration(calibration);
			// Ensure the state is restored after XStream object reconstruction
			fitConfig.getFitConfiguration().initialiseState();

			String prefix = String.format(format, i);

			// Save the settings
			String settingsFilename = saveRunSettings(prefix, imageFilename, fitConfig);

			// Run the fit engine
			if (settings.runPeakFit)
			{
				ResultsSettings resultsSettings = createResultsSettings(fitConfig, prefix);
				try
				{
					PeakFit peakFit = new PeakFit(fitConfig, resultsSettings);
					peakFit.setSilent(true);
					peakFit.run(imp, false);

					IJ.log(String.format("%s : %s : Size %d : Time = %s", imageFilename, settingsFilename,
							peakFit.getSize(), Utils.timeToString(peakFit.getTime())));
				}
				catch (Exception e)
				{
					// Ignore this as we assume this is from incorrect fit configuration
				}
			}
		}

		IJ.showStatus("");
	}

	private String getBaseName(String name)
	{
		int dot = name.lastIndexOf('.');
		return (dot == -1) ? name : name.substring(0, dot);
	}

	private String saveRunSettings(String prefix, String imageFilename, FitEngineConfiguration fitConfig)
	{
		BatchRun batchRun = new BatchRun(imageFilename, fitConfig);
		FileOutputStream fs = null;
		try
		{
			String settingsFilename = prefix + ".xml";
			fs = new FileOutputStream(settingsFilename);
			xs.toXML(batchRun, fs);
			return settingsFilename;
		}
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		catch (XStreamException e)
		{
			e.printStackTrace();
		}
		finally
		{
			if (fs != null)
			{
				try
				{
					fs.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
			}
		}
		return null;
	}

	private ResultsSettings createResultsSettings(FitEngineConfiguration fitConfig, String prefix)
	{
		ResultsSettings.Builder resultsSettings = ResultsSettings.newBuilder();
		resultsSettings.setShowDeviations(fitConfig.getFitConfiguration().isComputeDeviations());
		resultsSettings.getResultsFileSettingsBuilder().setResultsFilename(prefix + ".xls");
		resultsSettings.getResultsFileSettingsBuilder().setFileFormat(ResultsFileFormat.TEXT);
		return resultsSettings.build();
	}

	/**
	 * Ask for parameters
	 * 
	 * @return True if not cancelled
	 */
	private boolean showDialog()
	{
		ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
		gd.addHelp(About.HELP_URL);

		gd.addFilenameField("Config_filename", configFilename);
		if (Utils.isShowGenericDialog())
		{
			configFilenameText = (TextField) gd.getStringFields().get(0);

			// Add a button to create a configuration file
			gd.addAndGetButton("Create config file", new ActionListener()
			{
				public void actionPerformed(ActionEvent e)
				{
					Document doc = getDefaultSettingsXmlDocument();
					if (doc == null)
						return;

					try
					{
						// Look for nodes that are part of the fit configuration
						XPathFactory factory = XPathFactory.newInstance();
						XPath xpath = factory.newXPath();
						XPathExpression expr = xpath.compile("//gdsc.smlm.engine.FitEngineConfiguration//*");

						// For each node, add the name and value to the BatchParameters
						BatchSettings batchSettings = new BatchSettings();
						batchSettings.resultsDirectory = System.getProperty("java.io.tmpdir");
						batchSettings.images.add("/path/to/image.tif");

						Object result = expr.evaluate(doc, XPathConstants.NODESET);
						NodeList nodes = (NodeList) result;
						for (int i = 0; i < nodes.getLength(); i++)
						{
							Node node = nodes.item(i);
							if (node.getChildNodes().getLength() == 1) // Only nodes with a single text entry
							{
								batchSettings.parameters
										.add(new ParameterSettings(node.getNodeName(), node.getTextContent()));
							}
						}

						// Save the settings file
						String[] path = Utils.decodePath(configFilenameText.getText());
						OpenDialog chooser = new OpenDialog("Settings_file", path[0], path[1]);
						if (chooser.getFileName() != null)
						{
							String newFilename = chooser.getDirectory() + chooser.getFileName();
							newFilename = Utils.replaceExtension(newFilename, ".xml");
							FileOutputStream fs = null;
							try
							{
								fs = new FileOutputStream(newFilename);
								xs.toXML(batchSettings, fs);
							}
							finally
							{
								if (fs != null)
								{
									fs.close();
								}
							}

							// Update dialog filename
							configFilenameText.setText(newFilename);
						}
					}
					catch (Exception ex)
					{
						ex.printStackTrace();
					}
				}
			});
		}

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		configFilename = gd.getNextString().trim();

		return true;
	}

	private XStream createXStream()
	{
		XStream xs = new XStream(new DomDriver());
		xs.alias("gdsc.fitting.batchSettings", BatchSettings.class);
		xs.alias("parameter", ParameterSettings.class);
		xs.alias("gdsc.fitting.batchRun", BatchRun.class);

		xs.omitField(FitConfiguration.class, "flags");
		xs.omitField(FitConfiguration.class, "signalThreshold");
		//xs.omitField(FitConfiguration.class, "noise");
		xs.omitField(FitConfiguration.class, "enableValidation");

		return xs;
	}
}

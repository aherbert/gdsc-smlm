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

import ij.IJ;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import uk.ac.sussex.gdsc.core.ij.Utils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsFileFormat;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.settings.BatchRun;
import uk.ac.sussex.gdsc.smlm.ij.settings.BatchSettings;
import uk.ac.sussex.gdsc.smlm.ij.settings.ParameterSettings;

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

    private final XStream xs;
    private TextField configFilenameText;

    /**
     * Default constructor
     */
    public BatchPeakFit()
    {
        xs = createXStream();
    }

    /** {@inheritDoc} */
    @Override
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
     *            the configuration filename
     */
    private void runBatch(String configurationFilename)
    {
        final BatchSettings settings = loadSettings(configurationFilename);
        if (settings == null || settings.parameters.isEmpty())
        {
            IJ.log("No settings for the fitting engine");
            return;
        }

        if (!new File(settings.resultsDirectory).exists())
            if (!new File(settings.resultsDirectory).mkdirs())
            {
                IJ.log("Unable to create the results directory: " + settings.resultsDirectory);
                return;
            }

        final Document doc = getDefaultSettingsXmlDocument();
        if (doc == null)
            return;

        // Create XML for each variation
        final ArrayList<String> xmlSettings = new ArrayList<>();
        setParameters(settings.parameters, 0, doc, xmlSettings);

        // Run all the variants on the input images
        for (final String imageFilename : settings.images)
        {
            final ImagePlus imp = IJ.openImage(imageFilename);
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
            final String configXml = xs.toXML(new FitEngineConfiguration());
            doc = loadDocument(configXml);
        }
        catch (final XStreamException ex)
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
            final ParameterSettings param = parameters.get(i);
            final NodeList nodes = doc.getElementsByTagName(param.name);
            if (nodes.getLength() == 1)
            {
                // For each value, set the parameter and move to the next
                final String[] values = param.value.split(",");
                for (final String value : values)
                {
                    final Node node = nodes.item(0);
                    node.setTextContent(value);
                    setParameters(parameters, i + 1, doc, xmlSettings);
                }
            }
            else
                // Just move to the next parameter
                setParameters(parameters, i + 1, doc, xmlSettings);
        }
        else
        {
            // Add the final XML to the parameters to run
            final TransformerFactory tf = TransformerFactory.newInstance();
            Transformer transformer;
            try
            {
                transformer = tf.newTransformer();
                transformer.setOutputProperty(OutputKeys.OMIT_XML_DECLARATION, "yes");
                final StringWriter writer = new StringWriter();
                transformer.transform(new DOMSource(doc), new StreamResult(writer));
                xmlSettings.add(writer.getBuffer().toString());
            }
            catch (final TransformerConfigurationException e)
            {
                e.printStackTrace();
            }
            catch (final TransformerException e)
            {
                e.printStackTrace();
            }
        }
    }

    private BatchSettings loadSettings(String configurationFilename)
    {
        BatchSettings settings = null;
        try (FileInputStream fs = new FileInputStream(configurationFilename))
        {
            settings = (BatchSettings) xs.fromXML(fs);
        }
        catch (final ClassCastException ex)
        {
            //ex.printStackTrace();
        }
        catch (final FileNotFoundException ex)
        {
            ex.printStackTrace();
        }
        catch (final XStreamException ex)
        {
            ex.printStackTrace();
        }
        catch (final IOException e)
        {
            e.printStackTrace();
        }
        return settings;
    }

    private static Document loadDocument(String xml)
    {
        final DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder;
        try
        {
            docBuilder = docFactory.newDocumentBuilder();
            return docBuilder.parse(new ByteArrayInputStream(xml.getBytes()));
        }
        catch (final Exception e)
        {
            e.printStackTrace();
        }
        return null;
    }

    private void processImage(BatchSettings settings, ImagePlus imp, ArrayList<String> xmlSettings)
    {
        final String imageFilename = imp.getOriginalFileInfo().directory + imp.getOriginalFileInfo().fileName;
        final String basename = getBaseName(imp.getOriginalFileInfo().fileName);
        final String format = settings.resultsDirectory + File.separatorChar + basename + ".%05d";

        final String statusSuffix = String.format(" / %d: %s", xmlSettings.size(), imp.getOriginalFileInfo().fileName);

        final Calibration calibration = settings.getCalibration();
        for (int i = 0; i < xmlSettings.size(); i++)
        {
            IJ.showStatus((i + 1) + statusSuffix);

            // Create the configuration
            FitEngineConfiguration fitConfig = null;
            try
            {
                fitConfig = (FitEngineConfiguration) xs.fromXML(xmlSettings.get(i));
            }
            catch (final XStreamException e)
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

            final String prefix = String.format(format, i);

            // Save the settings
            final String settingsFilename = saveRunSettings(prefix, imageFilename, fitConfig);

            // Run the fit engine
            if (settings.runPeakFit)
            {
                final ResultsSettings resultsSettings = createResultsSettings(fitConfig, prefix);
                try
                {
                    final PeakFit peakFit = new PeakFit(fitConfig, resultsSettings);
                    peakFit.setSilent(true);
                    peakFit.run(imp, false);

                    IJ.log(String.format("%s : %s : Size %d : Time = %s", imageFilename, settingsFilename,
                            peakFit.getSize(), Utils.timeToString(peakFit.getTime())));
                }
                catch (final Exception e)
                {
                    // Ignore this as we assume this is from incorrect fit configuration
                }
            }
        }

        IJ.showStatus("");
    }

    private static String getBaseName(String name)
    {
        final int dot = name.lastIndexOf('.');
        return (dot == -1) ? name : name.substring(0, dot);
    }

    private String saveRunSettings(String prefix, String imageFilename, FitEngineConfiguration fitConfig)
    {
        final BatchRun batchRun = new BatchRun(imageFilename, fitConfig);
        final String settingsFilename = prefix + ".xml";
        try (FileOutputStream fs = new FileOutputStream(settingsFilename))
        {
            xs.toXML(batchRun, fs);
            return settingsFilename;
        }
        catch (final FileNotFoundException e)
        {
            e.printStackTrace();
        }
        catch (final XStreamException e)
        {
            e.printStackTrace();
        }
        catch (final IOException e)
        {
            e.printStackTrace();
        }
        return null;
    }

    private static ResultsSettings createResultsSettings(FitEngineConfiguration fitConfig, String prefix)
    {
        final ResultsSettings.Builder resultsSettings = ResultsSettings.newBuilder();
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
        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        gd.addHelp(About.HELP_URL);

        gd.addFilenameField("Config_filename", configFilename);
        if (Utils.isShowGenericDialog())
        {
            configFilenameText = (TextField) gd.getStringFields().get(0);

            // Add a button to create a configuration file
            gd.addAndGetButton("Create config file", new ActionListener()
            {
                @Override
                public void actionPerformed(ActionEvent e)
                {
                    final Document doc = getDefaultSettingsXmlDocument();
                    if (doc == null)
                        return;

                    try
                    {
                        // Look for nodes that are part of the fit configuration
                        final XPathFactory factory = XPathFactory.newInstance();
                        final XPath xpath = factory.newXPath();
                        // TODO: Check this still works after the package was refactored 
                        final XPathExpression expr = xpath
                                .compile("//uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration//*");

                        // For each node, add the name and value to the BatchParameters
                        final BatchSettings batchSettings = new BatchSettings();
                        batchSettings.resultsDirectory = System.getProperty("java.io.tmpdir");
                        batchSettings.images.add("/path/to/image.tif");

                        final Object result = expr.evaluate(doc, XPathConstants.NODESET);
                        final NodeList nodes = (NodeList) result;
                        for (int i = 0; i < nodes.getLength(); i++)
                        {
                            final Node node = nodes.item(i);
                            if (node.getChildNodes().getLength() == 1)
                                batchSettings.parameters
                                        .add(new ParameterSettings(node.getNodeName(), node.getTextContent()));
                        }

                        // Save the settings file
                        final String[] path = Utils.decodePath(configFilenameText.getText());
                        final OpenDialog chooser = new OpenDialog("Settings_file", path[0], path[1]);
                        if (chooser.getFileName() != null)
                        {
                            String newFilename = chooser.getDirectory() + chooser.getFileName();
                            newFilename = Utils.replaceExtension(newFilename, ".xml");
                            try (FileOutputStream fs = new FileOutputStream(newFilename))
                            {
                                xs.toXML(batchSettings, fs);
                            }

                            // Update dialog filename
                            configFilenameText.setText(newFilename);
                        }
                    }
                    catch (final Exception ex)
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

    private static XStream createXStream()
    {
        final XStream xs = new XStream(new DomDriver());
        XStream.setupDefaultSecurity(xs); // to be removed after 1.5
        // TODO: Check this still works after the package was refactored 
        xs.allowTypesByWildcard(new String[] { "uk.ac.sussex.gdsc.smlm.**" });
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

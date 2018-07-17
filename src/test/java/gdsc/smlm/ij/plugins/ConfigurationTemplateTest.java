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
package gdsc.smlm.ij.plugins;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import gdsc.smlm.data.config.FitProtos.DataFilterType;
import gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import gdsc.smlm.engine.FitConfiguration;
import gdsc.smlm.engine.FitEngineConfiguration;
import gdsc.smlm.ij.plugins.ConfigurationTemplate.TemplateResource;
import gdsc.smlm.ij.settings.SettingsManager;
import gdsc.smlm.results.filter.MultiFilter2;
import gdsc.test.TestSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;

@SuppressWarnings({ "javadoc" })
public class ConfigurationTemplateTest
{
	@Test
	public void canCreateResourceTemplate() throws IOException
	{
		// This is not really a test.
		// It creates some templates so that they can be put in the resources folder.
		File tmpFile = File.createTempFile("template", ".json.txt");
		tmpFile.deleteOnExit();

		FitEngineConfiguration config = new FitEngineConfiguration();
		FitConfiguration fitConfig = config.getFitConfiguration();
		fitConfig.setMaxIterations(150);
		fitConfig.setMinWidthFactor(0.2);
		fitConfig.setMinWidthFactor(5);
		fitConfig.setSmartFilter(true);
		fitConfig.setDirectFilter(new MultiFilter2(0, 22, 0.56, 2.55, 3.3, 0, 31, 0, 0));
		config.setSearch(0.607);
		config.setBorder(1);
		config.setFitting(3);
		config.setFailuresLimit(3);
		config.setIncludeNeighbours(true);
		config.setNeighbourHeightThreshold(0.3);
		config.setResidualsThreshold(0.4);
		config.setDataFilterType(DataFilterType.SINGLE);
		config.setDataFilter(DataFilterMethod.GAUSSIAN, 0.425, false, 0);

		//System.out.println(tmpFile.getPath());
		TemplateSettings.Builder builder = TemplateSettings.newBuilder();
		builder.setFitEngineSettings(config.getFitEngineSettings());
		Assert.assertTrue(SettingsManager.toJSON(builder.build(), tmpFile, SettingsManager.FLAG_JSON_WHITESPACE));
	}

	@Test
	public void canLoadResourceTemplates()
	{
		ConfigurationTemplate.clearTemplates();

		String[] before = ConfigurationTemplate.getTemplateNames(false);
		TemplateResource[] templates = ConfigurationTemplate.listTemplateResources();
		ConfigurationTemplate.loadTemplateResources(templates);
		String[] after = ConfigurationTemplate.getTemplateNames(false);

		checkLoaded("canLoadResourceTemplates", templates, before, after);
	}

	private static void checkLoaded(String test, TemplateResource[] templates, String[] before, String[] after)
	{
		// Subtract the before from the after
		HashSet<String> set = new HashSet<>(Arrays.asList(after));
		set.removeAll(Arrays.asList(before));

		Assert.assertEquals("Loaded incorrect number", templates.length, set.size());

		// Check all have been loaded
		for (TemplateResource template : templates)
		{
			if (set.contains(template.name))
			{
				TestSettings.info(test + " loaded: " + template);
				continue;
			}
			Assert.assertTrue(test + " could not load: " + template, false);
		}
	}

	@Test
	public void canListTemplateNameInOrder()
	{
		ConfigurationTemplate.clearTemplates();

		String[] names = new String[10];
		for (int i = 0; i < names.length; i++)
		{
			String name = "Test" + i;
			names[i] = name;
			ConfigurationTemplate.saveTemplate(name, TemplateSettings.getDefaultInstance(), null);
		}

		Assert.assertArrayEquals(names, ConfigurationTemplate.getTemplateNames());
	}

	@Test
	public void canLoadTemplateImageFromFile() throws IOException
	{
		ConfigurationTemplate.clearTemplates();

		Assert.assertEquals(0, ConfigurationTemplate.getTemplateNamesWithImage().length);

		// Create a dummy image
		int size = 20;
		float[] pixels = new float[size * size];
		RandomGenerator r = TestSettings.getRandomGenerator();
		for (int i = pixels.length; i-- > 0;)
			pixels[i] = r.nextFloat();
		ImagePlus imp = new ImagePlus("test", new FloatProcessor(size, size, pixels));
		File tmpFile = File.createTempFile("tmp", ".tif");
		tmpFile.deleteOnExit();
		IJ.save(imp, tmpFile.getPath());

		String name = "canLoadTemplateImageFromFile";
		File file = new File(Utils.replaceExtension(tmpFile.getPath(), ".xml"));
		ConfigurationTemplate.saveTemplate(name, TemplateSettings.getDefaultInstance(), file);

		Assert.assertEquals(1, ConfigurationTemplate.getTemplateNamesWithImage().length);

		ImagePlus imp2 = ConfigurationTemplate.getTemplateImage(name);

		Assert.assertNotNull(imp2);
		float[] data = (float[]) imp2.getProcessor().toFloat(0, null).getPixels();

		Assert.assertArrayEquals(pixels, data, 0);
	}

	@Test
	public void canLoadTemplateImageFromResources()
	{
		// This test requires that the system resources does have at least one template with an image

		ConfigurationTemplate.clearTemplates();
		ConfigurationTemplate.loadTemplateResources(ConfigurationTemplate.listTemplateResources());

		String[] names = ConfigurationTemplate.getTemplateNamesWithImage();
		if (names.length == 0)
			return;

		ImagePlus imp2 = ConfigurationTemplate.getTemplateImage(names[0]);

		Assert.assertNotNull(imp2);
		float[] data = (float[]) imp2.getProcessor().toFloat(0, null).getPixels();

		// Check data is not zero
		for (int i = 0; i < data.length; i++)
		{
			if (data[i] != 0)
			{
				// Check another pixel is different
				for (int j = i + 1; j < data.length; j++)
				{
					if (data[i] != data[j])
						return;
				}
			}
		}

		Assert.fail("Data is invalid");
	}
}

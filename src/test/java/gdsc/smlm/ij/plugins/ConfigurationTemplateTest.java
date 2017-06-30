package gdsc.smlm.ij.plugins;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.junit.Assert;
import org.junit.Test;

import gdsc.core.ij.Utils;
import gdsc.smlm.ij.plugins.ConfigurationTemplate.TemplateResource;
import gdsc.smlm.ij.settings.GlobalSettings;
import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;

@SuppressWarnings("deprecation")
public class ConfigurationTemplateTest
{
	@Test
	public void canLoadResourceTemplates()
	{
		canLoadResourceTemplates(true, true);
	}

	@Test
	public void canLoadAllMandatoryResourceTemplates()
	{
		canLoadResourceTemplates(true, false);
	}

	@Test
	public void canLoadAllOptionalResourceTemplates()
	{
		canLoadResourceTemplates(false, true);
	}

	private void canLoadResourceTemplates(boolean loadMandatory, boolean loadOptional)
	{
		ConfigurationTemplate.clearTemplates();

		// Load others
		if (loadMandatory ^ loadOptional)
		{
			TemplateResource[] templates = ConfigurationTemplate.listTemplates(!loadMandatory, !loadOptional);
			ConfigurationTemplate.loadTemplates(templates);
		}

		String[] before = ConfigurationTemplate.getTemplateNames(false);
		TemplateResource[] templates = ConfigurationTemplate.listTemplates(loadMandatory, loadOptional);
		ConfigurationTemplate.loadTemplates(templates);
		String[] after = ConfigurationTemplate.getTemplateNames(false);

		checkLoaded(String.format("Mandatory=%b, Optional=%b", loadMandatory, loadOptional), templates, before, after);
	}

	private void checkLoaded(String test, TemplateResource[] templates, String[] before, String[] after)
	{
		// Subtract the before from the after
		HashSet<String> set = new HashSet<String>(Arrays.asList(after));
		set.removeAll(Arrays.asList(before));

		Assert.assertEquals("Loaded incorrect number", set.size(), templates.length);

		// Check all have been loaded
		for (TemplateResource template : templates)
		{
			if (set.contains(template.name))
			{
				System.out.println(test + " loaded: " + template);
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
			ConfigurationTemplate.saveTemplate(name, new GlobalSettings(), null);
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
		RandomGenerator r = new Well19937c();
		for (int i = pixels.length; i-- > 0;)
			pixels[i] = r.nextFloat();
		ImagePlus imp = new ImagePlus("test", new FloatProcessor(size, size, pixels));
		File tmpFile = File.createTempFile("tmp", ".tif");
		tmpFile.deleteOnExit();
		IJ.save(imp, tmpFile.getPath());

		String name = "canLoadTemplateImageFromFile";
		File file = new File(Utils.replaceExtension(tmpFile.getPath(), ".xml"));
		ConfigurationTemplate.saveTemplate(name, new GlobalSettings(), file);

		Assert.assertEquals(1, ConfigurationTemplate.getTemplateNamesWithImage().length);

		ImagePlus imp2 = ConfigurationTemplate.getTemplateImage(name);

		Assert.assertNotNull(imp2);
		float[] data = (float[]) imp2.getProcessor().toFloat(0, null).getPixels();

		Assert.assertArrayEquals(pixels, data, 0);
	}

	@Test
	public void canLoadTemplateImageFromResources() throws IOException
	{
		// This test requires that the system resources does have at least one template with an image

		ConfigurationTemplate.clearTemplates();
		ConfigurationTemplate.loadTemplates(ConfigurationTemplate.listTemplates(true, true));

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

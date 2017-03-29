package gdsc.smlm.ij.plugins;

import java.util.Arrays;
import java.util.HashSet;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.ij.plugins.ConfigurationTemplate.TemplateResource;

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

		Assert.assertEquals(set.size(), templates.length);
		
		// Check all have been loaded
		for (TemplateResource template : templates)
		{
			if (set.contains(template.name))
			{
				System.out.println(test + " loaded: " + template.name);
				continue;
			}
			Assert.assertTrue(test + " could not load: " + template.name, false);
		}
	}
}

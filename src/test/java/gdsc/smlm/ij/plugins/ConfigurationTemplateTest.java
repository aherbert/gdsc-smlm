package gdsc.smlm.ij.plugins;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.ij.plugins.ConfigurationTemplate.TemplateResource;

public class ConfigurationTemplateTest
{
	@Test
	public void canLoadResourceTemplates()
	{
		TemplateResource[] templates = ConfigurationTemplate.listTemplates(true, true);
		ConfigurationTemplate.loadTemplates(templates);
		String[] after = ConfigurationTemplate.getTemplateNames(false);

		// Check all have been loaded
		outer: for (TemplateResource template : templates)
		{
			for (String name : after)
			{
				if (template.name.equals(name))
				{
					System.out.println("Loaded: " + template.name);
					continue outer;
				}
			}
			Assert.assertTrue("Could not load: " + template.name, false);
		}
	}
	
	// This test is only valid when executed alone since the loading is into static memory 
	// and the other test may already have loaded the template. 
	//@Test
	public void canLoadAllOptionalResourceTemplates()
	{
		String[] before = ConfigurationTemplate.getTemplateNames(false);
		TemplateResource[] templates = ConfigurationTemplate.listTemplates(false, true);
		ConfigurationTemplate.loadTemplates(templates);
		String[] after = ConfigurationTemplate.getTemplateNames(false);
		Assert.assertEquals(after.length, before.length + templates.length);
	}
}

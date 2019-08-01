/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import uk.ac.sussex.gdsc.core.utils.FileUtils;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.DataFilterType;
import uk.ac.sussex.gdsc.smlm.data.config.TemplateProtos.TemplateSettings;
import uk.ac.sussex.gdsc.smlm.engine.FitConfiguration;
import uk.ac.sussex.gdsc.smlm.engine.FitEngineConfiguration;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ConfigurationTemplate.TemplateResource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.filter.MultiFilter2;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.rng.RngUtils;
import uk.ac.sussex.gdsc.test.utils.TestLogUtils;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;

import org.apache.commons.rng.UniformRandomProvider;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class ConfigurationTemplateTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(ConfigurationTemplateTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  public void canCreateResourceTemplate() throws IOException {
    // This is not really a test.
    // It creates some templates so that they can be put in the resources folder.
    final File tmpFile = File.createTempFile("template", ".json.txt");
    tmpFile.deleteOnExit();

    final FitEngineConfiguration config = new FitEngineConfiguration();
    final FitConfiguration fitConfig = config.getFitConfiguration();
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

    // System.out.println(tmpFile.getPath());
    final TemplateSettings.Builder builder = TemplateSettings.newBuilder();
    builder.setFitEngineSettings(config.getFitEngineSettings());
    Assertions.assertTrue(
        SettingsManager.toJson(builder.build(), tmpFile, SettingsManager.FLAG_JSON_WHITESPACE));
  }

  @Test
  public void canLoadResourceTemplates() {
    ConfigurationTemplate.clearTemplates();

    final String[] before = ConfigurationTemplate.getTemplateNames(false);
    final TemplateResource[] templates = ConfigurationTemplate.listTemplateResources();
    ConfigurationTemplate.loadTemplateResources(templates);
    final String[] after = ConfigurationTemplate.getTemplateNames(false);

    checkLoaded("canLoadResourceTemplates", templates, before, after);
  }

  private static void checkLoaded(String test, TemplateResource[] templates, String[] before,
      String[] after) {
    // Subtract the before from the after
    final HashSet<String> set = new HashSet<>(Arrays.asList(after));
    set.removeAll(Arrays.asList(before));

    Assertions.assertEquals(templates.length, set.size(), "Loaded incorrect number");

    // Check all have been loaded
    for (final TemplateResource template : templates) {
      if (set.contains(template.name)) {
        logger.log(TestLogUtils.getRecord(Level.INFO, test + " loaded: " + template));
        continue;
      }
      Assertions.fail(test + " could not load: " + template);
    }
  }

  @Test
  public void canListTemplateNameInOrder() {
    ConfigurationTemplate.clearTemplates();

    final String[] names = new String[10];
    for (int i = 0; i < names.length; i++) {
      final String name = "Test" + i;
      names[i] = name;
      ConfigurationTemplate.saveTemplate(name, TemplateSettings.getDefaultInstance(), null);
    }

    Assertions.assertArrayEquals(names, ConfigurationTemplate.getTemplateNames());
  }

  @SeededTest
  public void canLoadTemplateImageFromFile(RandomSeed seed) throws IOException {
    ConfigurationTemplate.clearTemplates();

    Assertions.assertEquals(0, ConfigurationTemplate.getTemplateNamesWithImage().length);

    // Create a dummy image
    final int size = 20;
    final float[] pixels = new float[size * size];
    final UniformRandomProvider r = RngUtils.create(seed.getSeed());
    for (int i = pixels.length; i-- > 0;) {
      pixels[i] = r.nextFloat();
    }
    final ImagePlus imp = new ImagePlus("test", new FloatProcessor(size, size, pixels));
    final File tmpFile = File.createTempFile("tmp", ".tif");
    tmpFile.deleteOnExit();
    IJ.save(imp, tmpFile.getPath());

    final String name = "canLoadTemplateImageFromFile";
    final File file = new File(FileUtils.replaceExtension(tmpFile.getPath(), ".xml"));
    ConfigurationTemplate.saveTemplate(name, TemplateSettings.getDefaultInstance(), file);

    Assertions.assertEquals(1, ConfigurationTemplate.getTemplateNamesWithImage().length);

    final ImagePlus imp2 = ConfigurationTemplate.getTemplateImage(name);

    Assertions.assertNotNull(imp2);
    final float[] data = (float[]) imp2.getProcessor().toFloat(0, null).getPixels();

    Assertions.assertArrayEquals(pixels, data);
  }

  @Test
  public void canLoadTemplateImageFromResources() {
    // This test requires that the system resources does have at least one template with an image

    ConfigurationTemplate.clearTemplates();
    ConfigurationTemplate.loadTemplateResources(ConfigurationTemplate.listTemplateResources());

    final String[] names = ConfigurationTemplate.getTemplateNamesWithImage();
    if (names.length == 0) {
      return;
    }

    final ImagePlus imp2 = ConfigurationTemplate.getTemplateImage(names[0]);

    Assertions.assertNotNull(imp2);
    final float[] data = (float[]) imp2.getProcessor().toFloat(0, null).getPixels();

    // Check data is not zero
    for (int i = 0; i < data.length; i++) {
      if (data[i] != 0) {
        // Check another pixel is different
        for (int j = i + 1; j < data.length; j++) {
          if (data[i] != data[j]) {
            return;
          }
        }
      }
    }

    Assertions.fail("Data is invalid");
  }
}

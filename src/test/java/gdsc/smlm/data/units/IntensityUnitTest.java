package gdsc.smlm.data.units;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.units.IntensityUnit;
import gdsc.smlm.data.units.TypeConverter;

@SuppressWarnings("unchecked")
public class IntensityUnitTest
{
	@Test
	public void canConvert()
	{
		double offset = 120;
		double countPerPhoton = 45.5;
		for (int photon = 1; photon < 100; photon++)
		{
			//@formatter:off
    		check(offset, countPerPhoton, 
    			new ExpectedUnit<IntensityUnit>(IntensityUnit.COUNT, offset + photon * countPerPhoton), 
    			new ExpectedUnit<IntensityUnit>(IntensityUnit.PHOTON, photon)
    			);
    		check(0, countPerPhoton, 
        			new ExpectedUnit<IntensityUnit>(IntensityUnit.COUNT, photon * countPerPhoton), 
        			new ExpectedUnit<IntensityUnit>(IntensityUnit.PHOTON, photon)
        			);
    		check(countPerPhoton, 
        			new ExpectedUnit<IntensityUnit>(IntensityUnit.COUNT, photon * countPerPhoton), 
        			new ExpectedUnit<IntensityUnit>(IntensityUnit.PHOTON, photon)
        			);
    		//@formatter:on
		}
	}

	private void check(double offset, double countPerPhoton, ExpectedUnit<IntensityUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<IntensityUnit> c;
		for (int i = 0; i < n; i++)
		{
			IntensityUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				IntensityUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, offset, countPerPhoton);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}

	private void check(double countPerPhoton, ExpectedUnit<IntensityUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<IntensityUnit> c;
		for (int i = 0; i < n; i++)
		{
			IntensityUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				IntensityUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, countPerPhoton);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}

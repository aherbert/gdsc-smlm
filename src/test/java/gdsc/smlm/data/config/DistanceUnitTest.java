package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.SMLMSettings.DistanceUnit;
import gdsc.smlm.data.utils.TypeConverter;

@SuppressWarnings("unchecked")
public class DistanceUnitTest
{
	@Test
	public void canConvert()
	{
		double nmPerPixel = 104.5;
		for (int pixel = 1; pixel < 10; pixel++)
		{
    		//@formatter:off
    		check(nmPerPixel, 
    			new ExpectedUnit<DistanceUnit>(DistanceUnit.PIXEL, pixel), 
    			new ExpectedUnit<DistanceUnit>(DistanceUnit.UM, pixel * nmPerPixel / 1e3),
    			new ExpectedUnit<DistanceUnit>(DistanceUnit.NM, pixel * nmPerPixel)
    			);
    		//@formatter:on
		}
	}

	private void check(double nmPerPixel, ExpectedUnit<DistanceUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<DistanceUnit> c;
		for (int i = 0; i < n; i++)
		{
			DistanceUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				DistanceUnit u2 = expectedUnits[j].u;
				c =  UnitConverterFactory.createConverter(u1, u2, nmPerPixel);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}

package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.smlm.data.config.SMLMSettings.AngleUnit;
import gdsc.smlm.data.utils.TypeConverter;

@SuppressWarnings("unchecked")
public class AngleUnitTest
{
	@Test
	public void canConvert()
	{
		double degToRad = Math.PI / 180.0;
		for (int a = -360; a <= 360; a++)
		{
			//@formatter:off
    		check(degToRad, 
    			new ExpectedUnit<AngleUnit>(AngleUnit.DEGREE, a), 
    			new ExpectedUnit<AngleUnit>(AngleUnit.RADIAN, a * degToRad)
    			);
    		//@formatter:on
		}
	}

	private void check(double degToRad, ExpectedUnit<AngleUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<AngleUnit> c;
		for (int i = 0; i < n; i++)
		{
			AngleUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				AngleUnit u2 = expectedUnits[j].u;
				c =  UnitConverterFactory.createConverter(u1, u2);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}

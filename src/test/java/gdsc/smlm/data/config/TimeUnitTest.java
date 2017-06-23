package gdsc.smlm.data.config;

import org.junit.Assert;
import org.junit.Test;

import gdsc.core.data.utils.TypeConverter;
import gdsc.smlm.data.config.UnitConfig.TimeUnit;

@SuppressWarnings("unchecked")
public class TimeUnitTest
{
	@Test
	public void canConvert()
	{
		double msPerFrame = 35;
		for (int frame = 1; frame < 10; frame++)
		{
			//@formatter:off
    		check(msPerFrame, 
    			new ExpectedUnit<TimeUnit>(TimeUnit.FRAME, frame), 
    			new ExpectedUnit<TimeUnit>(TimeUnit.SECOND, frame * msPerFrame / 1e3),
    			new ExpectedUnit<TimeUnit>(TimeUnit.MILLISECOND, frame * msPerFrame)
    			);
    		//@formatter:on
		}
	}

	private void check(double msPerFrame, ExpectedUnit<TimeUnit>... expectedUnits)
	{
		int n = expectedUnits.length;
		TypeConverter<TimeUnit> c;
		for (int i = 0; i < n; i++)
		{
			TimeUnit u1 = expectedUnits[i].u;
			double v1 = expectedUnits[i].value;
			for (int j = 0; j < n; j++)
			{
				TimeUnit u2 = expectedUnits[j].u;
				c = UnitConverterFactory.createConverter(u1, u2, msPerFrame);
				double o = c.convert(v1);
				Assert.assertEquals(u1 + " to " + u2, expectedUnits[j].value, o, 1e-5);
			}
		}
	}
}

package gdsc.smlm.results;

import org.junit.Test;

@SuppressWarnings("deprecation")
public class CalibrationTest
{
	double[] test_a = { 100, 130, 160 };
	double[] test_s = { 80, 100, 140 };
	double[] test_N = { 1, 10, 30, 100, 1000 };
	double[] test_b2 = { 0, 1, 2, 4, 8 };
	int minPoints = 3, maxPoints = 20;

	@Test
	public void canGet()
	{
		Calibration c = new Calibration();
		c.getNmPerPixel();
		c.getGain();
		c.getExposureTime();
		c.getReadNoise();
		c.getBias();
		c.isEmCCD();
		c.getAmplification();		
	}

	@Test
	public void canGetWithNoException()
	{
		Calibration c = new Calibration(true);
		c.setNmPerPixel(1);
		c.setGain(2);
		c.setExposureTime(3);
		c.setReadNoise(4);
		c.setBias(5);
		c.setEmCCD(true);
		c.setAmplification(6);		
		
		c.getNmPerPixel();
		c.getGain();
		c.getExposureTime();
		c.getReadNoise();
		c.getBias();
		c.isEmCCD();
		c.getAmplification();		
	}
	
	@Test(expected=IllegalStateException.class)
	public void getNmPerPixelThrowsException()
	{
		new Calibration(true).getNmPerPixel();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getGainThrowsException()
	{
		new Calibration(true).getGain();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getExposureTimeThrowsException()
	{
		new Calibration(true).getExposureTime();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getReadNoiseThrowsException()
	{
		new Calibration(true).getReadNoise();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getBiasThrowsException()
	{
		new Calibration(true).getBias();
	}
	
	@Test(expected=IllegalStateException.class)
	public void isEmCCDThrowsException()
	{
		new Calibration(true).isEmCCD();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getAmplificationThrowsException()
	{
		new Calibration(true).getAmplification();
	}

	@Test(expected=IllegalStateException.class)
	public void getNmPerPixelThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setNmPerPixel(0);
		c.getNmPerPixel();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getGainThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setGain(0);
		c.getGain();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getExposureTimeThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setExposureTime(0);
		c.getExposureTime();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getReadNoiseThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setReadNoise(-1);
		c.getReadNoise();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getBiasThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setBias(-1);
		c.getBias();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getAmplificationThrowsExceptionWhenInvalid()
	{
		Calibration c = new Calibration(true);
		c.setAmplification(0);
		c.getAmplification();
	}
	

	@Test(expected=IllegalStateException.class)
	public void getNmPerPixelThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setNmPerPixel(1);
		c.clearHasNmPerPixel();
		c.getNmPerPixel();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getGainThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setGain(1);
		c.clearHasGain();
		c.getGain();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getExposureTimeThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setExposureTime(1);
		c.clearHasExposureTime();
		c.getExposureTime();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getReadNoiseThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setReadNoise(1);
		c.clearHasReadNoise();
		c.getReadNoise();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getBiasThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setBias(1);
		c.clearHasBias();
		c.getBias();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getEmCCDThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setEmCCD(true);
		c.clearHasEMCCD();
		c.isEmCCD();
	}
	
	@Test(expected=IllegalStateException.class)
	public void getAmplificationThrowsExceptionAfterClear()
	{
		Calibration c = new Calibration(true);
		c.setAmplification(1);
		c.clearHasAmplification();
		c.getAmplification();
	}
	
	@Test(expected=IllegalStateException.class)
	public void clearDoesNotResetFieldMissingFlag()
	{
		Calibration c = new Calibration(true);
		c.setAmplification(0);
		
		c.clear();
		c.getAmplification();
	}	
}

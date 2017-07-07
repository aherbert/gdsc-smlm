package gdsc.smlm.ij.settings;

/**
 * Contain the settings for the gdsc.fitting package
 * 
 * @deprecated The settings have been moved to the gdsc.smlm.data.config package.
 */
@Deprecated
public class GlobalSettings
{
	private OPTICSSettings opticsSettings = null;

	/**
	 * @return the opticsSettings
	 */
	public OPTICSSettings getOPTICSSettings()
	{
		if (opticsSettings == null)
			opticsSettings = new OPTICSSettings();
		return opticsSettings;
	}

	/**
	 * Check the setting is not currently null. If the setting is null then a call to the get() method will initialise a
	 * default object.
	 * 
	 * @return true if the setting is not null
	 */
	public boolean isOPTICSSettings()
	{
		return (opticsSettings != null);
	}
}

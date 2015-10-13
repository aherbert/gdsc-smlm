/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2013 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

import gdsc.smlm.ij.plugins.SMLMTools;
import ij.plugin.PlugIn;

/**
 * Default ImageJ plugin (no Java package) to run the gdsc.smlm.ij.plugins.SMLMTools plugin.
 * <p>
 * This class allows the project to be run in debug mode from an IDE (e.g. Eclipse). The Maven output directory will be
 * target/classes. Create a symbolic link to that directory from the project root and name it plugins. Optionally create
 * a link to the macros directory to allow the toolset to be loaded:
 * 
 * <pre>
 * ${smlm_home}/plugins -> ${smlm_home}/target/classes
 * ${smlm_home}/macros -> ${smlm_home}/target/classes/macros
 * </pre>
 * 
 * Set the project to run ij.ImageJ as the main class and use the smlm directory as the ImageJ path:
 * 
 * <pre>
 * ij.ImageJ -ijpath ${smlm_home}
 * </pre>
 * 
 * ImageJ will load this class from the plugins directory. This class can call all other plugins.
 */
public class SMLM_Plugins implements PlugIn
{
	/*
	 * (non-Javadoc)
	 * 
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */
	public void run(String arg)
	{
		// Create the SMLM Tools plugin.
		new SMLMTools();
	}
}

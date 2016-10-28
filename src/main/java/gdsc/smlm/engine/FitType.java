package gdsc.smlm.engine;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2015 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Define the type of fit that was performed
 */
public class FitType implements Cloneable
{
	public static final int MULTI = 1;
	public static final int MULTI_OK = 2;
	public static final int DOUBLET = 4;
	public static final int DOUBLET_OK = 8;
	public static final int MULTI_DOUBLET = 16;
	public static final int MULTI_DOUBLET_OK = 32;
	public static final int OK = 64;

	public static final int NO_OF_FLAGS = 7;

	private int flags;

	public FitType()
	{
	}

	public FitType(int flags)
	{
		this.flags = flags;
	}

	public int getFlags()
	{
		return flags;
	}

	public void setFlag(int flag, boolean enabled)
	{
		if (enabled)
			flags |= flag;
		else
			// TODO - check this
			flags &= flag;
	}

	public boolean getFlag(int flag)
	{
		return (flags & flag) == flag;
	}

	public void setMulti(boolean enabled)
	{
		setFlag(MULTI, enabled);
	}

	public void setMultiOK(boolean enabled)
	{
		setFlag(MULTI_OK, enabled);
	}

	public void setMultiDoublet(boolean enabled)
	{
		setFlag(MULTI_DOUBLET, enabled);
	}

	public void setMultiDoubletOK(boolean enabled)
	{
		setFlag(MULTI_DOUBLET_OK, enabled);
	}

	public void setDoublet(boolean enabled)
	{
		setFlag(DOUBLET, enabled);
	}

	public void setDoubletOK(boolean enabled)
	{
		setFlag(DOUBLET_OK, enabled);
	}

	public void setOK(boolean enabled)
	{
		setFlag(OK, enabled);
	}

	public boolean getMulti()
	{
		return getFlag(MULTI);
	}

	public boolean getMultiOK()
	{
		return getFlag(MULTI_OK);
	}

	public boolean getMultiDoublet()
	{
		return getFlag(MULTI_DOUBLET);
	}

	public boolean getMultiDoubletOK()
	{
		return getFlag(MULTI_DOUBLET_OK);
	}

	public boolean getDoublet()
	{
		return getFlag(DOUBLET);
	}

	public boolean getDoubletOK()
	{
		return getFlag(DOUBLET_OK);
	}

	public boolean getOK()
	{
		return getFlag(OK);
	}

	@Override
	public String toString()
	{
		if (flags == 0)
			return "None";
		StringBuilder sb = new StringBuilder();
		append(sb, getMulti(), "Multi");
		append(sb, getMultiOK(), "MultiOK");
		append(sb, getMultiDoublet(), "MultiDoublet");
		append(sb, getMultiDoubletOK(), "MultiDoubletOK");
		append(sb, getDoublet(), "Doublet");
		append(sb, getDoubletOK(), "DoubletOK");
		append(sb, getOK(), "OK");

		if (sb.length() == 0)
			return "None";
		return sb.toString();
	}

	private void append(StringBuilder sb, boolean enabled, String name)
	{
		if (enabled)
		{
			if (sb.length() != 0)
				sb.append("; ");
			sb.append(name);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public FitType clone()
	{
		try
		{
			return (FitType) super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			return null;
		}
	}
}
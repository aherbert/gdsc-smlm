package gdsc.smlm.results;

/*----------------------------------------------------------------------------- 
 * GDSC SMLM Software
 * 
 * Copyright (C) 2016 Alex Herbert
 * Genome Damage and Stability Centre
 * University of Sussex, UK
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *---------------------------------------------------------------------------*/

/**
 * Specifies a peak fitting result with assignable attributes. If the attributes are not set then the default values are
 * returned.
 */
public class AttributePeakResult extends PeakResult
{
	// Provide assignable attributes for all the attribute methods in the super class

	// State flags
	private static int FIELD_ID = 0x00000001;
	private static int FIELD_END_FRAME = 0x00000002;
	private static int FIELD_PRECISION = 0x00000004;
	private int fields = 0;

	private int id;
	private int endFrame;
	private float precision;

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#hasId()
	 */
	@Override
	public boolean hasId()
	{
		return ((fields & FIELD_ID) == FIELD_ID);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#hasEndFrame()
	 */
	@Override
	public boolean hasEndFrame()
	{
		return ((fields & FIELD_END_FRAME) == FIELD_END_FRAME);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#hasPrecision()
	 */
	@Override
	public boolean hasPrecision()
	{
		return ((fields & FIELD_PRECISION) == FIELD_PRECISION);
	}

	private void setHasId()
	{
		fields |= FIELD_ID;
	}

	private void setHasEndFrame()
	{
		fields |= FIELD_END_FRAME;
	}

	private void setHasPrecision()
	{
		fields |= FIELD_PRECISION;
	}

	public void clearHasId()
	{
		fields = fields & ~FIELD_ID;
	}

	public void clearHasEndFrame()
	{
		fields = fields & ~FIELD_END_FRAME;
	}

	public void clearHasPrecision()
	{
		fields = fields & ~FIELD_PRECISION;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#getId()
	 */
	@Override
	public int getId()
	{
		if (hasId())
			return id;
		return super.getId();
	}

	public void setId(int id)
	{
		// Allow ID to be anything, including zero
		//if (id != super.getId())
		//{
		setHasId();
		this.id = id;
		//}
		//else
		//	clearHasId();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#getEndFrame()
	 */
	@Override
	public int getEndFrame()
	{
		if (hasEndFrame())
			return endFrame;
		return super.getFrame();
	}

	public void setEndFrame(int endFrame)
	{
		// End frame must be after the start frame.
		if (endFrame > super.getFrame())
		{
			setHasEndFrame();
			this.endFrame = endFrame;
		}
		else
			clearHasEndFrame();
	}

	@Override
	public void setFrame(int frame)
	{
		// Set the new start frame
		super.setFrame(frame);
		// Validate the current end frame
		setEndFrame(endFrame);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.results.PeakResult#getPrecision()
	 */
	@Override
	public double getPrecision()
	{
		if (hasPrecision())
			return precision;
		return super.getPrecision();
	}

	public void setPrecision(double precision)
	{
		if (precision >= 0)
		{
			setHasPrecision();
			this.precision = (float) precision;
		}
		else
			clearHasPrecision();
	}

	/**
	 * Instantiates a new attribute peak result.
	 *
	 * @param startFrame
	 *            the start frame
	 * @param origX
	 *            the orig X
	 * @param origY
	 *            the orig Y
	 * @param origValue
	 *            the orig value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param params
	 *            the params
	 * @param paramsStdDev
	 *            the params std dev
	 */
	public AttributePeakResult(int startFrame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsStdDev)
	{
		super(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
	}

	/**
	 * Instantiates a new id peak result.
	 *
	 * @param x
	 *            the x
	 * @param y
	 *            the y
	 * @param sd
	 *            the sd
	 * @param signal
	 *            the signal
	 * @param id
	 *            the id
	 */
	public AttributePeakResult(float x, float y, float sd, float signal)
	{
		super(x, y, sd, signal);
	}
}

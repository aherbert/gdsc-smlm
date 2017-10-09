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
	private static final int FIELD_ID = 0x00000001;
	private static final int FIELD_END_FRAME = 0x00000002;
	private static final int FIELD_PRECISION = 0x00000004;
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
	 * @param frame
	 *            the frame
	 * @param origX
	 *            the original X position
	 * @param origY
	 *            the original Y position
	 * @param origValue
	 *            the original value
	 * @param error
	 *            the error
	 * @param noise
	 *            the noise
	 * @param params
	 *            the params (must not be null and must have at least {@value #STANDARD_PARAMETERS} parameters)
	 * @param paramsStdDev
	 *            the params standard deviations (if not null must match the length of the {@link #params} array)
	 * @throws IllegalArgumentException
	 *             the illegal argument exception if the parameters are invalid
	 */
	public AttributePeakResult(int startFrame, int origX, int origY, float origValue, double error, float noise,
			float[] params, float[] paramsStdDev) throws IllegalArgumentException
	{
		super(startFrame, origX, origY, origValue, error, noise, params, paramsStdDev);
	}

	/**
	 * Instantiates a new attribute peak result.
	 *
	 * @param frame
	 *            the frame
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 */
	public AttributePeakResult(int frame, float x, float y, float intensity)
	{
		super(frame, x, y, intensity);
	}

	/**
	 * Instantiates a new attribute peak result.
	 *
	 * @param x
	 *            the x position
	 * @param y
	 *            the y position
	 * @param intensity
	 *            the intensity
	 */
	public AttributePeakResult(float x, float y, float intensity)
	{
		super(x, y, intensity);
	}

	/**
	 * Instantiates a new attribute peak result copying all the attributes from the result.
	 * This is a deep copy of all the result data.
	 *
	 * @param peakResult
	 *            the peak result
	 */
	public AttributePeakResult(PeakResult peakResult)
	{
		super(peakResult.getFrame(), peakResult.origX, peakResult.origY, peakResult.origValue, peakResult.error,
				peakResult.noise, peakResult.params.clone(),
				(peakResult.paramStdDevs == null) ? null : peakResult.paramStdDevs.clone());
		if (peakResult.hasId())
			setId(peakResult.getId());
		if (peakResult.hasEndFrame())
			setEndFrame(peakResult.getEndFrame());
		if (peakResult.hasPrecision())
			setPrecision(peakResult.getPrecision());
	}
}

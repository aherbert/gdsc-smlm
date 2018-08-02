/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2018 Alex Herbert
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package uk.ac.sussex.gdsc.smlm.results;

import java.awt.Rectangle;
import java.util.Collection;

import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationProtos.Calibration;
import uk.ac.sussex.gdsc.smlm.data.config.PSFProtos.PSF;
import uk.ac.sussex.gdsc.smlm.results.predicates.PeakResultPredicate;

/**
 * Wraps peak results in memory and prevents modification of the results size.
 * <p>
 * Any method that modifies the size of the results set will throw a data exception.
 */
public class ImmutableMemoryPeakResults extends MemoryPeakResults
{
    private boolean built = false;

    /**
     * Instantiates a new immutable memory peak results with the original results store.
     *
     * @param results
     *            the results
     */
    public ImmutableMemoryPeakResults(MemoryPeakResults results)
    {
        this(results, false);
    }

    /**
     * Instantiates a new immutable memory peak results with an optional copy of the results store.
     *
     * @param results
     *            the results
     * @param copy
     *            Set to true to copy the original results store
     */
    public ImmutableMemoryPeakResults(MemoryPeakResults results, boolean copy)
    {
        super((PeakResultStoreList) ((copy) ? results.results.copy() : results.results));
        copySettings(results);
        built = true;
    }

    @Override
    PeakResult getfX(int index)
    {
        return new ImmutablePeakResult(super.getfX(index));
    }

    @Override
    public void setSource(ImageSource source)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setSource(source);
    }

    @Override
    public void setBounds(Rectangle bounds)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setBounds(bounds);
    }

    @Override
    public Rectangle getBounds()
    {
        // Prevent modification
        return new Rectangle(super.getBounds());
    }

    @Override
    public void setCalibration(Calibration calibration)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setCalibration(calibration);
    }

    @Override
    public void setPSF(PSF psf)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setPSF(psf);
    }

    @Override
    public void setConfiguration(String configuration)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setConfiguration(configuration);
    }

    @Override
    public void setName(String name)
    {
        if (built)
            throw new DataException("This results set is immutable");
        super.setName(name);
    }

    @Override
    public void add(PeakResult result)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void addAll(Collection<PeakResult> results)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void addAll(PeakResult[] results)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void addAll(PeakResultStore results)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void add(MemoryPeakResults results)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void removeNullResults()
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public boolean removeIf(PeakResultPredicate filter)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void begin()
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void add(int peak, int origX, int origY, float origValue, double chiSquared, float noise,
            float meanIntensity, float[] params, float[] paramsStdDev)
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public void end()
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public boolean isActive()
    {
        throw new DataException("This results set is immutable");
    }

    @Override
    public PeakResultView getPeakResultView()
    {
        return new DynamicPeakResultView(new ImmutablePeakResultStore(results));
    }

    @Override
    public PeakResultView getSnapshotPeakResultView()
    {
        return new CachedPeakResultView(new ImmutablePeakResultStore(results));
    }
}

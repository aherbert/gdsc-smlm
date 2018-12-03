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
package uk.ac.sussex.gdsc.smlm.engine;

/**
 * Define the type of fit that was performed.
 */
public class FitType implements Cloneable {
  /** Flag used for a single candidate fit with neighbours. */
  public static final int MULTI = 1;

  /** Flag used for a single candidate fit with neighbours was OK. */
  public static final int MULTI_OK = 2;

  /** Flag used for a double candidate fit. */
  public static final int DOUBLET = 4;

  /** Flag used for a double candidate fit was OK. */
  public static final int DOUBLET_OK = 8;

  /** Flag used for a double candidate fit with neighbours. */
  public static final int MULTI_DOUBLET = 16;

  /** Flag used for a double candidate fit with neighbours was OK. */
  public static final int MULTI_DOUBLET_OK = 32;

  /** Flag used when a fit was OK. */
  public static final int OK = 64;

  /** The number of flags. */
  public static final int NO_OF_FLAGS = 7;

  private int flags;

  /**
   * Instantiates a new fit type.
   */
  public FitType() {}

  /**
   * Instantiates a new fit type.
   *
   * @param flags the flags
   */
  public FitType(int flags) {
    this.flags = flags;
  }

  /**
   * Gets the flags.
   *
   * @return the flags
   */
  public int getFlags() {
    return flags;
  }

  /**
   * Sets the flag to enabled/disabled.
   *
   * @param flag the flag
   * @param enabled Set to true to enable
   */
  public void setFlag(int flag, boolean enabled) {
    if (enabled) {
      flags |= flag;
    } else {
      // TODO - check this
      flags &= flag;
    }
  }

  /**
   * Gets whether the flag is set to enabled/disabled.
   *
   * @param flag the flag
   * @return True if the flag is set to enabled
   */
  public boolean getFlag(int flag) {
    return (flags & flag) == flag;
  }

  /**
   * Sets the {@link #MULTI} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setMulti(boolean enabled) {
    setFlag(MULTI, enabled);
  }

  /**
   * Sets the {@link #MULTI_OK} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setMultiOK(boolean enabled) {
    setFlag(MULTI_OK, enabled);
  }

  /**
   * Sets the {@link #MULTI_DOUBLET} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setMultiDoublet(boolean enabled) {
    setFlag(MULTI_DOUBLET, enabled);
  }

  /**
   * Sets the {@link #MULTI_DOUBLET_OK} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setMultiDoubletOK(boolean enabled) {
    setFlag(MULTI_DOUBLET_OK, enabled);
  }

  /**
   * Sets the {@link #DOUBLET} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setDoublet(boolean enabled) {
    setFlag(DOUBLET, enabled);
  }

  /**
   * Sets the {@link #DOUBLET_OK} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setDoubletOK(boolean enabled) {
    setFlag(DOUBLET_OK, enabled);
  }

  /**
   * Sets the {@link #OK} flag to enabled/disabled.
   *
   * @param enabled Set to true to enable
   */
  public void setOK(boolean enabled) {
    setFlag(OK, enabled);
  }

  /**
   * Gets whether the {@link #MULTI} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getMulti() {
    return getFlag(MULTI);
  }

  /**
   * Gets whether the {@link #MULTI_OK} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getMultiOK() {
    return getFlag(MULTI_OK);
  }

  /**
   * Gets whether the {@link #MULTI_DOUBLET} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getMultiDoublet() {
    return getFlag(MULTI_DOUBLET);
  }

  /**
   * Gets whether the {@link #MULTI_DOUBLET_OK} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getMultiDoubletOK() {
    return getFlag(MULTI_DOUBLET_OK);
  }

  /**
   * Gets whether the {@link #DOUBLET} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getDoublet() {
    return getFlag(DOUBLET);
  }

  /**
   * Gets whether the {@link #DOUBLET_OK} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getDoubletOK() {
    return getFlag(DOUBLET_OK);
  }

  /**
   * Gets whether the {@link #OK} flag is set to enabled/disabled.
   *
   * @return Set to true if enabled
   */
  public boolean getOK() {
    return getFlag(OK);
  }

  @Override
  public String toString() {
    if (flags == 0) {
      return "None";
    }
    final StringBuilder sb = new StringBuilder();
    append(sb, getMulti(), "Multi");
    append(sb, getMultiOK(), "MultiOK");
    append(sb, getMultiDoublet(), "MultiDoublet");
    append(sb, getMultiDoubletOK(), "MultiDoubletOK");
    append(sb, getDoublet(), "Doublet");
    append(sb, getDoubletOK(), "DoubletOK");
    append(sb, getOK(), "OK");

    if (sb.length() == 0) {
      return "None";
    }
    return sb.toString();
  }

  private static void append(StringBuilder sb, boolean enabled, String name) {
    if (enabled) {
      if (sb.length() != 0) {
        sb.append("; ");
      }
      sb.append(name);
    }
  }

  /** {@inheritDoc} */
  @Override
  public FitType clone() {
    try {
      return (FitType) super.clone();
    } catch (final CloneNotSupportedException e) {
      return null;
    }
  }
}

/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Package
 *
 * Software for single molecule localisation microscopy (SMLM) in ImageJ
 * %%
 * Copyright (C) 2011 - 2023 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import ij3d.ContentNode;
import org.scijava.java3d.View;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Tuple3d;

/**
 * Node to allow display of ItemGroup content.
 */
public class ItemGroupNode extends ContentNode {
  private final ItemGroup itemGroup;

  /** The min. */
  protected Point3f min;

  /** The max. */
  protected Point3f max;

  /** The center. */
  protected Point3f center;

  /**
   * Instantiates a new item group node.
   *
   * @param itemGroup the item group
   */
  public ItemGroupNode(final ItemGroup itemGroup) {
    this.itemGroup = itemGroup;
    calculateMinMaxCenterPoint();
    addChild(itemGroup);
  }

  /**
   * Gets the item group.
   *
   * @return the item group
   */
  public ItemGroup getItemGroup() {
    return itemGroup;
  }

  @Override
  public void getMin(final Tuple3d min) {
    min.set(this.min);
  }

  @Override
  public void getMax(final Tuple3d max) {
    max.set(this.max);
  }

  @Override
  public void getCenter(final Tuple3d center) {
    center.set(this.center);
  }

  @Override
  public void channelsUpdated(final boolean[] channels) {
    // do nothing
  }

  @Override
  public void lutUpdated(final int[] red, final int[] green, final int[] blue, final int[] alpha) {
    // do nothing
  }

  @Override
  public void colorUpdated(final Color3f color) {
    itemGroup.setColor(color);
  }

  @Override
  public void eyePtChanged(final View view) {
    // do nothing
  }

  @Override
  public float getVolume() {
    return 0;
  }

  @Override
  public void shadeUpdated(final boolean shaded) {
    itemGroup.setShaded(shaded);
  }

  @Override
  public void thresholdUpdated(final int threshold) {
    // do nothing
  }

  @Override
  public void transparencyUpdated(final float transparency) {
    itemGroup.setTransparency(transparency);
  }

  private void calculateMinMaxCenterPoint() {
    min = new Point3f();
    max = new Point3f();
    center = new Point3f();
    if (itemGroup.size() == 0) {
      return;
    }

    final Point3f[] points = itemGroup.getPoints();
    CustomContentHelper.calculateMinMaxCenterPoint(min, max, center, points);
  }

  @Override
  public void restoreDisplayedData(final String path, final String name) {
    // do nothing
  }

  @Override
  public void swapDisplayedData(final String path, final String name) {
    // do nothing
  }

  @Override
  public void clearDisplayedData() {
    // do nothing
  }
}

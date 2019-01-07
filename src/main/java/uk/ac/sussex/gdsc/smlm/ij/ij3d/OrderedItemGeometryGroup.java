/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2019 Alex Herbert
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

import org.scijava.java3d.Appearance;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.Group;
import org.scijava.java3d.OrderedGroup;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;

/**
 * This class represents a list as a number of repeated shapes in the universe. The shape is defined
 * using a geometry array.
 *
 * <p>The order is fixed but can be user defined.
 *
 * @author Alex Herbert
 */
public class OrderedItemGeometryGroup extends ItemGeometryGroup implements UpdateableItemShape {
  /** The ordered group. */
  protected OrderedGroup orderedGroup;

  /**
   * Instantiates a new ordered item geometry group.
   *
   * @param points the points
   */
  public OrderedItemGeometryGroup(final Point3f[] points) {
    super(points, null, null, null, null, null);
  }

  /**
   * Instantiates a new ordered item geometry group. The geometry is scaled and translated for each
   * point. The default appearance is cloned per item and optionally can be updated with per-item
   * material colour and transparency. The transparency in the appearance is blended with per-item
   * alphas to create a per-item transparency.
   *
   * @param points the points.
   * @param ga the geometry array. If null then a default will be used.
   * @param appearance the default appearance of the shape. PolygonAttributes, Material and
   *        TransparencyAttributes are used.
   * @param sizes the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or
   *        points.length.
   * @param colors the per-item colors. Can be null.
   * @param alphas the per-item alphas. Can be null.
   */
  public OrderedItemGeometryGroup(final Point3f[] points, GeometryArray ga, Appearance appearance,
      Point3f[] sizes, Color3f[] colors, float[] alphas) {
    super(points, ga, appearance, sizes, colors, alphas);
  }

  /**
   * Gets the parent group to which all the shapes should be added.
   *
   * @return the parent group
   */
  @Override
  protected Group getParentGroup() {
    // Force the results to be ordered
    orderedGroup = new OrderedGroup();
    // orderedGroup.setCapability(OrderedGroup.ALLOW_CHILD_INDEX_ORDER_READ); // On by default
    orderedGroup.setCapability(OrderedGroup.ALLOW_CHILD_INDEX_ORDER_WRITE);
    addChild(orderedGroup);
    return orderedGroup;
  }

  /** {@inheritDoc} */
  @Override
  public void reorder(int[] indices) {
    reorderFast(indices);
  }

  /** {@inheritDoc} */
  @Override
  public void reorderFast(int[] indices) {
    orderedGroup.setChildIndexOrder(indices);
  }
}

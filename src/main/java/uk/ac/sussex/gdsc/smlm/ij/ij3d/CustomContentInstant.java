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
package uk.ac.sussex.gdsc.smlm.ij.ij3d;

import gnu.trove.map.hash.TIntObjectHashMap;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.io.OpenDialog;

import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.Group;
import org.scijava.java3d.Node;
import org.scijava.java3d.OrderedGroup;
import org.scijava.java3d.Switch;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransformGroup;
import org.scijava.java3d.View;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Matrix3f;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Vector3d;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Enumeration;

import customnode.CustomMesh;
import customnode.CustomMeshNode;
import ij3d.Content;
import ij3d.ContentInstant;
import ij3d.ContentNode;
import ij3d.UniverseSettings;
import ij3d.pointlist.PointListDialog;
import ij3d.pointlist.PointListPanel;
import ij3d.pointlist.PointListShape;
import ij3d.shapes.BoundingBox;
import ij3d.shapes.CoordinateSystem;
import isosurface.MeshGroup;
import orthoslice.MultiOrthoGroup;
import orthoslice.OrthoGroup;
import surfaceplot.SurfacePlotGroup;
import vib.PointList;
import voltex.VoltexGroup;

/**
 * Extend the ContentInstant class to avoid using an OrderedPath.
 */
public class CustomContentInstant extends ContentInstant {
  // Duplicate anything private in the super class

  // visibility flags
  private boolean locked = false;
  private boolean visible = true;
  private boolean bbVisible = false;
  private boolean coordVisible = UniverseSettings.showLocalCoordinateSystemsByDefault;
  private boolean showPL = false;

  // entries
  private ContentNode contentNode = null;

  // point list
  private PointListShape plShape = null;
  private PointListDialog plDialog = null;
  private PointListPanel plPanel = null;
  private PointList points;

  // scene graph entries
  // private final OrderedGroup ordered;
  private final Group ordered;

  private boolean available = true;

  private int customBefore = 0;
  private TIntObjectHashMap<Switch> switchMap;

  // Copy the entire contents of the super class

  /**
   * Instantiates a new custom content instant.
   *
   * @param name the name
   */
  public CustomContentInstant(final String name) {
    // Default is to be ordered
    this(name, true);
  }

  /**
   * Instantiates a new custom content instant.
   *
   * @param name the name
   * @param isOrdered Set to true to use an {@link OrderedGroup} for the children, the default is
   *        {@link Group}.
   */
  public CustomContentInstant(final String name, boolean isOrdered) {
    super(name);

    setCapability(BranchGroup.ALLOW_DETACH);
    setCapability(Node.ENABLE_PICK_REPORTING);

    // create transformation for pickeing
    localTranslate = new TransformGroup();
    localTranslate.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
    localTranslate.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    addChild(localTranslate);
    localRotate = new TransformGroup();
    localRotate.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
    localRotate.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    localTranslate.addChild(localRotate);

    ordered = (isOrdered) ? new OrderedGroup() : new Group();
    ordered.setCapability(Group.ALLOW_CHILDREN_WRITE);
    ordered.setCapability(Group.ALLOW_CHILDREN_EXTEND);
    for (int i = 0; i < 5; i++) {
      final Switch s = new Switch();
      s.setCapability(Switch.ALLOW_SWITCH_WRITE);
      s.setCapability(Switch.ALLOW_SWITCH_READ);
      s.setCapability(Group.ALLOW_CHILDREN_WRITE);
      s.setCapability(Group.ALLOW_CHILDREN_EXTEND);
      ordered.addChild(s);
    }

    localRotate.addChild(ordered);

    // create the point list
    points = new PointList();
    plShape = new PointListShape(points);
    plShape.setPickable(true);
    plPanel = new PointListPanel(name, points);
  }

  @Override
  public void displayAs(final int type) {
    if (image == null) {
      return;
    }
    // create content node and add it to the switch
    switch (type) {
      case VOLUME:
        contentNode = new VoltexGroup(this);
        break;
      case ORTHO:
        contentNode = new OrthoGroup(this);
        break;
      case SURFACE:
        contentNode = new MeshGroup(this);
        break;
      case SURFACE_PLOT2D:
        contentNode = new SurfacePlotGroup(this);
        break;
      case MULTIORTHO:
        contentNode = new MultiOrthoGroup(this);
        break;
      default:
        throw new IllegalArgumentException(
            "Specified type is neither VOLUME, ORTHO," + "SURFACE or SURFACEPLOT2D");
    }
    display(contentNode);
    // update type
    this.type = type;
  }

  public static int getDefaultThreshold(final ImagePlus imp, final int type) {
    if (type != SURFACE) {
      return 0;
    }
    final ImageStack stack = imp.getStack();
    final int d = imp.getStackSize();
    // compute stack histogram
    final int[] h = stack.getProcessor(1).getHistogram();
    for (int z = 1; z < d; z++) {
      final int[] tmp = stack.getProcessor(z + 1).getHistogram();
      for (int i = 0; i < h.length; i++) {
        h[i] += tmp[i];
      }

    }
    return imp.getProcessor().getAutoThreshold(h);
  }

  public static int getDefaultResamplingFactor(final ImagePlus imp, final int type) {
    final int w = imp.getWidth(), h = imp.getHeight();
    final int d = imp.getStackSize();
    final int max = Math.max(w, Math.max(h, d));
    switch (type) {
      case SURFACE:
        return (int) Math.ceil(max / 128f);
      case VOLUME:
        return (int) Math.ceil(max / 256f);
      case ORTHO:
        return (int) Math.ceil(max / 256f);
      case SURFACE_PLOT2D:
        return (int) Math.ceil(max / 128f);
    }
    return 1;
  }

  @Override
  public void display(final ContentNode node) {
    // remove everything if possible
    for (@SuppressWarnings("rawtypes")
    final Enumeration e = ordered.getAllChildren(); e.hasMoreElements();) {
      final Switch s = (Switch) e.nextElement();
      s.removeAllChildren();
    }
    // remove custom switch objects
    while (ordered.numChildren() > 5) {
      ordered.removeChild(ordered.numChildren() - 1);
    }
    customBefore = 0;
    if (switchMap != null) {
      switchMap.clear();
    }

    // create content node and add it to the switch
    contentNode = node;
    ((Switch) ordered.getChild(CO)).addChild(contentNode);

    // create the bounding box and add it to the switch
    final Point3d min = new Point3d();
    contentNode.getMin(min);
    final Point3d max = new Point3d();
    contentNode.getMax(max);
    BoundingBox bb = new BoundingBox(min, max);
    bb.setPickable(false);
    ((Switch) ordered.getChild(BS)).addChild(bb);
    bb = new BoundingBox(min, max, new Color3f(0, 1, 0));
    bb.setPickable(false);
    ((Switch) ordered.getChild(BB)).addChild(bb);

    // create coordinate system and add it to the switch
    final float cl = (float) Math.abs(max.x - min.x) / 5f;
    final CoordinateSystem cs = new CoordinateSystem(cl, new Color3f(0, 1, 0));
    cs.setPickable(false);
    ((Switch) ordered.getChild(CS)).addChild(cs);

    // create point list and add it to the switch
    ((Switch) ordered.getChild(PL)).addChild(plShape);

    // adjust the landmark point size properly
    plShape.setRadius((float) min.distance(max) / 100f);

    // initialize child mask of the switch
    setSwitch(BS, selected);
    setSwitch(CS, coordVisible);
    setSwitch(CO, visible);
    setSwitch(PL, showPL);

    // update type
    this.type = CUSTOM;
  }

  /**
   * Adds a custom switchable item to the content. The content can be optionally displayed. <p> The
   * before flag species if the switch should be inserted into the group before the standard
   * content, or added after. This is relevant if using an ordered group and transparent objects.
   * Any transparent object must be drawn after non-transparent objects.
   *
   * @param node the node
   * @param before the before flag
   * @return the switch number
   */
  public int addCustomSwitch(Node node, boolean before) {
    if (switchMap == null) {
      switchMap = new TIntObjectHashMap<>();
    }
    final int index = switchMap.size();
    final Switch s = new Switch();
    switchMap.put(index, s);
    s.setCapability(Switch.ALLOW_SWITCH_WRITE);
    s.setCapability(Switch.ALLOW_SWITCH_READ);
    s.setCapability(Group.ALLOW_CHILDREN_WRITE);
    s.setCapability(Group.ALLOW_CHILDREN_EXTEND);
    s.addChild(node);
    // Only a branch group can be added to a live scene
    final BranchGroup bg = new BranchGroup();
    bg.addChild(s);
    if (before) {
      ordered.insertChild(bg, 0);
      customBefore++;
    } else {
      ordered.addChild(bg);
    }
    return index; // Account for the standard switches
  }

  /**
   * Sets the custom switch content.
   *
   * @param which the which
   * @param on the on
   */
  public void setCustomSwitch(int which, final boolean on) {
    if (switchMap == null) {
      return;
    }
    final Switch s = switchMap.get(which);
    if (s == null) {
      return;
    }
    s.setWhichChild(on ? Switch.CHILD_ALL : Switch.CHILD_NONE);
  }

  private void setSwitch(int which, final boolean on) {
    // Account for custom switch data before the standard data
    which += customBefore;
    ((Switch) ordered.getChild(which)).setWhichChild(on ? Switch.CHILD_ALL : Switch.CHILD_NONE);
  }

  /*
   * ************************************************************ swapping
   *
   ***********************************************************/
  @Override
  public void swapDisplayedData() {
    if (!available) {
      return;
    }
    contentNode.swapDisplayedData(getDisplayedDataSwapfile(), getName());
    available = false;
  }

  @Override
  public void restoreDisplayedData() {
    System.out.println("restoreDisplayedData " + getName());
    if (available) {
      System.out.println("not restoring because it is not swapped");
      return;
    }
    contentNode.restoreDisplayedData(getDisplayedDataSwapfile(), getName());
    available = true;
  }

  @Override
  public void clearDisplayedData() {
    if (!available) {
      return;
    }
    contentNode.clearDisplayedData();
    available = false;
  }

  @Override
  public boolean isAvailable() {
    return available;
  }

  private String displayedDataSwapfile = null;

  private String getDisplayedDataSwapfile() {
    if (displayedDataSwapfile != null) {
      return displayedDataSwapfile;
    }
    File tmp = new File(System.getProperty("java.io.tmpdir"), "3D_Viewer");
    if (!tmp.exists()) {
      tmp.mkdirs();
    }
    tmp = new File(tmp, "displayed");
    if (!tmp.exists()) {
      tmp.mkdirs();
    }
    displayedDataSwapfile = new File(tmp, getName()).getAbsolutePath();
    return displayedDataSwapfile;
  }

  /*
   * ************************************************************ setters - visibility flags
   *
   ***********************************************************/

  @Override
  public void setVisible(final boolean b) {
    visible = b;
    setSwitch(CO, b);
    setSwitch(CS, b & coordVisible);
    // whichChild.set(BB, b && bbVisible);
    // only if hiding, hide the point list
    if (!b) {
      showPointList(false);
    }
  }

  @Override
  public void showBoundingBox(final boolean b) {
    bbVisible = b;
    setSwitch(BB, b);
  }

  @Override
  public void showCoordinateSystem(final boolean b) {
    coordVisible = b;
    setSwitch(CS, b);
  }

  @Override
  public void setSelected(final boolean selected) {
    this.selected = selected;
    final boolean sb = selected && UniverseSettings.showSelectionBox;
    setSwitch(BS, sb);
  }

  /*
   * ************************************************************ point list
   *
   ***********************************************************/

  @Override
  public void setPointListDialog(final PointListDialog p) {
    this.plDialog = p;
  }

  @Override
  public void showPointList(final boolean b) {
    if (plShape == null) {
      return;
    }

    setSwitch(PL, b);
    showPL = b;
    if (b && plDialog != null) {
      plDialog.addPointList(getName(), plPanel);
    } else if (!b && plDialog != null) {
      plDialog.removePointList(plPanel);
    }
  }

  @Override
  public void loadPointList() {
    final PointList points = PointList.load(image);
    if (points != null) {
      setPointList(points);
    }
  }

  @Override
  public void setPointList(final PointList points) {
    this.points = points;
    plPanel.setPointList(points);
    plShape.setPointList(points);
  }

  @Override
  public void savePointList() {
    String dir = OpenDialog.getDefaultDirectory();
    String n = this.getName();
    if (image != null) {
      final FileInfo fi = image.getFileInfo();
      dir = fi.directory;
      n = fi.fileName;
    }
    points.save(dir, n);
  }

  @Override
  public void savePointList(final PrintStream out) throws IOException {
    points.save(out, false);
  }

  /**
   * Adds the point list point.
   *
   * @param p the point
   * @deprecated To be removed
   */
  @Deprecated
  @Override
  public void addPointListPoint(final Point3d p) {
    points.add(p.x, p.y, p.z);
    if (plDialog != null) {
      plDialog.update();
    }
  }

  /**
   * Sets the list point point.
   *
   * @param i the index
   * @param p the point
   * @deprecated To be removed
   */
  @Override
  @Deprecated
  public void setListPointPos(final int i, final Point3d p) {
    points.placePoint(points.get(i), p.x, p.y, p.z);
  }

  @Override
  public float getLandmarkPointSize() {
    return plShape.getRadius();
  }

  @Override
  public void setLandmarkPointSize(final float r) {
    plShape.setRadius(r);
  }

  @Override
  public Color3f getLandmarkColor() {
    return plShape.getColor();
  }

  @Override
  public void setLandmarkColor(final Color3f color) {
    plShape.setColor(color);
  }

  @Override
  public PointList getPointList() {
    return points;
  }

  /**
   * Delete point list point.
   *
   * @param i the index
   * @deprecated To be removed
   */
  @Override
  @Deprecated
  public void deletePointListPoint(final int i) {
    points.remove(i);
    if (plDialog != null) {
      plDialog.update();
    }
  }

  /*
   * ************************************************************ setters - transform
   *
   **************************************************************/
  @Override
  public void toggleLock() {
    locked = !locked;
  }

  @Override
  public void setLocked(final boolean b) {
    locked = b;
  }

  @Override
  public void applyTransform(final double[] matrix) {
    applyTransform(new Transform3D(matrix));
  }

  @Override
  public void applyTransform(final Transform3D transform) {
    final Transform3D t1 = new Transform3D();
    localTranslate.getTransform(t1);
    final Transform3D t2 = new Transform3D();
    localRotate.getTransform(t2);
    t1.mul(t2);

    t1.mul(transform, t1);
    setTransform(t1);
  }

  @Override
  public void setTransform(final double[] matrix) {
    if (contentNode == null) {
      return;
    }
    setTransform(new Transform3D(matrix));
  }

  @Override
  public void setTransform(final Transform3D transform) {
    if (contentNode == null) {
      return;
    }
    final Transform3D t = new Transform3D();
    final Point3d c = new Point3d();
    contentNode.getCenter(c);

    final Matrix3f m = new Matrix3f();
    transform.getRotationScale(m);
    t.setRotationScale(m);
    // One might thing a rotation matrix has no translational
    // component, however, if the rotation is composed of
    // translation - rotation - backtranslation, it has indeed.
    final Vector3d v = new Vector3d();
    v.x = -m.m00 * c.x - m.m01 * c.y - m.m02 * c.z + c.x;
    v.y = -m.m10 * c.x - m.m11 * c.y - m.m12 * c.z + c.y;
    v.z = -m.m20 * c.x - m.m21 * c.y - m.m22 * c.z + c.z;
    t.setTranslation(v);
    localRotate.setTransform(t);

    final Vector3d v2 = new Vector3d();
    transform.get(v2);
    v2.sub(v);
    t.set(v2);
    localTranslate.setTransform(t);
  }

  /*
   * ************************************************************ setters - attributes
   *
   ***********************************************************/

  @Override
  public void setLUT(final int[] rLUT, final int[] gLUT, final int[] bLUT, final int[] aLUT) {
    this.rLUT = rLUT;
    this.gLUT = gLUT;
    this.bLUT = bLUT;
    this.aLUT = aLUT;
    if (contentNode != null) {
      contentNode.lutUpdated(rLUT, gLUT, bLUT, aLUT);
    }
  }

  @Override
  public void setChannels(final boolean[] channels) {
    final boolean channelsChanged = channels[0] != this.channels[0]
        || channels[1] != this.channels[1] || channels[2] != this.channels[2];
    if (!channelsChanged) {
      return;
    }
    this.channels = channels;
    if (contentNode != null) {
      contentNode.channelsUpdated(channels);
    }
  }

  @Override
  public void setThreshold(final int th) {
    if (th != threshold) {
      this.threshold = th;
      if (contentNode != null) {
        contentNode.thresholdUpdated(threshold);
      }
    }
  }

  @Override
  public void setShaded(final boolean b) {
    if (b != shaded) {
      this.shaded = b;
      if (contentNode != null) {
        contentNode.shadeUpdated(shaded);
      }
    }
  }

  @Override
  public boolean isShaded() {
    return shaded;
  }

  @Override
  public void setSaturatedVolumeRendering(final boolean b) {
    if (contentNode != null && type == VOLUME) {
      ((VoltexGroup) contentNode).getRenderer().getVolume().setSaturatedVolumeRendering(b);
    }
  }

  @Override
  public boolean isSaturatedVolumeRendering() {
    return contentNode != null && type == VOLUME
        && ((VoltexGroup) contentNode).getRenderer().getVolume().isSaturatedVolumeRendering();
  }

  @Override
  public void applySurfaceColors(final ImagePlus imp) {
    if (contentNode == null) {
      return;
    }
    CustomMesh mesh = null;
    switch (type) {
      case SURFACE:
        mesh = ((MeshGroup) contentNode).getMesh();
        break;
      case CUSTOM:
        mesh = ((CustomMeshNode) contentNode).getMesh();
        break;
    }
    if (mesh == null) {
      return;
    }
    mesh.loadSurfaceColorsFromImage(imp);
  }

  @Override
  public void setColor(final Color3f color) {
    if ((this.color == null && color == null)
        || (this.color != null && color != null && this.color.equals(color))) {
      return;
    }
    this.color = color;
    plShape.setColor(color);
    if (contentNode != null) {
      contentNode.colorUpdated(this.color);
    }
  }

  @Override
  public synchronized void setTransparency(float transparency) {
    transparency = transparency < 0 ? 0 : transparency;
    transparency = transparency > 1 ? 1 : transparency;
    if (Math.abs(transparency - this.transparency) < 0.01) {
      return;
    }
    this.transparency = transparency;
    if (contentNode != null) {
      contentNode.transparencyUpdated(this.transparency);
    }
  }

  /*
   * ************************************************************ UniverseListener interface
   *
   *************************************************************/
  @Override
  public void transformationStarted(final View view) {
    // Do nothing
  }

  @Override
  public void contentAdded(final Content c) {
    // Do nothing
  }

  @Override
  public void contentRemoved(final Content c) {
    if (plDialog != null) {
      plDialog.removePointList(plPanel);
    }
  }

  @Override
  public void canvasResized() {
    // Do nothing
  }

  @Override
  public void contentSelected(final Content c) {
    // Do nothing
  }

  @Override
  public void contentChanged(final Content c) {
    // Do nothing
  }

  @Override
  public void universeClosed() {
    if (plDialog != null) {
      plDialog.removePointList(plPanel);
    }
  }

  @Override
  public void transformationUpdated(final View view) {
    eyePtChanged(view);
  }

  @Override
  public void transformationFinished(final View view) {
    eyePtChanged(view);
  }

  @Override
  public void eyePtChanged(final View view) {
    if (contentNode != null) {
      contentNode.eyePtChanged(view);
    }
  }

  /*
   * ************************************************************* getters
   *
   **************************************************************/
  @Override
  public int getType() {
    return type;
  }

  @Override
  public ContentNode getContent() {
    return contentNode;
  }

  @Override
  public ImagePlus getImage() {
    return image;
  }

  @Override
  public boolean[] getChannels() {
    return channels;
  }

  @Override
  public void getRedLUT(final int[] l) {
    System.arraycopy(rLUT, 0, l, 0, rLUT.length);
  }

  @Override
  public void getGreenLUT(final int[] l) {
    System.arraycopy(gLUT, 0, l, 0, gLUT.length);
  }

  @Override
  public void getBlueLUT(final int[] l) {
    System.arraycopy(bLUT, 0, l, 0, bLUT.length);
  }

  @Override
  public void getAlphaLUT(final int[] l) {
    System.arraycopy(aLUT, 0, l, 0, aLUT.length);
  }

  @Override
  public Color3f getColor() {
    return color;
  }

  @Override
  public int getThreshold() {
    return threshold;
  }

  @Override
  public float getTransparency() {
    return transparency;
  }

  @Override
  public int getResamplingFactor() {
    return resamplingF;
  }

  @Override
  public TransformGroup getLocalRotate() {
    return localRotate;
  }

  @Override
  public TransformGroup getLocalTranslate() {
    return localTranslate;
  }

  @Override
  public void getLocalRotate(final Transform3D t) {
    localRotate.getTransform(t);
  }

  @Override
  public void getLocalTranslate(final Transform3D t) {
    localTranslate.getTransform(t);
  }

  @Override
  public boolean isLocked() {
    return locked;
  }

  @Override
  public boolean isVisible() {
    return visible;
  }

  @Override
  public boolean hasCoord() {
    return coordVisible;
  }

  @Override
  public boolean hasBoundingBox() {
    return bbVisible;
  }

  @Override
  public boolean isPLVisible() {
    return showPL;
  }
}

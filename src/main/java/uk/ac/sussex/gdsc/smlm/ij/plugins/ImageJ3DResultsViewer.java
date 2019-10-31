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

package uk.ac.sussex.gdsc.smlm.ij.plugins;

import uk.ac.sussex.gdsc.core.annotation.Nullable;
import uk.ac.sussex.gdsc.core.data.DataException;
import uk.ac.sussex.gdsc.core.data.utils.Rounder;
import uk.ac.sussex.gdsc.core.data.utils.RounderUtils;
import uk.ac.sussex.gdsc.core.data.utils.TypeConverter;
import uk.ac.sussex.gdsc.core.ij.ImageJTrackProgress;
import uk.ac.sussex.gdsc.core.ij.ImageJUtils;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog;
import uk.ac.sussex.gdsc.core.ij.gui.ExtendedGenericDialog.OptionListener;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper;
import uk.ac.sussex.gdsc.core.ij.process.LutHelper.LutColour;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicate;
import uk.ac.sussex.gdsc.core.ij.roi.CoordinatePredicateUtils;
import uk.ac.sussex.gdsc.core.logging.Ticker;
import uk.ac.sussex.gdsc.core.utils.MathUtils;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.core.utils.SortUtils;
import uk.ac.sussex.gdsc.core.utils.TextUtils;
import uk.ac.sussex.gdsc.core.utils.TurboList;
import uk.ac.sussex.gdsc.core.utils.rng.SplitMix;
import uk.ac.sussex.gdsc.smlm.data.NamedObject;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtos.PrecisionMethod;
import uk.ac.sussex.gdsc.smlm.data.config.FitProtosHelper;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettings.Builder;
import uk.ac.sussex.gdsc.smlm.data.config.GUIProtos.ImageJ3DResultsViewerSettingsOrBuilder;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsSettings;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.ij.gui.ListSelectionModelHelper;
import uk.ac.sussex.gdsc.smlm.ij.gui.PeakResultTableModel;
import uk.ac.sussex.gdsc.smlm.ij.gui.PeakResultTableModelFrame;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.CustomContent;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.CustomContentHelper;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.CustomContentInstant;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemGeometryGroup;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemGroup;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemGroupNode;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemPointMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemShape;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ItemTriangleMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.OrderedItemGeometryGroup;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.ReferenceItemMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.Shape3DHelper;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.Shape3DHelper.Rendering;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.TransparentItemPointMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.TransparentItemShape;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.TransparentItemTriangleMesh;
import uk.ac.sussex.gdsc.smlm.ij.ij3d.UpdateableItemShape;
import uk.ac.sussex.gdsc.smlm.ij.plugins.ResultsManager.InputSource;
import uk.ac.sussex.gdsc.smlm.ij.settings.SettingsManager;
import uk.ac.sussex.gdsc.smlm.results.MemoryPeakResults;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultsDigest;
import uk.ac.sussex.gdsc.smlm.results.procedures.PeakResultProcedureX;
import uk.ac.sussex.gdsc.smlm.results.procedures.PrecisionResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.RawResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.StandardResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyResultProcedure;
import uk.ac.sussex.gdsc.smlm.results.procedures.XyzResultProcedure;

import customnode.CustomLineMesh;
import customnode.CustomMesh;
import customnode.CustomMeshNode;
import customnode.CustomPointMesh;

import gnu.trove.map.hash.TObjectIntHashMap;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GUI;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.process.LUT;

import ij3d.Content;
import ij3d.ContentInstant;
import ij3d.ContentNode;
import ij3d.DefaultUniverse;
import ij3d.Image3DMenubar;
import ij3d.Image3DUniverse;
import ij3d.ImageCanvas3D;
import ij3d.ImageWindow3D;
import ij3d.UniverseListener;
import ij3d.UniverseSettings;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.lang3.tuple.Triple;
import org.scijava.java3d.Appearance;
import org.scijava.java3d.BranchGroup;
import org.scijava.java3d.Canvas3D;
import org.scijava.java3d.ColoringAttributes;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.LineAttributes;
import org.scijava.java3d.PickInfo;
import org.scijava.java3d.PickInfo.IntersectionInfo;
import org.scijava.java3d.PointAttributes;
import org.scijava.java3d.PolygonAttributes;
import org.scijava.java3d.SceneGraphPath;
import org.scijava.java3d.Shape3D;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TransformGroup;
import org.scijava.java3d.TransparencyAttributes;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.View;
import org.scijava.java3d.utils.pickfast.PickCanvas;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point2d;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3d;

import java.awt.Color;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

import javax.swing.DefaultListSelectionModel;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.WindowConstants;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * Draws a localisation results set using an ImageJ 3D image.
 *
 * @see <A href="https://imagej.net/3D_Viewer">https://imagej.net/3D_Viewer</a>
 */
public class ImageJ3DResultsViewer implements PlugIn {
  private static final String TITLE = "ImageJ 3D Results Viewer";
  private static final int NO_ENTRY = -1;

  // To debug this from Eclipse relies on being able to find the native
  // runtime libraries for Open GL. See the README in the eclipse project folder.

  private static final String[] SIZE_MODE = SettingsManager.getNames((Object[]) SizeMode.values());
  private static final String[] RENDERING = SettingsManager.getNames((Object[]) Rendering.values());
  private static final String[] DEPTH_MODE =
      SettingsManager.getNames((Object[]) DepthMode.values());
  private static final String[] TRANSPARENCY_MODE =
      SettingsManager.getNames((Object[]) TransparencyMode.values());
  private static final String[] SORT_MODE = SettingsManager.getNames((Object[]) SortMode.values());

  /** The executor service for message digests. */
  private static final ExecutorService executorService = Executors.newCachedThreadPool();

  /** The identity transform. */
  private static final Transform3D IDENTITY = new Transform3D();

  // No need to store this in settings as when the plugin is first run there are no windows
  private static AtomicReference<String> lastWindow = new AtomicReference<>("");

  private static final Map<PeakResultsDigest,
      Triple<PeakResultTableModel, ListSelectionModel, PeakResultTableModelFrame>> resultsTables =
          new ConcurrentHashMap<>();

  private static final AtomicReference<ResultsTableSettings> resultsTableSettings =
      new AtomicReference<>();
  private static final AtomicBoolean addToSelection = new AtomicBoolean();
  private static final AtomicReference<Color3f> highlightColor = new AtomicReference<>();

  private Image3DUniverse univ;
  private JMenuItem resetRotation;
  private JMenuItem resetTranslation;
  private JMenuItem resetZoom;
  private JMenuItem resetAll;
  private JMenuItem changeColour;
  private JMenuItem changePointSize;
  private JMenuItem increasePointSize;
  private JMenuItem decreasePointSize;
  private JMenuItem resetSelectedView;
  private JMenuItem findEyePoint;
  private JMenuItem sortBackToFront;
  private JMenuItem sortFrontToBack;
  private JMenuItem colourSurface;
  private JMenuItem toggleTransparent;
  private JMenuItem toggleShaded;
  private JMenuItem updateSettings;
  private JMenuItem cropResults;
  private JCheckBoxMenuItem toggleDynamicTransparency;

  /** A map between the name of the AWT colour and the colour as a SciJava colour. */
  private static HashMap<String, Color3f> colours;

  static {
    colours = new HashMap<>();
    final Field[] fields = Color.class.getFields();
    for (final Field field : fields) {
      if (Modifier.isStatic(field.getModifiers()) && field.getType() == Color.class) {
        try {
          final Color c = (Color) field.get(null);
          colours.put(field.getName().toLowerCase(Locale.US), new Color3f(c));
        } catch (final IllegalArgumentException | IllegalAccessException ex) {
          // Ignore
        }
      }
    }
  }

  //@formatter:off
  private enum SizeMode implements NamedObject
  {
    FIXED_SIZE { @Override
    public String getName() { return "Fixed"; }},
    XY_PRECISION { @Override
    public String getName() { return "XY Precision"; }},
    XYZ_DEVIATIONS { @Override
    public String getName() { return "XYZ Deviations"; }},
        ;

    static final SizeMode[] values = SizeMode.values();

    @Override
    public String getShortName()
    {
      return getName();
    }

    public static SizeMode forNumber(int number)
    {
      if (number < 0 || number >= values.length) {
        throw new IllegalArgumentException();
      }
      return values[number];
    }
  }

  private enum DepthMode implements NamedObject
  {
    NONE { @Override
    public String getName() { return "None"; }},
    INTENSITY { @Override
    public String getName() { return "Intensity"; }},
    DITHER { @Override
    public String getName() { return "Dither"; }},
        ;

    static final DepthMode[] values = DepthMode.values();

    @Override
    public String getShortName()
    {
      return getName();
    }

    public static DepthMode forNumber(int number)
    {
      if (number < 0 || number >= values.length) {
        throw new IllegalArgumentException();
      }
      return values[number];
    }
  }

  private enum TransparencyMode implements NamedObject
  {
    NONE { @Override
    public String getName() { return "None"; }},
    SIZE { @Override
    public String getName() { return "Size"; }},
    INTENSITY { @Override
    public String getName() { return "Intensity"; }},
    XY_PRECISION { @Override
    public String getName() { return "XY Precision"; }},
    XYZ_DEVIATIONS { @Override
    public String getName() { return "XYZ Deviations"; }},
        ;

    static final TransparencyMode[] values = TransparencyMode.values();

    @Override
    public String getShortName()
    {
      return getName();
    }

    public static TransparencyMode forNumber(int number)
    {
      if (number < 0 || number >= values.length) {
        throw new IllegalArgumentException();
      }
      return values[number];
    }
  }

  private enum SortMode implements NamedObject
  {
    NONE { @Override
    public String getName() { return "None"; }
    @Override
    public String getDescription() { return ""; }},
    XYZ { @Override
    public String getName() { return "XYZ"; }
    @Override
    public String getDescription() { return "Sort using XYZ. The order is defined by the direction with the major component used first, e.g. 1,2,3 for zyx ascending, -3,-2,-1 for xyz descending."; }},
    OTHOGRAPHIC { @Override
    public String getName() { return "Othographic"; }
    @Override
    public String getDescription() { return "Project all points to the plane defined by the direction. Rank by distance to the plane."; }},
    PERSPECTIVE { @Override
    public String getName() { return "Perspective"; }
    @Override
    public String getDescription() { return "Rank by distance to the eye position for true depth perspective rendering."; }},
        ;

    static final SortMode[] values = SortMode.values();

    @Override
    public String getShortName()
    {
      return getName();
    }

    public static SortMode forNumber(int number)
    {
      if (number < 0 || number >= values.length) {
        throw new IllegalArgumentException();
      }
      return values[number];
    }

    public abstract String getDescription();

    public String getDetails()
    {
      return getName() + ": " + getDescription();
    }
  }
  //@formatter:on

  private static class ResultsMetaData implements ListSelectionListener {
    PeakResultTableModel peakResultTableModel;
    ListSelectionModel listSelectionModel;

    Color3f highlightColor;
    final ImageJ3DResultsViewerSettings settings;

    /** The results when the object mesh was constructed. */
    final MemoryPeakResults results;

    /**
     * The actual positions of each item in the object mesh. The coordinates may differ from the
     * results as 2D datasets can be dithered.
     */
    final TurboList<Point3f> points;

    /** The sizes of each rendered point. */
    final Point3f[] sizes;

    /** The rendering mode. */
    final Rendering rendering;

    /** Used to test if this is the same results set. */
    PeakResultsDigest digest;

    CustomContentInstant contentInstance;
    CustomMesh outline;
    TurboList<PeakResult> selected = new TurboList<>();
    TurboList<TransformGroup> selectedNode = new TurboList<>();

    ResultsMetaData(ImageJ3DResultsViewerSettings settings, MemoryPeakResults results,
        TurboList<Point3f> points, Point3f[] sizes) {
      this.settings = settings;
      this.results = results;
      this.points = points;
      this.sizes = sizes;
      Rendering outlineRendering = Rendering.forNumber(settings.getRendering());
      if (outlineRendering.isHighResolution()) {
        // Don't draw a mesh with too much detail.
        // Note the Icosahedron does not envelope the shape but the next division does.
        outlineRendering = Rendering.LOW_RES_SPHERE;
      }
      this.rendering = outlineRendering;
      highlightColourUpdated();
    }

    void createClickSelectionNode(CustomContentInstant contentInstance) {
      // Note: Allow multiple items to be picked.
      // Maintain a list of the results that are picked (the table model).
      // At each new click, check the list does not contain the points
      // and add it.
      // Selection is handled by a selection model.
      // For multiple items we add new switches. Each new selected point uses the
      // first non-visible switch for display (or creates a new one).
      // If a point is removed then set the switch off.

      this.contentInstance = contentInstance;
      outline = createOutline();
    }

    /**
     * Create the point outline for the rendering.
     *
     * @return the custom mesh
     */
    private CustomMesh createOutline() {
      if (settings.getRendering() == 0) {
        final TurboList<Point3f> points = new TurboList<>(1);
        points.add(new Point3f());
        final ItemPointMesh mesh = new ItemPointMesh(points, highlightColor, 0);
        mesh.setPointSize((float) settings.getPixelSize());
        mesh.getAppearance().getPolygonAttributes().setPolygonMode(PolygonAttributes.POLYGON_LINE);
        return mesh;
      }

      List<Point3f> pointOutline;
      if (rendering.is2D()) {
        // Handle all the 2D objects to create an outline.
        pointOutline = Shape3DHelper.createLocalisationObjectOutline(rendering);

        final CustomLineMesh mesh =
            new CustomLineMesh(pointOutline, CustomLineMesh.CONTINUOUS, highlightColor, 0);
        mesh.setAntiAliasing(true);
        mesh.setPattern(LineAttributes.PATTERN_SOLID);
        // mesh.setLineWidth(0.5f); Hard to see
        return mesh;
      }

      // 3D objects can use the same rendering but then post-process to line polygon.
      // Note using a line mesh would work but does not cull the backface
      final Shape3D shape = Shape3DHelper.createShape(rendering, 0);

      final Appearance appearance = shape.getAppearance();

      final PolygonAttributes pa = appearance.getPolygonAttributes();
      pa.setCullFace(PolygonAttributes.CULL_BACK);
      pa.setBackFaceNormalFlip(false);
      pa.setPolygonMode(PolygonAttributes.POLYGON_LINE);
      // final ColoringAttributes ca = appearance.getColoringAttributes();
      // ca.setShadeModel(ColoringAttributes.SHADE_FLAT);
      appearance.setMaterial(null);
      final LineAttributes la = new LineAttributes();
      la.setLineWidth(0.5f);
      appearance.setLineAttributes(la);

      return new ItemMesh(new Point3f(), (GeometryArray) shape.getGeometry(), appearance, null,
          highlightColor, 0f);
    }

    void select(PeakResult result) {
      select(results.indexOf(result));
    }

    void select(int index) {
      if (index < 0 || index >= points.size()) {
        return;
      }

      final PeakResult r = results.get(index);

      // Find in the list of selected
      int switchIndex = findSelected(r);

      if (switchIndex == -1) {
        switchIndex = addToSelected(r);
      }

      // Position correctly
      final Transform3D t = new Transform3D();

      final Vector3d centre = new Vector3d(points.get(index));
      t.setTranslation(centre);

      // Enlarge the outline
      if (rendering != Rendering.POINT) {
        final Vector3d scale = new Vector3d((sizes.length == 0) ? sizes[0] : sizes[index]);
        scale.scale(1.1);
        t.setScale(scale);
      }

      selectedNode.getf(switchIndex).setTransform(t);
      contentInstance.setCustomSwitch(switchIndex, true);
    }

    private int findSelected(PeakResult result) {
      for (int i = 0; i < selected.size(); i++) {
        final PeakResult r2 = selected.getf(i);
        if (r2 != null && result.equals(r2)) {
          return i;
        }
      }
      return -1;
    }

    private int addToSelected(PeakResult result) {
      // Find the first available index
      int switchIndex = findEmpty();
      if (switchIndex == -1) {
        final TransformGroup tg = new TransformGroup();
        tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        tg.addChild(new Shape3D(outline.getGeometry(), outline.getAppearance()));
        tg.setPickable(false);
        // Add the outline before to support transparency.
        switchIndex = contentInstance.addCustomSwitch(tg, true);
        selected.add(result);
        selectedNode.add(tg);
      } else {
        selected.setf(switchIndex, result);
      }
      return switchIndex;
    }

    private int findEmpty() {
      for (int i = 0; i < selected.size(); i++) {
        final PeakResult r2 = selected.getf(i);
        if (r2 == null) {
          return i;
        }
      }
      return -1;
    }

    public void clearSelected() {
      for (int i = 0; i < selected.size(); i++) {
        contentInstance.setCustomSwitch(i, false);
        selected.setf(i, null);
      }
    }

    void setPointSize(float size) {
      if (rendering == Rendering.POINT) {
        ((ItemPointMesh) outline).setPointSize(size);
      }
    }

    float getPointSize() {
      if (rendering == Rendering.POINT) {
        ((ItemPointMesh) outline).getPointSize();
      }
      return 0;
    }

    void highlightColourUpdated() {
      this.highlightColor = ImageJ3DResultsViewer.highlightColor.get();
      if (outline != null) {
        outline.setColor(highlightColor);
      }
    }

    void addSelectionModel(
        Triple<PeakResultTableModel, ListSelectionModel, PeakResultTableModelFrame> triplet) {
      this.peakResultTableModel = triplet.getLeft();
      this.listSelectionModel = triplet.getMiddle();
      listSelectionModel.addListSelectionListener(this);
      updateSelection();
    }

    @Override
    public void valueChanged(ListSelectionEvent event) {
      if (event.getValueIsAdjusting()) {
        return;
      }
      updateSelection();
    }

    private void updateSelection() {
      final int[] indices = ListSelectionModelHelper.getSelectedIndices(listSelectionModel);

      if (indices.length == 0) {
        clearSelected();
        return;
      }

      // If the selection is output to a table then it may have been sorted
      // and we map the index from the table to the data model
      final PeakResultTableModelFrame table = findTable(this);
      if (table != null) {
        table.convertRowIndexToModel(indices);
      }

      // Try to to preserve those currently selected

      // Find all those currently selected (old selection)
      final TObjectIntHashMap<PeakResult> oldSelection =
          new TObjectIntHashMap<>(selected.size(), 0.5f, NO_ENTRY);
      for (int i = 0; i < selected.size(); i++) {
        final PeakResult r = selected.getf(i);
        if (r != null) {
          oldSelection.put(r, i);
        }
      }

      if (oldSelection.isEmpty()) {
        // Just select the new indices
        for (int i = 0; i < indices.length; i++) {
          select(peakResultTableModel.get(indices[i]));
        }
        return;
      }

      // Process the new selection, checking if already selected.
      final TurboList<PeakResult> toSelect = new TurboList<>(indices.length);
      for (int i = 0; i < indices.length; i++) {
        final PeakResult r = peakResultTableModel.get(indices[i]);
        // Check if already selected
        if (oldSelection.remove(r) == NO_ENTRY) {
          // Do this later to save space
          toSelect.addf(r);
        }
      }

      // Remove the old selection no longer required
      oldSelection.forEachEntry((result, index) -> {
        contentInstance.setCustomSwitch(index, false);
        selected.setf(index, null);
        return true;
      });

      // Select the new additions
      for (int i = toSelect.size(); i-- > 0;) {
        select(toSelect.getf(i));
      }
    }

    void removeSelectionModel() {
      if (listSelectionModel != null) {
        listSelectionModel.removeListSelectionListener(this);
        listSelectionModel = null;
      }
    }

    void removeFromSelectionModel(PeakResult result) {
      if (peakResultTableModel == null) {
        return;
      }
      // Find in the model
      int index = peakResultTableModel.indexOf(result);
      if (index != -1) {
        // If the selection is output to a table then it may have been sorted
        // and we map the index from the data to the table
        final PeakResultTableModelFrame table = findTable(this);
        if (table != null) {
          index = table.convertRowIndexToView(index);
        }
        listSelectionModel.removeSelectionInterval(index, index);
      }
    }

    void addToSelectionModel(PeakResult result) {
      if (peakResultTableModel == null) {
        return;
      }
      // Find in the model
      int index = peakResultTableModel.indexOf(result);

      if (index == -1) {
        // Not currently in the table so add it
        index = peakResultTableModel.getRowCount();
        peakResultTableModel.add(this, result);
      }

      // If the selection is output to a table then it may have been sorted
      // and we map the index from the data to the table
      PeakResultTableModelFrame table;
      if (resultsTableSettings.get().getShowTable()) {
        table = createTable(results, this);
      } else {
        table = findTable(this);
      }

      if (table != null) {
        index = table.convertRowIndexToView(index);
      }

      // Highlight the localisation using an outline.
      // Add to or replace previous selection.
      if (addToSelection.get()) {
        listSelectionModel.addSelectionInterval(index, index);
      } else {
        listSelectionModel.setSelectionInterval(index, index);
      }
    }

    private static PeakResultTableModelFrame findTable(ResultsMetaData data) {
      // There is a single TableModel and SelectionModel for each unique results set.
      // This may be displayed in a window.
      final Triple<PeakResultTableModel, ListSelectionModel, PeakResultTableModelFrame> t =
          resultsTables.get(data.digest);
      final PeakResultTableModelFrame table = t.getRight();
      if (table != null && table.isVisible()) {
        return table;
      }
      return null;
    }

    private static PeakResultTableModelFrame createTable(MemoryPeakResults results,
        ResultsMetaData data) {
      // There is a single TableModel and SelectionModel for each unique results set.
      // This is displayed in a window. Show the window if it is not visible.
      // Set the window to have dispose on close (to save memory).

      // Clicks just select from the selection model, and add results to the table model.

      final Triple<PeakResultTableModel, ListSelectionModel, PeakResultTableModelFrame> triplet =
          resultsTables.get(data.digest);

      PeakResultTableModelFrame table = triplet.getRight();
      if (table != null) {
        if (table.isVisible()) {
          return table;
        }

        // Just in case the listeners are still active
        table.cleanUp();
      }

      // No table or not visible so create a new one
      table = new PeakResultTableModelFrame(triplet.getLeft(), triplet.getMiddle());
      table.setTitle(TITLE + " " + results.getName());
      table.setReadOnly(false);
      // Ensure cleanup
      table.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
      final PeakResultTableModelFrame finalTable = table;
      table.addWindowListener(new WindowAdapter() {
        @Override
        public void windowClosed(WindowEvent event) {
          // We must unmap the selection since we use the selection model
          // across all active views of the same dataset.
          final int[] indices = ListSelectionModelHelper.getSelectedIndices(triplet.getMiddle());
          finalTable.convertRowIndexToModel(indices);
          finalTable.cleanUp(); // Remove listeners
          ListSelectionModelHelper.setSelectedIndices(triplet.getMiddle(), indices);
        }
      });
      table.setVisible(true);
      resultsTables.put(data.digest, Triple.of(triplet.getLeft(), triplet.getMiddle(), table));

      return table;
    }
  }

  private static class CustomSortObject {
    final float f1;
    final float f2;
    final float f3;
    final int index;

    CustomSortObject(int index, float f1, float f2, float f3) {
      this.f1 = f1;
      this.f2 = f2;
      this.f3 = f3;
      this.index = index;
    }

    static int compare(CustomSortObject r1, CustomSortObject r2) {
      if (r1.f1 < r2.f1) {
        return -1;
      }
      if (r1.f1 > r2.f1) {
        return 1;
      }
      if (r1.f2 < r2.f2) {
        return -1;
      }
      if (r1.f2 > r2.f2) {
        return 1;
      }
      if (r1.f3 < r2.f3) {
        return -1;
      }
      if (r1.f3 > r2.f3) {
        return 1;
      }
      return 0;
    }
  }

  private static class MouseListenerWrapper implements MouseListener {
    static final int MOUSE_CLICKED = 0x01;
    static final int MOUSE_PRESSED = 0x02;
    static final int MOUSE_RELEASED = 0x04;
    static final int MOUSE_ENTERED = 0x08;
    static final int MOUSE_EXITED = 0x10;

    final MouseListener listener;
    final int flags;

    MouseListenerWrapper(MouseListener listener, int flags) {
      this.listener = listener;
      this.flags = flags;
    }

    private boolean run(MouseEvent event, int flag) {
      return (flags & flag) != 0 && !event.isConsumed();
    }

    @Override
    public void mouseClicked(MouseEvent event) {
      if (run(event, MOUSE_CLICKED)) {
        listener.mouseClicked(event);
      }
    }

    @Override
    public void mousePressed(MouseEvent event) {
      if (run(event, MOUSE_PRESSED)) {
        listener.mousePressed(event);
      }
    }

    @Override
    public void mouseReleased(MouseEvent event) {
      if (run(event, MOUSE_RELEASED)) {
        listener.mouseReleased(event);
      }
    }

    @Override
    public void mouseEntered(MouseEvent event) {
      if (run(event, MOUSE_ENTERED)) {
        listener.mouseEntered(event);
      }
    }

    @Override
    public void mouseExited(MouseEvent event) {
      if (run(event, MOUSE_EXITED)) {
        listener.mouseExited(event);
      }
    }
  }

  @SuppressWarnings("unused")
  private static class MouseMotionListenerWrapper implements MouseMotionListener {
    static final int MOUSE_DRAGGED = 0x01;
    static final int MOUSE_MOVED = 0x02;

    final MouseMotionListener listener;
    final int flags;

    MouseMotionListenerWrapper(MouseMotionListener listener, int flags) {
      this.listener = listener;
      this.flags = flags;
    }

    private boolean run(MouseEvent event, int flag) {
      return (flags & flag) != 0 && !event.isConsumed();
    }

    @Override
    public void mouseDragged(MouseEvent event) {
      if (run(event, MOUSE_DRAGGED)) {
        listener.mouseDragged(event);
      }
    }

    @Override
    public void mouseMoved(MouseEvent event) {
      if (run(event, MOUSE_MOVED)) {
        listener.mouseMoved(event);
      }
    }
  }

  private static class LocalUniverseListener implements UniverseListener {
    @Override
    public void transformationStarted(View view) {
      // Ignore
    }

    @Override
    public void transformationUpdated(View view) {
      // Ignore.
      // This is called when the zoom is adjusted. We can update clipping if required.
    }

    @Override
    public void transformationFinished(View view) {
      // Ignore
    }

    @Override
    public void contentAdded(Content content) {
      // Ignore
    }

    @Override
    public void contentRemoved(Content content) {
      // Unregister from the selection model
      if (content != null && content.getUserData() instanceof ResultsMetaData) {
        final ResultsMetaData data = (ResultsMetaData) content.getUserData();
        data.removeSelectionModel();
      }
    }

    @Override
    public void contentChanged(Content content) {
      // Ignore
    }

    @Override
    public void contentSelected(Content content) {
      // Ignore
    }

    @Override
    public void canvasResized() {
      // Ignore
    }

    @Override
    public void universeClosed() {
      // Ignore
    }
  }

  @Override
  public void run(String arg) {
    // For testing
    // if (true || Utils.isExtraOptions())
    // {
    // new ImageJ3DResultsViewerTest().run(arg);
    // return;
    // }

    SmlmUsageTracker.recordPlugin(this.getClass(), arg);

    if (ImageJ3DViewerUtils.JAVA_3D_VERSION == null) {
      IJ.error(TITLE, "Java 3D is not available");
      return;
    }

    if (MemoryPeakResults.isMemoryEmpty()) {
      IJ.error(TITLE, "There are no fitting results in memory");
      return;
    }

    final ImageJ3DResultsViewerSettings.Builder settings =
        SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();

    addToSelection.set(settings.getAddToSelection());

    // Get a list of the window titles available. Allow the user to select
    // an existing window or a new one.
    final String title = TITLE;
    final List<Image3DUniverse> univList = new TurboList<>();
    final List<String> titleList = new TurboList<>();
    titleList.add("New window");
    buildWindowList(title, univList, titleList);
    final String[] titles = titleList.toArray(new String[titleList.size()]);

    final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
    gd.addMessage("Select a dataset to display");
    ResultsManager.addInput(gd, settings.getInputOption(), InputSource.MEMORY);
    // The behaviour is to allow the settings to store if the user prefers a new window
    // or to reuse an existing window. If 'new window' then a new window should always
    // be selected. Otherwise open in the same window as last time. If there was no last
    // window then the settings will carried over from the last ImageJ session.
    final String window = (settings.getNewWindow()) ? "" : lastWindow.get();
    gd.addChoice("Window", titles, window);
    gd.addSlider("Transparancy", 0, 0.9, settings.getTransparency(), new OptionListener<Double>() {
      @Override
      public boolean collectOptions(Double value) {
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Transparancy options", null);
        egd.addCheckbox("Support_dynamic_transparency", settings.getSupportDynamicTransparency());
        egd.addCheckbox("Enable_dynamic_transparency", settings.getEnableDynamicTransparency());
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        settings.setSupportDynamicTransparency(egd.getNextBoolean());
        settings.setEnableDynamicTransparency(egd.getNextBoolean());
        return true;
      }
    });
    gd.addChoice("Colour", LutHelper.getLutNames(), settings.getLut());
    gd.addChoice("Rendering", RENDERING, settings.getRendering(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setRendering(value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Drawing mode options", null);
        final int rendering = settings.getRendering();
        if (rendering != 0) {
          return false;
        }
        egd.addNumericField("Pixel_size", settings.getPixelSize(), 2, 6, "px");
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        settings.setPixelSize(egd.getNextNumber());
        return true;
      }
    });
    gd.addCheckbox("Shaded", settings.getShaded());
    gd.addChoice("Size_mode", SIZE_MODE, settings.getSizeMode(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setSizeMode(value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Size mode options", null);
        final SizeMode mode = SizeMode.forNumber(settings.getSizeMode());
        if (mode == SizeMode.FIXED_SIZE) {
          egd.addNumericField("Size", settings.getSize(), 2, 6, "nm");
        } else {
          // Other modes do not require options
          return false;
        }
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        settings.setSize(egd.getNextNumber());
        return true;
      }
    });
    gd.addChoice("Sort_mode", SORT_MODE, settings.getSortMode(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setSortMode(value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Sort mode options", null);
        final SortMode mode = SortMode.forNumber(settings.getSortMode());
        if (mode == SortMode.NONE) {
          return false;
        }
        egd.addMessage(
            TextUtils.wrap("Note: The sort mode is used to correctly render transparent objects. "
                + "For non-transparent objects faster rendering is achieved with a reverse "
                + "sort to put close objects at the front.", 80));
        egd.addMessage(TextUtils.wrap(mode.getDetails(), 80));
        egd.addMessage("Define the direction of the view");
        egd.addNumericField("Direction_x", settings.getSortDirectionX(), 3, 10, "");
        egd.addNumericField("Direction_y", settings.getSortDirectionY(), 3, 10, "");
        egd.addNumericField("Direction_z", settings.getSortDirectionZ(), 3, 10, "");
        if (mode == SortMode.PERSPECTIVE) {
          egd.addMessage("Define the view eye position");
          egd.addNumericField("Eye_x", settings.getSortEyeX(), 3, 10, "nm");
          egd.addNumericField("Eye_y", settings.getSortEyeY(), 3, 10, "nm");
          egd.addNumericField("Eye_z", settings.getSortEyeZ(), 3, 10, "nm");
        }
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        settings.setSortDirectionX(egd.getNextNumber());
        settings.setSortDirectionY(egd.getNextNumber());
        settings.setSortDirectionZ(egd.getNextNumber());
        if (mode == SortMode.PERSPECTIVE) {
          settings.setSortEyeX(egd.getNextNumber());
          settings.setSortEyeY(egd.getNextNumber());
          settings.setSortEyeZ(egd.getNextNumber());
        }
        return true;
      }
    });
    gd.addChoice("Transparency_mode", TRANSPARENCY_MODE, settings.getTransparencyMode(),
        new OptionListener<Integer>() {
          @Override
          public boolean collectOptions(Integer value) {
            settings.setTransparencyMode(value);
            return collectOptions(false);
          }

          @Override
          public boolean collectOptions() {
            return collectOptions(true);
          }

          private boolean collectOptions(boolean silent) {
            final ExtendedGenericDialog egd =
                new ExtendedGenericDialog("Transparency mode options", null);
            final TransparencyMode mode =
                TransparencyMode.forNumber(settings.getTransparencyMode());
            if (mode == TransparencyMode.NONE) {
              return false;
            }
            egd.addSlider("Min_transparancy", 0, 0.95, settings.getMinTransparency());
            egd.addSlider("Max_transparancy", 0, 0.95, settings.getMaxTransparency());
            egd.setSilent(silent);
            egd.showDialog(true, gd);
            if (egd.wasCanceled()) {
              return false;
            }
            settings.setMinTransparency(egd.getNextNumber());
            settings.setMaxTransparency(egd.getNextNumber());
            return true;
          }
        });
    gd.addMessage("2D options");
    gd.addChoice("Depth_mode", DEPTH_MODE, settings.getDepthMode(), new OptionListener<Integer>() {
      @Override
      public boolean collectOptions(Integer value) {
        settings.setDepthMode(value);
        return collectOptions(false);
      }

      @Override
      public boolean collectOptions() {
        return collectOptions(true);
      }

      private boolean collectOptions(boolean silent) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Depth mode options", null);
        final DepthMode mode = DepthMode.forNumber(settings.getDepthMode());
        if (mode == DepthMode.NONE) {
          return false;
        }
        egd.addNumericField("Depth_range", settings.getDepthRange(), 2, 6, "nm");
        if (mode == DepthMode.DITHER) {
          egd.addNumericField("Dither_seed", settings.getDitherSeed(), 0);
        }
        egd.setSilent(silent);
        egd.showDialog(true, gd);
        if (egd.wasCanceled()) {
          return false;
        }
        settings.setDepthRange(egd.getNextNumber());
        if (mode == DepthMode.DITHER) {
          settings.setDitherSeed((int) egd.getNextNumber());
        }
        return true;
      }
    });

    gd.showDialog();
    if (gd.wasCanceled()) {
      return;
    }
    final String name = ResultsManager.getInputSource(gd);
    final int windowChoice = gd.getNextChoiceIndex();
    lastWindow.set(titles[windowChoice]);
    settings.setInputOption(name);
    settings.setTransparency(gd.getNextNumber());
    settings.setLut(gd.getNextChoiceIndex());
    settings.setRendering(gd.getNextChoiceIndex());
    settings.setShaded(gd.getNextBoolean());
    settings.setSizeMode(gd.getNextChoiceIndex());
    settings.setSortMode(gd.getNextChoiceIndex());
    settings.setTransparencyMode(gd.getNextChoiceIndex());
    settings.setDepthMode(gd.getNextChoiceIndex());
    gd.collectOptions();

    if (windowChoice == 0) {
      // Store if the user chose a new window when they had a choice of an existing window
      if (titleList.size() > 1) {
        settings.setNewWindow(true);
        // Otherwise they had no choice so leave the preferences as they are.
      }
    } else {
      // This was not a new window
      settings.setNewWindow(false);
    }

    SettingsManager.writeSettings(settings);
    MemoryPeakResults results = ResultsManager.loadInputResults(name, false, null, null);
    if (MemoryPeakResults.isEmpty(results)) {
      IJ.error(TITLE, "No results could be loaded");
      return;
    }

    // Determine if the drawing mode is supported and compute the point size
    final Point3f[] sphereSize = createSphereSize(results, settings);
    if (sphereSize == null) {
      return;
    }

    // Cache the table settings
    resultsTableSettings.set(settings.getResultsTableSettings());

    // Create a 3D viewer.
    if (windowChoice == 0) {
      univ = createImage3DUniverse(title, titleList);
    } else {
      univ = univList.get(windowChoice - 1); // Ignore the new window
    }

    lastWindow.set(univ.getWindow().getTitle());

    results = results.copy();

    // Commence a digest
    final Future<PeakResultsDigest> futureDigest =
        PeakResultsDigest.digestLater(executorService, results.toArray());

    final TurboList<Point3f> points = getPoints(results, settings);

    final ResultsMetaData data = new ResultsMetaData(settings.build(), results, points, sphereSize);

    sort(data);

    final float[] alpha = createAlpha(results, settings, sphereSize);

    final float transparency = getTransparency(settings);

    final Color3f[] colors = createColour(results, settings);

    ContentNode contentNode;

    // Build to support true transparency (depends on settings).
    // Currently this is not supported for PointArrays as they require colouring
    // in the coordinate data.
    IJ.showStatus("Creating 3D geometry ...");
    if (settings.getSupportDynamicTransparency()) {
      final ItemGeometryGroup pointGroup =
          createItemGroup(settings, sphereSize, points, alpha, transparency, colors);
      if (pointGroup == null) {
        IJ.showStatus("");
        return;
      }

      if (settings.getEnableDynamicTransparency()) {
        final long total = points.size() + getTotalTransparentObjects(univ, name);
        activateDynamicTransparency(univ, total, settings.getEnableDynamicTransparency());
      } else {
        activateDynamicTransparency(univ, 0, settings.getEnableDynamicTransparency());
      }

      contentNode = new ItemGroupNode(pointGroup);
    } else {
      final ItemMesh mesh = createItemMesh(settings, points, sphereSize, transparency, alpha);
      if (mesh == null) {
        IJ.showStatus("");
        return;
      }

      setColour(mesh, colors);

      contentNode = new CustomMeshNode(mesh);
    }

    IJ.showStatus("Creating 3D content ...");

    // Use custom content to support adding new switchable nodes
    final CustomContent content =
        new CustomContent(name, !settings.getSupportDynamicTransparency());
    final CustomContentInstant contentInstant = (CustomContentInstant) content.getCurrent();
    contentInstant.setTransparency((float) settings.getTransparency());
    contentInstant.setShaded(settings.getShaded());
    contentInstant.showCoordinateSystem(UniverseSettings.showLocalCoordinateSystemsByDefault);
    contentInstant.display(contentNode);

    createHighlightColour(settings.getHighlightColour());

    content.setUserData(data);
    // Prevent relative rotation
    content.setLocked(true);

    // Set up the click selection node
    data.createClickSelectionNode(contentInstant);

    // Set up the results selection model
    data.digest = PeakResultsDigest.waitForDigest(futureDigest, -1);
    if (data.digest == null) {
      IJ.error(TITLE, "Failed to identify repeat results set");
      IJ.showStatus("");
      return;
    }
    Triple<PeakResultTableModel, ListSelectionModel, PeakResultTableModelFrame> triplet =
        resultsTables.get(data.digest);
    if (triplet == null) {
      triplet = Triple.of(new PeakResultTableModel(results, false,
          // Note the settings do not matter until the table is set live
          resultsTableSettings.get()), new DefaultListSelectionModel(), null);
      triplet.getLeft().setCheckDuplicates(true);
      resultsTables.put(data.digest, triplet);
    }

    // Preserve orientation on the content
    final boolean auto = univ.getAutoAdjustView();
    final Content oldContent = univ.getContent(name);
    if (oldContent == null) {
      univ.setAutoAdjustView(true);
    } else {
      univ.removeContent(name);
      univ.setAutoAdjustView(false);
    }

    IJ.showStatus("Drawing 3D content ... ");
    final StopWatch sw = StopWatch.createStarted();
    final Future<Content> future = univ.addContentLater(content);
    Content added = null;
    for (;;) {
      try {
        // Wait for 1 second
        for (int i = 0; i < 20; i++) {
          Thread.sleep(50);
          if (future.isDone()) {
            // Only get the result when finished, so avoiding a blocking wait
            added = future.get();
            break;
          }
        }
        if (added != null) {
          break;
        }

        final long seconds = sw.getTime(TimeUnit.SECONDS);
        if (seconds % 20 == 0) {
          final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE, null);
          egd.addMessage("Current wait time is " + sw.toString());
          egd.setOKLabel("Wait");
          egd.showDialog();
          if (egd.wasCanceled()) {
            future.cancel(true);
            break;
          }
        }

        IJ.showStatus("Drawing 3D content ... " + seconds);
      } catch (final InterruptedException ex) {
        Thread.currentThread().interrupt();
      } catch (final ExecutionException ex) {
        break;
      }
    }
    univ.setAutoAdjustView(auto);

    // Initialise the selection model
    if (added != null) {
      data.addSelectionModel(triplet);
    }

    IJ.showStatus("");
  }

  private static Point3f[] createSphereSize(MemoryPeakResults results, Builder settings) {
    // Special case for point rendering
    if (settings.getRendering() == 0) {
      final float size = getFixedSize(settings.getPixelSize());
      return new Point3f[] {new Point3f(size, size, size)};
    }

    // Store XYZ size for each localisation
    final SizeMode mode = SizeMode.forNumber(settings.getSizeMode());
    switch (mode) {
      case FIXED_SIZE:
        final float size = getFixedSize(settings.getSize());
        final Point3f[] sizes = new Point3f[results.size()];
        Arrays.fill(sizes, new Point3f(size, size, size));
        return sizes;
      case XYZ_DEVIATIONS:
        return createSphereSizeFromDeviations(results);
      case XY_PRECISION:
        return createSphereSizeFromPrecision(results);
      default:
        throw new IllegalStateException("Unknown drawing mode: " + mode);
    }
  }

  private static float getFixedSize(double size) {
    return (size > 0) ? (float) size : 1f;
  }

  @Nullable
  private static Point3f[] createSphereSizeFromDeviations(MemoryPeakResults results) {
    if (!results.hasDeviations()) {
      IJ.error(TITLE, "The results have no deviations");
      return null;
    }
    // Currently the rendering is in nm
    final TypeConverter<DistanceUnit> dc = results.getDistanceConverter(DistanceUnit.NM);

    final Point3f[] size = new Point3f[results.size()];
    final boolean failed = results.forEach(new PeakResultProcedureX() {
      int index;

      @Override
      public boolean execute(PeakResult peakResult) {
        final float x = peakResult.getParameterDeviation(PeakResult.X);
        final float y = peakResult.getParameterDeviation(PeakResult.Y);
        float z = peakResult.getParameterDeviation(PeakResult.Z);
        // Check x & y are not zero.
        // This should be OK as 2D fitting should provide these.
        if (x == 0 || y == 0) {
          return true;
        }
        if (z == 0) {
          z = (float) Math.sqrt((x * x + y * y) / 2); // Mean variance
        }
        // z = (x + y) / 2; // Mean Std Dev
        size[index++] = new Point3f(dc.convert(x), dc.convert(y), dc.convert(z));
        return false;
      }
    });
    return (failed) ? null : size;
  }

  @Nullable
  private static Point3f[] createSphereSizeFromPrecision(MemoryPeakResults results) {
    final PrecisionResultProcedure p = new PrecisionResultProcedure(results);
    try {
      final PrecisionMethod m = p.getPrecision();
      IJ.log("Using precision method " + FitProtosHelper.getName(m));
      final Point3f[] size = new Point3f[results.size()];
      for (int i = 0, j = 0; i < p.precisions.length; i++) {
        // Precision is in NM which matches the rendering
        final float v = (float) p.precisions[i];
        size[j++] = new Point3f(v, v, v);
      }
      return size;
    } catch (final DataException ex) {
      IJ.error(TITLE, "The results have no precision: " + ex.getMessage());
      return null;
    }
  }

  @Nullable
  private static float[] createAlpha(MemoryPeakResults results, Builder settings,
      Point3f[] sphereSize) {
    if (settings.getTransparencyMode() == 0) {
      return null;
    }

    final double min = MathUtils.clip(0, 1, settings.getMinTransparency());
    final double max = MathUtils.clip(0, 1, settings.getMaxTransparency());
    if (min == max) {
      // No per item transparency
      ImageJUtils.log("No per-item transparency as min == max");
      return null;
    }
    // Convert to alpha
    final double minA = 1 - max;
    final double maxA = 1 - min;

    final SizeMode sizeMode = SizeMode.forNumber(settings.getSizeMode());
    final TransparencyMode mode = TransparencyMode.forNumber(settings.getTransparencyMode());
    switch (mode) {
      case INTENSITY:
        return createAlphaFromIntensity(results, minA, maxA);
      case XYZ_DEVIATIONS:
        if (sizeMode != SizeMode.XYZ_DEVIATIONS) {
          sphereSize = createSphereSizeFromDeviations(results);
        }
        return createAlphaFromSize(minA, maxA, sphereSize);
      case XY_PRECISION:
        if (sizeMode != SizeMode.XY_PRECISION) {
          sphereSize = createSphereSizeFromPrecision(results);
        }
        return createAlphaFromSize(minA, maxA, sphereSize);
      case SIZE:
        return createAlphaFromSize(results, settings, minA, maxA, sphereSize);
      default:
        throw new IllegalStateException("Unknown transparency mode: " + mode);
    }
  }

  @Nullable
  private static float[] createAlphaFromIntensity(MemoryPeakResults results, double minA,
      double maxA) {
    final RawResultProcedure p = new RawResultProcedure(results);
    p.getI();

    final float[] intensity = p.intensity;
    final float[] limits = MathUtils.limits(intensity);
    final float min = limits[0];
    final float max = limits[1];
    if (min == max) {
      ImageJUtils.log("No per-item transparency as intensity is fixed");
      return null;
    }

    final double range = (maxA - minA) / (max - min);
    final float[] alpha = new float[intensity.length];
    for (int i = 0; i < alpha.length; i++) {
      // Lowest intensity has lowest alpha (more transparent)
      alpha[i] = (float) (minA + range * (intensity[i] - min));
    }
    return alpha;
  }

  @Nullable
  private static float[] createAlphaFromSize(MemoryPeakResults results, Builder settings,
      double minA, double maxA, Point3f[] sphereSize) {
    final SizeMode sizeMode = SizeMode.forNumber(settings.getSizeMode());
    if (sizeMode == SizeMode.FIXED_SIZE) {
      ImageJUtils.log("No per-item transparency as size is fixed");
      return null;
    }

    if (settings.getRendering() == 0) {
      // No size was created for fixed point rendering so create it now
      switch (sizeMode) {
        case XYZ_DEVIATIONS:
          sphereSize = createSphereSizeFromDeviations(results);
          break;
        case XY_PRECISION:
          sphereSize = createSphereSizeFromPrecision(results);
          break;
        default:
          throw new IllegalStateException("Unknown drawing mode: " + sizeMode);
      }
    }

    return createAlphaFromSize(minA, maxA, sphereSize);
  }

  @Nullable
  private static float[] createAlphaFromSize(double minA, double maxA, Point3f[] sphereSize) {
    if (sphereSize == null) {
      return null;
    }

    final double[] d = new double[sphereSize.length];
    for (int i = 0; i < d.length; i++) {
      final Point3f p = sphereSize[i];
      if (p.x == p.y && p.y == p.z) {
        d[i] = 3.0 * p.x * p.x;
      } else {
        // Use the squared distance. This is the equivalent of the area of the shape projected to
        // 2D.
        d[i] = (double) p.x * p.x + p.y * p.y + p.z * p.z;
      }

      // Use the average radius. This is the equivalent of the mean radius of an enclosing
      // ellipsoid.
      // d[i] = Math.sqrt(d[i] / 3);
    }
    final double[] limits = MathUtils.limits(d);
    final double min = limits[0];
    final double max = limits[1];
    if (min == max) {
      ImageJUtils.log("No per-item transparency as size is fixed");
      return null;
    }

    final double range = (maxA - minA) / (max - min);
    final float[] alpha = new float[d.length];
    for (int i = 0; i < alpha.length; i++) {
      // Largest distance has lowest alpha (more transparent)
      alpha[i] = (float) (minA + range * (max - d[i]));
    }
    return alpha;
  }

  private static float getTransparency(ImageJ3DResultsViewerSettingsOrBuilder settings) {
    float transparency = MathUtils.clip(0, 1, (float) settings.getTransparency());
    // Do not allow fully transparent objects
    if (transparency == 1) {
      transparency = 0;
    }
    return transparency;
  }

  /**
   * Gets the points.
   *
   * @param results the results
   * @param settings the settings
   * @return the points
   */
  static TurboList<Point3f> getPoints(MemoryPeakResults results,
      ImageJ3DResultsViewerSettingsOrBuilder settings) {
    final TurboList<Point3f> points = new TurboList<>(results.size());
    if (results.is3D()) {
      results.forEach(DistanceUnit.NM, (XyzResultProcedure) (x, y, z) -> {
        points.addf(new Point3f(x, y, z));
      });
    } else {
      results.forEach(DistanceUnit.NM, (XyResultProcedure) (x, y) -> {
        points.addf(new Point3f(x, y, 0));
      });

      final double range = settings.getDepthRange();
      if (range > 0 && results.size() > 1) {
        final DepthMode mode = DepthMode.forNumber(settings.getDepthMode());
        final double min = -settings.getDepthRange() / 2;
        switch (mode) {
          case DITHER:
            final SplitMix r = SplitMix.new64(settings.getDitherSeed());
            for (int i = points.size(); i-- > 0;) {
              points.getf(i).z += (min + r.nextDouble() * range);
            }
            break;
          case INTENSITY:
            // Rank by intensity, highest first
            final StandardResultProcedure p = new StandardResultProcedure(results);
            p.getI();
            final int[] indices = SimpleArrayUtils.natural(results.size());
            SortUtils.sortIndices(indices, p.intensity, true);
            final double inc = range / indices.length;
            for (int i = 0; i < indices.length; i++) {
              // The standard rendering has +z going away so put the highest rank at min
              points.getf(indices[i]).z += (min + i * inc);
            }
            break;
          case NONE:
            break;
          default:
            throw new IllegalStateException("Unknown depth mode: " + mode);
        }
      }
    }
    return points;
  }

  private static void sort(ResultsMetaData data) {
    final SortMode mode = SortMode.forNumber(data.settings.getSortMode());
    switch (mode) {
      case NONE:
        return;
      case PERSPECTIVE:
        sortPerspective(data);
        break;
      case OTHOGRAPHIC:
        sortOrthographic(data);
        break;
      case XYZ:
        sortXyz(data);
        break;
      default:
        throw new IllegalStateException("Unknown sort mode " + mode);
    }
  }

  private static void sortPerspective(ResultsMetaData data) {
    final ImageJ3DResultsViewerSettingsOrBuilder settings = data.settings;
    final Vector3d direction = getViewDirection(data.settings);
    if (direction == null) {
      throw new IllegalStateException("The view direction is not valid");
    }
    final Point3d eye =
        new Point3d(settings.getSortEyeX(), settings.getSortEyeY(), settings.getSortEyeZ());

    final double[] d = getDistance(data.points, direction, eye);

    final int[] indices = SimpleArrayUtils.natural(d.length);
    SortUtils.sortIndices(indices, d, true);

    reorder(indices, data);
  }

  private static double[] getDistance(TurboList<Point3f> points, Vector3d direction, Point3d eye) {
    final double[] d = new double[points.size()];
    for (int i = 0; i < d.length; i++) {
      final Point3f p = points.getf(i);

      final Vector3d v2 = new Vector3d(p.x - eye.x, p.y - eye.y, p.z - eye.z);

      // Compute distance of all points from the eye.
      d[i] = v2.length();

      // We need to know if point is in-front of the eye or behind.
      // Compute dot product (if positive then this is an acute angle).
      // We use a descending sort so acute angles (in front of the eye)
      // should be higher (ranked first) and obtuse angles (behind the
      // eye) ranked later.
      if (v2.dot(direction) < 0) {
        d[i] = -d[i];
      }
    }
    return d;
  }

  private static Vector3d getViewDirection(ImageJ3DResultsViewerSettingsOrBuilder settings) {
    final Vector3d dir = new Vector3d(settings.getSortDirectionX(), settings.getSortDirectionY(),
        settings.getSortDirectionZ());
    final double l1 = dir.lengthSquared();
    if (!Double.isFinite(l1)) {
      return null;
    }
    return dir;
  }

  private static void reorder(int[] indices, ResultsMetaData data) {
    final MemoryPeakResults results = data.results;
    final TurboList<Point3f> points = data.points;

    final PeakResult[] originalPeakResults = results.toArray();
    final Point3f[] originalPoints = points.toArray(new Point3f[points.size()]);

    // We need another array to store the output
    final PeakResult[] peakResults = new PeakResult[originalPeakResults.length];

    // Rewrite order
    for (int i = 0; i < indices.length; i++) {
      final int index = indices[i];
      points.setf(i, originalPoints[index]);
      peakResults[i] = originalPeakResults[index];
    }

    // Bulk update the results
    results.setSortAfterEnd(false);
    results.begin();
    results.addAll(peakResults);
    results.end();

    final Point3f[] sizes = data.sizes;
    if (sizes.length == indices.length) {
      final Point3f[] originalSizes = sizes.clone();
      // Rewrite order
      for (int i = 0; i < indices.length; i++) {
        final int index = indices[i];
        sizes[i] = originalSizes[index];
      }
    }
  }

  private static void sortOrthographic(ResultsMetaData data) {
    Vector3d direction = getViewDirection(data.settings);
    if (direction == null) {
      // Default
      direction = new Vector3d(0, 0, -1);
    }
    direction.normalize();
    final double a = direction.x;
    final double b = direction.y;
    final double c = direction.z;

    final TurboList<Point3f> points = data.points;

    final double[] d = new double[points.size()];
    for (int i = 0; i < d.length; i++) {
      final Point3f p = points.getf(i);

      // Compute signed distance of all points from the plane
      // defined by normal v and point (0,0,0)
      d[i] = a * p.x + b * p.y + c * p.z;
    }

    final int[] indices = SimpleArrayUtils.natural(d.length);
    SortUtils.sortIndices(indices, d, true);

    reorder(indices, data);
  }

  private static void sortXyz(ResultsMetaData data) {
    Vector3d direction = getViewDirection(data.settings);
    if (direction == null) {
      // Default to z axis (away), then y then x in ascending order
      direction = new Vector3d(-1, -2, 3);
    }
    direction.normalize();

    // Use the vector lengths in each dimension to set the order
    int[] indices = SimpleArrayUtils.natural(3);
    final double[] values = new double[] {direction.x, direction.y, direction.z};
    final double[] absValues =
        new double[] {Math.abs(direction.x), Math.abs(direction.y), Math.abs(direction.z)};
    SortUtils.sortIndices(indices, absValues, true);

    final int ix = search(indices, 0);
    final int iy = search(indices, 1);
    final int iz = search(indices, 2);

    // Use the vector signs to set the direction
    final int sx = (values[0] <= 0) ? 1 : -1;
    final int sy = (values[1] <= 0) ? 1 : -1;
    final int sz = (values[2] <= 0) ? 1 : -1;

    // Sort using the points since these have dithered positions for 2D results.
    final TurboList<Point3f> points = data.points;
    final CustomSortObject[] toSort = new CustomSortObject[points.size()];
    final float[] f = new float[3];
    for (int i = 0; i < toSort.length; i++) {
      final Point3f p = points.getf(i);
      f[ix] = sx * p.x;
      f[iy] = sy * p.y;
      f[iz] = sz * p.z;
      toSort[i] = new CustomSortObject(i, f[0], f[1], f[2]);
    }

    Arrays.sort(toSort, CustomSortObject::compare);

    indices = new int[toSort.length];
    for (int i = 0; i < toSort.length; i++) {
      indices[i] = toSort[i].index;
    }

    reorder(indices, data);
  }

  private static int search(int[] values, int key) {
    for (int i = 0; i < values.length; i++) {
      if (values[i] == key) {
        return i;
      }
    }
    return -1;
  }

  private static void updateAppearance(CustomMesh mesh,
      final ImageJ3DResultsViewerSettingsOrBuilder settings) {
    mesh.setShaded(settings.getShaded());

    final Appearance appearance = mesh.getAppearance();
    final PolygonAttributes pa = appearance.getPolygonAttributes();

    // For all 3D polygons we want to support a true face orientation so transparency works
    final Rendering r = Rendering.forNumber(settings.getRendering());
    if (r.is2D()) {
      pa.setCullFace(PolygonAttributes.CULL_NONE);
      pa.setBackFaceNormalFlip(true);
    } else {
      pa.setCullFace(PolygonAttributes.CULL_BACK);
      pa.setBackFaceNormalFlip(false);
    }

    // TransparencyAttributes ta = appearance.getTransparencyAttributes();
    // ta.setSrcBlendFunction(TransparencyAttributes.BLEND_SRC_ALPHA);
    // ta.setDstBlendFunction(TransparencyAttributes.BLEND_ONE);
    // ta.setDstBlendFunction(TransparencyAttributes.BLEND_ONE_MINUS_SRC_ALPHA); // Default

    ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.FASTEST);
    // ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.SCREEN_DOOR);
    // ItemTriangleMesh.setTransparencyMode(TransparencyAttributes.BLENDED);

    final ColoringAttributes ca = appearance.getColoringAttributes();
    if (r.isHighResolution() || r.is2D()) {
      // Smooth across vertices. Required to show 2D surfaces smoothly
      ca.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
    } else {
      // Faster polygon rendering with flat shading
      ca.setShadeModel(ColoringAttributes.SHADE_FLAT);
    }
  }

  private static void createHighlightColour(String highlightColour) {
    highlightColor.set(null);
    for (int i = 0; i < highlightColour.length(); i++) {
      if (Character.isDigit(highlightColour.charAt(i))) {
        // Try and extract RGB
        final String[] split = highlightColour.split("[\\s,:]+");
        if (split.length > 2) {
          try {
            // RGB
            final int red = Integer.parseInt(split[0]);
            final int green = Integer.parseInt(split[1]);
            final int blue = Integer.parseInt(split[2]);
            highlightColor.set(new Color3f(new Color(red, green, blue)));
            return;
          } catch (final NumberFormatException ex) {
            // Ignore
          }
        }
      }
    }
    highlightColor.set(colours.get(highlightColour.toLowerCase(Locale.US)));
  }

  @SuppressWarnings("null")
  private static Color3f[] createColour(MemoryPeakResults results,
      ImageJ3DResultsViewerSettingsOrBuilder settings) {
    // Colour by z
    final LUT lut = LutHelper.createLut(LutColour.forNumber(settings.getLut()), false);

    StandardResultProcedure procedure = null;
    float range = 0;
    float[] limits = null;
    if (results.is3D()) {
      procedure = new StandardResultProcedure(results);
      procedure.getZ();
      limits = MathUtils.limits(procedure.z);
      range = limits[1] - limits[0];
    }

    if (range == 0) {
      return new Color3f[] {new Color3f(new Color(lut.getRGB(255)))};
    }

    // Create 256 Colors
    final float scale = 255f / range;
    final Color3f[] colors = new Color3f[256];
    for (int i = 0; i < 256; i++) {
      final Color c = new Color(lut.getRGB(i));
      colors[i] = new Color3f(c);
    }

    final float minimum = limits[0];
    final Color3f[] allColors = new Color3f[results.size()];
    for (int i = 0, size = results.size(); i < size; i++) {
      float value = procedure.z[i];
      value = value - minimum;
      if (value < 0f) {
        value = 0f;
      }
      int ivalue = (int) ((value * scale) + 0.5f);
      if (ivalue > 255) {
        ivalue = 255;
      }
      allColors[i] = colors[ivalue];
    }
    return allColors;
  }

  private static void changeColour(ItemShape itemShape, MemoryPeakResults results,
      ImageJ3DResultsViewerSettingsOrBuilder settings) {
    setColour(itemShape, createColour(results, settings));
  }

  private static void setColour(ItemShape itemShape, Color3f[] colours) {
    if (colours.length == 1) {
      itemShape.setItemColor(colours[0]);
    } else {
      itemShape.setItemColor(colours);
    }
  }

  /**
   * Builds the window list of all visible windows starting with the title prefix.
   *
   * @param titlePrefix the title prefix
   * @param univList the univ list
   * @param titleList the title list
   */
  private static void buildWindowList(String titlePrefix, List<Image3DUniverse> univList,
      List<String> titleList) {
    for (final Image3DUniverse univ : Image3DUniverse.universes) {
      final ImageWindow3D w = univ.getWindow();
      if (w != null && w.isVisible() && w.getTitle().startsWith(titlePrefix)) {
        univList.add(univ);
        titleList.add(w.getTitle());
      }
    }
  }

  /**
   * Creates the image 3D universe with a unique name.
   *
   * @param title the title
   * @param titleList the title list (of titles to ignore)
   * @return the image 3D universe
   */
  private Image3DUniverse createImage3DUniverse(String title, List<String> titleList) {
    // Get a unique name by appending numbers to the end
    String title2 = title;
    int counter = 2;
    while (titleList.contains(title2)) {
      title2 = title + " " + (counter++);
    }

    final Image3DUniverse universe = new Image3DUniverse();

    universe.addUniverseListener(new LocalUniverseListener());

    universe.setShowBoundingBoxUponSelection(false);
    universe.showAttribute(DefaultUniverse.ATTRIBUTE_SCALEBAR, false);

    // Capture a canvas mouse click/region and identify the coordinates.
    final ImageCanvas3D canvas = (ImageCanvas3D) universe.getCanvas();
    final BranchGroup scene = universe.getScene();

    final MouseListener mouseListener = new MouseAdapter() {
      @Override
      public void mouseClicked(final MouseEvent event) {
        if (!consumeEvent(event)) {
          return;
        }

        // Consume the event
        event.consume();

        // This finds the vertex indices of the rendered object.
        final Pair<Content, IntersectionInfo> pair =
            getPickedContent(canvas, scene, event.getX(), event.getY());
        if (pair == null) {
          universe.select(null); // Do the same as the mouseClicked in Image3DUniverse
          return;
        }

        // Only process content added from localisations
        final Content c = pair.getKey();
        if (!(c.getUserData() instanceof ResultsMetaData)) {
          universe.select(c); // Do the same as the mouseClicked in Image3DUniverse
          return;
        }

        final ResultsMetaData data = (ResultsMetaData) c.getUserData();

        final MemoryPeakResults results = data.results;

        // Look up the localisation from the clicked vertex
        final ContentInstant content = c.getCurrent();
        int index = -1;
        if (content.getContent() instanceof CustomMeshNode) {
          final CustomMeshNode node = (CustomMeshNode) content.getContent();
          final CustomMesh mesh = node.getMesh();
          int vertexCount;
          final GeometryArray ga = (GeometryArray) mesh.getGeometry();
          // Default to the number of vertices
          vertexCount = ga.getValidVertexCount();

          final int countPerLocalisation = vertexCount / results.size();

          // Determine the localisation
          final int vertexIndex = pair.getValue().getVertexIndices()[0];
          index = vertexIndex / countPerLocalisation;
        } else if (content.getContent() instanceof ItemGroupNode) {
          // All shapes have the index as the user data
          final Object o = pair.getValue().getGeometry().getUserData();
          if (o instanceof Integer) {
            index = (Integer) pair.getValue().getGeometry().getUserData();
          }
        }
        if (index == -1) {
          return;
        }

        final PeakResult result = results.get(index);

        if (event.getClickCount() > 1) {
          // Centre on the localisation
          final Point3d coordinate = new Point3d();
          coordinate.set(data.points.get(index));

          // Handle the local transform of the content ...
          final Transform3D vWorldToLocal = getVworldToLocal(content);
          vWorldToLocal.transform(coordinate);

          universe.centerAt(coordinate);
        } else if (event.isShiftDown()) {
          // Ctrl+Shift held down to remove selected
          data.removeFromSelectionModel(result);
        } else {
          // Ctrl held down to set selection
          data.addToSelectionModel(result);
        }
      }

      private boolean consumeEvent(final MouseEvent event) {
        // Consume left-mouse clicks with the Ctrl or Alt key down.
        // Single clicks only if showing the results table.
        // Double clicks for centring the universe.

        if (event.isConsumed() || event.getButton() != MouseEvent.BUTTON1
            || !(event.isControlDown())) {
          return false;
        }
        if (event.getClickCount() == 1) {
          return true;
        }
        return (event.getClickCount() == 2);
      }
    };

    // 0 = ImageCanvas3D
    // 1 = DefaultUniverse
    // 2 = Image3DUniverse
    final MouseListener[] l = canvas.getMouseListeners();
    for (int i = 0; i < l.length; i++) {
      if (l[i].getClass().getName().contains("Image3DUniverse")) {
        // We want to be before the Image3DUniverse to allow consuming the click event.
        // Only allow the click event.
        // This disables the right-click pop-up menu.
        // It doesn't have anything of use for localisations anyway.
        canvas.removeMouseListener(l[i]);
        canvas.addMouseListener(mouseListener);
        canvas.addMouseListener(new MouseListenerWrapper(l[i], MouseListenerWrapper.MOUSE_CLICKED
        // |MouseListenerWrapper.MOUSE_PRESSED|MouseListenerWrapper.MOUSE_RELEASED
        ));
      }
    }

    // 0 = ImageCanvas3D
    // 1 = DefaultUniverse
    // 2 = Image3DUniverse
    // 3 = EventCatcher (from scijava)
    final MouseMotionListener[] ml = canvas.getMouseMotionListeners();
    for (int i = 0; i < ml.length; i++) {
      if (ml[i].getClass().getName().contains("Image3DUniverse")) {
        // Ignore this as it just shows the name in the IJ status bar
        canvas.removeMouseMotionListener(ml[i]);
      }
    }

    // Finally display the window
    universe.show();
    final ImageWindow3D window = universe.getWindow();
    GUI.center(window);
    window.setTitle(title2);

    WindowManager.addWindow(window);
    window.addWindowListener(new WindowAdapter() {
      @Override
      public void windowClosing(WindowEvent event) {
        WindowManager.removeWindow(window);
      }
    });

    // Add a new menu for SMLM functionality
    createSmlmMenuBar(universe);

    return universe;
  }

  private static long getTotalTransparentObjects(Image3DUniverse univ, String ignoreName) {
    long size = 0;
    for (final Iterator<Content> it = univ.contents(); it.hasNext();) {
      final Content c = it.next();
      if (ignoreName.equals(c.getName())) {
        continue;
      }
      if (!(c.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }
      final ContentInstant content = c.getCurrent();
      if (content.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) content.getContent();
        final ItemGroup g = node.getItemGroup();
        if (g instanceof ItemGeometryGroup) {
          size += g.size();
        }
      }
    }
    return size;
  }

  private static void activateDynamicTransparency(Image3DUniverse univ, long size, boolean enable) {
    if (size == 0) {
      IJ.log("No transparent objects");
      enable = false;
    }

    final View view = univ.getViewer().getView();
    final boolean isEnabled =
        view.getTransparencySortingPolicy() == View.TRANSPARENCY_SORT_GEOMETRY;
    if (enable == isEnabled) {
      return;
    }

    if (enable) {
      if (size > 20000L) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
        egd.addMessage("The results contain " + size
            + " transparent objects.\nDynamic transparency may take a long time to render.");
        egd.setOKLabel("Dynamic");
        egd.setCancelLabel("Standard");
        egd.showDialog();
        if (egd.wasCanceled()) {
          return;
        }
      }

      view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);
      IJ.log("Enabled dynamic transparency");
      // Q. I am not sure if this is required if objects are sorted?
      // view.setDepthBufferFreezeTransparent(false)
    } else {
      view.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_NONE);
      IJ.log("Disabled dynamic transparency");
    }

    final JMenuBar menubar = univ.getMenuBar();
    final JMenu menu = menubar.getMenu(menubar.getMenuCount() - 1);
    for (int i = 0; i < menu.getItemCount(); i++) {
      final JMenuItem item = menu.getItem(i);
      if (item == null) {
        continue;
      }
      if (item.getText().equals("Toggle dynamic transparency")) {
        ((JCheckBoxMenuItem) item).setSelected(enable);
      }
    }
  }

  /**
   * Get the Content and closest intersection point at the specified canvas position.
   *
   * <p>Adapted from Picker.getPickedContent(...).
   *
   * @param canvas the canvas
   * @param scene the scene
   * @param x the x
   * @param y the y
   * @return the Content and closest intersection point
   */
  private static Pair<Content, IntersectionInfo> getPickedContent(Canvas3D canvas,
      BranchGroup scene, final int x, final int y) {
    final PickCanvas pickCanvas = new PickCanvas(canvas, scene);
    pickCanvas.setMode(PickInfo.PICK_GEOMETRY);
    pickCanvas.setFlags(PickInfo.SCENEGRAPHPATH | PickInfo.CLOSEST_GEOM_INFO
    // | PickInfo.CLOSEST_INTERSECTION_POINT
    );
    pickCanvas.setTolerance(3);
    pickCanvas.setShapeLocation(x, y);
    try {
      final PickInfo[] result = pickCanvas.pickAllSorted();
      if (result == null) {
        return null;
      }
      for (int i = 0; i < result.length; i++) {
        final SceneGraphPath path = result[i].getSceneGraphPath();
        Content content = null;
        for (int j = path.nodeCount(); j-- > 0;) {
          if (path.getNode(j) instanceof Content) {
            content = (Content) path.getNode(j);
          }
        }
        if (content == null) {
          continue;
        }
        return Pair.of(content, result[i].getIntersectionInfos()[0]);
      }
      return null;
    } catch (final Exception ex) {
      return null;
    }
  }

  /**
   * Creates the SMLM menu bar.
   *
   * @param univ the universe
   */
  private void createSmlmMenuBar(Image3DUniverse univ) {
    final Image3DMenubar menubar = (Image3DMenubar) univ.getMenuBar();

    final JMenu menu = new JMenu("GDSC SMLM");
    menu.setMnemonic(KeyEvent.VK_G);

    resetRotation = new JMenuItem("Reset global rotation", KeyEvent.VK_R);
    resetRotation.addActionListener(this::menuActionPerformed);
    menu.add(resetRotation);

    resetTranslation = new JMenuItem("Reset global translation", KeyEvent.VK_T);
    resetTranslation.addActionListener(this::menuActionPerformed);
    menu.add(resetTranslation);

    resetZoom = new JMenuItem("Reset global zoom", KeyEvent.VK_Z);
    resetZoom.addActionListener(this::menuActionPerformed);
    menu.add(resetZoom);

    menu.addSeparator();

    resetAll = new JMenuItem("Reset all transformations", KeyEvent.VK_A);
    resetAll.addActionListener(this::menuActionPerformed);
    menu.add(resetAll);

    resetSelectedView = new JMenuItem("Reset selected transformation", KeyEvent.VK_S);
    resetSelectedView.addActionListener(this::menuActionPerformed);
    menu.add(resetSelectedView);

    findEyePoint = new JMenuItem("Find eye point", KeyEvent.VK_E);
    findEyePoint.addActionListener(this::menuActionPerformed);
    menu.add(findEyePoint);

    sortBackToFront = new JMenuItem("Sort Back-to-Front", KeyEvent.VK_B);
    sortBackToFront.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed B"));
    sortBackToFront.addActionListener(this::menuActionPerformed);
    menu.add(sortBackToFront);

    sortFrontToBack = new JMenuItem("Sort Front-to-Back", KeyEvent.VK_F);
    sortFrontToBack.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed R"));
    sortFrontToBack.addActionListener(this::menuActionPerformed);
    menu.add(sortFrontToBack);

    menu.addSeparator();

    changeColour = new JMenuItem("Change colour", KeyEvent.VK_O);
    changeColour.addActionListener(this::menuActionPerformed);
    menu.add(changeColour);

    changePointSize = new JMenuItem("Change point size", KeyEvent.VK_H);
    changePointSize.addActionListener(this::menuActionPerformed);
    menu.add(changePointSize);

    increasePointSize = new JMenuItem("Increase point size", KeyEvent.VK_I);
    increasePointSize.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed PERIOD"));
    increasePointSize.addActionListener(this::menuActionPerformed);
    menu.add(increasePointSize);

    decreasePointSize = new JMenuItem("Decrease point size", KeyEvent.VK_D);
    decreasePointSize.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed COMMA"));
    decreasePointSize.addActionListener(this::menuActionPerformed);
    menu.add(decreasePointSize);

    toggleTransparent = new JMenuItem("Toggle transparent", KeyEvent.VK_P);
    toggleTransparent.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed E"));
    toggleTransparent.addActionListener(this::menuActionPerformed);
    menu.add(toggleTransparent);

    toggleShaded = new JMenuItem("Toggle shaded", KeyEvent.VK_S);
    toggleShaded.addActionListener(this::menuActionPerformed);
    menu.add(toggleShaded);

    toggleDynamicTransparency = new JCheckBoxMenuItem("Toggle dynamic transparency");
    toggleDynamicTransparency.setMnemonic(KeyEvent.VK_D);
    toggleDynamicTransparency.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed D"));
    toggleDynamicTransparency.setSelected(
        univ.getViewer().getView().getTransparencySortingPolicy() != View.TRANSPARENCY_SORT_NONE);
    toggleDynamicTransparency.addActionListener(this::menuActionPerformed);
    menu.add(toggleDynamicTransparency);

    colourSurface = new JMenuItem("Colour surface from 2D image", KeyEvent.VK_I);
    colourSurface.addActionListener(this::menuActionPerformed);
    menu.add(colourSurface);

    menu.addSeparator();

    cropResults = new JMenuItem("Crop results", KeyEvent.VK_C);
    cropResults.setAccelerator(KeyStroke.getKeyStroke("ctrl pressed X"));
    cropResults.addActionListener(this::menuActionPerformed);
    menu.add(cropResults);

    menu.addSeparator();

    updateSettings = new JMenuItem("Update settings", KeyEvent.VK_U);
    updateSettings.addActionListener(this::menuActionPerformed);
    menu.add(updateSettings);

    menubar.add(menu);
    // Add back so it is redrawn
    univ.setMenubar(menubar);
    // 4.0.3 method
    // univ.refreshShortcuts();
  }

  private interface ContentAction {
    /**
     * Run the action.
     *
     * @param content The content
     * @return negative for error. No further content can be processed.
     */
    public int run(Content content);

    public void finish();
  }

  private abstract static class BaseContentAction implements ContentAction {
    @Override
    public void finish() {
      // Ignore
    }
  }

  private static class ChangeColourContentAction extends BaseContentAction {
    ImageJ3DResultsViewerSettings.Builder settings;

    @Override
    public int run(Content content) {
      if (!(content.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }

      final ResultsMetaData data = (ResultsMetaData) content.getUserData();

      final MemoryPeakResults results = data.results;

      // Change the colour
      if (settings == null) {
        // Use the latest settings
        settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        // Transparency can be set interactively using: Edit > Change Transparency
        // so not option here.
        gd.addChoice("Colour", LutHelper.getLutNames(), settings.getLut());
        gd.showDialog();
        if (gd.wasCanceled()) {
          return -1;
        }
        settings.setLut(gd.getNextChoiceIndex());
        SettingsManager.writeSettings(settings);
      }

      final ContentInstant contentInstant = content.getCurrent();
      ItemShape itemShape = null;
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        if (mesh instanceof ItemShape) {
          itemShape = (ItemShape) mesh;
        }
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        itemShape = node.getItemGroup();
      }
      if (itemShape != null) {
        changeColour(itemShape, results, settings);
      }
      return 0;
    }
  }

  private static class ChangePointSizeContentAction extends BaseContentAction {
    float pointSize = -1;

    @Override
    public int run(Content content) {
      final ContentInstant contentInstant = content.getCurrent();
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        if (mesh instanceof ItemMesh) {
          final ItemMesh t = (ItemMesh) mesh;
          if (t.isPointArray()) {
            // Change the point size
            if (!getSettings()) {
              return -1;
            }
            t.getAppearance().getPointAttributes().setPointSize(pointSize);
          }
        } else if (mesh instanceof CustomPointMesh) {
          // Change the point size
          if (!getSettings()) {
            return -1;
          }
          ((CustomPointMesh) mesh).setPointSize(pointSize);
        }
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        if (!getSettings()) {
          return -1;
        }
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup g = node.getItemGroup();
        g.setPointSize(pointSize);
      }

      if (content.getUserData() instanceof ResultsMetaData) {
        final ResultsMetaData data = (ResultsMetaData) content.getUserData();
        data.setPointSize(pointSize);
      }

      return 0;
    }

    private boolean getSettings() {
      if (pointSize == -1) {
        // Use the latest settings
        final ImageJ3DResultsViewerSettings.Builder settings =
            SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
        final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
        gd.addNumericField("Pixel_size", settings.getPixelSize(), 2, 6, "px");
        gd.showDialog();
        if (gd.wasCanceled()) {
          return false;
        }
        settings.setPixelSize(gd.getNextNumber());
        SettingsManager.writeSettings(settings);
        pointSize = (float) settings.getPixelSize();
      }
      return true;
    }
  }

  private static class UpdatePointSizeContentAction extends BaseContentAction {
    static final UpdatePointSizeContentAction INCREASE = new UpdatePointSizeContentAction(true);
    static final UpdatePointSizeContentAction DECREASE = new UpdatePointSizeContentAction(false);

    private final boolean increase;

    private UpdatePointSizeContentAction(boolean increase) {
      this.increase = increase;
    }

    @Override
    public int run(Content content) {
      final ContentInstant contentInstant = content.getCurrent();
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        if (mesh instanceof ItemMesh) {
          final ItemMesh t = (ItemMesh) mesh;
          if (t.isPointArray()) {
            final PointAttributes pa = t.getAppearance().getPointAttributes();
            pa.setPointSize(updatePointSize(pa.getPointSize()));
          }
        } else if (mesh instanceof CustomPointMesh) {
          CustomPointMesh cmesh = ((CustomPointMesh) mesh);
          cmesh.setPointSize(updatePointSize(cmesh.getPointSize()));
        }
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup g = node.getItemGroup();
        g.setPointSize(updatePointSize(g.getPointSize()));
      }

      if (content.getUserData() instanceof ResultsMetaData) {
        final ResultsMetaData data = (ResultsMetaData) content.getUserData();
        data.setPointSize(updatePointSize(data.getPointSize()));
      }

      return 0;
    }

    private float updatePointSize(float pointSize) {
      return increase ? pointSize * 1.5f : pointSize / 1.5f;
    }
  }

  private static class ResetViewContentAction extends BaseContentAction {
    final boolean error;

    public ResetViewContentAction(boolean error) {
      this.error = error;
    }

    @Override
    public int run(Content content) {
      if (content.isLocked()) {
        if (error) {
          IJ.error(TITLE, content.getName() + " is locked");
          return -1;
        }
        return 0;
      }
      final Transform3D transform = new Transform3D();
      content.setTransform(transform);
      return 0;
    }
  }

  private static class ToggleTransparentAction extends BaseContentAction {
    private static class TransparencyData {
      float transparency;
      float[] alpha;

      public void save(CustomMesh mesh) {
        transparency = mesh.getTransparency();
        mesh.setTransparency(0);
        if (mesh instanceof ItemMesh) {
          final ItemMesh itemMesh = (ItemMesh) mesh;
          if (itemMesh.hasColor4()) {
            save((TransparentItemShape) mesh);
          }
        } else if (mesh instanceof TransparentItemShape) {
          save((TransparentItemShape) mesh);
        }
      }

      private void save(TransparentItemShape shape) {
        final int size = shape.size();
        if (alpha == null || alpha.length != size) {
          alpha = new float[size];
        }
        shape.getItemAlpha(alpha);
        shape.setItemAlpha(1);
      }

      public void save(ItemGroup group) {
        if (group instanceof ItemGeometryGroup) {
          save((ItemGeometryGroup) group);
        } else {
          transparency = group.getTransparency();
          group.setTransparency(0f);
        }
      }

      public void save(ItemGeometryGroup group) {
        final int size = group.size();
        if (alpha == null || alpha.length != size) {
          alpha = new float[size];
        }
        group.getItemAlpha(alpha);
        transparency = group.getTransparency();
        group.setItemAlpha(1f, 0f);
      }

      public void restore(CustomMesh mesh) {
        mesh.setTransparency(transparency);
        if (mesh instanceof ItemMesh) {
          final ItemMesh t = (ItemMesh) mesh;
          if (t.hasColor4()) {
            restore((TransparentItemShape) mesh);
          }
        } else if (mesh instanceof TransparentItemShape) {
          restore((TransparentItemShape) mesh);
        }
      }

      private void restore(TransparentItemShape shape) {
        if (alpha != null && alpha.length == shape.size()) {
          shape.setItemAlpha(alpha);
        }
      }

      public void restore(ItemGroup group) {
        if (group instanceof ItemGeometryGroup) {
          restore((ItemGeometryGroup) group);
        } else {
          group.setTransparency(transparency);
        }
      }

      public void restore(ItemGeometryGroup group) {
        final int size = group.size();
        if (alpha != null && alpha.length == size) {
          group.setItemAlpha(alpha, transparency);
        }
      }
    }

    @Override
    public int run(Content content) {
      if (!(content.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }
      final ContentInstant contentInstant = content.getCurrent();
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();

        // Polygons can just switch the transparency mode
        final TransparencyAttributes ta = mesh.getAppearance().getTransparencyAttributes();
        boolean off;
        if (ta.getTransparencyMode() == TransparencyAttributes.NONE) {
          ta.setTransparencyMode(ItemTriangleMesh.getTransparencyMode());
          off = false;
        } else {
          ta.setTransparencyMode(TransparencyAttributes.NONE);
          off = true;
        }

        // The point mesh does not support the transparency mode switching off.
        // So switch the actual transparency.
        if (mesh instanceof CustomPointMesh
            || (mesh instanceof ItemMesh && ((ItemMesh) mesh).isPointArray())) {
          TransparencyData data;
          if (mesh.getUserData() instanceof TransparencyData) {
            data = (TransparencyData) mesh.getUserData();
          } else {
            data = new TransparencyData();
            mesh.setUserData(data);
          }

          if (off) {
            data.save(mesh);
          } else {
            data.restore(mesh);
          }
        }
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup group = node.getItemGroup();

        TransparencyData data;
        if (group.getUserData() instanceof TransparencyData) {
          data = (TransparencyData) group.getUserData();
        } else {
          data = new TransparencyData();
          group.setUserData(data);
        }

        // All shapes have their own transparency to switch
        if (group.isTransparent()) {
          data.save(group);
        } else {
          data.restore(group);
        }
      }

      return 0;
    }
  }

  private static class ToggleShadedAction extends BaseContentAction {
    @Override
    public int run(Content content) {
      if (!(content.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }
      final ContentInstant contentInstant = content.getCurrent();
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        mesh.setShaded(!mesh.isShaded());
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup g = node.getItemGroup();
        g.setShaded(!g.isShaded());
      }

      return 0;
    }
  }

  private class FindEyePointContentAction extends BaseContentAction {
    ImageJ3DResultsViewerSettings.Builder settings;
    final Point3d eyePtInVWorld = new Point3d();
    final Point3d dir0InVWorld = new Point3d();
    final Point3d dir1InVWorld = new Point3d(0, 0, -1);

    Point3d eye;
    Vector3d direction;

    FindEyePointContentAction() {
      this.settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
      if (!settings.getSaveEyePoint()) {
        settings = null;
      }

      final Transform3D ipToVWorld = new Transform3D();
      univ.getCanvas().getImagePlateToVworld(ipToVWorld);

      univ.getCanvas().getCenterEyeInImagePlate(eyePtInVWorld);
      ipToVWorld.transform(eyePtInVWorld);

      // Work out where the camera is looking in the virtual world
      final Transform3D cameraToVWorld = new Transform3D();
      univ.getVworldToCameraInverse(cameraToVWorld);
      cameraToVWorld.transform(dir0InVWorld);
      cameraToVWorld.transform(dir1InVWorld);
    }

    @Override
    public int run(Content content) {
      final Transform3D vWorldToLocal = getVworldToLocal(content.getCurrent());

      final boolean identity = vWorldToLocal.equals(IDENTITY);

      eye = new Point3d(eyePtInVWorld);

      final Point3d dir0InLocalWorld = new Point3d(dir0InVWorld);
      final Point3d dir1InLocalWorld = new Point3d(dir1InVWorld);

      if (!identity) {
        // Since we require the eye-point relative to the local world
        // the matrix is inverted
        vWorldToLocal.invert();
        vWorldToLocal.transform(eye);
        vWorldToLocal.transform(dir0InLocalWorld);
        vWorldToLocal.transform(dir1InLocalWorld);
      }

      direction = new Vector3d();
      direction.sub(dir1InLocalWorld, dir0InLocalWorld);

      // Print the eye coords and direction in the virtual world.
      // This can be used for a custom sort.
      final Rounder rounder = RounderUtils.create(4);
      ImageJUtils.log("%s : Eye point = (%s,%s,%s) : Direction = (%s,%s,%s)", content.getName(),
          rounder.round(eye.x), rounder.round(eye.y), rounder.round(eye.z),
          rounder.round(direction.x), rounder.round(direction.y), rounder.round(direction.z));

      if (settings != null) {
        settings.setSortEyeX(eye.x);
        settings.setSortEyeY(eye.y);
        settings.setSortEyeZ(eye.z);
        settings.setSortDirectionX(direction.x);
        settings.setSortDirectionY(direction.y);
        settings.setSortDirectionZ(direction.z);
      }

      return 0;
    }

    @Override
    public void finish() {
      if (settings != null) {
        SettingsManager.writeSettings(settings);
      }
    }
  }

  private class SortContentAction extends FindEyePointContentAction {
    final boolean reverse;

    SortContentAction(boolean reverse) {
      this.reverse = reverse;
    }

    @Override
    public int run(Content content) {
      final int result = super.run(content);
      if (result != 0) {
        return result;
      }

      final ContentInstant contentInstant = content.getCurrent();
      if (content.getUserData() == null) {
        return 0;
      }
      UpdateableItemShape updateable = null;
      boolean reorderData = false;
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        if (!(mesh instanceof UpdateableItemShape)) {
          return 0;
        }
        updateable = (UpdateableItemShape) mesh;
        reorderData = true;
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup g = node.getItemGroup();
        if (!(g instanceof UpdateableItemShape)) {
          return 0;
        }
        updateable = (UpdateableItemShape) g;
      }

      final ResultsMetaData data = (ResultsMetaData) content.getUserData();

      if (reverse) {
        direction.negate();
      }

      final double[] d = getDistance(data.points, direction, eye);

      final int[] indices = SimpleArrayUtils.natural(d.length);
      SortUtils.sortIndices(indices, d, true);

      if (updateable != null) {
        // Switch to fast mode when not debugging
        updateable.reorderFast(indices);
      }

      // We reorder the data that is used to create colours and clicked point size.
      // This is not needed for the ItemGeometryNode as it uses the indices directly.
      // Do this second as points is updated in-line so may break reordering the mesh
      // if it has a reference to the points list (e.g. ItemPointMesh initially uses
      // the points list but will create a new internal list when it is re-ordered).
      if (reorderData) {
        reorder(indices, data);
      }

      return 0;
    }
  }

  private static class ColourSurfaceContentAction extends BaseContentAction {
    private static final AtomicReference<Options> optionsRef = new AtomicReference<>(new Options());

    private static class Options {
      String title = "";
      boolean resetTransparency = true;
    }

    @Override
    public int run(Content content) {
      if (!(content.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }

      final Options options = optionsRef.get();

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      final String[] list = ImageJUtils.getImageList(ImageJUtils.SINGLE);
      if (list.length == 0) {
        return -1;
      }
      gd.addChoice("Image", list, options.title);
      gd.addCheckbox("Reset_transparency", options.resetTransparency);
      gd.showDialog();
      if (gd.wasCanceled()) {
        return -1;
      }
      options.title = gd.getNextChoice();
      options.resetTransparency = gd.getNextBoolean();
      ImagePlus imp = WindowManager.getImage(options.title);
      optionsRef.set(options);
      if (imp == null) {
        return -1;
      }

      final ContentInstant contentInstant = content.getCurrent();
      ItemShape itemShape = null;
      if (contentInstant.getContent() instanceof CustomMeshNode) {
        final CustomMeshNode node = (CustomMeshNode) contentInstant.getContent();
        final CustomMesh mesh = node.getMesh();
        if (!(mesh instanceof ItemShape)) {
          return 0;
        }
        if (options.resetTransparency) {
          mesh.setTransparency(0);
        }
        itemShape = (ItemShape) mesh;
      } else if (contentInstant.getContent() instanceof ItemGroupNode) {
        final ItemGroupNode node = (ItemGroupNode) contentInstant.getContent();
        final ItemGroup g = node.getItemGroup();
        if (options.resetTransparency) {
          g.setTransparency(0);
        }
        itemShape = g;
      }
      if (itemShape != null) {
        CustomContentHelper.loadSurfaceColorsFromImage2D(itemShape, imp);
      }
      return 0;
    }
  }

  private class CropResultsAction extends BaseContentAction {
    CoordinatePredicate shape;
    private final Point2d p2d = new Point2d();
    ImageJ3DResultsViewerSettings.Builder settings;

    @Override
    public int run(Content content) {
      if (!(content.getUserData() instanceof ResultsMetaData)) {
        return 0;
      }

      final Canvas3D canvas = univ.getCanvas();
      if (shape == null) {
        final ImageCanvas3D canvas3D = (ImageCanvas3D) canvas;
        Roi roi = canvas3D.getRoi();
        if (roi == null) {
          roi = new Roi(0, 0, canvas.getWidth(), canvas.getHeight());
        }
        shape = CoordinatePredicateUtils.createContainsPredicate(roi);
        if (shape == null) {
          return -1;
        }
        settings = SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();
      }

      final ContentInstant contentInstant = content.getCurrent();

      // Adapted from CustomTriangleMesh.retain
      final Transform3D virtualWorldToIp = new Transform3D();
      canvas.getImagePlateToVworld(virtualWorldToIp);
      virtualWorldToIp.invert();

      final Transform3D vWorldToLocal = getVworldToLocal(contentInstant);
      if (!virtualWorldToIp.equals(IDENTITY)) {
        vWorldToLocal.invert(); // To make localToVworld
        virtualWorldToIp.mul(vWorldToLocal); // To make localToIP
      }

      final ResultsMetaData data = (ResultsMetaData) content.getUserData();

      // Transform all points to the image plate and test if they are in the ROI
      final MemoryPeakResults results = data.results;
      final TurboList<Point3f> points = data.points;
      final MemoryPeakResults newResults = new MemoryPeakResults();
      newResults.copySettings(results);
      // Get the output name
      String outputName;
      if (settings.getNameOption() == CropResults.NAME_OPTION_NAME) {
        final ExtendedGenericDialog egd = new ExtendedGenericDialog("Crop results");
        final String name =
            (TextUtils.isNullOrEmpty(settings.getOutputName())) ? (results.getName() + " Cropped")
                : settings.getOutputName();
        egd.addStringField("Output_name", name, MathUtils.clip(60, 120, name.length()));
        egd.showDialog();
        if (egd.wasCanceled()) {
          return -1;
        }
        settings.setOutputName(egd.getNextString());
        outputName = settings.getOutputName();
        if (TextUtils.isNullOrEmpty(outputName)) {
          IJ.error(TITLE, "No output name");
          return -1;
        }
      } else if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX) {
        final String suffix = settings.getNameSuffix();
        if (TextUtils.isNullOrEmpty(suffix)) {
          IJ.error(TITLE, "No output suffix");
          return -1;
        }
        outputName = results.getName() + suffix;
      } else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE) {
        outputName = results.getName();
        final String suffix = settings.getNameSuffix();
        if (!TextUtils.isNullOrEmpty(suffix)) {
          outputName += suffix;
        }
        final int counter = settings.getNameCounter();
        outputName += counter;
        settings.setNameCounter(counter + 1); // Increment for next time
      } else {
        IJ.error(TITLE, "No output name");
        return -1;
      }
      newResults.setName(outputName);
      newResults.begin();
      ImageJUtils.showStatus("Cropping " + results.getName());
      final Ticker ticker = ImageJUtils.createTicker(results.size(), 0);
      final PeakResult[] allResults = results.toArray();
      for (int i = 0, size = results.size(); i < size; i++) {
        final Point3d locInImagePlate = new Point3d(points.getf(i));
        virtualWorldToIp.transform(locInImagePlate);
        canvas.getPixelLocationFromImagePlate(locInImagePlate, p2d);
        if (shape.test(p2d.x, p2d.y)) {
          newResults.add(allResults[i]);
        }
        ticker.tick();
      }

      final int size = newResults.size();
      ImageJUtils.finished("Cropped " + TextUtils.pleural(size, "result"));
      newResults.end();
      if (size != 0) {
        MemoryPeakResults.addResults(newResults);
      }

      return 0;
    }

    @Override
    public void finish() {
      if (settings != null) {
        SettingsManager.writeSettings(settings, 0);
      }
    }
  }

  private static class UpdateHighlightColourAction extends BaseContentAction {
    @Override
    public int run(Content content) {
      if (content.getUserData() instanceof ResultsMetaData) {
        final ResultsMetaData data = (ResultsMetaData) content.getUserData();
        data.highlightColourUpdated();
      }
      return 0;
    }
  }

  private static Transform3D getVworldToLocal(ContentInstant content) {
    final Transform3D vWorldToLocal = new Transform3D();
    final Transform3D translate = new Transform3D();
    final Transform3D rotate = new Transform3D();
    content.getLocalTranslate(translate);
    content.getLocalRotate(rotate);
    vWorldToLocal.mul(translate, rotate);
    return vWorldToLocal;
  }

  private void menuActionPerformed(ActionEvent event) {
    final Object src = event.getSource();

    ContentAction action = null;

    // Universe actions
    // Adapted from univ.resetView()
    if (src == resetRotation) {
      univ.fireTransformationStarted();
      // rotate so that y shows downwards
      final Transform3D t = new Transform3D();
      final AxisAngle4d aa = new AxisAngle4d(1, 0, 0, Math.PI);
      t.set(aa);
      univ.getRotationTG().setTransform(t);
      univ.fireTransformationUpdated();
      univ.fireTransformationFinished();
      return;
    }
    if (src == resetTranslation) {
      univ.fireTransformationStarted();
      final Transform3D t = new Transform3D();
      univ.getTranslateTG().setTransform(t);
      univ.recalculateGlobalMinMax();
      univ.getViewPlatformTransformer().centerAt(univ.getGlobalCenterPoint());
      univ.fireTransformationUpdated();
      univ.fireTransformationFinished();
      return;
    }
    if (src == resetZoom) {
      univ.fireTransformationStarted();
      final Transform3D t = new Transform3D();
      univ.getZoomTG().setTransform(t);
      final Point3d max = new Point3d();
      final Point3d min = new Point3d();
      univ.getGlobalMaxPoint(max);
      univ.getGlobalMinPoint(min);
      final float range = (float) (max.x - min.x);
      final double d = (range) / Math.tan(Math.PI / 8);
      univ.getViewPlatformTransformer().zoomTo(d);
      univ.fireTransformationUpdated();
      univ.fireTransformationFinished();
      return;
    }
    if (src == updateSettings) {
      final ImageJ3DResultsViewerSettings.Builder settings =
          SettingsManager.readImageJ3DResultsViewerSettings(0).toBuilder();

      final ExtendedGenericDialog gd = new ExtendedGenericDialog(TITLE);
      final ResultsSettings.Builder s = ResultsSettings.newBuilder();
      ResultsTableSettings localResultsTableSettings = resultsTableSettings.get();
      s.setResultsTableSettings(localResultsTableSettings);
      gd.addMessage("Click on the image to view localisation data.\nCtrl/Alt key must be pressed.");
      final TextField[] tf = new TextField[1];
      gd.addStringField("Highlight_colour", settings.getHighlightColour(),
          new OptionListener<String>() {
            @Override
            public boolean collectOptions(String value) {
              createHighlightColour(value);
              int red;
              int green;
              int blue;
              final Color3f color = highlightColor.get();
              if (color == null) {
                red = blue = 0;
                green = 255;
              } else {
                red = (int) (color.x * 255);
                green = (int) (color.y * 255);
                blue = (int) (color.z * 255);
              }
              final ExtendedGenericDialog egd = new ExtendedGenericDialog("Highlight colour", null);
              egd.addSlider("Red", 0, 255, red);
              egd.addSlider("Green", 0, 255, green);
              egd.addSlider("Blue", 0, 255, blue);
              egd.showDialog(true, gd);
              if (egd.wasCanceled()) {
                return false;
              }
              red = (int) egd.getNextNumber();
              green = (int) egd.getNextNumber();
              blue = (int) egd.getNextNumber();
              final Color c = new Color(red, green, blue);
              final String cvalue = c.getRed() + "," + c.getGreen() + "," + c.getBlue();
              tf[0].setText(cvalue);
              return true;
            }

            @Override
            public boolean collectOptions() {
              return false;
            }
          });
      tf[0] = gd.getLastTextField();
      gd.addCheckbox("Add_to_selection", settings.getAddToSelection());
      ResultsManager.addTableResultsOptions(gd, s, ResultsManager.FLAG_NO_SECTION_HEADER);
      gd.addMessage("Allow the 'Find Eye Point' command to save to settings");
      gd.addCheckbox("Save_eye_point", settings.getSaveEyePoint());
      // Same as CropResults
      gd.addChoice("Crop_name_option", CropResults.NAME_OPTIONS, settings.getNameOption(),
          new OptionListener<Integer>() {
            @Override
            public boolean collectOptions(Integer value) {
              settings.setNameOption(value);
              final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
              if (settings.getNameOption() == CropResults.NAME_OPTION_NAME) {
                return false;
              } else if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX) {
                final String name = (TextUtils.isNullOrEmpty(settings.getNameSuffix())) ? " Cropped"
                    : settings.getNameSuffix();
                egd.addStringField("Name_suffix", name, MathUtils.clip(20, 60, name.length()));
              } else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE) {
                final String name = settings.getNameSuffix();
                egd.addStringField("Name_suffix", name, MathUtils.clip(20, 60, name.length()));
                int counter = settings.getNameCounter();
                if (counter < 1) {
                  counter = 1;
                }
                egd.addNumericField("Name_counter", counter, 0);
              } else {
                throw new IllegalStateException("Unknown name option: " + settings.getNameOption());
              }
              egd.showDialog(true, gd);
              if (egd.wasCanceled()) {
                return false;
              }
              if (settings.getNameOption() == CropResults.NAME_OPTION_SUFFIX) {
                settings.setNameSuffix(egd.getNextString());
              } else if (settings.getNameOption() == CropResults.NAME_OPTION_SEQUENCE) {
                settings.setNameSuffix(egd.getNextString());
                settings.setNameCounter(Math.max(1, (int) egd.getNextNumber()));
              }

              return true;
            }

            @Override
            public boolean collectOptions() {
              return false;
            }
          });
      gd.addCheckbox("Update_existing_tables", localResultsTableSettings.getUpdateExistingTables());
      gd.showDialog();
      if (gd.wasCanceled()) {
        return;
      }
      settings.setHighlightColour(gd.getNextString());
      boolean add = gd.getNextBoolean();
      addToSelection.set(add);
      settings.setAddToSelection(add);
      final ResultsTableSettings.Builder resultsTableSettingsBuilder =
          s.getResultsTableSettingsBuilder();
      resultsTableSettingsBuilder.setShowTable(gd.getNextBoolean());
      settings.setSaveEyePoint(gd.getNextBoolean());
      settings.setNameOption(gd.getNextChoiceIndex());
      resultsTableSettingsBuilder.setUpdateExistingTables(gd.getNextBoolean());

      createHighlightColour(settings.getHighlightColour());

      // Save updated settings
      localResultsTableSettings = resultsTableSettingsBuilder.build();
      settings.setResultsTableSettings(localResultsTableSettings);
      SettingsManager.writeSettings(settings);

      // Update the table settings for all the selection models
      if (resultsTableSettingsBuilder.getUpdateExistingTables()) {
        for (final Triple<PeakResultTableModel, ?, ?> t : resultsTables.values()) {
          t.getLeft().setTableSettings(localResultsTableSettings);
        }
      }

      action = new UpdateHighlightColourAction();
    }
    if (src == toggleDynamicTransparency) {
      final long total = getTotalTransparentObjects(univ, "");
      final View view = univ.getViewer().getView();
      final boolean activate = view.getTransparencySortingPolicy() == View.TRANSPARENCY_SORT_NONE;
      activateDynamicTransparency(univ, total, activate);
      return;
    }

    // Actions to perform on content
    if (src == changeColour) {
      action = new ChangeColourContentAction();
    } else if (src == resetAll) {
      univ.resetView();
      univ.select(null);
      action = new ResetViewContentAction(false);
    } else if (src == resetSelectedView) {
      action = new ResetViewContentAction(true);
    } else if (src == findEyePoint) {
      action = new FindEyePointContentAction();
    } else if (src == sortBackToFront) {
      action = new SortContentAction(false);
    } else if (src == sortFrontToBack) {
      action = new SortContentAction(true);
    } else if (src == colourSurface) {
      action = new ColourSurfaceContentAction();
    } else if (src == toggleTransparent) {
      action = new ToggleTransparentAction();
    } else if (src == toggleShaded) {
      action = new ToggleShadedAction();
    } else if (src == changePointSize) {
      action = new ChangePointSizeContentAction();
    } else if (src == increasePointSize) {
      action = UpdatePointSizeContentAction.INCREASE;
    } else if (src == decreasePointSize) {
      action = UpdatePointSizeContentAction.DECREASE;
    } else if (src == cropResults) {
      action = new CropResultsAction();
    }
    if (action == null) {
      return;
    }

    if (univ.getSelected() != null) {
      action.run(univ.getSelected());
    } else {
      for (final Iterator<Content> it = univ.contents(); it.hasNext();) {
        if (action.run(it.next()) < 0) {
          break;
        }
      }
    }

    action.finish();
  }

  private static ItemGeometryGroup createItemGroup(
      final ImageJ3DResultsViewerSettings.Builder settings, final Point3f[] sphereSize,
      final TurboList<Point3f> points, float[] alpha, float transparency, Color3f[] colors) {
    final Rendering rendering = Rendering.forNumber(settings.getRendering());
    // All objects have colour using the appearance not per vertex colours.
    // The exception is points which do not support colour from appearance.
    final int colorDepth = (rendering == Rendering.POINT) ? 4 : 0;
    final Shape3D shape = Shape3DHelper.createShape(rendering, colorDepth);

    // Use max so that points get a value of 1
    final int triangles = Math.max(Shape3DHelper.getNumberOfTriangles(rendering), 1);

    final GeometryArray ga = (GeometryArray) shape.getGeometry();

    final long size = (long) points.size() * triangles;
    if (size > 10000000L) {
      final String name = (rendering == Rendering.POINT) ? "points" : "triangles";
      final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
      egd.addMessage("The results will generate a large dataset of " + size + " " + name
          + ".\nThis may take a long time to render and may run out of memory.");
      egd.setOKLabel("Continue");
      egd.showDialog();
      if (egd.wasCanceled()) {
        return null;
      }
    }

    final Appearance appearance = shape.getAppearance();
    final TransparencyAttributes ta = new TransparencyAttributes();
    ta.setTransparency(transparency);
    ta.setTransparencyMode(
        (transparency == 0) ? TransparencyAttributes.NONE : TransparencyAttributes.FASTEST);
    appearance.setTransparencyAttributes(ta);
    if (rendering == Rendering.POINT) {
      appearance.getPointAttributes().setPointSize(sphereSize[0].x);
    }
    if (settings.getSupportDynamicTransparency()) {
      return new ItemGeometryGroup(points.toArray(new Point3f[points.size()]), ga, appearance,
          sphereSize, colors, alpha);
    }
    return new OrderedItemGeometryGroup(points.toArray(new Point3f[points.size()]), ga, appearance,
        sphereSize, colors, alpha);
  }

  @SuppressWarnings("unused")
  private static Shape3D createShape(Builder settings) {
    final TurboList<Point3f> points = new TurboList<>(1);
    points.addf(new Point3f());

    // We try and match the geometry and appearance of the standard mesh.
    // Do this by creating a mesh with a single point and get the Geometry and Appearance.
    GeometryArray ga;
    CustomMesh mesh;

    final float transparency = getTransparency(settings);

    // Support drawing as points ...
    if (settings.getRendering() == 0) {
      mesh = new TransparentItemPointMesh(points, null, transparency);
      ((ItemPointMesh) mesh).setPointSize((float) settings.getPixelSize());

      updateAppearance(mesh, settings);

      // Assume the TransparentItemPointMesh sets COLOR_4
      ga = (GeometryArray) mesh.getGeometry();
    } else {
      final Rendering r = Rendering.forNumber(settings.getRendering());

      final List<Point3f> point = Shape3DHelper.createLocalisationObject(r);
      final Point3f[] vertices = point.toArray(new Point3f[1]);

      // Correct the direction
      ItemTriangleMesh.checkFacets(vertices);

      final double creaseAngle = (r.isHighResolution()) ? 44 : 0;
      mesh = new ItemTriangleMesh(vertices, points.toArray(new Point3f[1]), null, null,
          transparency, creaseAngle, null);

      updateAppearance(mesh, settings);

      final int nVertices = vertices.length;
      ga = new TriangleArray(nVertices, GeometryArray.COORDINATES | GeometryArray.NORMALS);

      // Copy the coords and normals. We don't require the vertex colours.
      final float[] coords = new float[nVertices * 3];
      final float[] normals = new float[nVertices * 3];
      final GeometryArray gaToCopy = (GeometryArray) mesh.getGeometry();
      gaToCopy.getCoordinates(0, coords);
      gaToCopy.getNormals(0, normals);
      ga.setCoordinates(0, coords);
      ga.setNormals(0, normals);
      ga.setValidVertexCount(nVertices);
    }

    return new Shape3D(ga, mesh.getAppearance());
  }

  private static ItemMesh createItemMesh(final ImageJ3DResultsViewerSettingsOrBuilder settings,
      TurboList<Point3f> points, final Point3f[] sphereSize, float transparency, float[] alpha) {
    final Rendering rendering = Rendering.forNumber(settings.getRendering());
    final int colorDepth = (alpha != null) ? 4 : 3;
    final Shape3D shape = Shape3DHelper.createShape(rendering, colorDepth);

    final GeometryArray ga = (GeometryArray) shape.getGeometry();
    final Appearance appearance = shape.getAppearance();

    // Estimate the largest array required for the data.
    // The mesh is created by reference using an array for coords, normals and colors.

    final int singlePointVertexSize = ga.getValidVertexCount();
    int singlePointIndexSize = 0;

    final int stride = Math.max(3, colorDepth);
    if (ga instanceof IndexedGeometryArray) {
      // Indexed arrays may have much larger index array than the vertex array
      singlePointIndexSize = ((IndexedGeometryArray) ga).getIndexCount();
    }

    final int singlePointSize = Math.max(singlePointIndexSize, stride * singlePointVertexSize);

    final long arraySize = (long) points.size() * singlePointSize;
    if (arraySize > CustomContentHelper.MAX_ARRAY_SIZE) {
      final double capacity = (double) arraySize / CustomContentHelper.MAX_ARRAY_SIZE;
      //@formatter:off
      IJ.error(TITLE,
          TextUtils.wrap(String.format(
              "The results will generate data of %d values. " +
              "This is amount of data is not supported (%.2fx capacity). " +
              "Please choose a different dataset with fewer points or " +
              "different rendering model.",
              arraySize, capacity), 80));
      //@formatter:on
      return null;
    }

    // Support drawing as points ...
    if (settings.getRendering() == 0) {
      final ItemMesh mesh = new ReferenceItemMesh(points.toArray(new Point3f[points.size()]), ga,
          appearance, null, null, transparency);
      if (alpha != null) {
        mesh.setItemAlpha(alpha);
      }
      mesh.getAppearance().getPointAttributes().setPointSize(sphereSize[0].x);
      return mesh;
    }

    final int triangles = Shape3DHelper.getNumberOfTriangles(rendering);
    final long size = (long) points.size() * triangles;
    if (size > 10000000L) {
      final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
      egd.addMessage("The results will generate a large mesh of " + size
          + " triangles.\nThis may take a long time to render and may run out of memory.");
      egd.setOKLabel("Continue");
      egd.showDialog();
      if (egd.wasCanceled()) {
        return null;
      }
    }

    IJ.showStatus("Creating 3D mesh ...");
    final ItemMesh mesh = new ReferenceItemMesh(points.toArray(new Point3f[points.size()]), ga,
        appearance, sphereSize, null, transparency);
    if (alpha != null) {
      mesh.setItemAlpha(alpha);
    }
    return mesh;
  }

  @SuppressWarnings("unused")
  private static CustomMesh createMesh(final ImageJ3DResultsViewerSettingsOrBuilder settings,
      TurboList<Point3f> points, final Point3f[] sphereSize, float transparency, float[] alpha) {
    int stride = 3 + 3; // Coordinates + color
    if (alpha != null) {
      stride++; // add color alpha
    }

    // Support drawing as points ...
    if (settings.getRendering() == 0) {
      final long arraySize = (long) points.size() * stride;
      if (arraySize > CustomContentHelper.MAX_ARRAY_SIZE) {
        final double capacity = (double) arraySize / CustomContentHelper.MAX_ARRAY_SIZE;
        //@formatter:off
        IJ.error(TITLE,
            TextUtils.wrap(String.format(
                "The results will generate data of %d values. " +
                "This is amount of data is not supported (%.2fx capacity). " +
                "Please choose a different dataset with fewer points.",
                arraySize, capacity), 80));
        //@formatter:on
        return null;
      }

      CustomPointMesh mesh;
      if (alpha != null) {
        final TransparentItemPointMesh mesh2 =
            new TransparentItemPointMesh(points, null, transparency);
        mesh = mesh2;
        mesh2.setItemAlpha(alpha);
      } else {
        mesh = new ItemPointMesh(points, null, transparency);
      }
      mesh.setPointSize(sphereSize[0].x);
      return mesh;
    }

    final Rendering r = Rendering.forNumber(settings.getRendering());

    // Repeated mesh creation is much faster as the normals are cached.
    // There does not appear to be a difference in the speed the image responds
    // to user interaction between indexed or standard triangles.

    // Currently the RepeatedIndexedTriangleMesh computes the normals a different way to
    // the super class to preserve the orientation of the normals. So if the coordinates
    // are modified through the mesh then the appearance will change. For now just use
    // the RepeatedTriangleMesh.

    // TODO - check this. It may not be true if the shading mode is flat...

    // Also the IndexedTriangleMesh has one normal per vertex and this causes a colour fall-off
    // on the triangle plane towards the edges. The TriangleMesh colours the entire surface
    // of each triangle the same which looks 'normal'.

    final List<Point3f> point = Shape3DHelper.createLocalisationObject(r);

    stride += 3; // + normals

    final int singlePointSize = point.size();
    final long size = (long) points.size() * singlePointSize;
    final long arraySize = size * stride;
    if (arraySize > CustomContentHelper.MAX_ARRAY_SIZE) {
      final double capacity = (double) arraySize / CustomContentHelper.MAX_ARRAY_SIZE;
      //@formatter:off
      IJ.error(TITLE,
          TextUtils.wrap(String.format(
              "The results will generate data of %d values. " +
              "This is amount of data is not supported (%.2fx capacity). " +
              "Please choose a different rendering model with fewer vertices.",
              arraySize, capacity), 80));
      //@formatter:on
      return null;
    }
    if (size > 10000000L) {
      final ExtendedGenericDialog egd = new ExtendedGenericDialog(TITLE);
      egd.addMessage("The results will generate a large mesh of " + size
          + " vertices.\nThis may take a long time to render and may run out of memory.");
      egd.setOKLabel("Continue");
      egd.showDialog();
      if (egd.wasCanceled()) {
        return null;
      }
    }

    IJ.showStatus("Creating 3D mesh ...");
    final double creaseAngle = (r.isHighResolution()) ? 44 : 0;
    final ImageJTrackProgress progress = null; // Used for debugging construction time
    if (alpha != null) {
      final TransparentItemTriangleMesh mesh = new TransparentItemTriangleMesh(
          point.toArray(new Point3f[singlePointSize]), points.toArray(new Point3f[points.size()]),
          sphereSize, null, transparency, creaseAngle, progress);
      mesh.setItemAlpha(alpha);
      return mesh;
    }

    return new ItemTriangleMesh(point.toArray(new Point3f[singlePointSize]),
        points.toArray(new Point3f[points.size()]), sphereSize, null, transparency, creaseAngle,
        progress);
  }
}

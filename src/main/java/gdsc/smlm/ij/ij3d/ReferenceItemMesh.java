package gdsc.smlm.ij.ij3d;

import java.util.Arrays;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.GeometryStripArray;
import org.scijava.java3d.GeometryUpdater;
import org.scijava.java3d.IndexedGeometryArray;
import org.scijava.java3d.IndexedGeometryStripArray;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Color4f;
import org.scijava.vecmath.Point3f;

/**
 * Use a mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 * <p>
 * Note: TransparentItemShape is only supported if {@link #hasColor4()} is true, i.e. the input geometry has per vertex
 * colours with alpha.
 * <p>
 * The created geometry array is created by reference.
 */
public class ReferenceItemMesh extends ItemMesh
{
	/**
	 * Instantiates a new item mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 * <p>
	 * The input geometry array vertex format is checked and unsupported formats throw an exception. Currently this only
	 * supports coordinates, normals and color. Strip and Fan arrays are supported.
	 *
	 * @param point
	 *            the point
	 * @param ga
	 *            the geometry array. If null then a default will be used. Assumed to be centred on the origin.
	 * @param appearance
	 *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
	 * @param sizes
	 *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
	 * @param color
	 *            the color
	 * @param transparency
	 *            the transparency
	 */
	public ReferenceItemMesh(Point3f point, GeometryArray ga, Appearance appearance, Point3f[] sizes,
			final Color3f color, final float transparency)
	{
		super(point, ga, appearance, sizes, color, transparency);
	}

	/**
	 * Instantiates a new item mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 * <p>
	 * The input geometry array vertex format is checked and unsupported formats throw an exception. Currently this only
	 * supports coordinates, normals and color. Strip and Fan arrays are supported.
	 *
	 * @param points
	 *            the points
	 * @param ga
	 *            the geometry array. If null then a default will be used. Assumed to be centred on the origin.
	 * @param appearance
	 *            the default appearance of the shape. PolygonAttributes, Material and TransparencyAttributes are used.
	 * @param sizes
	 *            the sizes of each point. Can be null (no scaling); length=1 (fixed scaling); or points.length.
	 * @param color
	 *            the color
	 * @param transparency
	 *            the transparency
	 */
	public ReferenceItemMesh(Point3f[] points, GeometryArray ga, Appearance appearance, Point3f[] sizes,
			final Color3f color, final float transparency)
	{
		super(points, ga, appearance, sizes, color, transparency);
	}

	protected GeometryArray createGeometry(float[] coords, GeometryArray sourceGA)
	{
		final GeometryArray ga = createGeometryArray(sourceGA, GeometryArray.BY_REFERENCE);

		ga.setCoordRefFloat(coords);
		ga.setCapability(GeometryArray.ALLOW_REF_DATA_READ);
		ga.setCapability(GeometryArray.ALLOW_REF_DATA_WRITE);

		// Fan are extensions of GeometryStripArray so do not need extra code.

		// Handle indexed array
		if (isIndexGeometryArray())
		{
			IndexedGeometryArray sourceIGA = (IndexedGeometryArray) sourceGA;
			IndexedGeometryArray iga = (IndexedGeometryArray) ga;
			int objectIndexCount = sourceIGA.getValidIndexCount();
			int[] objectIndices = new int[objectIndexCount];
			int[] allIndices = new int[objectIndices.length * points.length];
			sourceIGA.getCoordinateIndices(0, objectIndices);
			duplicateIndices(objectIndices, allIndices);
			iga.setCoordinateIndices(0, allIndices);

			// Check if we need the color and normal indices 
			if ((vertexFormat & GeometryArray.USE_COORD_INDEX_ONLY) == 0)
			{
				// Duplicate the other indices
				if (hasNormals())
				{
					sourceIGA.getNormalIndices(0, objectIndices);
					duplicateIndices(objectIndices, allIndices);
					iga.setNormalIndices(0, allIndices);
				}

				if (hasColor())
				{
					sourceIGA.getColorIndices(0, objectIndices);
					duplicateIndices(objectIndices, allIndices);
					iga.setColorIndices(0, allIndices);
				}
			}
		}
		
		// Handle normals
		if (hasNormals())
		{
			float[] objectNormals = new float[vertexCount * 3];
			sourceGA.getNormals(0, objectNormals);
			float[] allNormals = new float[objectNormals.length * points.length];
			duplicate(objectNormals, 0, objectNormals.length, points.length, allNormals, 0);
			ga.setNormalRefFloat(allNormals);
		}
		
		// Handle colors
		if (hasColor())
		{
			colorUpdater = ArrayColorUpdater.create(vertexCount, hasColor4());
			ga.setColorRefFloat(new float[colorUpdater.size() * size()]);
		}		

		return ga;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.UpdatedableItemMesh#reorderFast(int[])
	 */
	public void reorderFast(int[] indices) throws IllegalArgumentException
	{
		changed = true;

		final int oldSize = size();
		final int size = (indices == null) ? 0 : Math.min(oldSize, indices.length);

		if (size == 0)
		{
			points = new Point3f[0];
			sizes = new Point3f[0];
			this.setGeometry(null);
			return;
		}

		// From here on we assume the current geometry will not be null
		// as this only happens when the original size is zero. Size has 
		// been checked at this point to be the smaller of new and old. 
		GeometryArray ga = (GeometryArray) getGeometry();

		points = reorder(points, indices);
		// Sizes could be null or a single size
		if (sizes != null && sizes.length == points.length)
			sizes = reorder(sizes, indices);

		// Reorder all things in the geometry: coordinates and colour.
		// The normals, indices, strip counts are are unchanged.
		//int objectSize = vertexCount;

		int n = vertexCount * 3;
		final float[] oldCoords = ga.getCoordRefFloat();
		final float[] coords = new float[size * n];
		for (int i = 0; i < size; i++)
		{
			int j = indices[i];
			int ii = i * n;
			int jj = j * n;
			System.arraycopy(oldCoords, jj, coords, ii, n);
		}

		final float[] colors;
		if (hasColor())
		{
			n = colorUpdater.size();
			float[] oldColors = ga.getColorRefFloat();
			colors = new float[size * n];
			for (int i = 0; i < size; i++)
			{
				int j = indices[i];
				int ii = i * n;
				int jj = j * n;
				System.arraycopy(oldColors, jj, colors, ii, n);
			}
		}
		else
		{
			colors = null;
		}

		ga.updateData(new GeometryUpdater()
		{
			public void updateData(Geometry geometry)
			{
				GeometryArray ga = (GeometryArray) geometry;
				// We re-use the geometry and just truncate the vertex count
				ga.setCoordRefFloat(coords);
				if (colors != null)
					ga.setColorRefFloat(colors);

				if (size != oldSize)
				{
					if (isIndexGeometryArray())
					{
						if (isStripGeometryArray())
						{
							int[] indices = new int[indexCount * oldSize];
							((IndexedGeometryStripArray) ga).getStripIndexCounts(indices);
							indices = Arrays.copyOf(indices, indexCount * size);
							((IndexedGeometryStripArray) ga).setStripIndexCounts(indices);
						}
						else
						{
							((IndexedGeometryArray) ga).setValidIndexCount(size * indexCount);
						}
					}
					else if (isStripGeometryArray())
					{
						int[] indices = new int[vertexCount * oldSize];
						((GeometryStripArray) ga).getStripVertexCounts(indices);
						indices = Arrays.copyOf(indices, vertexCount * size);
						((GeometryStripArray) ga).setStripVertexCounts(indices);
					}
					else
					{
						ga.setValidVertexCount(size * vertexCount);
					}
				}
			}
		});
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemShape#setItemColor(org.scijava.vecmath.Color3f)
	 */
	public void setItemColor(Color3f color)
	{
		if (color == null)
			color = DEFAULT_COLOR;
		this.color = color;
		if (!hasColor())
		{
			if (isColorByMaterial)
				getAppearance().getMaterial().setDiffuseColor(color);
			else
				getAppearance().getColoringAttributes().setColor(color);
			return;
		}
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		final float[] colors = ga.getColorRefFloat();
		if (hasColor3())
		{
			float[] tmp = new float[3];
			color.get(tmp);
			duplicate(tmp, 0, 3, colors.length / 3, colors, 0);
		}
		else
		{
			// Preserve alpha
			for (int i = 0; i < colors.length; i += 4)
			{
				colors[i] = color.x;
				colors[i + 1] = color.y;
				colors[i + 2] = color.z;
			}
		}
		ga.setColorRefFloat(colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.ItemMesh#setItemColor(org.scijava.vecmath.Color3f[])
	 */
	public void setItemColor(Color3f[] color) throws IllegalArgumentException
	{
		if (!hasColor())
		{
			setItemColor(color[0]);
			return;
		}
		this.color = null;
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int n = colorUpdater.size();
		final float[] colors = ga.getColorRefFloat();
		if (hasColor3())
		{
			for (int i = 0; i < color.length; i++)
			{
				System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
			}
		}
		else
		{
			// Preserve alpha
			for (int i = 0; i < color.length; i++)
			{
				int offset = i * n;
				colorUpdater.getColors(color[i], colors[offset + 3]);
				System.arraycopy(colorUpdater.pointColor, 0, colors, offset, n);
			}
		}
		ga.setColorRefFloat(colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemColor4(org.scijava.vecmath.Color4f[])
	 */
	public void setItemColor4(Color4f[] color) throws IllegalArgumentException
	{
		if (!hasColor4())
			throw new IllegalArgumentException("Per-item alpha not supported");

		this.color = null;
		int size = size();
		if (color.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int n = colorUpdater.size();
		final float[] colors = ga.getColorRefFloat();
		for (int i = 0; i < color.length; i++)
		{
			System.arraycopy(colorUpdater.getColors(color[i]), 0, colors, i * n, n);
		}
		ga.setColorRefFloat(colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemAlpha(float[])
	 */
	public void setItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		if (!hasColor4())
			throw new IllegalArgumentException("Per-item alpha not supported");

		int size = size();
		if (alpha.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int n = colorUpdater.size();
		// Preserve color
		final float[] colors = ga.getColorRefFloat();
		for (int i = 0; i < size; i++)
		{
			int offset = i * n;
			for (int j = 3; j < n; j += 4)
				colors[j + offset] = alpha[i];
		}
		ga.setColorRefFloat(colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#setItemAlpha(float)
	 */
	public void setItemAlpha(float alpha) throws IllegalArgumentException
	{
		if (!hasColor4())
			throw new IllegalArgumentException("Per-item alpha not supported");

		int size = size();
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int n = colorUpdater.size();
		// Preserve color
		final float[] colors = ga.getColorRefFloat();
		for (int i = 0; i < size; i++)
		{
			int offset = i * n;
			for (int j = 3; j < n; j += 4)
				colors[j + offset] = alpha;
		}
		ga.setColorRefFloat(colors);
		changed = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see gdsc.smlm.ij.ij3d.TransparentItemShape#getItemAlpha(float[])
	 */
	public void getItemAlpha(float[] alpha) throws IllegalArgumentException
	{
		if (!hasColor4())
			throw new IllegalArgumentException("Per-item alpha not supported");

		int size = size();
		if (alpha.length != size)
			throw new IllegalArgumentException("list of size " + size + " expected");
		final GeometryArray ga = (GeometryArray) getGeometry();
		if (ga == null)
			return;
		int n = colorUpdater.size();
		final float[] colors = ga.getColorRefFloat();
		for (int i = 0; i < size; i++)
		{
			// Get only alpha
			alpha[i] = colors[i * n + 3];
		}
	}
}

package gdsc.smlm.ij.ij3d;

import java.util.Arrays;
import java.util.List;

import org.scijava.java3d.Appearance;
import org.scijava.java3d.Geometry;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.TriangleArray;
import org.scijava.java3d.utils.geometry.GeometryInfo;
import org.scijava.java3d.utils.geometry.NormalGenerator;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import customnode.CustomTriangleMesh;
import gdsc.core.logging.NullTrackProgress;
import gdsc.core.logging.Ticker;
import gdsc.core.logging.TrackProgress;

/**
 * Use a triangle mesh object to represent a set of points. The object is duplicated, scaled and translated for
 * each point.
 */
public class RepeatedTriangleMesh extends CustomTriangleMesh
{
	protected Point3f[] objectVertices;
	protected Vector3f[] objectNormals;
	protected Point3f[] points;
	protected Point3f[] sizes;
	private boolean dirty = false;

	/**
	 * Instantiates a new repeated indexed triangle mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 *
	 * @param objectVertices
	 *            the vertices of the object for a single point
	 * @param points
	 *            the points
	 * @param sizes
	 *            the size of each point
	 * @param color
	 *            the color
	 * @param transp
	 *            the transparency
	 */
	public RepeatedTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes, Color3f color,
			float transp)
	{
		this(objectVertices, points, sizes, color, transp, -1, NullTrackProgress.INSTANCE);
	}

	/**
	 * Instantiates a new repeated indexed triangle mesh.
	 * <p>
	 * This will repeat the object for each input point. The object
	 * is assumed to be centred on the origin. It will be scaled and
	 * translated for each input point.
	 * <p>
	 * The crease angle is used to collapse facets normals at a vertex into a single normal for smoothing shading. Set
	 * to 0 to draw the polygon with no shading.
	 *
	 * @param objectVertices
	 *            the vertices of the object for a single point
	 * @param points
	 *            the points
	 * @param sizes
	 *            the size of each point
	 * @param color
	 *            the color
	 * @param transp
	 *            the transparency
	 * @param creaseAngle
	 *            the crease angle (in degrees). Set to negative to ignore. The default is 44.
	 * @param progress
	 *            the progress
	 */
	public RepeatedTriangleMesh(Point3f[] objectVertices, Point3f[] points, Point3f[] sizes, Color3f color,
			float transp, double creaseAngle, TrackProgress progress)
	{
		// Create empty 
		super(null, color, transp);

		if (sizes != null && points.length != sizes.length)
			throw new IllegalArgumentException("Points and sizes must be the same length");

		progress = NullTrackProgress.createIfNull(progress);

		if (progress.isStatus())
			progress.status("Standardising object vertices");
		this.objectVertices = objectVertices;

		checkFacets(objectVertices);

		objectNormals = getNormals(objectVertices, creaseAngle);

		// Now build the actual vertices by repeating the points.
		this.points = points;
		this.sizes = sizes;

		Point3f[] vertices = new Point3f[objectVertices.length * points.length];

		final int n = objectVertices.length;
		if (progress.isStatus())
			progress.status("Computing vertices");
		Ticker ticker = Ticker.createStarted(progress, n, false);
		boolean sameSize = false;
		if (sizes == null || (sameSize = sameSize(sizes)))
		{
			if (sameSize)
			{
				// Scale the input object
				final Point3f s = sizes[0];
				final float sx = s.x;
				final float sy = s.y;
				final float sz = s.z;
				for (int j = 0; j < n; j++)
				{
					objectVertices[j].x *= sx;
					objectVertices[j].y *= sy;
					objectVertices[j].z *= sz;
				}
			}

			// Translate
			for (int i = 0, k = 0; i < points.length; i++)
			{
				final Point3f p = points[i];
				final float dx = p.x;
				final float dy = p.y;
				final float dz = p.z;
				for (int j = 0; j < n; j++)
				{
					vertices[k++] = new Point3f(objectVertices[j].x + dx, objectVertices[j].y + dy,
							objectVertices[j].z + dz);
				}
				ticker.tick();
			}
		}
		else
		{
			// Translate and scale
			for (int i = 0, k = 0; i < points.length; i++)
			{
				final Point3f p = points[i];
				final float dx = p.x;
				final float dy = p.y;
				final float dz = p.z;
				final Point3f s = sizes[i];
				final float sx = s.x;
				final float sy = s.y;
				final float sz = s.z;
				for (int j = 0; j < n; j++)
				{
					vertices[k++] = new Point3f(objectVertices[j].x * sx + dx, objectVertices[j].y * sy + dy,
							objectVertices[j].z * sz + dz);
				}
				ticker.tick();
			}
		}

		ticker.stop();

		this.mesh = Arrays.asList(vertices);

		// Update the geometry
		if (progress.isStatus())
			progress.status("Creating geometry");
		this.setGeometry(createGeometry());
	}

	/**
	 * Check if all the points are the same size.
	 *
	 * @param sizes
	 *            the sizes (must be an array of at least 1)
	 * @return true, if successful
	 */
	public static boolean sameSize(Point3f[] sizes)
	{
		if (sizes.length == 1)
			return true;
		final Point3f s = sizes[0];
		for (int j = 1; j < sizes.length; j++)
		{
			if (!sizes[j].equals(s))
				return false;
		}
		return true;
	}

	@Override
	public void setCoordinate(final int i, final Point3f p)
	{
		dirty = true;
		super.setCoordinate(i, p);
	}

	@Override
	public void setCoordinates(final int[] indices, final Point3f p)
	{
		dirty = true;
		super.setCoordinates(indices, p);
	}

	@Override
	public void addTriangle(Point3f p1, Point3f p2, Point3f p3)
	{
		dirty = true;
		super.addTriangle(p1, p2, p3);
	}

	@Override
	public void addTriangles(Point3f[] v)
	{
		dirty = true;
		super.addTriangles(v);
	}

	@Override
	public void removeTriangle(int index)
	{
		dirty = true;
		super.removeTriangle(index);
	}

	@Override
	public void removeTriangles(int[] indices)
	{
		dirty = true;
		super.removeTriangles(indices);
	}

	@Override
	protected void addVertices(Point3f[] v)
	{
		dirty = true;
		super.addVertices(v);
	}

	@Override
	protected void addVerticesToGeometryArray(Point3f[] v)
	{
		dirty = true;
		super.addVerticesToGeometryArray(v);
	}

	@Override
	protected void addVerticesToGeometryStripArray(Point3f[] v)
	{
		dirty = true;
		super.addVerticesToGeometryStripArray(v);
	}

	@Override
	protected void removeVertices(int[] indices)
	{
		dirty = true;
		super.removeVertices(indices);
	}

	@Override
	protected GeometryArray createGeometry()
	{
		if (mesh == null || mesh.size() < 3)
			return null;
		final List<Point3f> tri = mesh;
		final int nValid = tri.size();
		final int nAll = 2 * nValid;

		final Point3f[] coords = new Point3f[nValid];
		tri.toArray(coords);

		final Color3f colors[] = new Color3f[nValid];
		if (null == color)
		{
			// Vertex-wise colors are not stored
			// so they have to be retrieved from the geometry:
			final GeometryArray gaOld = (GeometryArray) getGeometry();
			if (null != gaOld)
			{
				gaOld.getColors(0, colors);
			}
			else
			{
				Arrays.fill(colors, new Color3f(DEFAULT_COLOR));
			}
		}
		else
		{
			Arrays.fill(colors, color);
		}

		final GeometryArray ta = new TriangleArray(nAll,
				GeometryArray.COORDINATES | GeometryArray.COLOR_3 | GeometryArray.NORMALS);

		ta.setCoordinates(0, coords);
		ta.setColors(0, colors);

		// generate normals
		final GeometryArray result;

		if (dirty)
		{
			final GeometryInfo gi = new GeometryInfo(ta);
			final NormalGenerator ng = new NormalGenerator();
			ng.generateNormals(gi);
			//			for (int i = 0; i < nValid; i++)
			//			{
			//				System.out.printf("[%d] %s : %s  %b\n", i, gi.getNormals()[gi.getNormalIndices()[i]], objectNormals[i],
			//						gi.getNormals()[gi.getNormalIndices()[i]].equals(objectNormals[i]));
			//			}
			result = gi.getGeometryArray();
		}
		else
		{
			// Use the same normals for each repeated object
			final Vector3f[] normals = new Vector3f[nValid];

			// Binary fill
			int fill = objectNormals.length;
			System.arraycopy(objectNormals, 0, normals, 0, fill);
			for (int i = 2; i < points.length; i *= 2)
			{
				System.arraycopy(normals, 0, normals, fill, fill);
				fill *= 2;
			}
			// Final fill
			System.arraycopy(normals, 0, normals, fill, normals.length - fill);

			//			final Vector3f[] normals2 = normals.clone();
			//			for (int i = 0, k = 0; i < points.length; i++)
			//			{
			//				for (int j = 0; j < objectNormals.length; j++)
			//				{
			//					normals[k++] = objectNormals[j];
			//				}
			//			}
			//			assert Arrays.equals(normals, normals2) : "Not the same fill";

			ta.setNormals(0, normals);

			//final GeometryInfo gi = new GeometryInfo(ta);
			//result = gi.getGeometryArray();
			result = ta;
		}

		result.setCapability(GeometryArray.ALLOW_NORMAL_WRITE);
		result.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
		result.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
		result.setCapability(GeometryArray.ALLOW_COUNT_WRITE);
		result.setCapability(GeometryArray.ALLOW_COUNT_READ);
		result.setCapability(GeometryArray.ALLOW_FORMAT_READ);
		result.setCapability(Geometry.ALLOW_INTERSECT);
		result.setValidVertexCount(nValid);

		return result;
	}

	@Override
	protected Appearance createAppearance()
	{
		return super.createAppearance();

		//		final Appearance appearance = new Appearance();
		//		appearance.setCapability(Appearance.ALLOW_TRANSPARENCY_ATTRIBUTES_READ);
		//
		//		final PolygonAttributes polyAttrib = new PolygonAttributes();
		//		polyAttrib.setCapability(PolygonAttributes.ALLOW_MODE_WRITE);
		//		if (this.shaded)
		//			polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_FILL);
		//		else
		//			polyAttrib.setPolygonMode(PolygonAttributes.POLYGON_LINE);
		//		polyAttrib.setCullFace(PolygonAttributes.CULL_NONE);
		//
		//		// This is what makes the polygons look the same on both sides!
		//		//polyAttrib.setBackFaceNormalFlip(true);
		//
		//		appearance.setPolygonAttributes(polyAttrib);
		//
		//		final ColoringAttributes colorAttrib = new ColoringAttributes();
		//		colorAttrib.setShadeModel(ColoringAttributes.SHADE_GOURAUD);
		//		if (null != color) // is null when colors are vertex-wise
		//			colorAttrib.setColor(color);
		//		appearance.setColoringAttributes(colorAttrib);
		//
		//		final TransparencyAttributes tr = new TransparencyAttributes();
		//		final int mode = TransparencyAttributes.FASTEST;
		//		tr.setCapability(TransparencyAttributes.ALLOW_VALUE_WRITE);
		//		tr.setCapability(TransparencyAttributes.ALLOW_MODE_WRITE);
		//		tr.setTransparencyMode(mode);
		//		tr.setTransparency(transparency);
		//		appearance.setTransparencyAttributes(tr);
		//
		//		final Material material = new Material();
		//		material.setCapability(Material.ALLOW_COMPONENT_WRITE);
		//		material.setAmbientColor(0.1f, 0.1f, 0.1f);
		//		material.setSpecularColor(0.1f, 0.1f, 0.1f);
		//		material.setDiffuseColor(0.1f, 0.1f, 0.1f);
		//		appearance.setMaterial(material);
		//		return appearance;
	}

	/**
	 * Check the facet normals point out from the centre 0,0,0. If the normal points inwards then the vertices will be
	 * swapped.
	 *
	 * @param vertices
	 *            the vertices
	 * @return the count of the number swapped
	 */
	public static int checkFacets(Point3f[] vertices)
	{
		int count = 0;
		int nVertices = vertices.length;
		final Vector3f v1 = new Vector3f(), v2 = new Vector3f();
		for (int i = 0; i < nVertices; i += 3)
		{
			// Use the same order as that used to compute facet normals in 
			// org.scijava.java3d.utils.geometry.NormalGenerator
			v1.sub(vertices[i + 2], vertices[i + 1]);
			v2.sub(vertices[i], vertices[i + 1]);
			v1.cross(v1, v2);
			v1.normalize();

			// Project point (x,y,z) to plane with normal (a,b,c) and point (d,e,f)
			// t = (ad - ax + be - by + cd - cz) / (a^2 + b^2 + c^2)
			// projected point = (x+ta,y+tb,z+tc)

			// Project 0,0,0 to the facet
			double a = v1.x;
			double b = v1.y;
			double c = v1.z;
			double d = vertices[i].x;
			double e = vertices[i].y;
			double f = vertices[i].z;
			double t = a * d + b * e + c * f;
			if (t < 0)
			{
				count++;
				swap(vertices, i + 2, i);
			}
		}
		//System.out.printf("Swapped %d\n", count);
		return count;
	}

	private static void swap(Point3f[] vertices, int i, int j)
	{
		Point3f tmp = vertices[i];
		vertices[i] = vertices[j];
		vertices[j] = tmp;
	}

	/**
	 * Gets the normals assuming triangle vertices.
	 *
	 * @param vertices
	 *            the vertices
	 * @param creaseAngle
	 *            the crease angle (in degrees)
	 * @return the normals
	 */
	public static Vector3f[] getNormals(Point3f[] vertices, double creaseAngle)
	{
		int nVertices = vertices.length;
		Vector3f[] normals = new Vector3f[nVertices];

		final GeometryArray ta = new TriangleArray(nVertices, GeometryArray.COORDINATES | GeometryArray.NORMALS);
		ta.setCoordinates(0, vertices);
		final GeometryInfo gi = new GeometryInfo(ta);
		final NormalGenerator ng = new NormalGenerator();
		if (creaseAngle >= 0 && creaseAngle <= 180)
			ng.setCreaseAngle(creaseAngle * Math.PI / 180.0);
		ng.generateNormals(gi);
		Vector3f[] n = gi.getNormals();
		int[] indices = gi.getNormalIndices();
		for (int i = 0; i < nVertices; i++)
		{
			normals[i] = n[indices[i]];
		}

		return normals;
	}
}

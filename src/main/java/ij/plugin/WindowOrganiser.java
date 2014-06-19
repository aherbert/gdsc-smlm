package ij.plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Dimension;

/**
 * Extend the standard ImageJ window organiser plugin and make the methods public
 */
public class WindowOrganiser extends ij.plugin.WindowOrganizer
{
	private static final int XSTART = 4, YSTART = 80, XOFFSET = 8, YOFFSET = 24, MAXSTEP = 200, GAP = 2;
	private int titlebarHeight = IJ.isMacintosh() ? 40 : 20;

	@Override
	public void tileWindows(int[] wList)
	{
		try
		{
			super.tileWindows(wList);
		}
		catch (Throwable t)
		{
			// Some versions of ImageJ / Java do not allow this so call the duplicated function
			copyOfTileWindows(wList);
		}
	}

	@Override
	public void cascadeWindows(int[] wList)
	{
		try
		{
			super.cascadeWindows(wList);
		}
		catch (Exception e)
		{
			// Some versions of ImageJ / Java do not allow this so call the duplicated function
			copyOfCascadeWindows(wList);
		}
	}

	void copyOfTileWindows(int[] wList)
	{
		Dimension screen = IJ.getScreenSize();
		int minWidth = Integer.MAX_VALUE;
		int minHeight = Integer.MAX_VALUE;
		double totalWidth = 0;
		double totalHeight = 0;
		for (int i = 0; i < wList.length; i++)
		{
			ImageWindow win = getWindow(wList[i]);
			if (win == null)
				continue;
			Dimension d = win.getSize();
			int w = d.width;
			int h = d.height + titlebarHeight;
			if (w < minWidth)
				minWidth = w;
			if (h < minHeight)
				minHeight = h;
			totalWidth += w;
			totalHeight += h;
		}
		int nPics = wList.length;
		double averageWidth = totalWidth / nPics;
		double averageHeight = totalHeight / nPics;
		int tileWidth = (int) averageWidth;
		int tileHeight = (int) averageHeight;
		//IJ.write("tileWidth, tileHeight: "+tileWidth+" "+tileHeight);
		int hspace = screen.width - 2 * GAP;
		if (tileWidth > hspace)
			tileWidth = hspace;
		int vspace = screen.height - YSTART;
		if (tileHeight > vspace)
			tileHeight = vspace;
		int hloc, vloc;
		boolean theyFit;
		do
		{
			hloc = XSTART;
			vloc = YSTART;
			theyFit = true;
			int i = 0;
			do
			{
				i++;
				if (hloc + tileWidth > screen.width)
				{
					hloc = XSTART;
					vloc = vloc + tileHeight;
					if (vloc + tileHeight > screen.height)
						theyFit = false;
				}
				hloc = hloc + tileWidth + GAP;
			} while (theyFit && (i < nPics));
			if (!theyFit)
			{
				tileWidth = (int) (tileWidth * 0.98 + 0.5);
				tileHeight = (int) (tileHeight * 0.98 + 0.5);
			}
		} while (!theyFit);
		hloc = XSTART;
		vloc = YSTART;

		for (int i = 0; i < nPics; i++)
		{
			if (hloc + tileWidth > screen.width)
			{
				hloc = XSTART;
				vloc = vloc + tileHeight;
			}
			ImageWindow win = getWindow(wList[i]);
			if (win != null)
			{
				win.setLocation(hloc, vloc);
				//IJ.write(i+" "+w+" "+tileWidth+" "+mag+" "+IJ.d2s(zoomFactor,2)+" "+zoomCount);
				ImageCanvas canvas = win.getCanvas();
				while (win.getSize().width * 0.85 >= tileWidth && canvas.getMagnification() > 0.03125)
					canvas.zoomOut(0, 0);
				win.toFront();
			}
			hloc += tileWidth + GAP;
		}
	}

	ImageWindow getWindow(int id)
	{
		ImageWindow win = null;
		ImagePlus imp = WindowManager.getImage(id);
		if (imp != null)
			win = imp.getWindow();
		return win;
	}

	void copyOfCascadeWindows(int[] wList)
	{
		Dimension screen = IJ.getScreenSize();
		int x = XSTART;
		int y = YSTART;
		int xstep = 0;
		int xstart = XSTART;
		for (int i = 0; i < wList.length; i++)
		{
			ImageWindow win = getWindow(wList[i]);
			if (win == null)
				continue;
			Dimension d = win.getSize();
			if (i == 0)
			{
				xstep = (int) (d.width * 0.8);
				if (xstep > MAXSTEP)
					xstep = MAXSTEP;
			}
			if (y + d.height * 0.67 > screen.height)
			{
				xstart += xstep;
				if (xstart + d.width * 0.67 > screen.width)
					xstart = XSTART + XOFFSET;
				x = xstart;
				y = YSTART;
			}
			win.setLocation(x, y);
			win.toFront();
			x += XOFFSET;
			y += YOFFSET;
		}
	}
}

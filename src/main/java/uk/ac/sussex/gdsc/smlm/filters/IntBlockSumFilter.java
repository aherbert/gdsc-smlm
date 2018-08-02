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
package uk.ac.sussex.gdsc.smlm.filters;

/**
 * Computes the sum of the neighbourhood for each point within the array using a rolling block algorithm.
 * <p>
 * rollingBlock algorithm uses two consecutive 1D passes using (2n+1) strips. Each pass is computed using a rolling
 * total thus each pixel sum can be computed using a single addition and subtraction of the end pixels of the strip.
 * Speed ~ Order(1).
 * <p>
 * Note: Due to lack of small dimension checking the routines will fail if maxx or maxy are less than 2. All routines
 * are OK for 3x3 images and larger.
 * <p>
 * This algorithm does not mirror edge pixels in contrast to the BlockSumFilter.
 */
public class IntBlockSumFilter extends BaseFilter
{
    private int[] buffer = null;
    private int[] intRowBuffer = null;

    /**
     * Compute the filter within a 2n+1 size block around each point.
     * Only pixels with a full block are processed. Pixels within border regions
     * are unchanged.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     * @param n
     *            The block size
     */
    public void rollingBlockFilterInternal(int[] data, final int maxx, final int maxy, final int n)
    {
        if (n == 1)
            rollingBlockFilter3x3Internal(data, maxx, maxy);
        else
            rollingBlockFilterNxNInternal(data, maxx, maxy, n);
    }

    /**
     * Compute the filter within a 2n+1 size block around each point.
     * Only pixels with a full block are processed. Pixels within border regions
     * are unchanged.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     * @param n
     *            The block size
     */
    void rollingBlockFilterNxNInternal(int[] data, final int maxx, final int maxy, final int n)
    {
        final int blockSize = 2 * n + 1;
        if (maxx < blockSize || maxy < blockSize)
            return;

        final int[] wdata = initialise(data, maxx, maxy, n, true);

        // X-direction
        for (int y = 0; y < maxy; y++)
        {
            // Initialise the rolling sum
            int sum = 0;

            int endIndex = y * maxx;
            int x = 0;
            while (x < blockSize)
            {
                sum += wdata[endIndex];
                endIndex++;
                x++;
            }

            // Rolling sum over the X-direction
            int startIndex = y * maxx;
            int centreIndex = startIndex + n;

            buffer[centreIndex] = sum;

            while (x < maxx)
            {
                centreIndex++;

                sum += wdata[endIndex] - wdata[startIndex];

                buffer[centreIndex] = sum;

                x++;
                startIndex++;
                endIndex++;
            }
        }

        // Y-direction.
        // Only sweep over the interior
        for (int x = n; x < maxx - n; x++)
        {
            // Initialise the rolling sum
            int sum = 0;

            int endIndex = x;
            int y = 0;
            while (y < blockSize)
            {
                sum += buffer[endIndex];
                endIndex += maxx;
                y++;
            }

            // Rolling sum over the Y-direction
            int startIndex = x;
            int centreIndex = startIndex + n * maxx;

            data[centreIndex] = sum;

            while (y < maxy)
            {
                centreIndex += maxx;

                sum += buffer[endIndex] - buffer[startIndex];

                data[centreIndex] = sum;

                y++;
                startIndex += maxx;
                endIndex += maxx;
            }
        }
    }

    /**
     * Compute the filter within a 3x3 size block around each point.
     * Only pixels with a full block are processed. Pixels within border regions
     * are unchanged.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     */
    void rollingBlockFilter3x3Internal(int[] data, final int maxx, final int maxy)
    {
        if (maxx < 3 || maxy < 3)
            return;

        final int[] wdata = initialise(data, maxx, maxy, 1, true);

        // X-direction
        for (int y = 0; y < maxy; y++)
        {
            // Initialise the rolling sum
            int startIndex = y * maxx;
            int centreIndex = startIndex + 1;
            int endIndex = centreIndex + 1;
            int sum = wdata[startIndex] + wdata[centreIndex] + wdata[endIndex];

            // Rolling sum over the X-direction
            buffer[centreIndex++] = sum;

            for (int x = 0; x < maxx - 3; x++)
            {
                sum += wdata[++endIndex] - wdata[startIndex++];
                buffer[centreIndex++] = sum;
            }
        }

        // Y-direction.
        // Only sweep over the interior
        for (int x = 1; x < maxx - 1; x++)
        {
            // Initialise the rolling sum
            int startIndex = x;
            int centreIndex = startIndex + maxx;
            int endIndex = centreIndex + maxx;
            int sum = buffer[startIndex] + buffer[centreIndex] + buffer[endIndex];

            // Rolling sum over the Y-direction
            data[centreIndex] = sum;

            for (int y = 0; y < maxy - 3; y++)
            {
                centreIndex += maxx;
                endIndex += maxx;
                sum += buffer[endIndex] - buffer[startIndex];
                data[centreIndex] = sum;
                startIndex += maxx;
            }
        }
    }

    /**
     * Initialise for filtering.
     *
     * @param data
     *            the data
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     * @param n
     *            The block size
     * @param internal
     *            the internal flag
     * @return the working data (which may be weighted)
     */
    private int[] initialise(int[] data, final int maxx, final int maxy, final int n, boolean internal)
    {
        final int size = data.length;
        createIntBuffer(size);
        return data;
    }

    /**
     * Creates the int buffer.
     *
     * @param size
     *            the size
     * @return the int buffer
     */
    private int[] createIntBuffer(int size)
    {
        if (buffer == null || buffer.length < size)
            buffer = new int[size];
        return buffer;
    }

    /**
     * Compute the filter within a 2n+1 size block around each point.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     * @param n
     *            The block size
     */
    public void rollingBlockFilter(int[] data, final int maxx, final int maxy, final int n)
    {
        if (n == 1)
            rollingBlockFilter3x3(data, maxx, maxy);
        else
            rollingBlockFilterNxN(data, maxx, maxy, n);
    }

    /**
     * Compute the filter within a 2n+1 size block around each point.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     * @param n
     *            The block size
     */
    void rollingBlockFilterNxN(int[] data, final int maxx, final int maxy, final int n)
    {
        final int[] wdata = initialise(data, maxx, maxy, n, false);

        // NOTE:
        // To increase speed when sweeping the arrays and allow for reusing code:
        //   buffer is XY ordinal => x * maxy + y
        //   data is YX ordinal    => y * maxx + x

        // X-direction
        int width = maxx;
        int height = maxy;
        int[] inData = wdata;
        int[] outData = buffer;

        int[] row = intRowBuffer(width + 2 * n);
        for (int y = 0; y < height; y++)
        {
            extractRow(inData, y, width, n, row);

            // Initialise rolling sum
            int sum = n * row[0] + row[n];
            int endIndex = n + 1;
            for (int i = 0; i < n; i++)
                sum += row[endIndex++];

            int centreIndex = y;
            outData[centreIndex] = sum;

            // Rolling sum over the X-direction
            for (int x = 0; x < width - 1; x++)
            {
                sum += row[endIndex++] - row[x];
                centreIndex += height;
                outData[centreIndex] = sum;
            }
        }

        // Y-direction.
        width = maxy;
        height = maxx;
        inData = buffer;
        outData = data;

        row = intRowBuffer(width + 2 * n);
        for (int y = 0; y < height; y++)
        {
            extractRow(inData, y, width, n, row);

            // Initialise rolling sum
            int sum = n * row[0] + row[n];
            int endIndex = n + 1;
            for (int i = 0; i < n; i++)
                sum += row[endIndex++];

            int centreIndex = y;
            outData[centreIndex] = sum;

            // Rolling sum over the X-direction
            for (int x = 0; x < width - 1; x++)
            {
                sum += row[endIndex++] - row[x];
                centreIndex += height;
                outData[centreIndex] = sum;
            }
        }
    }

    /**
     * Compute the filter within a 3x3 size block around each point.
     * <p>
     * Note: the input data is destructively modified
     *
     * @param data
     *            The input/output data (packed in YX order)
     * @param maxx
     *            The width of the data
     * @param maxy
     *            The height of the data
     */
    void rollingBlockFilter3x3(int[] data, final int maxx, final int maxy)
    {
        final int[] wdata = initialise(data, maxx, maxy, 1, false);

        // NOTE:
        // To increase speed when sweeping the arrays and allow for reusing code:
        //   buffer is XY ordinal => x * maxy + y
        //   data is YX ordinal    => y * maxx + x

        // X-direction
        int width = maxx;
        int height = maxy;
        int[] inData = wdata;
        int[] outData = buffer;

        int[] row = intRowBuffer(width + 2);
        for (int y = 0; y < height; y++)
        {
            extractRow1(inData, y, width, row);

            // Initialise rolling sum
            int sum = row[0] + row[1] + row[2];
            int endIndex = 3;

            int centreIndex = y;
            outData[centreIndex] = sum;

            // Rolling sum over the X-direction
            for (int x = 0; x < width - 1; x++)
            {
                sum += row[endIndex++] - row[x];
                centreIndex += height;
                outData[centreIndex] = sum;
            }
        }

        // Y-direction.
        width = maxy;
        height = maxx;
        inData = buffer;
        outData = data;

        row = intRowBuffer(width + 2);
        for (int y = 0; y < height; y++)
        {
            extractRow1(inData, y, width, row);

            // Initialise rolling sum
            int sum = row[0] + row[1] + row[2];
            int endIndex = 3;

            int centreIndex = y;
            outData[centreIndex] = sum;

            // Rolling sum over the X-direction
            for (int x = 0; x < width - 1; x++)
            {
                sum += row[endIndex++] - row[x];
                centreIndex += height;
                outData[centreIndex] = sum;
            }
        }
    }

    private int[] intRowBuffer(int size)
    {
        if (intRowBuffer == null || intRowBuffer.length < size)
            intRowBuffer = new int[size];
        return intRowBuffer;
    }

    private static void extractRow(int[] inData, int y, int width, final int n, int[] row)
    {
        final int index = y * width;

        // Pad ends
        for (int i = 0; i < n; i++)
        {
            row[i] = 0;
            row[i + n + width] = 0;
            // This will mirror the edge pixel.
            //row[i] = inData[index];
            //row[i + n + width] = inData[index + width - 1];
        }

        // Fill in data
        System.arraycopy(inData, index, row, n, width);
    }

    private static void extractRow1(int[] inData, int y, int width, int[] row)
    {
        final int index = y * width;

        // Pad ends
        row[0] = 0;
        row[1 + width] = 0;
        // This will mirror the edge pixel.
        //row[0] = inData[index];
        //row[1 + width] = inData[index + width - 1];

        // Fill in data
        System.arraycopy(inData, index, row, 1, width);
    }

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Object#clone()
     */
    @Override
    public IntBlockSumFilter clone()
    {
        final IntBlockSumFilter o = (IntBlockSumFilter) super.clone();
        o.buffer = null;
        o.intRowBuffer = null;
        return o;
    }
}

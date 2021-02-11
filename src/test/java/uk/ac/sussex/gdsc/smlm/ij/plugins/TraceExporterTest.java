/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2020 Alex Herbert
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

import java.io.File;
import java.io.IOException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.Matrix;
import us.hebi.matlab.mat.types.Source;
import us.hebi.matlab.mat.types.Sources;

@SuppressWarnings({"javadoc"})
class TraceExporterTest {
  @Test
  void canReadWriteMatFile() throws IOException {
    // Write a double matrix
    final int rows = 4;
    final int cols = 5;
    @SuppressWarnings("resource")
    final Matrix out = Mat5.newMatrix(rows, cols);

    // row, col
    final int row = 3;
    final int col = 4;
    out.setDouble(row, col, 5.0);
    Assertions.assertEquals(5.0, out.getDouble(row, col));

    // From AbstractArray.getColumnMajorIndex(row, col)
    // Column major index: row + rows * col
    Assertions.assertEquals(5.0, out.getDouble(row + rows * col), "column major index");

    // Fill using index
    for (int i = 0, size = rows * cols; i < size; i++) {
      out.setDouble(i, i + 1);
    }

    final File file = File.createTempFile("double", ".mat");
    file.deleteOnExit();

    final String name = "tracks";
    @SuppressWarnings("resource")
    final MatFile matFile = Mat5.newMatFile().addArray(name, out);
    Mat5.writeToFile(matFile, file);

    try (Source source = Sources.openFile(file)) {
      @SuppressWarnings("resource")
      final MatFile mat = Mat5.newReader(source).readMat();
      @SuppressWarnings("resource")
      final Matrix in = mat.getMatrix(name);
      Assertions.assertEquals(rows, in.getNumRows());
      Assertions.assertEquals(cols, in.getNumCols());
      for (int i = 0, size = rows * cols; i < size; i++) {
        Assertions.assertEquals(i + 1, in.getDouble(i));
      }
    }
  }
}

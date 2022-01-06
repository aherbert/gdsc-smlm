/*-
 * #%L
 * Genome Damage and Stability Centre SMLM ImageJ Plugins
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2022 Alex Herbert
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
import us.hebi.matlab.mat.types.Cell;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.Matrix;
import us.hebi.matlab.mat.types.Source;
import us.hebi.matlab.mat.types.Sources;

@SuppressWarnings({"javadoc"})
class TraceExporterTest {
  @SuppressWarnings("resource")
  @Test
  void canReadWriteMatFile() throws IOException {
    // Write a double matrix
    final int rows = 4;
    final int cols = 5;
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
    final MatFile matFile = Mat5.newMatFile().addArray(name, out);
    Mat5.writeToFile(matFile, file);

    try (Source source = Sources.openFile(file)) {
      final MatFile mat = Mat5.newReader(source).readMat();
      final Matrix in = mat.getMatrix(name);
      Assertions.assertEquals(rows, in.getNumRows());
      Assertions.assertEquals(cols, in.getNumCols());
      for (int i = 0, size = rows * cols; i < size; i++) {
        Assertions.assertEquals(i + 1, in.getDouble(i));
      }
    }
  }

  @SuppressWarnings("resource")
  @Test
  void canReadWriteMatCellFile() throws IOException {
    // Create a cell: cell(1,2)
    final int crows = 1;
    final int ccols = 2;
    final Cell cell = Mat5.newCell(crows, ccols);

    // Write a matrix to two cells
    final int rows1 = 2;
    final int cols1 = 3;
    final Matrix m1 = Mat5.newMatrix(rows1, cols1);
    for (int i = 0, size = rows1 * cols1; i < size; i++) {
      m1.setDouble(i, i + 1);
    }
    final int rows2 = 4;
    final int cols2 = 3;
    final Matrix m2 = Mat5.newMatrix(rows2, cols2);
    for (int i = 0, size = rows2 * cols2; i < size; i++) {
      m2.setDouble(i, i + 10);
    }
    // zero-indexed not 1-indexed as per matlab
    cell.set(0, 0, m1);
    cell.set(0, 1, m2);

    final File file = File.createTempFile("double", ".mat");
    file.deleteOnExit();

    final String name = "tracks";
    final MatFile matFile = Mat5.newMatFile().addArray(name, cell);
    Mat5.writeToFile(matFile, file);

    try (Source source = Sources.openFile(file)) {
      final MatFile mat = Mat5.newReader(source).readMat();
      final Cell in = mat.getCell(name);
      Assertions.assertEquals(crows, in.getNumRows());
      Assertions.assertEquals(ccols, in.getNumCols());
      final Matrix m1b = in.getMatrix(0, 0);
      for (int i = 0, size = rows1 * cols1; i < size; i++) {
        Assertions.assertEquals(i + 1, m1b.getDouble(i));
      }
      final Matrix m2b = in.getMatrix(0, 1);
      for (int i = 0, size = rows2 * cols2; i < size; i++) {
        Assertions.assertEquals(i + 10, m2b.getDouble(i));
      }
    }
  }
}

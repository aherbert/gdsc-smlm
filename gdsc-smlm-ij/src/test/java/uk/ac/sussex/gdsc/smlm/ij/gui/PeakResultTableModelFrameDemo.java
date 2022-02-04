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

package uk.ac.sussex.gdsc.smlm.ij.gui;

import java.awt.EventQueue;
import javax.swing.DefaultListSelectionModel;
import javax.swing.ListSelectionModel;
import uk.ac.sussex.gdsc.core.utils.rng.SplitMix;
import uk.ac.sussex.gdsc.smlm.data.config.CalibrationWriter;
import uk.ac.sussex.gdsc.smlm.data.config.ResultsProtos.ResultsTableSettings;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.DistanceUnit;
import uk.ac.sussex.gdsc.smlm.data.config.UnitProtos.IntensityUnit;
import uk.ac.sussex.gdsc.smlm.results.ArrayPeakResultStore;
import uk.ac.sussex.gdsc.smlm.results.PeakResult;
import uk.ac.sussex.gdsc.smlm.results.PeakResultStoreList;

/**
 * Test application to show the PeakResultsTableModelFrame.
 */
public class PeakResultTableModelFrameDemo {
  /**
   * Launch the application.
   *
   * @param args the arguments
   */
  public static void main(String[] args) {
    final SplitMix r = SplitMix.new64(System.currentTimeMillis());
    final int n = 20;

    final ListSelectionModel selectionModel = new DefaultListSelectionModel();

    EventQueue.invokeLater((Runnable) () -> {
      try {
        final PeakResultStoreList store = new ArrayPeakResultStore(10);
        for (int i = n; i-- > 0;) {
          store.add(new PeakResult(r.nextInt(), r.nextInt(), r.nextInt(), r.nextFloat(),
              r.nextDouble(), r.nextFloat(), r.nextFloat(), PeakResult.createParams(r.nextFloat(),
                  r.nextFloat(), r.nextFloat(), r.nextFloat(), r.nextFloat()),
              null));
        }

        final CalibrationWriter cw = new CalibrationWriter();
        cw.setNmPerPixel(100);
        cw.setCountPerPhoton(10);
        cw.setDistanceUnit(DistanceUnit.PIXEL);
        cw.setIntensityUnit(IntensityUnit.COUNT);

        final ResultsTableSettings.Builder tableSettings = ResultsTableSettings.newBuilder();
        tableSettings.setDistanceUnit(DistanceUnit.NM);
        tableSettings.setIntensityUnit(IntensityUnit.PHOTON);
        tableSettings.setShowFittingData(true);
        tableSettings.setShowNoiseData(true);
        tableSettings.setShowPrecision(true);
        tableSettings.setRoundingPrecision(4);

        final PeakResultTableModel model =
            new PeakResultTableModel(store, cw.getCalibration(), null, tableSettings.build());

        final PeakResultTableModelFrame d =
            new PeakResultTableModelFrame(model, null, selectionModel);
        d.setTitle("D");
        d.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
        d.setVisible(true);

        // Selecting in one list activates the other list

        final PeakResultTableModelFrame d2 =
            new PeakResultTableModelFrame(model, null, selectionModel);
        d2.setTitle("D2");
        // Since we have the same selection model we need the same row sorter,
        // otherwise the selection is scrambled by sorting.
        // The alternative would be to get the source for the selection event (the table)
        // and get the row sorter to do the mapping.
        // However this breaks deletion of data as the row sorter double processes the deletion.
        // Basically only one table can use the same selection model when sorting is desired.
        // d2.table.setRowSorter(d.table.getRowSorter())
        d2.setDefaultCloseOperation(javax.swing.JFrame.EXIT_ON_CLOSE);
        d2.setVisible(true);
      } catch (final Exception ex) {
        ex.printStackTrace();
      }
    });
  }
}

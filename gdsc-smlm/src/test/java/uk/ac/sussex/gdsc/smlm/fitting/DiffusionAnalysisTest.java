/*-
 * #%L
 * Genome Damage and Stability Centre SMLM Package
 *
 * Software for single molecule localisation microscopy (SMLM)
 * %%
 * Copyright (C) 2011 - 2025 Alex Herbert
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

package uk.ac.sussex.gdsc.smlm.fitting;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import uk.ac.sussex.gdsc.test.api.Predicates;
import uk.ac.sussex.gdsc.test.api.TestAssertions;

@SuppressWarnings({"javadoc"})
class DiffusionAnalysisTest {
  @ParameterizedTest
  @CsvSource({
    // Data obtained from the python code:
    // https://gitlab.com/tjian-darzacq-lab/Spot-On-cli
    // fastspt.fastspt.C_AbsorBoundAUTO
    "0.35171539143446967, 0.01, 0.75, 0.0001, 1.0",
    "0.4409547738693467, 0.01, 1.5, 1.0, 0.9711316588468534",
    "0.7424623115577889, 0.01, 1.5, 1.0, 0.042506726296412944",
    "-0.7047738693467337, 0.02, 1.5, 1.0, 0.1789001958617121",
    "-0.2242462311557789, 0.03, 0.75, 1.0, 0.447549370861025",
    "0.3071608040201005, 0.03, 0.75, 1.5, 0.16233614943364597",
    "0.32223618090452255, 0.03, 0.75, 0.25, 0.3333969348168449",
    "0.28831658291457285, 0.03, 0.75, 0.025, 0.9747891295238313",
    // Boundary case.
    // Evaluates as -5.009431107187975e-17 so we clip to zero.
    "0.375, 0.03, 0.75, 1.5, 0",
  })
  void testWithinBound(double z, double dt, double dz, double d, double p) {
    TestAssertions.assertTest(p, DiffusionAnalysis.withinBound(z, dt, dz, d),
        Predicates.doublesAreRelativelyClose(1e-10));
  }

  @ParameterizedTest
  @CsvSource({
    // Data obtained from the python code:
    // https://gitlab.com/tjian-darzacq-lab/Spot-On-cli
    // Using integration of: fastspt.fastspt.C_AbsorBoundAUTO
    // Integration performed using scipy.integrate.quad
    // using default epsabs=1.49e-08, epsrel=1.49e-08.
    "0.03, 0.75, 0.025, 0.9175948369017197, 1e-7",
    "0.03, 0.75, 0.25, 0.7394119937446663, 1e-9",
    "0.01, 0.75, 2.5, 0.5244804117348121, 1e-8",
    "0.01, 0.75, 10, 0.14020993007941418, 1e-8",
    "0.02, 1.25, 5, 0.43129485401041096, 1e-8",
    // Cases requiring more than n=200 evaluations are not as accurate
    "0.01, 1.75, 0.025, 0.9796100101011112, 1e-5",
    "0.001, 1.75, 0.02, 0.994232839949668, 1e-3",
    // "Fixed" molecules
    "0.01, 0.75, 0.0001, 0.9969909888877455, 1e-5",
  })
  void testRemaining(double dt, double dz, double d, double p, double eps) {
    TestAssertions.assertTest(p, DiffusionAnalysis.remaining(dt, dz, d),
        Predicates.doublesAreRelativelyClose(eps));
  }
}

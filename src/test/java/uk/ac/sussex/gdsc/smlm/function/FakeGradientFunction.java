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

package uk.ac.sussex.gdsc.smlm.function;

import uk.ac.sussex.gdsc.core.data.NotImplementedException;
import uk.ac.sussex.gdsc.core.utils.PseudoRandomSequence;
import uk.ac.sussex.gdsc.core.utils.SimpleArrayUtils;
import uk.ac.sussex.gdsc.test.rng.RngUtils;

@SuppressWarnings({"javadoc"})
public class FakeGradientFunction
    implements ExtendedGradient2Function, Gradient2Function, Gradient1Function, NonLinearFunction {
  private final int maxx;
  private final int size;
  private final int nparams;
  private final PseudoRandomSequence rng;

  /**
   * Instantiates a new fake gradient function.
   *
   * @param maxx the maxx
   * @param nparams the nparams
   */
  public FakeGradientFunction(int maxx, int nparams) {
    this(maxx, nparams, 1000, 30051977, 10.0);
  }

  /**
   * Instantiates a new fake gradient function.
   *
   * @param maxx the maxx
   * @param nparams the nparams
   * @param scale the scale
   */
  public FakeGradientFunction(int maxx, int nparams, double scale) {
    this(maxx, nparams, 1000, 30051977, scale);
  }

  /**
   * Instantiates a new fake gradient function.
   *
   * @param maxx the maxx
   * @param nparams the nparams
   * @param randomSize the random size
   * @param randomSeed the random seed
   * @param scale the scale
   */
  public FakeGradientFunction(int maxx, int nparams, int randomSize, int randomSeed, double scale) {
    this.maxx = maxx;
    this.size = maxx * maxx;
    this.nparams = nparams;
    this.rng = new PseudoRandomSequence(randomSize, RngUtils.create(randomSeed), scale);
  }

  @Override
  public int size() {
    return size;
  }

  @Override
  public void initialise(double[] a) {
    int seed = 0;
    for (int i = a.length; i-- > 0;) {
      seed += Double.hashCode(a[i]);
    }
    // logger.fine(FunctionUtils.getSupplier("seed = %d", seed);
    rng.setSeed(seed);
  }

  @Override
  public void initialise0(double[] a) {
    initialise(a);
  }

  @Override
  public void initialise1(double[] a) {
    initialise(a);
  }

  @Override
  public void initialise2(double[] a) {
    initialise(a);
  }

  @Override
  public void initialiseExtended2(double[] a) {
    initialise(a);
  }

  @Override
  public int[] gradientIndices() {
    return SimpleArrayUtils.natural(nparams);
  }

  @Override
  public int getNumberOfGradients() {
    return nparams;
  }

  @Override
  public void forEach(ValueProcedure procedure) {
    // Simulate a 2D forEach
    for (int y = 0; y < maxx; y++) {
      for (int x = 0; x < maxx; x++) {
        procedure.execute(rng.nextDouble());
      }
    }
  }

  @Override
  public void forEach(Gradient1Procedure procedure) {
    final double[] gradient = new double[nparams];

    // Simulate a 2D forEach
    for (int y = 0; y < maxx; y++) {
      for (int x = 0; x < maxx; x++) {
        for (int j = nparams; j-- > 0;) {
          gradient[j] = rng.nextDouble() * x + y;
        }
        // System.out.println(Arrays.toString(dyDa));
        procedure.execute(rng.nextDouble(), gradient);
      }
    }
  }

  @Override
  public void forEach(Gradient2Procedure procedure) {
    final double[] gradient1 = new double[nparams];
    final double[] gradient2 = new double[nparams];

    // Simulate a 2D forEach
    for (int y = 0; y < maxx; y++) {
      for (int x = 0; x < maxx; x++) {
        for (int j = nparams; j-- > 0;) {
          gradient1[j] = rng.nextDouble() * x + y;
          gradient2[j] = rng.nextDouble() * x + y;
        }
        // System.out.println(Arrays.toString(dyDa));
        procedure.execute(rng.nextDouble(), gradient1, gradient2);
      }
    }
  }

  @Override
  public void forEach(ExtendedGradient2Procedure procedure) {
    final double[] gradient1 = new double[nparams];
    final double[] gradient2 = new double[nparams * nparams];

    // Simulate a 2D forEach
    for (int y = 0; y < maxx; y++) {
      for (int x = 0; x < maxx; x++) {
        for (int j = nparams; j-- > 0;) {
          gradient1[j] = rng.nextDouble() * x + y;
        }
        for (int j = gradient2.length; j-- > 0;) {
          gradient2[j] = rng.nextDouble() * x + y;
        }

        // System.out.println(Arrays.toString(dyDa));
        procedure.executeExtended(rng.nextDouble(), gradient1, gradient2);
      }
    }
  }

  @Override
  public double eval(int index, double[] dyda) {
    // Unpack the predictor to the 2D coordinates
    final int y = index / maxx;
    final int x = index % maxx;
    for (int j = nparams; j-- > 0;) {
      dyda[j] = rng.nextDouble() * x + y;
    }
    // System.out.println(Arrays.toString(dyDa));
    return rng.nextDouble();
  }

  @Override
  public double eval(int x) {
    return rng.nextDouble();
  }

  @Override
  public double evalw(int x, double[] dyda, double[] weight) {
    throw new NotImplementedException();
  }

  @Override
  public double evalw(int x, double[] weight) {
    throw new NotImplementedException();
  }

  @Override
  public boolean canComputeWeights() {
    return false;
  }
}

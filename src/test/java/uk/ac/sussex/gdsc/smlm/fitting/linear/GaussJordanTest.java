package uk.ac.sussex.gdsc.smlm.fitting.linear;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

@SuppressWarnings({"javadoc"})
public class GaussJordanTest {
  private static Logger logger;

  @BeforeAll
  public static void beforeAll() {
    logger = Logger.getLogger(GaussJordanTest.class.getName());
  }

  @AfterAll
  public static void afterAll() {
    logger = null;
  }

  @Test
  public void canSolveLinearEquation() {
    final GaussJordan solver = new GaussJordan();

    // Solves (one) linear equation, a x = b, for x[n]

    // Example taken from http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html
    final float[][] a =
        new float[][] {new float[] {1, 2, 3}, new float[] {4, 5, 6}, new float[] {7, 8, 10}};
    final float[] b = new float[] {3, 3, 4};
    final float[] expecteds = new float[] {-2, 1, 1};

    final boolean result = solver.solve(a, b);

    Assertions.assertTrue(result);
    Assertions.assertArrayEquals(expecteds, b, 1e-4f);

    if (logger.isLoggable(Level.INFO)) {
      logger.info(() -> String.format("x = %s", Arrays.toString(b)));
      for (int i = 0; i < b.length; i++) {
        final int ii = i;
        logger.info(() -> String.format("a[%d] = %s", ii, Arrays.toString(a[ii])));
      }
    }
  }
}

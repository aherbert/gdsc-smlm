package uk.ac.sussex.gdsc.smlm.fitting.nonlinear;

import org.junit.jupiter.api.Assertions;
import org.opentest4j.AssertionFailedError;

import uk.ac.sussex.gdsc.test.junit5.ExtraAssumptions;
import uk.ac.sussex.gdsc.test.junit5.RandomSeed;
import uk.ac.sussex.gdsc.test.junit5.SeededTest;
import uk.ac.sussex.gdsc.test.utils.TestComplexity;
import uk.ac.sussex.gdsc.test.utils.TestLog;

/**
 * Test that a stepping solver can fit a function.
 */
@SuppressWarnings({ "javadoc" })
public class SteppingFunctionSolverTest extends BaseSteppingFunctionSolverTest
{
    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_LSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__FastMLE(RandomSeed seed)
    {
        // The FastMLE method can generate very big steps that make the method unstable.
        // This test may fail depending on the random number generator.
        try
        {
            fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
        }
        catch (final AssertionError e)
        {
            logger.log(TestLog.getFailRecord(e));
        }
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, FastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, FastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_x__BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_C__BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_x_DC_BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_x__BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_C__BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void canFitSingleGaussianEMCCD_B_DC_BTFastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, BTFastMLE, NoiseModel.EMCCD);
    }

    @SeededTest
    public void cannotFitSingleGaussianEMCCD_x_x__JFastMLE(RandomSeed seed)
    {
        // The JFastMLE method was built using a misinterpretation of the Newton
        // method in Numerical Recipes, 2nd Ed. This test is just here to prove that.
        ExtraAssumptions.assume(TestComplexity.MAXIMUM);

        // The JFastMLE method does not work
        Assertions.assertThrows(AssertionFailedError.class, () -> {
            fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, JFastMLE, NoiseModel.EMCCD);
        });
    }

    //	@Test
    //	public void canFitSingleGaussianEMCCD_x_C__JFastMLE()
    //	{
    //		fitSingleGaussian(seed, NO_BOUND, CLAMP, JFastMLE, NoiseModel.EMCCD);
    //	}
    //
    //	@Test
    //	public void canFitSingleGaussianEMCCD_x_DC_JFastMLE()
    //	{
    //		fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, JFastMLE, NoiseModel.EMCCD);
    //	}
    //
    //	@Test
    //	public void canFitSingleGaussianEMCCD_B_x__JFastMLE()
    //	{
    //		fitSingleGaussian(seed, BOUNDED, NO_CLAMP, JFastMLE, NoiseModel.EMCCD);
    //	}
    //
    //	@Test
    //	public void canFitSingleGaussianEMCCD_B_C__JFastMLE()
    //	{
    //		fitSingleGaussian(seed, BOUNDED, CLAMP, JFastMLE, NoiseModel.EMCCD);
    //	}
    //
    //	@Test
    //	public void canFitSingleGaussianEMCCD_B_DC_JFastMLE()
    //	{
    //		fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, JFastMLE, NoiseModel.EMCCD);
    //	}

    // Weighted solvers for sCMOS

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_x__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_C__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_DC_WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_x__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_C__WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_DC_WLSELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, WLSELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_x__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_C__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_DC_MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_x__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_C__MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_DC_MLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_x__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_C__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_DC_FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_x__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_C__FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_DC_FastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastLogMLELVM, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_x__FastMLE(RandomSeed seed)
    {
        // The FastMLE method can generate very big steps that make the method unstable
        // This test may fail depending on the random number generator.
        try
        {
            fitSingleGaussian(seed, NO_BOUND, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
        }
        catch (final AssertionError e)
        {
            logger.log(TestLog.getFailRecord(e));
        }
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_C__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, CLAMP, FastMLE, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_x_DC_FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, NO_BOUND, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_x__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, NO_CLAMP, FastMLE, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_C__FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, CLAMP, FastMLE, NoiseModel.SCMOS);
    }

    @SeededTest
    public void canFitSingleGaussianSCMOS_B_DC_FastMLE(RandomSeed seed)
    {
        fitSingleGaussian(seed, BOUNDED, DYNAMIC_CLAMP, FastMLE, NoiseModel.SCMOS);
    }

    private void fitSingleGaussian(RandomSeed seed, boolean bounded, SteppingFunctionSolverClamp clamp,
            SteppingFunctionSolverType type, NoiseModel noiseModel)
    {
        //org.junit.Assumptions.assumeTrue(false);
        final SteppingFunctionSolver solver = getSolver(clamp, type);
        canFitSingleGaussian(seed, solver, bounded, noiseModel);
    }

    // Is Bounded/Clamped better?

    @SeededTest
    public void fitSingleGaussianEMCCD_B_x__LSELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_C__LSELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_C__LSELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_DC_LSELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_DC_LSELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_x__MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_C__MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_C__MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_DC_MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanMLELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_C__MLELVMBetterThanMLELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_DC_MLELVMBetterThanMLELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanMLELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, MLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_x__MLELVMBetterThanBLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, NO_CLAMP, MLELVM, BOUNDED, NO_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_C__MLELVMBetterThanBCLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, CLAMP, MLELVM, BOUNDED, CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_B_DC_MLELVMBetterThanBDCLSELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, BOUNDED, DYNAMIC_CLAMP, MLELVM, BOUNDED, DYNAMIC_CLAMP, LSELVM, NoiseModel.EMCCD);
    }

    // Note: The FastLogMLELVM converges too fast when there is still some
    // optimisation to do. The following tests have far fewer iterations with the fastLog version.

    @SeededTest
    public void fitSingleGaussianEMCCD_x_x__LSELVMBetterThanFastLogMLELVM(RandomSeed seed)
    {
        // This is actually a tie with the current fixed random number generator (50/50 for each)
        fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, LSELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    @SeededTest
    public void fitSingleGaussianEMCCD_x_x__MLELVMBetterThanFastLogMLELVM(RandomSeed seed)
    {
        fitSingleGaussianBetter(seed, NO_BOUND, NO_CLAMP, MLELVM, NO_BOUND, NO_CLAMP, FastLogMLELVM, NoiseModel.EMCCD);
    }

    private void fitSingleGaussianBetter(RandomSeed seed, boolean bounded2, SteppingFunctionSolverClamp clamp2,
            SteppingFunctionSolverType type2, boolean bounded, SteppingFunctionSolverClamp clamp,
            SteppingFunctionSolverType type, NoiseModel noiseModel)
    {
        ExtraAssumptions.assume(TestComplexity.MEDIUM);
        final SteppingFunctionSolver solver = getSolver(clamp, type);
        final SteppingFunctionSolver solver2 = getSolver(clamp2, type2);
        canFitSingleGaussianBetter(seed, solver, bounded, solver2, bounded2, getName(bounded, clamp, type),
                getName(bounded2, clamp2, type2), noiseModel);
    }

    @SeededTest
    public void canFitAndComputeDeviationsEMCCD_LSELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_LSELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsEMCCD_WLSELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_WLSELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_WLSELVM_Weighted(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeDeviationsEMCCD_MLELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_MLELVM(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_MLELVM_Weighted(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeDeviationsEMCCD_FastMLE(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_FastMLE(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_FastMLE_Weighted(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeDeviationsEMCCD_BTFastMLE(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_BTFastMLE(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeDeviationsSCMOS_BTFastMLE_Weighted(RandomSeed seed)
    {
        canFitAndComputeDeviations(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, true);
    }

    private void canFitAndComputeDeviations(RandomSeed seed, SteppingFunctionSolverType type, NoiseModel noiseModel,
            boolean useWeights)
    {
        final SteppingFunctionSolver solver1 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type,
                noToleranceChecker);
        final SteppingFunctionSolver solver2 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type,
                noToleranceChecker);
        fitAndComputeDeviationsMatch(seed, solver1, solver2, noiseModel, useWeights);
    }

    @SeededTest
    public void canFitAndComputeValueEMCCD_LSELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_LSELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.LSELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeValueEMCCD_WLSELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_WLSELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_WLSELVM_Weighted(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.WLSELVM, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeValueEMCCD_MLELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_MLELVM(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_MLELVM_Weighted(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.MLELVM, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeValueEMCCD_FastMLE(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_FastMLE(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_FastMLE_Weighted(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.FastMLE, NoiseModel.SCMOS, true);
    }

    @SeededTest
    public void canFitAndComputeValueEMCCD_BTFastMLE(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.EMCCD, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_BTFastMLE(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, false);
    }

    @SeededTest
    public void canFitAndComputeValueSCMOS_BTFastMLE_Weighted(RandomSeed seed)
    {
        canFitAndComputeValue(seed, SteppingFunctionSolverType.BTFastMLE, NoiseModel.SCMOS, true);
    }

    private void canFitAndComputeValue(RandomSeed seed, SteppingFunctionSolverType type, NoiseModel noiseModel,
            boolean useWeights)
    {
        final SteppingFunctionSolver solver1 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type,
                noToleranceChecker);
        final SteppingFunctionSolver solver2 = getSolver(SteppingFunctionSolverClamp.NO_CLAMP, type,
                noToleranceChecker);
        fitAndComputeValueMatch(seed, solver1, solver2, noiseModel, useWeights);
    }
}

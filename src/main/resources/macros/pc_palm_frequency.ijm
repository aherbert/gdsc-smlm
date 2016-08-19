// This is a test macro to run the PC-PALM simulation and then analyse and fit
// the results.
// It is used to test the estimation of the cluster parameters for different 
// cluster simulations.
// See Sengupta, et al (2013). Quantifying spatial resolution in 
// point-localisation superresolution images using pair correlation analysis. 
// Nature Protocols 8, pp345-354.

molecules = 3000;
precision = 20;
blinking = 3;

nm_per_pixel = 2;
correlationInterval = 15;
samples = 16;

cluster_number = 3;
cluster_variation = 0;
cluster_radius = 200;

correlation_distance = 400;

// None
// Circles
// Non-overlapping circles
// Circles Mask
distribution = "Poisson";
cluster_simulation = "Non-overlapping circles";
options = " blinking_distribution="+distribution+" cluster_simulation=["+cluster_simulation+"]";

curve_file = "/tmp/pcpalm.b"+blinking+"_p"+precision+"_c"+cluster_number+
    "_r"+cluster_radius+"_sim"+cluster_simulation+"_n"+molecules;


r=1;
datasets="";

// Run once to clear memory
run("PC-PALM Molecules", "run_mode=Simulation molecules=10 " +
    " simulation_size=16 blinking_rate=1 average_precision=20" + 
    " image_size=1024 roi_size=4 nm_per_pixel_limit=0 clear_results blinking_distribution=None");

for (j=0; j<10; j++) {

run("PC-PALM Molecules", "run_mode=Simulation molecules="+molecules+
    " simulation_size=16 blinking_rate=" + blinking + " average_precision=" + precision + 
    " image_size=1024 roi_size=4 nm_per_pixel_limit=0 " + options +
    " cluster_number="+cluster_number+" cluster_variation="+cluster_variation+
    " cluster_radius="+cluster_radius);

selectWindow("Molecule Simulation Binary Image (low res)");
getSelectionBounds(x, y, size, height);
w = getWidth() - size*0.75;
h = getHeight() - size*0.75;

i=0;
for (x=0; x<w && i<samples; x+=size) {
for (y=0; y<h && i<samples; y+=size, i++) {
    selectWindow("Molecule Simulation Binary Image (low res)");
    makeRectangle(x, y, size, size);
    run("PC-PALM Analysis", "correlation_distance="+correlation_distance+
        " blinking_rate=" + blinking + " nm_per_pixel=" + nm_per_pixel + 
        " show_error_bars apply_window correlation_interval="+ correlationInterval);
    datasets += "r_" + r + "=[" + 2*r + "*:] ";
    r++;
}}

} // End j loop

run("PC-PALM Fitting", "input=[Select PC-PALM Analysis results] " + datasets + 
    " estimated_precision=" + precision + 
    " blinking_rate=" + blinking + " show_error_bars fit_clustered_models" + 
    " fit_above_estimated_precision=1 fitting_tolerance=50 gr_random_threshold=1.50" +
    " save_correlation_curve output_correlation_file=[" + curve_file+"]"); 


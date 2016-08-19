// This is a test macro to run the PC-PALM simulation and then analyse and fit
// the results. 
// It is used to test the estimation of the blinking rate for different 
// blinking distributions.

molecules = 10000;

for (ii=0; ii<5; ii++) {
//molecules = 100 * (ii+1);
blinking = ii + 2;
for (j=0; j<5; j++) {

precision = 20;
//blinking = 1.5;
//blinking = 3;

nm_per_pixel = 2;
correlationInterval = 15;
samples = 12;
options = " blinking_distribution=Poisson";

run("PC-PALM Molecules", "run_mode=Simulation molecules="+molecules+
    " simulation_size=16 blinking_rate=" + blinking + " average_precision=" + precision + 
    " image_scale=1024 roi_size=4 nm_per_pixel_limit=0 clear_results" + options);

selectWindow("Molecule Simulation Binary Image (low res)");
getSelectionBounds(x, y, size, height);
w = getWidth() - size*0.75;
h = getHeight() - size*0.75;
r=1;
datasets="";
i=0;

for (x=0; x<w && i<samples; x+=size) {
for (y=0; y<h && i<samples; y+=size, i++) {
    selectWindow("Molecule Simulation Binary Image (low res)");
    makeRectangle(x, y, size, size);
    run("PC-PALM Analysis", "correlation_distance=300 blinking_rate=" + blinking + " nm_per_pixel=" + nm_per_pixel + 
        " show_error_bars  correlation_interval="+ correlationInterval);
    datasets += "r_" + r + "=[" + 2*r + "*:] ";
    r++;
}}
run("PC-PALM Fitting", 
    "input=[Select PC-PALM Analysis results] " + datasets + 
    " estimated_precision=" + precision + 
    " blinking_rate=" + blinking + 
    " show_error_bars fit_above_estimated_precision=1 fitting_tolerance=0 gr_random_threshold=1.50"); 

}
}

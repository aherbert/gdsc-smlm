// This is a test macro to run the PC-PALM simulation and cluster the results. 

// It is used to test the estimation of the flourophore photoefficiency 
// parameters and oligomer size as per the method of Puchnar:  
// Puchnar, et al (2013). Counting molecules in single organelles with 
// superresolution microscopy allows tracking of the endosome maturation 
// trajectory. PNAS. doi:10.1073/pnas.1309676110

// Binomial distribution parameters: (n,p)
for (p=30; p<=70; p+=10) {
for (n=2; n<=4; n++) {

// Molecules (will be incremented)
m = 1000;

precision = 25;
radius = 50;
max_n=1;

for (i=1; i<10; i++)
{
run("PC-PALM Molecules", "run_mode=Simulation image_size=1024 roi_size=4 nm_per_pixel_limit=0" +
    " clear_results molecules="+m+" simulation_size=16 blinking_rate=" + n +
    " blinking_distribution=Binomial average_precision="+precision+
    " cluster_simulation=None p-value="+ p);
run("Select None");
if (max_n!=0) {
    max_n=n;
}
run("PC-PALM Clusters", "distance="+radius+" algorithm=Pairwise max_n=" + max_n);
m+=1000;
}

for (i=0; i<10; i++)
{
run("PC-PALM Molecules", "run_mode=Simulation image_size=1024 roi_size=4 nm_per_pixel_limit=0" +
    " clear_results molecules="+m+" simulation_size=16 blinking_rate=" + n +
    " blinking_distribution=Binomial average_precision="+precision+
    " cluster_simulation=None p-value="+ p);
run("Select None");
if (max_n!=0) {
    max_n=n;
}
run("PC-PALM Clusters", "distance="+radius+" algorithm=Pairwise max_n=" + max_n);
m*=1.25;
}

}}

Moving Average Estimation of irregularly sampled data using Allan Variance (AVAR)
==========================================================================

This particular MATLAB directory contains:
* AVAR2 function: takes a temproal signal along with a list of window lengths and outputs the AVAR of the data for each window length
* Demo scripts (AVAR_irregular_sampling_demo_[...].m)
* other files which are used in demo scripts

After cloning the directory,
* AVAR_irregular_sampling_demo_simple.m script calculates AVAR of a deterministic signal currupted with noise. You can change the parameter dynamics at the very bottom of the script in the get_truth_at function
* AVAR_irregular_sampling_demo_random_walk.m script calculates AVAR of a Random Walk signal currupted with noise
* AVAR_irregular_sampling_demo_flicker.m script calculates AVAR of a Flicker signal currupted with noise
* AVAR_irregular_sampling_demo_jerath_signal.m script calculates AVAR of any arbirtary signal (default is jerath_signal.m which is a bias instability signal) currupted with noise

You can change noise characteristic from default noise (white Gaussian) to uniform white or flicker noise by changing noise type at line 22.


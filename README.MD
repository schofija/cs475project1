Project 1 for CS475

This is a basic Monte Carlo simulation-- see the details here:
https://web.engr.oregonstate.edu/~mjb/cs575/Projects/proj01.html

The included bash script compiles the program with a variety of thread counts/numtrials to
test the performance of the simulation when parallelized.

Results are redirected to a file named montecarlo.csv in the format of:

numthreads, numtrials, probability, maxperformance(megatrials/sec)

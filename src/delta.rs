/**
* Computes bounds on the value of delta by merging precomputed results
* 
* The program first reads a high-precision rational number on the standard
* input. It should be the exact result of the sum of the rk values with k < 10^6
* 
* Then it reads bounds on the partial sum. Each line is the sum between two indices
* with the following format:
* - "%d %d %la %la" (lower index) (upper index) (lower bound) (upper bound).
*
* The program then computes missing approximations till NB (2.16e8) and the remainder
* after NB. It finally sums everything to get an enclosure of delta.
* 
* NOTE: This runs after all the others I'm pretty sure
*/


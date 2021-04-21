## Distributed-Interference-Management-for-6G-inX-subnetworks


![6G in-X Use Cases](/images/6GinXUseCases.png)


## Usage

1. Run the file main_run.m with the deployment area dimension X, Y in meters and the type of SINR thresholding for channel switching decision
 * 'single' - for single SINR threshold. The value can be changed by specifying the index in the nCIndex variable (1-8 allowed)
 * 'anyOtherString' - SINR threshold for each bandwidth configuration as defined in envConstants
2. Modify envConstants.m to change simulation, deployment and propagation settings

## Distributed Algorithms Implemented
The algorithms are described in details in [1]
1. Random channel group selection
2. Greedy channel group selection 
3. Nearest neighbour conflict avaoidance (NNCA)
4. Greed channel group selection with repetition order adaptation
5. NNCA with repetition order adaptation 
6. Supervised Learning based channel group selection
7. Centralized graph coloring - presented as a benchmark

# References 

If you use this repository or any part of it, please cite the following articles

[1] Ramoni Adeogun, Gilberto Berardinelli, Ignacio R., Preben Mogensen, Distributed 


 
   

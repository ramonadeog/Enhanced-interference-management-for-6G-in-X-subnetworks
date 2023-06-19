# Enhanced interference management for 6G in-X subnetworks

Wireless systems are continuously faced with the demand to support higher reliability, lower latency, increased data rate and improved coverage. Consequently, beyond 5th Generation (5G) systems may need to support up to 10x lower latency and higher reliability than the 1 ms and 99.9999% limits in 5G. For example, industrial closed-loop control at the sensor-actuator level may require sub-milliseconds communication latencies with extremely high reliability in order to preserve stability of the control loop. Similar extreme connectivity requirements may also be demanded in other evolving life-critical wirelessapplications such as brake and ignition control in intra-vehicle communication, wireless heart pace-maker in intra-body networks and intra-avionics communication. Although use-cases and applications that may feature such extreme requirements are still evolving, recent visions on 6th Generation (6G) networks have identified independent and uncoordinated subnetworks (i.e., short range cells comprising of a controller acting as the access point for multiple devices) as potential solutions for supporting extreme  connectivity. Recently we introduced visions and design concepts for such 6G in-X subnetworks. The term in-X for insideeverything was introduced in [6] to highlight the emerging scenarios such as in-robots, in-vehicles, in-aircrafts, and inhuman bodies where these subnetworks are expected to be installed. These use cases can lead to situations with fastmoving subnetworks and hence highly dynamic interference conditions. Also, the lack of coordination may lead in some cases to high interference power translating to higher failure rates than tolerated. Thus, efficient and robust algorithms that
are capable of adapting utilization of the available multidimensional radio resources (such as frequency bands, time slots and transmit power) to dynamic interference conditions
under super-tight latency constraints are crucial for these systems. Since 6G subnetworks are expected to operate with very low power (e.g., -10dBm per channel), optimizing power usage may not yield any significant gains. We therefore focus on methods for dynamically managing the time-frequency resources.

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

## References 

If you use this repository or any part of it, please cite the following articles

[1] R. Adeogun, G. Berardinelli, I. Rodriguez and P. Mogensen, "Distributed Dynamic Channel Allocation in 6G in-X Subnetworks for Industrial Automation," 2020 IEEE Globecom Workshops (GC Wkshps, 2020, pp. 1-6, doi: 10.1109/GCWkshps50303.2020.9367532.

[2]. Adeogun, R. O., Berardinelli, G., & E. Mogensen, P. (2021). Learning to Dynamically Allocate Radio Resources in Mobile 6G in-X Subnetworks. Manuscript submitted for publication. In IEEE International Symposium on Personal, Indoor and Mobile Radio Communications (PIMRC)


 
   

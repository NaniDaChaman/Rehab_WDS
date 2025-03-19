2025-03-0919:51

Tags : #wds_modeling

Status :  ongoing

## Reading the ILP results section 
### Computational Results 
#### Models
#### Real world Networks 
#### Synthetic network
- 5 networks created 
- called randomly ==generated branched networks== 
- number of children nodes range  : 1-5
- elevation range (m): 100-300
- demand range (litres/sec): 0.01-5
- length of links : 500 -5000

#### potential solution 
- number of nodes and how they are connected : dnt know yet 
- elevation : randarray(N,range) => sort it in descanding order distribute 
- demand : randarray(N,range) => sort it in descanding order distribute 
- length of links : randarray(NE,range) => sort it in descanding order distribute 


## Reading optimal scheduling section 

### test network : 
- obtained from the paper [[Optimal_Water-Power_Flow_Problem_Formulation.pdf]] 
- its a looped network and they tested all their models on this network 
- For us its important that : we get optimal output and we get to test it on a real network 


## Creating Node Networks Questions
- [ ] How did the ILP model create water networks in the problem stated ?
- elevation, demand, length is known
	- [ ] graph gen is not known
	- [ ] Tank parameters and Pump parameters is not known
	- [ ] demand factor and time slices is not known 
- [ ] Are there ways to create water networks online ?
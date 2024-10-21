
## Introduction : 
- Why are we doing / interested in Water Distribution System Optimization in the first place ? 
- What made them choose a model like this ? 
- What's up with this shit ?

## Input of our Model : 


## The First Proposed Models : 
- The source of current study is is restricted to determining the pipe diameters in a single source acyclic (branched) network.
- Assumptions : 
1. ==Single source + gravity fed== - no pumps required anywhere else in the network !!!
2. No cycles - common in rural settings 
3. ==maintains pressure at a given head!!=
- Nomenclature :
1. C (S) = total pipe cost function of the diameter of pipe segment in each link ( $\sum_{i=1}^{NL}\sum_{j=1}^{NP}Li*xij*C(Ji)$)
2. Nl = number of links 
3. Di = Pipe diameter for link i (for one pipe per link problem)
4. Li = Length of link i 
5. C(Di) = cost/length of diameter i (function of diameter)
6. NP = Number of commercially available pipe diameters
7. xij = boolean variable if ith link uses jth daimeter  $\sum_{i=1}^{NL}\sum_{j=1}^{NP}xij=Nl$ [maybe this is our decison variable]
9. Pn = min pa to be maintained at node n 
10. Hr = head supplied by node R
10. En =Elevation of node n 
11. Sn = Set of pipes that connect n->R
12. HLij = Headloss in link i due to diameter j

### Input 
![[Pasted image 20240923094315.png]]

- MBR : the source node : consist of pressure supplied water (kept at an elevation higher than the other nodes)
- Nodes = water demands, min pressure requirments and elevation  < MBR's elevation
- Link = startnode, endnode, length
- commercial pipe diameter = cost per unit length and roughness

### Output 
- Length and diameter of pipe segments for each link 

## Objective 
- Total cost of installing the pipes. (minimize)
1. One pipe per link $xij\in I$
$\sum_{i=1}^{NL}\sum_{j=1}^{NP}Li*xij*C(Dij)$
2. Proposed Model $xij\in R$
$\sum_{i=1}^{NL}\sum_{j=1}^{NP}Li*xij*C(Dij)$
## Constraints
- Pressure at each node must exceed min pressure specified
1. Hazen Williams Equations to find out pressure at each node !!! : 
			$Hr - En -\sum_{i\in Sn}\sum_{j=1}^{Np}HLij*Li*xij$
			$HLij=\frac{10.68*Li*\frac{flowi}{roughnessi}^a}{diameter^b}$
			a=1.852
			b=4.87			
1. min pressure constraint : 
			$Pn <= Hr - En -\sum_{i\in Sn}\sum_{j=1}^{Np}HLij*xij$
- Water demand must be met for each node 
- Pipe diameters can only take values from the commercially provided pipe diameters (implicitly modelled in as j)
- One Pipe per link : $\sum_{j=1}^{j=NP}xij = 1$
- Proposed model link constraints : $\sum_{j=1}^{j=NP}Li*xij = Li => \sum_{j=1}^{j=NP}xij = 1$ 

## The second proposed model : 
- 
# Questions : 
- [ ] What is the Hazen-Williams equation to calculate pressure at each node ? [check here](https://en.wikipedia.org/wiki/Hazen%E2%80%93Williams_equation)
- [ ] What is convexity conditions for the 2 pip segment problem ? 
- [ ] Head loss in Link i due diameter j (HLij what is it and how is it calculated + what id there are multiple daimeters present in link j at different lengths?)
- HLij for length varied systems becomes Li*xij  or lij
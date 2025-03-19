2024-09-1008:50

Tags : #ilp #pipeselprob

Status : ongoing


## Todo 
- [x] Note down new stuff parameters being used
- [ ] Look at how primary and secondary networks are created, refer to the paper before 
- [x] in the input file some demand nodes are empty what does that mean ?
- Could be that they are intermediary nodes and have no demand
- could be that they do have demand 
- [x] FLp and FLp should be variables ?
- No
 - [x] Read inputs  from excel file 
 -  [ ] Why do all pipes have the same roughness = default roughness ?  
		 - Check the tutorial video it might have some answers !
- [ ] Add link and tank constraints
-  [x] Adding head loss constraint at every link
-  [ ] Adding effective head caused by a tank

# Sections
1. Introduction 
2. Initial Model
3. Pipe headloss improvement
4. Tank Cost Improvement
5. Tank Config Improvement
6. Edge-based model
7. Computational results
8. Conclusion

# Important Terms and their Definition 
- Head Loss
- Water Pressure Head
- Use of Tanks 
- Pipe diameter selection 

## Initial Model
### Inputs :
- Mis stuff : 

	1. water density : $\delta$, ==not in branched model==
	2. $g$,==not in branched model==

- Primary , Secondary Supply hours ,min/maximum head loss per km, max water speed 

	1. hours water supply in primary : $PH$
	2. hours water supply in secondary : $SH$
	3. The primary network is the main system of pipes and infrastructure that transports water from the source (such as a water treatment plant or reservoir) to various storage facilities
	4. The secondary network is the subsystem that distributes water from the primary network (or storage tanks) to individual demand nodes, such as households, businesses, and public facilities
	
- Links in a node : 

	1. Number of existing links : $NL$
	2. Length of the link : $L_i$ where $i \in {1,..,NL}$
	3. Flow in ith link, if in secondary network $FL_i^s$,==only had fi==
	4.  Flow in ith link, if in primary network $FL_i^p$,==only had fi==

- Diameters of commercially available pipes with their roughness.

	1. Number of pipes : $NP$
	2. Diameter of pipe : $Dj$
	3. Cost per unit length : $C_j$
	4. Roughnes of pipe : $R_j$
	
-  Nodes : elevation and location

	1. Number of nodes : $NN$,==N in the previous network==
	2. Minimum pressure requirement of a node : $PR_n$,==P in the previous network==
	3. Elevation of a node : $E_n$==E in the previous network==
	4. Water demand of node $DE_n$,==D in the previous network==
	5. Set of ancestor nodes : $A_n$, ==$S_n$ in the previous network==
	6. Set of descendand nodes : $D_n$, ==new to this network==
	7. parent nodes : $P_n$,==new to this network , but can be dervied from S_n==
	8. set of incoming links for node n : $I_n$,==new to this network==
	9. set of outgoing links : $O_n$, ==going out of this network==
	
- Source node : $S$

- Tanks : , ==new to this network==

	1. Number of rows in TCT : $NE$
	2. Base cost of kth row : $B_k$
	3. Unit cost of kth row : $UN_k$
	4. Upper limit : $UP_k$
	5. Lower Limit : $LO_k$
	6. Min tank height allowed : $T_{min}$
	7. Max tank height allowed : $T_{max}$
	
   - Tank Cost Table : 
   1. Tank1 -> base cost (cost associated with its building capacity), unit cost (cost associated per unit of tank capacity) additional cost associated with building that tank, capacity ranges (what is the min/max capacity for this row of the tank ) 
   2. The cost of a tank is a piecewise linear function
   
- Pumps : ,==new==

1. Capital cost/kW : $CP$
2. Energy cost/kWh : $EP$, ==where does it come in==
3. Discount Factor over its lifetime (==how the pipe's energy lowers per day ?==): $DF$
4. Inflation rate : $INFR$
5. Interest rate (==rise of value for money==): $INTR$
6. Pump efficiency : $\eta$
7. min pump power: $PP_{min}$
8. max pump power: $PP_{max}$
9. design life time (==could be a how long a pipe can run ?==) : $Y$
10. pipes that cannot have pumps (==maybe they don't have the right attatchment !!!!==)

- Valves : location, pressure rating
### Outputs : 
- length(l) and diameter (j) of pipe segment for each link(i) $l_{ij}$
	1. length in primary network : $l_{i,j}^p$
	2. if link(i) is part of the primary network : $f_i$
	
- total head across link i : $hl_i$
	1. effective head provided to link(i) by starting node (n) : $h_n$ 
	
- total demand being served by tank (k) at node (n) : $z_{n,k}$
	1. if tank(k) at node (n) is with : $e_{n,k}$ 

- if node (m) being served by a tank(n) : $s_{nm}$ 
	1. demand of each node : $d_{n}$
	2. tank height : $t_n$
	3. if source of water for link (i) is its immediate upstream node (n) : $es_{n,i}$
	
- link and power of pump : $p_i$
	1. if primary : $p_i^p$ 
	2.  if secondary : $p_i^s$ 
	3. head by a pump at link i : $ph_i$
	4. if pump is installed at link i : $pe_i$
### Objective : 
- Minimize capex cost (pipe+tank+pump) and total energy cost (pump)
- Capital cost of pipes = $\sum_{i=1}^{NL}\sum_{j=1}^{NP}C_j*D_j*l_{i,j}$
- Tank Capital cost = $\sum_{n=1}^{NN}\sum_{k=1}^{NE}(B_k+UN_k(d_n-LO_k))e_{n,k}$
- Pump Capital cost = $\sum_{i=1}^{NL}CP*p_i$
- Energy Pump cost = $EP*DF*(\sum_{i=1}^{NL}PH*p_i^p+\sum_{i=1}^{NL}SH*p_i^s)$ ==this will become a parameter of the lower level model==

### Constraints : 
- Pressure at each node must be at least the minimum pressure specified $h_n$<-$hl_i$<-$l_{ij},s_{nm}$
- Water demand must be met 
- Calculation of pressure loss across a link $\sum_{j=1}^{NP}(HL^p_{i,j}*l^p_{i,j}+HL^s_{i,j}*(l_{i,j}-l^p_{i,j})-ph_i)$

## Julia Implementation
- Cant find interest rates anywhere for pump opex calc
```julia

```

## Tanks as Consideration : 
- When we add tanks to a Water network : it divides our network into primary and secondary networks !!!!
- Primary Network : The primary network is the main system of pipes and infrastructure that transports water from the source (such as a water treatment plant or reservoir) to various storage facilities, such as tanks or elevated storage reservoirs (ESRs).
- Secondary Network : The secondary network is the subsystem that distributes water from the primary network (or storage tanks) to individual demand nodes, such as households, businesses, and public facilities.
- We assume that nodes and links that connect a source -> tank do not demand!!! which allows them to be modeled with a binary variable !!!
$h_n=h_{m,i}-hl_{i}$
$h_{m,i}=(t_m+E_m)*es_{m,i}+h_m*(1-es_{m,i})$
$es_{m,i}=s_{m,m}*(1-f_i)$

### tank hierarchy
- $s_{m,m}<=s_{n,n}, m\in D_n$
- $s_{n,m}<=s_{n,n},m\in D_n$
- $\sum_{p\in A_n} s_{p,n}=1$
- $d_n=\sum_{m\in D_n}s_{n,m}*DE$
- $f_i=s_{n,n},i\in I_n$
- In is the incoming node of node n
- Dn are descendants of node n 
- An is the ascendants of node n
- DE is the demand of all nodes 


## Inputs for a Tank : 
- We are given a tank cost table !!!
- Tank min size 
- Tank max size

### Parameters : 
- NE - number of rows in tank cost table 
- Bk - base cost of the kth row of the tank cost table 
- UNk- Unit cost per capacity of the kth row of the tank cost table
- UPk - Upper limit capacity of the kth row of the tank cost table 
- LOk- Lower limit capacity of the kth row of the tank cost table 
- Tmin 
- Tmax

### Cont vars : 
- tn :height of the effective tank 

## Binary Vars : 
- enk : if tank at node n is costed by the kth node of the
- fi : if the link at node i is part of the binary network or the secondary network 
- esni : if the source for water for link i is its immediate upstream node n, 0 otherwise
- snm : 1 if the tank at n node provides water to m node 

## Objective Function : 
$\sum_{i \in links}CP*p_i$
- $EP*DF*(\sum_{i \in links}PH*p_i)$


## Pumps constraint : 
- $p_i=\frac{(\rho*g*FL_i*ph_i)}{\eta}$
- $Pmin*pe_i<=p_i<=Pmax*pe_i$
- $hl'_i=hl_i-ph_i$
- 
## Terms we don't understand : 
- Headloss :
1. Its the heatloss created by friction when a fluid flows through pipes. 
2. Its also related to hazen williams equations : 
- Water Demand for nodes 
- Tank Cost Table : maybe find it in the previous paper or this papers previous parts !!!
- Efficency of a pump : 
- Flow of water : the amount of water (in volume) passing through a point (cross-section of area) per second.
1. It probably relates to water demands of various nodes as more the demand more needs to be the flow from a pipe
2. It also probably relates to the supply hours of a node 
- Primary and Secondary Networks : 
1. Primary Network is one where water is supplied from source to tanks (links and nodes which are part of that)
2. Secondary Network : is one where water is supplied to the demand nodes from the source and tanks 
3. What if a node is in middle of tank and pump ? 
$lij*fi = lijp$
$lij*(1-fi) =lijs$
for a given link i fi will determine if its part of the primary and secondary network but it can't be part of both !!!
4. Hyposthesis : we will limit our model to exclude cases like this 
## ILP examples for WDS systems Questions
- [ ] How did they manage to improve the water distribution system model ? (why are their changes an improvement ?)
- [ ] What is the papers' main claim ? 
- They have build a model that can perform pump ,pipes  ,valves and tank optimization as well !!!
- India has software which can do pipe selection optimization but they have done that better somehow !!!???
- [ ] Can their models guarantee liveness and other requirements ?
- [ ] How can we extend this paper ? 
- Include looped networks
- scheduling of water distribution !!!
- [ ] What is the initial model used for JalTantra  ?  
- [ ] How does the tank serve as a secondary source of water ? 
- In rural network its common that tanks will serve as a second source of water rather than a way to stop water.
- Water comes in at a high pressure from the nodes above !! -> incoming into tank -> water comes out of tank at a lower pressure !!! 
- ==How can we express that thru equations?==
- [ ] How do pumps and valves affect the Hazen williams equation!!!
- [ ] What is the decision variable of the pump !!!! -> power of the pump at link i !!!! (which then decides the head provided at link i!!!)
- [ ] How do we convert all our constraints and variables into a standard ILP form ? 
- [ ] What is the relation between power of a tank and the head created by it !!!
- [ ] How can location and power of each pump be used to serve the other pumps !!!
- [ ] How do you create a feasible problem case out of this !
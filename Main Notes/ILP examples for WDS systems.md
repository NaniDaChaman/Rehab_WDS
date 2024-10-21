2024-09-1008:50

Tags : #ilp #pipeselprob

Status : ongoing

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
- Primary , Secondary Supply hours ,min/maximum head loss per km, max water speed 
- Links in a node : length
- Diameters of commercially available pipes with their roughness.
-  Nodes : elevation and location
- Already existing pipes : nodeS->nodeE, length, diameter, parallel allowed, roughness
- Source node : 
- Tanks : 
   - Tank Cost Table : 
   1. Tank1 -> base cost (cost associated with its building capacity), unit cost (cost associated per unit of tank capacity) additional cost associated with building that tank, capacity ranges (what is the min/max capacity for this row of the tank ) 
   2. The cost of a tank is a piecewise linear function
- Pumps : 
1. min pump size : (what is the variable associated with it and how does it affect our model ?)
2. effeciency : (")
3. design life time (==could be a how long a pipe can run ?==)
4. discount/interest rate (==how the pipe's energy lowers per day ?==)
5. pipes that cannot have pumps (==maybe they don't have the right attatchment !!!!==)
- Valves : 
### Outputs : 
### Objective : 
### Constraints : 

## Tanks as Consideration : 
- When we add tanks to a Water network : it divides our network into primary and secondary networks !!!!
- Primary Network : The primary network is the main system of pipes and infrastructure that transports water from the source (such as a water treatment plant or reservoir) to various storage facilities, such as tanks or elevated storage reservoirs (ESRs).
- Secondary Network : The secondary network is the subsystem that distributes water from the primary network (or storage tanks) to individual demand nodes, such as households, businesses, and public facilities.
- We assume that nodes and links that connect a source -> tank do not demand!!! which allows them to be modeled with a binary variable !!!

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
$pipeCost+pumpCost+\sum_{k=1}$
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
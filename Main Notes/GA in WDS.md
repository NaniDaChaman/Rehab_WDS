2024-10-1817:58

Tags :

Status : 

# Genetic Algorithm 
- [ ] What are features of a Genetic Algorithm ? 
- Population : a set of individuals which die and are replaced by offspring
- Breeding : the offspring are formed by a combination of genes from their parents taken from the population by breeding
- Selection : this involves selecting parents for breeding from the population (all aren't chosen)
- [ ] How are parents chosen for the GA ? 
- [ ] What is the fitness function and how is stuff evaluated from it ? 

### Population : 
- [ ] How is the initial population created : 
- design parameters are expressed as numerical vales
- a binary string with each of these parameters is expressed as 
- A bunch of these are made randomly keeping in mind our convex set ??

## Breeding : 
- [ ] How do we select which induviduals should breed ? 
- Fitness function
- Roulete wheel 
- Tournament 

## CrossOver : 
- Take the design char of two parents and mix them to create two children
- Takes the best traits of our current generations and adds them to our next generation
### General Cross over
- Pool together the genetic characteristic of parents
- grow a feasible offspring
### Cross over for WDS 
1. One point Crossover : [1,1,0,1] x [1,0,1,0] , randomly choose a point :2 => [1,1,1,0],[1,0,0,1]
2. 2 point crossover : "",choose two points after and before which genetic material is switched , creates 2 points
3. k point crossover : , creates $2$ points

## Mutations : 
- Take the design char of a child and flip one of them based on a random draw with a prob of 1/l
- Adds variety to our next generation

## Elitism : 
- Take the best solutions and put it in the next generation 

## Applying GA to WDS

### Population Initialization :
1. A subset of binary string reps possible combo of pipe sizes
### Network Costs : 
1. Hydraulic analysis : for every string there will be heads associated with each node and discharges, these will be gathered using a ==hydraulic net solver==
2. Penalty Costs : When a pipe is not within the pressure and velocity limits we calculate a penalty cost for the string!
3. Overall net cost(nc) = net cost + penalty cost for each string 
4. fitness = f(nc) of the overall net cost ,1/nc usually
5. generation of new population : selection (init_pop)
6. Crossover : 
7. Mutation

## GA in WDS Questions
- [ ] WHat is a design parameter ? 
- It cnt be a decision variable  (binary string 1 is associated with weather that design parameter is present) which can't represent desc vars
- [ ] How can we represent a design decision ? 
2024-12-2210:46

Tags : #proposal

Status : ongoing


## Reading SMGA : 

## Introduction : 
- [ ] Why do networks require rehabilitation ?
- Integrate new technology
- Pipes may increase their roughness : fix that
- Nodal demand may have increased so they might wanna deal with that.
- [ ] What is multi objective opt and why is network rehab a multi obj opt
- An optimization problem which consist of two objectives. 
- Find a parametric solution S(x) for both of the optimization problem which relies on input from the other optimisation problem to find the best answer. 
- Find the best solution in the parametric solution : paretor optimum
- Network rehab : will have two opposed effects -> cost of the network which we wanna minimize, benefit that the net provides should be maximized.
- [ ] What is the SMGA algorithm  ?
- Its a genetic algorithm used for MOO problems
- it takes simple partial solutions and develops more complex solutions in every iteration
- it outputs a range of solution for the moo problem, which range from low to high on both objectives.
- [ ] What are the adv and disadv of single objective problems in context of wds rehab ? 
- [ ] How is SMGA tested ?


## Sections : 
- Intro
- Approach (talks about the paper sections I assume)
- Problem formulation (Objective function)
- SMGA (Describes the actual algo)
- SMGA parameters 
- Anytown (Describes the network + What should our software give as an output!+experiment)
- Conclusion

## Approach 
- [ ] Why is modelling pumps and tanks in a steady state simulator (used in GA) more difficult ? 
- Know what a GA is + what is exactly modelled
- Know what a SMGA is 
- Then you'd get to know what difficulty arises
- (Think about pumps and tanks difficulty )
- [ ] What is the difference between the two approaches of SMGA ? 

## Problem formulation : 
- Optimal Design and mode of operation for a wds with pumps and tanks 
- desc vars :- 
- two objectives : $max_if(i)=Benfit(i)$, $min_iF(i)=Cost(i)$ 
- [ ] How can benefit function be formulated and what are its parameter + vars ?
	- [ ] What is Nodal Pressure Shortfall ? 
	- [ ] What is Tank Flow Difference ? 
	- [ ] What is TLD ? 
	- demand balancing flow : its the trail in/out flows for a trail volume produced at the start of SMGA iterations
	- Tank demand flw : its the demand in/outflow for a trail volume produced at the start of SMGA
	- Nodal demand of a tank = $\sum_{out-nodes}D_n$
	- Nodal heads : simulated head f the tank for the nodal demand and min and max head are also found for each simulation period
	- Tank operatigng level diff =$\sum_{i=1}^4minh_i-trailminl_i$+$\sum_{i=1}^4maxh_i-trailmaxl_i$
	- ==why did we end up calculating in and out flprlows of trail (demand balancing flow) and one produced by demands==
- [ ] How can cost function be formulated and what are its parameter + vars ?



## Any town network 
- What is the problem 
	 1. Dec vars : pipes that can be installed or reinforced
	 2. new pumps installed or replaced in existing pumping station
	 3. new tanks installed
	 4. Mode of operation : when should pumps be on/off 
	 5. Objective function : The idea is that we have a projection of demands for a given number of years , that take into account pumping cost and capital expenditure 
- What are we given is the network ? 
	1. Water Reservoir and three pumps : ==what are the 5 point char for these pumps?==
	2. Flow rate/discharge rate = the volume it can pump/ transport across a pipe in sec
	3. Pump head : the pressure that the pump imparts 
	4. Efficiency : the transfer from electric energy to mechanical energy 
	5. These are fixed speed pumps so they can output more or less their head imparted is fixed and require a fixed amount of power per second to operate !
	6. Node demand for a given set of years (1985-2005), along with elevation
	7. 2 tanks (water level 225-250)
	8. Variation of water use/demand in table 3 ==what is a demand factor and how does it relate to node demand== , its like the std deviation from avg daily demand 
	9. Min pressure at each node 40 psi (peak flow) 20 psi (fire flow)
	 peak flow : the idea is that its those times of days where more than avg flow is required
	fire flow : those times where some nodes require more water than other nodes
	6. Cost of laying pipe, cleaning pipe per unit length
	7. Capex Cost of pumps on discharge rates $Q_r$ and head $H_r$
	8. Opex Cost is unit cost of energy 0.12 $/KWh on 12% interest  consideration of 20 years
	9. Tank cost table :  up/lw tank capacity ,cost per unit capacity, base cost 
	10. Problem modelling described before 
- How does our MOO help in finding solutions for this network ?
- 

## Problem Description : 
-  Problem formulation The problem is to find the optimal design and mode of operation of a pumped water distribution system with storage.
- 
## Talk w Nitya Questions : 
- [ ] WHat is a perto optimum ? 
- [ ] in 2 person games we had , our input variables being controlled by another player, after which we found a pareto optimum, here we don't have that case do we !!!

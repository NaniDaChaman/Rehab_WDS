2024-10-2109:02

Tags : #scheduling 

Status : ongoing

## Sections
- Intro : talks about why is pump optimization kind of important and diff approaches to pump opt
- Section 2 : General model for wds 
- Section 3 : OWF problem to min electricity op cost for fixed speed pumps 
- Section 4-5 : convex relaxation + penalty term + when are the two min same ? 
- Section 6 : benchmark
## Abstract : 
- [ ] What is the optimal scheduling problem for WDS 
- [ ] How can we express this problem as the optimal scheduling problem ? 

## Motivation for Optimal Pump scheduling 
1.  a recent survey on WDS optimization identifies pump scheduling and water quality as the two focus areas [5]. Recognizing that 4% of the total electricity consumption in the United States is attributed to water network operations [6], and that the electricity cost for pumping constitutes the largest expenditure for water utilities [7],
2. Pump at night, but in the ==smart city vision with dynamic electricity pricing and demand-response programs WDS schedules minimize operational costs incentivize flexible wds net==

## Water Network Modeling 
- G = (M,P) this our water network
## Nodes (M)
- $M=Mr+Mb+Mc$ and $Mr\cap Mc\cap Mb=\phi$
- Mr are the nodes hosting reservoirs
- Mb are nodes hosting tanks
- Mc are nodes hosting Consumers
- $d_m^t=$ Water flowing out of m at t, ==in our wds dm is the demand of each node (is a parameter) is it a variable here whereas it was a parameter of our model ==, its the negative of our demand for each node
- $m\in Mr \rightarrow d_m^t>=0$
-  $m\in Mb \rightarrow d_m^t \in R$
-  $m\in Mc \rightarrow d_m^t<=0$
- $d_m^t=0$ in junctions
## Pipes (P)
- $|P|$, all edges are directional (m,n)
- $(m,n)\in P \rightarrow (n,m) \not\in P$
- $P_a=$ edges which denote the ideal pump
- $\overline P=$ lossy pipes
- $d_{mn}^t =$ water flow from node m,n, $FL_i$ in our equations
- $d_{mn}^t >=0$ if water flows from  m to n ,rlse negative !!!
- $d_{m}^t = \sum_{k:(m,k)\in P}d_{m,k}^t - \sum_{k:(k,m)\in P}d_{k,m}^t$

## Pressure (h)
- $h_m^t=$ pressure at that node at time t
- In static analysis = $h_m^t=E_n+t_n$ 
- Minimum pressure is a common minimum value for all nodes m which our pressure should be above for water to flow $h_m^t>= \underline{h_m}$ ==same as our model==
- Darcy-Welsbach equation : gives us pipe loss due to friction ==in our equation we use the henry Williams equations==
	1. $h_m^t-h_n^t=c_{mn}*sign(d_{m,n}^t)*(d_{m,n}^t)^2$
	2. $c_{m,n}=\frac{l_{m,n}*f_{m,n}}{4*\pi^2*r_{m,n}^5*g}$
	3. $f_{m,n}$ is a constant
	4. Big M trick : 
		1. $-M*(1-x_{m,n}^t)<=d_{m,n}^t<=M*x_{m,n}^t$
		2. $-M*(1-x_{m,n}^t)<=h_m^t-h_n^t-c_{m,n}*(d_{m,n}^t)^2<=M*x_{m,n}^t$
		3. $-M*(1-x_{m,n}^t)<=h_m^t-h_n^t+c_{m,n}*(d_{m,n}^t)^2<=M*x_{m,n}^t$
		4. $x_{m,n}^t =1$ if water flows from m->n else 0

## Pumps
1. $P_a\in P$ , ideal pumps 
2. $P_a'^=P/P_a$ , lossy pipes
3. fixed speed pump : they operate at a fixed speed ???, they pump out water and the pressure given out by them constant
4. if a pump is running there is a fixed flow of water in the link !!!
5. reference to $pe_{m,n}$ is refers to the ideal segment of the pump
6. Let $x_{m,n}^t$ , indicate the time when the pump is on and when the pump is off 
7. Let $(n,m)\in P_a$ then $g_{n,m}*x_{n,m}^t=(h_n^t-h_m^t)$
8. Let $(n,m) \in P_a$ then $dmin_{n,m}*x_{n,m}^t<=d_{n,m}^t*x_{n,m}^t<=dmax_{n,m}*x_{n,m}^t$ ==to be expanded further==
9. ==In our model we consider variable speed pumps !!! : the power of the pump changes the pressure that is applied to it ==
10. ==the power and head created by the pump stays constant for the scheme of the optimization==
11. ==the choice from our first model ($p_i,ph_i$) becomes the parameter of our second model($g_{n,m}$)== 

## Reservoir
1. Our decision is when should we draw water from this reservoir attached to node m
2. $\alpha_m^t$ is our decision variable
3. $0<=d_m^t<=M*\alpha_m^t$ , water shall flow if there is pressure 
4. $h_m^t<=hrev_m+M*(1-\alpha_m^t)$
5. ==This is not given by the upper level model so we'll probably ignore it !!!==
6. ==alt we can look to expand the upper level model so that they can make this decision==

## Tanks 
1. $l_m^t=l_m^{t-1}-\frac{d^t_m*\delta}{A_m}$, the level of the tank changes, according to the water that flows out of it and into the connecting node m in duration $\delta$
2. $lmin_m<=l_m^t<=lmax_m$, ==set by parameters of the upper level problem==
3. $l_m^0=l_m^T$
4. let $\alpha_m^t$ indicate if the tank is connected at time , by the outlet valve to node m
5. let $\beta_m^t$ indicate if the tank if filling at time t, by the inlet valve to node n
6. When $\beta_m^t=0$, $h_m^t<=l_m^t$, when $\beta_m^t$=1, $h_m^t>=l_m^t$
7.  ...
8. ==the ability of tanks to connect and disconnect was not present in our first model (where should it go ?)== , we treat the upper models pressure identified as the bounds of the lower model, and optimise it as a whole model
9. ==the ability of tanks to go to be filled or not (where should this go in our upper model)==, this is a empty node created for the tank ( connected to the real node n by a lossless pipe )

## Objective Function 

1. $\sum_{t=1}^T\sum_{m,n\in P_a}e*c_{m,n}*d^t_{m,n}$
2. we want to minimize the cost of electricity for the given control period T, this is opex of pumps!
		1. I think it fits in with the objective function of the upper level problem!
3. $c_{m,n}=\frac{\delta*\rho*g*g'}{\eta}$

## Problem Formulation
- [ ] What is our decision variables, objectives and constraints ? 
- [ ] What is the time duration $\delta$ ? 
- We discretize our continuous system into time chunks t=0,1,2,...T
- 0 represents the start of our system while T represents the end
- All of these chunks are of a constant duration $\delta$ 

## Testing our Algorithms 
- [ ] Are there sample inputs/outputs provided ? 
- [ ] Is the model available ? 

## Recovering Pressure Differences in P1 : 
1. Given ${dbar,d}$ $\forall m \in V, \forall t \in T$ : 
2. [Pumps eq 6](## Pumps) and [Pipes eq 1](## pipes) can give us the values of pressure differences !!! call this a vector of pressure diff $b(dbar^t,d^t)$
3. Recover Pressure from pdiff Create $A(d^t)$ which has possitive or neg signs on the nodes the pipes flow from 
	- $A(d^t)h^t=b(dbar^t,d^t)$ $\forall t \in T$
	- h

## Optimal Scheduling of Water Distribution Systems Questions
- [ ] Paper claim
- Develops a general model for WDS components
- OWF problem : minimize electricity cost for fixed speed pumping 
- Convex relaxation of the OWF problem and when does it solve the og OWF problem 
- [ ] What is the variable , fixed speed model in pumps 
- [ ] Can water flow backward in our case ? 
- Yes
- [ ] What is the difference between $f_{m,n}$ and $d_{m,n}$
- [ ] Understand the relation between pump pressure gain and water flow created in that pipe 
- [ ] Why is $hmt<=h'm$ in reservoirs and $hmt<=ltm$ in tanks instead of being equal ? 
- [ ] what does it mean for the tank to be filling as it seems that in our upper level model the tank is always connected to its parent node .
- [ ] What is the duration of T's being considered ? in this scenario 
- $\delta$=1hour 
- T=12
- electricity Prices were obtained from a day ahead average prices
- [ ] how is tank connections modeled here
- tank node $m \in M_b$ has no consumption demand $\forall m\in M_b \implies DE_m=0$
- We described $DE_n$ as the water per sec the node should get in the upper model, in the lower model it could be that it should at least get this much water (at every time frame) or (on avg)
- at every time frame => $d_m^t>=DE_m$ 
- [ ] Why is there not a mention of minimum pressure to be maintained in the paper itself ? 
- [ ] Is there a minimum water demand to be met mentioned by the paper 
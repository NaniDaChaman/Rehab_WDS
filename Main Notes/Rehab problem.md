2025-03-0517:46

Tags :

Status : 

## Task :
- [ ] Create head reference constraint 
- [ ] Create min pressure constraint
- [ ] Model the previously built network  
- [ ] Model benefit of upgrade ? 
- [ ] model elevation of nodes !!!
## Primary and Secondary Network 

### Primary Supply Hours 
- The primary network refers to the portion of the water distribution system that transports water from the source (such as a reservoir or treatment plant) to the storage tanks within the system.
- Primary supply hours refer to the time period during which the primary network is actively delivering water from the source (such as a reservoir or treatment plant) to the storage tanks.

### Secondary Supply Hours 
- The secondary network typically connects the storage tanks to the demand nodes (such as households or commercial establishments) that require water. This part of the network distributes water from the tanks to the end-users.
- Secondary supply hours refer to the time period during which the secondary network distributes water from the storage tanks to the demand nodes (e.g., households, commercial areas).

- Coordination => Effective coordination between primary and secondary supply hours is vital. The primary network must replenish the storage tanks during primary supply hours to ensure that the secondary network can meet the demands during its operational hours.

## Variable Introductions

### Tank variable
- h_bar due tank

## Objective Function 
## pipe cost 
- Cost per unit length : $C_j$
- Number of pipes : $NP$
- Capital cost of pipes = $\sum_{i=1}^{NL}\sum_{j=1}^{NP}C_j*D_j*l_{i,j}$
### tank cost 
Number of tanks in level : $NE$
Number of nodes : $NN$
Constants
	1. Base cost of kth row : $B_k$
	2. Unit cost of kth row : $UN_k$
	3. Upper limit : $UP_k$
	4. Lower Limit : $LO_k$
	5. Min tank height allowed : $T_{min}$
	6. Max tank height allowed : $T_{max}$
	-  Tank Capital cost = $\sum_{n=1}^{NN}\sum_{k=1}^{NE}(B_k+UN_k(dmax_n-LO_k))e_{n,k})$

### pump capex cost 
- $P_{i}=\frac{\rho*g*g'}{\eta}$ => power of the pump
- Pump Capital cost = $\sum_{i=1}^{NL}CP*P_{i}$

### pump opex cost 
-  $C=\frac{\delta*\rho*g*g'}{\eta}$ => electricity consumption coefficient
- $\pi_{t}$ => avg elec prices at time t
- $\sum_{t=1}^T\sum_{i=1}^{NL}\pi_{t}*C*dp_{i,t}$
### Length Constraints
- $\sum_{j=1}^{NP}l_{i,j}=L_i$
- $l_{i,j}*f_i=lp_{i,j}$
### head loss calculation
- $hl_{i,t}=\sum_{j=1}^{NP}{HL_{i,j,t}^p*l^p_{i,j}+HL_{i,j,t}^s*(l_{i,j}-l^p_{i,j})}-g*x_{i,t}$ 
- $HL^p_{i,j,t}=10.68*\frac{{\frac{FL^p_i}{R_j}}^{1.852}}{D_j^{4.87}}$
- $FL^s_i=FL^p_i\frac{PH}{SH}$

### min pressure requirments
- $h_{n,t}$ => head at a given node at a time
- $P_n$ => minimum pressure requirements
- $P_n<=h_{n,t}$

### head node calculation
- $h_{n,t}$ => head at a given node at t time
- $hl_{i,t}$ => headloss across link i at time t
- $he_{m,i,t}$=> effective head caused by node m across link i 
- $es_{m,i}$=>if source of water for link i is its immediate upstream node m
- $\alpha_{m,t}$=> tank at node m is connected at time t
- $f_i$=>if the link i is part of the primary and secondary network 
- $s_{n,m}$ => source for m is n
- $hbar_{m,t}$ => head caused by the tank 
#### equations : 
- $h_{n,t}=he_{m,i,t}-hl_{i,t}$
- $he_{m,i,t}=hbar_{m,t}*\alpha_{m,t} + h_{m,t}*(1-\alpha_{m,t})$
- $\alpha_{m,t}<=\sum_{i=1}^{NL} es_{m,i}$
- $es_{m,i}=s_{m,m}*(1-f_i)$
- $\forall n$ where $m \in P_n$ and $i \in O_m$

###  calculating head at tank and its operation
-  Max tank height allowed : $T_{max}$
-  $hbar_{m,t}$ => head caused by the tank 
- $tl_{m,t}$ => tank level at node m at time t
- $\beta_{m,t}$ => if tank is filling at node m at t

#### equations : 
- When $\beta_m^t=0$, $h_m^t<=l_m^t$, when $\beta_m^t$=1, $h_m^t>=l_m^t$
- $-M*(1-\beta_{m,t})+(Tmax+E_n)<=hbar_{t,m}<=(tl_{m,t}+E_n)+M*\beta_{m,t}$
- 
### Demand that tank serves  
- demand served by tank at time t  : $d_{n,t}$
-  Water demand of node : $DE_n$
-  $s_{n,m}$ => source for m is n
- $D_n$ => descands of node n
- $I_n$ => link with start node as n
- $FLp_{i,t},FLs_{i,t}$ => Flow in each link if its part of the primary and secondary network
-  $f_i$=>if the link i is part of the primary and secondary network 
- $\alpha_{n,t}$=> tank at node n is connected at time t
- $dmax_n$ =>max demand served by tank at node n
#### equations
- $\sum_{m\in D_n}{s_{n,m}*DE_n} - FLp_{i,t}*f_i - FLs_{i,t}*(1-f_i) =d_n$ for $i\in I_n$
- $-M*\alpha_{n,t}<=d_{n,t}<=M*\alpha_{n,t}$
- $-M*\beta_{n,t}<=d_{n,t}<=M*(1-\beta_{n,t})$
- $x$

### Change in tank level due to change in water flow in tank
- $tl_{m,t}$ => tank level at node m at time t
- demand served by tank at time t  : $d_{n,t}$
- $UP_k$ the upper limit of capacity of a tank k
- $Tmax$ : max tank level 
- $Tmin$ : minimum tank level allowed
- $\delta$: time period to control pump operations
- $y_{m,k,t}$ : demand served by node m having tank k at time t

#### Idea behind it : 
- The first and last tank level should be the same 
- the first tank level is given
- change is tank level is depended on the water flowing out of it, and the pressure outside the tank 
- if pressure is low enough tank will empty
- if pressure is high tank enough tank will fill
- we can disconnect it to do neither 
- this can help tank fill up when demand is low and empty when demand is high

#### equations 
- $l_m^t=l_m^{t-1}-\sum_{k=1}^{NE}\frac{Tmax*y_{m,k,t}*\delta}{UP_k}$
- $y_{m,k,t}=d_{m,t}*e_{m,k}$ Linearized 
- $Tmin<=l_m^t<=Tmax$
- $l_m^0=l_m^T=Tmin$
### Pump Constraints
- $ph_{i,t} =g*x_{i,t}$ , head at a pump head
- $dp_{i,t}$ => water flowing through pump when its on
-  $FLp_{i,t},FLs_{i,t}$ => Flow in each link if its part of the primary and secondary network
- $x_{i,t}$ => weather pump is on  or off  at link i i at time t
- $pe_i$ => pump is at link i or not 

### equations 
- $ph_{i,t}=g*x_{i,t}$
- $-M*(1-x_{i,t})<=FLp_{i,t}*f_i+FLs_{i,t}*(1-f_i) -dp_{i,t}<=M*(1-x_{i,t})$
- $-M*x_{i,t}<=dp_{i,t}<=M*x_{i,t}$
- $x_{i,t}<=pe_i$

## Error in our model
## Rehab problem Questions
 - [x] 1 tank equations give us an infeasible output ? 
 - Related to tank level problem as it happened when i introduced it b into the equations 
 - in feasible means that it is also overdetermined : ???
 - what are the parameters involved in the tank feasibility equations 
 $\sum_{k=1}^{NE}e_{n,k}<=1$
 - if one of my equations make it so that at least 2 of these variable=1 then our model is infeasible 
 - Value of d and y could be determined by the tank level equation => e
 - ==Our model equations for the model were wrong !!!==
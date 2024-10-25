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

## Water Network Modeling 
- G = (M,P) this our water network
## Nodes (M)
- $M=Mr+Mb+Mc$ and $Mr\cap Mc\cap Mb=\phi$
- Mr are the nodes hosting reservoirs
- Mb are nodes hosting tanks
- Mc are nodes hosting Consumers
- $d_m^t=$ Water flowing into m at t
- $m\in Mr \rightarrow d_m^t>=0$
-  $m\in Mb \rightarrow d_m^t=0$
-  $m\in Mc \rightarrow d_m^t<=0$
## Pipes (P)
- $|P| =|M|*|M|$, all edges are directional (m,n)
- $(m,n)\in P \rightarrow (n,m) \not\in P$
- $P_a=$ edges which denote the ideal pump
- $\overline P=$ lossy pipes
- $d_{mn}^t =$ water flow from node m,n
- $d_{mn}^t >=0$ if water flows from  m to n
- $d_{m}^t = \sum_{k:(m,k)\in P}d_{m,k}^t - \sum_{k:(k,m)\in P}d_{k,m}^t$

## Pressure (h)
- $h_m^t=$ pressure at that node at time t
- In static analysis = $h_m^t=E_n+t_n$ 
- Minimum pressure is a common minimum value for all nodes n which our pressure should be above for water to flow $h_m^t>= \underline{h_m}$
- Darcy-Wesbach equation : gives us pipe loss due to friction
	1. $h_m^t-h_n^t=c_{mn}*sign(d_{m,n}^t)*(d_{m,n}^t)^2$
	2. $c_{m,n}=\frac{l_{m,n}*f_{m,n}}{4*\pi^2*r_{m,n}^5*g}$
	3. $f_{m,n}$ is a constant
	4. Big M trick : 
		1. $-M*(1-x_{m,n}^t)<=d_{m,n}^t<=M*x_{m,n}^t$
		2. $-M*(1-x_{m,n}^t)<=h_m^t-h_n^t-c_{m,n}*(d_{m,n}^t)^2<=M*x_{m,n}^t$
		3. $-M*(1-x_{m,n}^t)<=h_m^t-h_n^t+c_{m,n}*(d_{m,n}^t)^2<=M*x_{m,n}^t$
		4. $x_{m,n}^t =1$ if water flows from m->n else 0
## Problem Formulation
- [ ] What is our decision variables, objectives and constraints ? 
- [ ] What is the time duration $\delta$ ? 
- We discretize our continuous system into time chunks t=0,1,2,...T
- 0 represents the start of our system while T represents the end
- All of these chunks are of a constant duration $\delta$ 


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
using JuMP, GLPK ,HiGHS
using ExcelReaders
using Base.Iterators
using DataStructures
using LinearAlgebra
using Dates
using Random

#-p[m,n]*x[m,n,t]=h[m,t]-h[n,t]
#youll need change in demand at every hour and average demand for each node to be able to do anything 

mutable struct demand_node
	node_id::Int
	elevation::Float64
	water_demand::Float64
	min_pressure::Float64
end

mutable struct reference_node
	node_id::Int
	elevation::Float64
	water_demand::Float64
	pressure::Float64
end

mutable struct pipe
	pipe_id::Int
	diameter::Float64
	roughness::Float64
	cost::Float64
	length::Float64
end

mutable struct link
	link_id::Int
	start_node::Int
	end_node::Int
	length::Float64
end

mutable struct tank
    tank_id::Int
    LO::Float64
    UP::Float64
    B::Float64
    UN::Float64
end

mutable struct time_period
    time_id::Int
    start_time::Time
    end_tme::Time
    demand_factor::Float64
    elec_price::Float64
end


function check_feasibility(model)
	optimize!(model)
	if termination_status(model)==INFEASIBLE::TerminationStatusCode
		println("model is INFEASIBLE")
		println(solution_summary(model,verbose=true))
	else
		println("Model status : $(termination_status(model))")
	end
end

#reading source nodes 
data = readxl("input/Sample_input_ESR_Pump.xls", "General!D15:D18")
ref_node=reference_node(data[1],data[3],0,data[4])
water_density=1
g=9.8
println(ref_node)

#reading general data 
data = readxl("input/Sample_input_ESR_Pump.xls", "General!D8:D18")
P_min::Float64=data[1]
R_def::Float64=data[2]
PH::Int=data[7]
println("Min pressure : $(P_min), Default Roughness : $(R_def), Primary Hours : $(PH)")


# reading nodes
data = readxl("input/Sample_input_ESR_Pump.xls", "General!A24:E32")
#println(data)
N=11
number_of_rows = length(data[:,1])
#println(number_of_rows)
def_link=link(0,0,0,0)
link_structure = fill(def_link,N, N)
demand_nodes=Array{demand_node}(undef,N)
for i = 1:number_of_rows 
    d1=demand_node(data[i,1],data[i,3],data[i,4],P_min)
    demand_nodes[d1.node_id]=d1
end 
demand_nodes[ref_node.node_id]=demand_node(ref_node.node_id,ref_node.elevation,0,0)
println(demand_nodes)


#readin pipes
data = readxl("input/Sample_input_ESR_Pump.xls", "General!A52:C64")
NP=13#number of pipes=number of rows
possible_pipes= Array{pipe}(undef,NP)
for i = 1:NP 
    if data[i,2]>0
        possible_pipes[i]=pipe(i,data[i,1],data[i,2],data[i,3],0)
    else 
        possible_pipes[i]=pipe(i,data[i,1],R_def,data[i,3],0)
    end
end
println(possible_pipes)


# reading links
data = readxl("input/Sample_input_ESR_Pump.xls", "General!A38:G46")
NL=9#number of links = number of rows
pipe_structure =zeros(N,N,NP) #tells us about the existing pipes
adj_mat=Matrix{Float64}(I,N,N)
link_list = Array{link}(undef,NL)
for i = 1:NL
	link1 = link(data[i,1],data[i,2],data[i,3],data[i,4])
    link_structure[Int(data[i,2]),Int(data[i,3])]=link1
    link_list[i]=link1
	diameter_present=data[i,5]
	if diameter_present>0
		for j =1:NP
			if possible_pipes[j].diameter==diameter_present
                println("link $(i) from node $([Int(data[i,2]),Int(data[i,3])]) : $(diameter_present)")
				pipe_structure[Int(data[i,2]),Int(data[i,3]),j] = data[i,4]
			end
		end
	end
	#println(diameter_present)
end 
println("link structure is  : $(link_structure[8,:])")
println("Link List is : $(link_list[5])")
#println("pipe structure is : $(pipe_structure[3,7,:])")



#tank general data 
data = readxl("input/Sample_input_ESR_Pump.xls", "General!D69:D75")
SH::Float64=data[2]
tank_cp::Float64=data[3]
T_max::Float64=data[4]
T_min::Float64=0
#tank data 
data = readxl("input/Sample_input_ESR_Pump.xls", "General!A80:D89")
NE=10
possible_tanks=Array{tank}(undef,NE)
for i =1:NE
    possible_tanks[i]=tank(i,data[i,1],data[i,2],data[i,3],data[i,4])
end
println("Possible tanks : $(possible_tanks)")

# pump general data
eff::Float64=0
CP::Float64=0
EP::Float64=0
DF::Float64=0
INFR::Float64=0
INTR::Float64=0
PP_min::Float64=-Inf
PP_max::Float64=Inf
Y::Float64=1

data = readxl("input/Sample_input_ESR_Pump.xls", "General!D93:D101")
PP_min=data[2]
eff=data[3]
CP=data[4]
EP=data[5]
Y=data[6]
DF=data[7]
INFR=data[8]

println("\nMin Pump Size : $(PP_min), Pump Efficiency : $(eff), Capital Cost : $(CP)\n
Opex : $(EP), Pump Lifetime : $(Y), Discount Rate : $(DF), Inflation Rate : $(INFR)")

data = readxl("input/Sample_input_ESR_Pump.xls", "General!A142:E149")
#println(data)
T=length(data[:,1])
time_slices=Array{time_period}(undef,T)
time_delta=3*60*60
delta=Time(3,0)
for i=1:T
    println(data[i,2],"\t",data[i,3])
    td1 = time_period(Int(data[i,1]),data[i,2],data[i,3],data[i,4],data[1,5])
    time_slices[i]=td1
end

# data = readxl("input/elc_prices.xlsx", "elec_prices!Q2:Q9")
# for i=1:T
#     time_slices[i].elec_price=data[i]
# end

owf_model=Model(HiGHS.Optimizer)
set_silent(owf_model)


# variables
# pressure at each and every node 
@variable(owf_model,h[1:N*T])
@variable(owf_model,x[1:NL*T],Bin)
@variable(owf_model,d_bar[1:NL*T])
@variable(owf_model,a[1:N*T],Bin)
@variable(owf_model,d[1:N*T],Bin)
@variable(owf_model,tl[1:N*T])
@variable(owf_model,b[1:N*T],Bin)
@variable(owf_model,tank_a[1:N*T],Bin)
@variable(owf_model,h_bar[1:N*T])
@variable(owf_model,e[1:N*NE],Bin)
@variable(owf_model,s[1:N*N],Bin)

#objective function
ph_star=[6.776263578034403e-21,
3.5970678814571146e-14,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0]
#electricity prices 
pie=[t.elec_price for t in time_slices]
#electricity consumption coefficient 
elc_c=(water_density*g/eff)*ph_star
#elc consumed per time
# I_elc=Matrix{Float64}(I,NL,NL)*pie
# elc_t=[]
# for t=1:T
#     global elc_t=[elc_t;I_elc]
# end
# elc_t=elc_t*elc_c
println((elc_c*pie')[2,:])
K=reshape(pie*elc_c',NL*T,1)
println(K[T+1:2T])
#K*.d_bar
@objective(owf_model,Min,sum(K.*d_bar))

# pump pressure -flow requirments 
H_r=ref_node.pressure
println("Source node pressure : $(H_r)")
P_n=zeros(Float64,N)
E_n=zeros(Float64,N)
for i=1:N
	if isassigned(demand_nodes,i)
		P_n[i]=demand_nodes[i].min_pressure
		E_n[i]=demand_nodes[i].elevation
	end
end

M_p=2*maximum(P_n)

for n=1:N
    @constraint(owf_model,P_n[n].<=h[(n-1)*T+1:n*T])
end

println("\nMin pressure feasibility\n")
check_feasibility(owf_model)

flow_p=[ 17.33333333333334,
12.0,
11.666666666666666,
22.66666666666666,
20.399999999999995,
 5.333333333333333,
 4.8,
14.000000000000002,
12.600000000000001]

flow_s=[  8.66666666666667,
6.0,
5.833333333333333,
11.33333333333333,
10.199999999999998,
2.6666666666666665,
2.4,
7.000000000000001,
6.300000000000001]

f=[1.0,
0.0,
1.0,
0.0,
1.0,
0.0,
1.0,
0.0,
1.0]

pe=[1
1
0.0
0.0
0.0
0.0
0.0
0.0
0.0]

df_t=[t.demand_factor for t in time_slices]
flow_st=reshape(df_t*flow_s',NL*T,1)
flow_pt=reshape(df_t*flow_p',NL*T,1)

#pressure at pumps requirments
M_d=2*maximum([flow_pt;flow_st])
for i =1:NL
    m=link_list[i].end_node
    n=link_list[i].start_node
    #@constraint(owf_model,pe[i]*(h[(n-1)*T+1:n*T]-h[(m-1)*T+1:m*T]).==ph_star[i]*x[(i-1)*T+1:i*T])
    @constraint(owf_model,-M_d*(ones(T)-x[(i-1)*T+1:i*T]).+-M_d*pe[i].<=flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i])-d_bar[(i-1)*T+1:i*T])
    @constraint(owf_model,M_d*(ones(T)-x[(i-1)*T+1:i*T]).+M_d*pe[i].>=flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i])-d_bar[(i-1)*T+1:i*T])
    #-M*(1-x[(i-1)*T+1:i*T]).<=d[(i-1)*T+1:i*T]-d_bar[(i-1)*T+1:i*T]
   @constraint(owf_model,zeros(T).<=d_bar[(i-1)*T+1:i*T])
   @constraint(owf_model,M_d.*x[(i-1)*T+1:i*T].>=d_bar[(i-1)*T+1:i*T]) #modeled in case we need it 
end


println("\nPump pressure feasibility\n")
check_feasibility(owf_model)

#headloss due to pipes
#cost vector ,roughness vector and Daimeter vectors also used in ilp's objective function
C=zeros(Float64,NP)
R=zeros(Float64,NP)
D=zeros(Float64,NP)

for i =1:NP
	C[i]=possible_pipes[i].cost
	R[i]=possible_pipes[i].roughness
	D[i]=possible_pipes[i].diameter
end

# creating headloss per unit length
num_pt=10.68*flow_pt.^1.852
num_st=10.68*flow_st.^1.852
denom=1 ./(R.*D.^4.87)
println("size of numerator : $(size(num_pt)) \t size of denominator : $(size(denom))")
HL_pt=num_pt*denom'
HL_st=num_st*denom'

#creating l_star for our model
L=zeros(Float64,NL)
for i =1:NL
	L[i]=link_list[i].length
end
l=[]
l_s=zeros(Float64,NL*NP)
for i=1:NL
	one_hot_vector = zeros(Int64, NP)
	one_hot_index = rand(1:NP)
	one_hot_vector[one_hot_index] = L[i]
	global l=[l;one_hot_vector]
end
for i =1:NL
	l_s[(i-1)*NP+1:NP*i]=f[i]*l[(i-1)*NP+1:NP*i]
end
#println(l_s)
for i=1:NL
	m=link_list[i].end_node
    n=link_list[i].start_node
	for t=1:T
		@constraint(owf_model,h[(n-1)*T+t]-h[(m-1)*T+t].==sum(HL_pt[T*(i-1)+t,:].*l[(i-1)*NP+1:NP*i]+HL_st[T*(i-1)+t,:].*(l[(i-1)*NP+1:NP*i]-l_s[(i-1)*NP+1:NP*i]))+ph_star[i]*x[(i-1)*T+t])
	end
end

println("\nlink pressure loss feasibility\n")
check_feasibility(owf_model)


#pressure due to reseviors 
M_h=2*H_r
for i=1:NL
	m=link_list[i].end_node
	n=link_list[i].start_node
	if ref_node.node_id==m
		@constraint(owf_model,h[(m-1)*T+1:m*T].<=H_r.+M_h*(ones(T)-a[(n-1)*T+1:n*T]))
		@constraint(owf_model,flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i]).<=M_d*a[(n-1)*T+1:n*T])
		#flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i]).<=M_d*a[]
		break
	end
end

println("\nresevior pressure-flow feasibility\n")
check_feasibility(owf_model)

#water flow and pressure due to tanks
#demand of water nodes 
DE=zeros(N)
for n = 1:N
	if isassigned(demand_nodes,n)
		DE[n]=demand_nodes[n].water_demand
	end
end
println(DE)

#descandants of all nodes 
function D_n(start_node,adj_mat)
	s=Stack{Int}()
	push!(s,start_node)
	visited=Set{Int}()
	while !isempty(s)
		c_n=pop!(s)
		#println("Current node : $(c_n)\n")
		push!(visited,c_n)
		adj_links=adj_mat[c_n,:]
		for i=1:N
			if adj_links[i].link_id!=0
				#print("Links : $(adj_links[i]) , end node : $(i)\t")
				n=adj_links[i].end_node# should be equal to I
				#println("i : $(i) == n : $(n)")
				if !(n in visited)
					push!(s,i)
					#path[c_n]=adj_mat[i,c_n]
				end
			end
		end
		#println("\nNodes Visited : $(visited)\n\n")
	end
	return visited
end
#parent of a node
function Parent(node,link_list)
	for i=1:NL
		lk=link_list[i]
		if lk.end_node==node
			return lk
		end
	end
end

d=[0.0,
0.0,
15.7,
12.1,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0,
0.0
]

d=reshape(df_t*d',N*T,1)
println(d[(3-1)*T+1:3*T])

t=[0.0,
0.0,
25.0,
25.0,
 0.0,
 0.0,
 0.0,
0.0,
0.0,
 0.0,
 0.0]

#e=[]
# s=[]

# for n=1:N
# 	# one_hot_vector = zeros(Int64, NE)
# 	# one_hot_index = rand(1:NE)
# 	# one_hot_vector[one_hot_index] = t[n]/maximum([t[n],1])
# 	# global e=[e;one_hot_vector]
# 	one_hot_vector = zeros(Int64, N)
# 	one_hot_index=n
# 	one_hot_vector[one_hot_index] = t[n]/maximum([t[n],1])
# 	global s=[s;one_hot_vector]
# end

# println(s)

# parameters for possible tanks
B=zeros(Float64,NE)
UN =zeros(Float64,NE)
LO =zeros(Float64,NE)
UP =zeros(Float64,NE) #not needed

for i =1:NE
	B[i]=possible_tanks[i].B
	UN[i]=possible_tanks[i].UN
	LO[i]=possible_tanks[i].LO
	UP[i]=possible_tanks[i].UP
end

# A=zeros(N)

# for k=1:N
# #A[i]*t[i].==sum(UP.*e[(i-1)*N+1:i*N])
# 	A[k]=sum((UP./T_max).*e[(k-1)*NE+1:k*NE])
# end

for m=1:N
	P=Parent(m,link_list)
	for t=1:(T-1)
		if P!=nothing
			i=P.link_id
			@constraint(owf_model,tl[(m-1)*T+t+1].==tl[(m-1)*T+t]-sum(T_max ./UP.*e[(m-1)*NE+1:m*NE])*(d[(m-1)*T+t]-flow_pt[(i-1)*T+t]*f[i]-flow_st[(i-1)*T+t]*(1-f[i]))*time_delta)
		end
	end
	@constraint(owf_model,tl[(m-1)*T+1:m*T].<=T_max)
	@constraint(owf_model,tl[(m-1)*T+1:m*T].>=T_min)
	@constraint(owf_model,tl[(m-1)*T+1].==tl[(m)*T])
	@constraint(owf_model,tl[(m-1)*T+1].==T_min)
end
#l[(m-1)*T+t+1].==l[(m-1)*T+t]-sum(TT_max*1./UP.*e[(i-1)*NE+1:i*NE])*d[(m-1)*T+t]*delta

println("\nTank Level feasibility\n")
check_feasibility(owf_model)
# println(owf_model)
# throw(DomainError("Stoped by user"))
#tank operation with head and flow variables



for m=1:N
	P=Parent(m,link_list)
	#need to use s variable here : tank_a[(n-1)*T+1:T].<=s[(n-1)*T+n]
	@constraint(owf_model,tank_a[(m-1)*T+1:T].<=s[(m-1)*T+m]) # this constraint is making the whole thing in feasible
#pressure constraints 1 : -M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T]
	@constraint(owf_model,-M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T])
#pressure constraint 2 : M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).>=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T]
	@constraint(owf_model,M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).>=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T])
#flow constraint , tank connection 1: -M_d*tank_a[(m-1)*T+1:m*T].<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
	if P!=nothing
		i=P.link_id	
		@constraint(owf_model,-M_d*tank_a[(m-1)*T+1:m*T].<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T])
	#flow constraint , tank connection 2: M_d*tank_a[(m-1)*T+1:m*T].>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(owf_model,M_d*tank_a[(m-1)*T+1:m*T].>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T])
	#flow constraint , tank filling 1 : M_d*(ones(T)-b[(m-1)*T+1:m*T]).>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(owf_model,M_d*(ones(T)-b[(m-1)*T+1:m*T]).>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T])
	#flow constraint , tank filling 2 : -M_d*(b[(m-1)*T+1:m*T]).<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(owf_model,-M_d*(b[(m-1)*T+1:m*T]).<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T])
	end
#pressure constraint , tank filling 1 : T_max.-M_h*(ones(T)-b[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]
	@constraint(owf_model,T_max.-M_h*(ones(T)-b[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T])
#pressure constraint , tank filling 2 : l[(m-1)*T+1:m*T]+M_h*(b[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]
	@constraint(owf_model,tl[(m-1)*T+1:m*T]+M_h*(b[(m-1)*T+1:m*T]).+E_n[m].<=h_bar[(m-1)*T+1:m*T])
end

println("\ntank operation feasibility\n")
check_feasibility(owf_model)
# println(owf_model)
# throw(DomainError("Stoped by user"))
using JuMP, GLPK ,HiGHS
using ExcelReaders
using Base.Iterators
using DataStructures
using LinearAlgebra
using Dates
using Random

#-p[m,n]*x[m,n,t]=h[m,t]-h[n,t]
#youll need change in demand at every hour and average demand for each node to be able to do anything 
h_star=nothing
hbar_star=nothing
hs_star=nothing
tl_star=nothing
l_star=nothing 
lp_star=nothing
f_star=nothing
function Parent(node,link_list)
	for i=1:NL
		lk=link_list[i]
		if lk.end_node==node
			return lk
		end
	end
end

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

function write_model(model,name)
    dir="output"
    fpath=joinpath(dir,name)
    write_to_file(model,fpath)
end

function check_feasibility(model)
	optimize!(model)
	if termination_status(model)==INFEASIBLE::TerminationStatusCode
		println("model is INFEASIBLE")
		println(solution_summary(model,verbose=true))
        write_model(model,"rehab_problem.lp")
        throw(DomainError("Stoped by user"))
    elseif (termination_status(model)==OPTIMAL::TerminationStatusCode)
        global h_star=value.(h_t)
        global hbar_star=value.(hbar_t)
        global hs_star=value.(hs_t)
        global tl_star=value.(tl_t)
        global l_star=value.(l)
        global lp_star=value.(l_p)
        global f_star=value.(f)
        global l_star=[l_star[(i-1)*NP+1:i*NP] for i=1:NL]
        global lp_star=[lp_star[(i-1)*NP+1:i*NP] for i=1:NL]
        global tl_star=[tl_star[(n-1)*T+1:n*T] for n=1:N]
        global h_star=[h_star[(n-1)*T+1:n*T] for n=1:N]
        println("Model status : $(termination_status(model))")
	else
		println("Model status : $(termination_status(model))")
        write_model(model,"rehab_problem.lp")
        throw(DomainError("Stoped by user"))
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

t=1
index_n=[((n-1)*T+t) for n=1:N]
index_m=[]
index_i=[]
for n=1:N
    P=Parent(n,link_list)
    if P!=nothing
        m=P.start_node
        global index_m=[index_m;m]
        i=P.link_id
        global index_i=[index_i;i]
    end
end
    
eff_index=[(index_m[n]-1)*NL+(index_i[n]-1)*T+t for n=1:NL]



rehab_problem=Model(HiGHS.Optimizer)
set_silent(rehab_problem)


# variables

# pipe variables 
#none from the second model / makes sense as pipe opt does not change with time
@variable(rehab_problem,l[1:NP*NL].>=0)#possible lengths of pipes
@variable(rehab_problem,l_p[1:NP*NL].>=0)#possible lengths of primary pipes
@variable(rehab_problem,f[1:NL],Bin)#links in primary
@variable(rehab_problem,hl_t[1:NL*T].>=0)#headloss across links for every slice of time

#node variables
@variable(rehab_problem,h_t[1:N*T].>=0) # head at each node at each given time slice
@variable(rehab_problem,hs_t[1:N*NL*T]) #effective head given to the outgoing link at each given time slice 
@variable(rehab_problem,d[1:N*T]) #water that flows out of tank at node n at time slice t

#tank variables 
@variable(rehab_problem,z[1:NE*N].>=0) # max demand being served by tank k at node N
@variable(rehab_problem,e[1:NE*N],Bin) #if tank k is at node n 
@variable(rehab_problem,s[1:N*N],Bin) #if node m is served by tank at node n
@variable(rehab_problem,tl_t[1:N*T].>=0)#water level at tank at node n at T
@variable(rehab_problem,t[1:N].>=0) #height of tank at node N
@variable(rehab_problem,es[1:NL*N],Bin) #if tank at node n serves its immediated downstream link I
@variable(rehab_problem,b_t[1:N*T],Bin) # if the tank is filling 1, or no 0 at t
@variable(rehab_problem,atank_t[1:N*T],Bin) # if the tank is connected 1, or no 0 at t
@variable(rehab_problem,hbar_t[1:N*T]) # effective head if tank is connected at t at node n
@variable(rehab_problem,dmax[1:N]) # max demand served by node n 

#pump variables
@variable(rehab_problem,p[1:NL].>=0) #power of pump at link I
@variable(rehab_problem,p_p[1:NL].>=0) #power of pump at primary networks link I
@variable(rehab_problem,p_s[1:NL].>=0) #power of pump at secondary networks link I

@variable(rehab_problem,pe[1:NL],Bin) #if pump is installed at link I
@variable(rehab_problem,ph_t[1:NL*T].>=0) #head provided by the pump at link i
@variable(rehab_problem,x_t[1:NL*T],Bin) #pump on of at every link <=pe[i]
@variable(rehab_problem,dp_t[1:NL*T].>=0) #flow in pump if pump is connected 
@variable(rehab_problem,y_t[1:N*NE*T]) #auxiliray variable to fill the tank

# resevior connection
@variable(rehab_problem,a[1:N*T],Bin)

C=zeros(Float64,NP)
R=zeros(Float64,NP)
D=zeros(Float64,NP)

for i =1:NP
	C[i]=possible_pipes[i].cost
	R[i]=possible_pipes[i].roughness
	D[i]=possible_pipes[i].diameter/1000
end

println("\n Cost Vector : $(C)\n")

I_c=Matrix{Float64}(I,NP,NP)
C_l=[]
for i=1:NL
    global C_l=[C_l;I_c]
end
C_l=C_l*C

println("\nCost per length : $(C_l[1:NP+NP])\n")
# the vector is formed as rows of all possible pipes(total = NP) for each link(total = NL)
println("Size of cost vector is $(size(C_l))")

#objective = Sum(C_l.*l)

#considering capex cost of pumps
# p[i] is the variable which determines the power of pump at each link
#we need a CP [1:NL] vector which is a parameter for capital cost for each pump at each link 
# objective = Sum(CP*p)

#considering opex of pumps
#p_p[i] at i th link the power if i is part of the primary networks
#objective = Sum(PH*EF*DF*p_p)

# forming Base cost, Unit cost, upper limit and lower limit tanks 
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

I_n=Matrix{Float64}(I,NE,NE)
T_n=[]
#UN_n=[]
#LO_n=[]
#UP_n=[]
for i=1:N
    global T_n=[T_n;I_n]
	#global UN_n=[UN_n;I_n]
	#global LO_n=[LO_n;I_n]
	#global UP_n=[UP_n;I_n]
end
B_n=T_n*B
UN_n=T_n*UN
LO_n=T_n*LO
#UP_n=T_n*UP
CT_n=LO_n.*UN_n + B_n


#pump and elec constants
#pump head
ph=14*ones(NL)
#electricity prices 
pie=[t.elec_price for t in time_slices]
#electricity consumption coefficient 
elc_c=(water_density*g/eff)*ph

println((elc_c*pie')[2,:])
K=reshape(pie*elc_c',NL*T,1)
println(K[T+1:2T])

#objective function
@objective(rehab_problem,Min,sum(C_l.*l)+sum(CP*p)+sum(UN_n.*z+CT_n.*e) + DF*sum(K.*dp_t)) # opex cost will be sum(K.*d_bar) and with 

# water flowing into each node
flow_p=zeros(N)
flow_s=zeros(N)
for i=1:N
	if isassigned(demand_nodes,i)
		flow_p[i]=demand_nodes[i].water_demand*PH/(PH+SH)
		flow_s[i]=demand_nodes[i].water_demand*SH/(PH+SH)
	else
		flow_p[i]=0.0
		flow_s[i]=0.0
	end
	if flow_p[i]==0
		adj_mat[i,:]=[if link.link_id > 0 1 else 0 end for link in link_structure[i,:]]
	end
end
flow_n=flow_p
for i =1:NL
	global flow_n=adj_mat*flow_n+flow_p
end
flow_p=flow_n

flow_n=flow_s
for i =1:NL
	global flow_n=adj_mat*flow_n+flow_s
end
flow_s=flow_n

println("Demand nodes : $(demand_nodes)")
println("Water flowing into each node at Primary networks : $(flow_p)\n Water Flowing into each node at Secondary Networks : $(flow_s)")


#println("Water flowing into each node : $(flow)")

L=zeros(Float64,NL)
FL_p=zeros(Float64,NL)
FL_s=zeros(Float64,NL)
for i =1:NL
	L[i]=link_list[i].length
	FL_p[i]=flow_p[Int(link_list[i].end_node)]/1000
	FL_s[i]=flow_s[Int(link_list[i].end_node)]/1000
end

df_t=[t.demand_factor for t in time_slices]
flow_st=reshape(df_t*FL_s',NL*T,1)
flow_pt=reshape(df_t*FL_p',NL*T,1)



function S_n(start_node,adj_mat,goal_node)
	path=Dict{Int,link}()
	s=Stack{Int}()
	push!(s,start_node)
	while !isempty(s)
		c_n=pop!(s)
		#println("Current node : $(c_n)\n")
		if c_n==goal_node
			return path
		end
		adj_links=adj_mat[:,c_n]
		for i=1:N
			if adj_links[i].link_id!=0
				#print("Links : $(adj_links[i]) , end node : $(i)\t")
				push!(s,i)
				path[c_n]=adj_mat[i,c_n]
			end
		end
		#println("\nCurrent Path : $(path)\n\n")
	end
	return Dict{Int,link}()
end

path=S_n(7,link_structure,8)
#println("\nPath from 8 to 7 is $(path)\n")
#navigating the path
function get_curr_path(start_node,end_node,path)
	curr_node=end_node
	ordered_path=[]
	while curr_node!=start_node
		cur_link=path[curr_node]
		pushfirst!(ordered_path,cur_link)
		curr_node=cur_link.start_node
		#println(curr_node)
	end
	return ordered_path
end

println("Path from 8->7 is : $(get_curr_path(8,7,path))")

for j=1:NL
	@constraint(rehab_problem,sum(l[1+(j-1)*NP:(j)*NP])==L[j])
end

println("\nLength of pipes constraints")
check_feasibility(rehab_problem)

#Headloss per unit length in each link
num_p=10.68*flow_pt.^1.852
num_s=10.68*flow_st.^1.852
denom=1 ./((R.^1.852).*(D.^4.87))
HL_pt=num_p*denom'
HL_st=num_s*denom'
println("Headloss for each possible assignment in the Primary network: $(HL_pt[1,:]) \n its size $(size(HL_pt))")

for i=1:NL
    for t= 1:T
	    @constraint(rehab_problem,hl_t[(i-1)*T+t]==sum(HL_pt[(i-1)*T+t,:].*l_p[1+(i-1)*NP:(i)*NP]+HL_st[(i-1)*T+t,:].*(l-l_p)[1+(i-1)*NP:(i)*NP])-ph[i]*x_t[(i-1)*T+t])
    end
	@constraint(rehab_problem,l_p[1+(i-1)*NP:(i)*NP].<=L[i]*f[i])
	@constraint(rehab_problem,l[1+(i-1)*NP:(i)*NP]-l_p[1+(i-1)*NP:(i)*NP].<=L[i]*(1-f[i]))
	@constraint(rehab_problem,l_p[1+(i-1)*NP:(i)*NP].<=l[1+(i-1)*NP:(i)*NP])
end

println("\nHeadloss across links constraints")
check_feasibility(rehab_problem)
#throw(DomainError)

# for the pressure at each node constraint
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

M_p=6*maximum(H_r)
println("Maximum pressure :$(M_p)")



for n=1:N
	P=Parent(n,link_list)
	if P !== nothing
		m=P.start_node
		i=P.link_id
		println(" starts node $(m) ,link : $(i) ends at $(n)")
        for t=1:T
		    @constraint(rehab_problem,h_t[(n-1)*T+t]==hs_t[(m-1)*NL+(i-1)*T+t]-hl_t[(i-1)*T+t])
            @constraint(rehab_problem,hbar_t[(m-1)*T+t]-M_p*(1-atank_t[(m-1)*T+t])<=hs_t[(m-1)*NL+(i-1)*T+t])
		    @constraint(rehab_problem,hs_t[(m-1)*NL+(i-1)*T+t]<=hbar_t[(m-1)*T+t]+M_p*(1-atank_t[(m-1)*T+t]))
		    @constraint(rehab_problem,h_t[(m-1)*T+t]-M_p*(atank_t[(m-1)*T+t])<=hs_t[(m-1)*NL+(i-1)*T+t])
		    #@constraint(rehab_problem,hs_t[(m-1)*NL+(i-1)*T+t]<=M_p*atank_t[(m-1)*T+t]+h_t[(m-1)*T+t]) #something wrong with the constraint!!!
            @constraint(rehab_problem,es[(m-1)*NL+i]<=s[(m-1)*N+m])
		    @constraint(rehab_problem,es[(m-1)*NL+i]<=1-f[i])
		    @constraint(rehab_problem,es[(m-1)*NL+i]>=s[(m-1)*N+m]-f[i])
            @constraint(rehab_problem,es[(m-1)*NL+i]>=atank_t[(m-1)*T+t])
        end
    end
    @constraint(rehab_problem,P_n[n]+E_n[n] .<=h_t[(n-1)*T+1:n*T]-tl_t[(n-1)*T+1:n*T])
end
ref_id=ref_node.node_id
@constraint(rehab_problem,H_r.==h_t[(ref_id-1)*T+1:ref_id*T])



println("\nnode_tank pressure")
#print(rehab_problem)
check_feasibility(rehab_problem)


throw(DomainError("Stoped by user"))
#println(rehab_problem)

DE=zeros(N)
for n = 1:N
	if isassigned(demand_nodes,n)
		DE[n]=demand_nodes[n].water_demand
	end
end

DE_t=reshape(df_t*DE',N*T,1)
#println(DE)
# tank heirarchy constraints
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

#function to get incoming node
#use Parent

#desc = D_n(10,link_structure)
#println(DE[collect(desc)])
for n=1:N	
	desc=D_n(n,link_structure)
	asc=push!(collect(keys(S_n(n,link_structure,ref_node.node_id))),ref_node.node_id)
	P=Parent(n,link_list)
	#O=setdiff!(union!(desc,asc),n)
	#O_index=collect([(o-1)*N+o for o in O])
	#println("Ascendands and descandants of $(n) are : $(O)")
	for m in setdiff!(desc,Set(n))
		asc=push!(collect(keys(S_n(m,link_structure,ref_node.node_id))),ref_node.node_id)
		# O=setdiff!(union!(desc,asc),n)
		# O_index=collect([(o-1)*N+o for o in O])
		# for o in O_index
		# 	#@constraint(rehab_problem,s[(n-1)*N+m]+s[o]<=1) giving infeasible 
		# end
		 # each node only has 1 tank located in a path
		@constraint(rehab_problem,s[(n-1)*N+n]>=s[(m-1)*N+m]) #if n can't serve itself then its descandants m can't serve itself, should be removed
		@constraint(rehab_problem,s[(n-1)*N+n]>=s[(n-1)*N+m])#if n can serve itself it can serve its descandants m
	end
	desc=collect(desc)
	index_desc=[(n-1)*N+m for m in desc]
	index_asc=[(m-1)*N+n for m in asc]
	if desc !== nothing
        for t=1:T
            
            index_t=[(m-1)*T+t for m in desc]
            # if P!=nothing
            #     i=P.link_id
            #     @constraint(rehab_problem,d[(n-1)*T+t]>=sum(DE_t[index_t].*s[index_desc])-flow_pt[(i-1)*T+t]*f[i]-flow_st[(i-1)*T+t]*(1-f[i]))
            # else
		        @constraint(rehab_problem,d[(n-1)*T+t]==sum(DE_t[index_t].*s[index_desc])) #total demand served by node n is the sum of demands of nodes that relies on it
            #end
        end
    end    
	if asc !==nothing
		@constraint(rehab_problem,sum(s[index_asc])==1) #only 1 Ascending node m can serve n
	end
	if P !==nothing
		i=P.link_id
		@constraint(rehab_problem,f[i]==s[(n-1)*N+n])
	end
end

println("Tank heirarchy feasibility\n")
check_feasibility(rehab_problem)

M_d=2*maximum([flow_pt;flow_st])
#tank level constraintS
for m=1:N
	P=Parent(m,link_list)
	for t=1:(T-1)
		if P!=nothing
			i=P.link_id
            y_index=[(m-1)*NE+(k-1)*T+t for k=1:NE]
			@constraint(rehab_problem,tl_t[(m-1)*T+t+1].==tl_t[(m-1)*T+t]-sum(T_max ./UP.*y_t[y_index]*time_delta))
            #multiplication of e and d variables
            #need to linearize
            for k=1:NE
                 @constraint(rehab_problem,y_t[(m-1)*NE+(k-1)*T+t]<=M_d*(1-e[(m-1)*NE+k])+d[(m-1)*T+t])
                 @constraint(rehab_problem,y_t[(m-1)*NE+(k-1)*T+t]>=-M_d*(1-e[(m-1)*NE+k])+d[(m-1)*T+t])
                 @constraint(rehab_problem,y_t[(m-1)*NE+(k-1)*T+t]<=M_d*e[(m-1)*NE+k])
                 @constraint(rehab_problem,y_t[(m-1)*NE+(k-1)*T+t]>=-M_d*e[(m-1)*NE+k])
            end
        end

	end
	@constraint(rehab_problem,tl_t[(m-1)*T+1:m*T].<=T_max)
	@constraint(rehab_problem,tl_t[(m-1)*T+1:m*T].>=T_min)
	@constraint(rehab_problem,tl_t[(m-1)*T+1].==tl_t[(m)*T])
	@constraint(rehab_problem,tl_t[(m-1)*T+1].==T_min)
end
#l[(m-1)*T+t+1].==l[(m-1)*T+t]-sum(TT_max*1./UP.*e[(i-1)*NE+1:i*NE])*d[(m-1)*T+t]*delta

println("\nTank Level feasibility\n")
check_feasibility(rehab_problem)
#throw(DomainError("Stoped by user"))

#tank operation constraints
for m=1:N
	P=Parent(m,link_list)
#pressure constraints 1 : -M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T]
	#@constraint(rehab_problem,-M_h*(ones(T)-atank_t[(m-1)*T+1:m*T]).<=hbar_t[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T])# dnt need
#pressure constraint 2 : M_h*(ones(T)-tank_a[(m-1)*T+1:m*T]).>=h_bar[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T]
	#@constraint(rehab_problem,M_h*(ones(T)-atank_t[(m-1)*T+1:m*T]).>=hbar_t[(m-1)*T+1:m*T]-h[(m-1)*T+1:m*T])
#flow constraint , tank connection 1: -M_d*tank_a[(m-1)*T+1:m*T].<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
	if P!=nothing
		i=P.link_id	
		@constraint(rehab_problem,-M_d*atank_t[(m-1)*T+1:m*T].<=d[(m-1)*T+1:m*T])
	#flow constraint , tank connection 2: M_d*tank_a[(m-1)*T+1:m*T].>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(rehab_problem,M_d*atank_t[(m-1)*T+1:m*T].>=d[(m-1)*T+1:m*T])
	#flow constraint , tank filling 1 : M_d*(ones(T)-b[(m-1)*T+1:m*T]).>=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(rehab_problem,M_d*(ones(T)-b_t[(m-1)*T+1:m*T]).>=d[(m-1)*T+1:m*T])
	#flow constraint , tank filling 2 : -M_d*(b[(m-1)*T+1:m*T]).<=flow_pt[(i-1)*T+1:i*T]*f[i]+flow_st[(i-1)*T+1:i*T]*(1-f[i])-d[(m-1)*T+1:m*T]
		@constraint(rehab_problem,-M_d*(b_t[(m-1)*T+1:m*T]).<=d[(m-1)*T+1:m*T])
	end
#pressure constraint , tank filling 1 : T_max.-M_h*(ones(T)-b[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]
	@constraint(rehab_problem,(T_max+E_n[m]).-M_p*(ones(T)-b_t[(m-1)*T+1:m*T]).<=hbar_t[(m-1)*T+1:m*T])
#pressure constraint , tank filling 2 : l[(m-1)*T+1:m*T]+M_h*(b[(m-1)*T+1:m*T]).<=h_bar[(m-1)*T+1:m*T]
	@constraint(rehab_problem,tl_t[(m-1)*T+1:m*T]+M_p*(b_t[(m-1)*T+1:m*T]).+E_n[m].>=hbar_t[(m-1)*T+1:m*T])
    #@constraint(rehab_problem,0 .<=hbar_t[(m-1)*T+1:m*T])
   # @constraint(rehab_problem,hbar_t[(m-1)*T+1:m*T].<=T_max.+E_n[m])
   #@constraint(rehab_problem,hbar_t[(m-1)*T+1:m*T]==tl_t[(m-1)*T+1:m*T].+E_n[m])
    #establishing dmax_t
    #@constraint(rehab_problem,d[(m-1)*T+1:m*T].<=dmax[m]) #since ou problem is a minimisation problem it will take the minimum of all these max values
   # @constraint(rehab_problem,dmax[m]<=sum(e[(m-1)*N:m*N]*))
    
end

println("\ntank operation feasibility\n")
check_feasibility(rehab_problem)

#pump operation
#pressure at pumps requirments

for i =1:NL
    m=link_list[i].end_node
    n=link_list[i].start_node
    @constraint(rehab_problem,ph_t[i].==ph*x_t[i])
    #@constraint(owf_model,pe[i]*(h[(n-1)*T+1:n*T]-h[(m-1)*T+1:m*T]).==ph_star[i]*x[(i-1)*T+1:i*T])
    @constraint(rehab_problem,-M_d*(ones(T)-x_t[(i-1)*T+1:i*T]).<=flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i])-dp_t[(i-1)*T+1:i*T])
    @constraint(rehab_problem,M_d*(ones(T)-x_t[(i-1)*T+1:i*T]).>=flow_pt[(i-1)*T+1:i*T].*f[i]+flow_st[(i-1)*T+1:i*T].*(1-f[i])-dp_t[(i-1)*T+1:i*T])
    @constraint(rehab_problem,-M_d*(x_t[(i-1)*T+1:i*T]).<=dp_t[(i-1)*T+1:i*T])
    @constraint(rehab_problem,M_d*(x_t[(i-1)*T+1:i*T]).>=dp_t[(i-1)*T+1:i*T])
    #-M*(1-x[(i-1)*T+1:i*T]).<=d[(i-1)*T+1:i*T]-d_bar[(i-1)*T+1:i*T]
    @constraint(rehab_problem,x_t[(i-1)*T+1:i*T].<=pe[i])
end

#problem with our model is that if x =0 dp_t is free to be a max value

println("\npump operation feasibility")
check_feasibility(rehab_problem)

#setting up znk
#z[(n-1)*NE+k]=e[(n-1)*NE+k]*d[n]
M_c=2*maximum([t.UP for t in possible_tanks])
println("maximum z value is $(M_c)")

for n =1:N
	for k=1:NE
		@constraint(rehab_problem,z[(n-1)*NE+k]<=M_c*e[(n-1)*NE+k])
		@constraint(rehab_problem,z[(n-1)*NE+k]<=dmax[n])
		@constraint(rehab_problem,z[(n-1)*NE+k]>=-M_c*(1-e[(n-1)*NE+k])+dmax[n])
	end
end

println("\nTank selection feasibility")
check_feasibility(rehab_problem)

#one possible tank for a node
#e[1+(k-1)*NE:(k)*NE]==1
for n =1:N
	@constraint(rehab_problem,sum(e[1+(n-1)*NE:(n)*NE])==1) # making model infeasible
end
println("\n1 Tank feasibility")
check_feasibility(rehab_problem)

t=1
index_n=[((n-1)*T+t) for n=1:N]
index_m=[]
index_i=[]
for n=1:N
    P=Parent(n,link_list)
    if P!=nothing
        m=P.start_node
        global index_m=[index_m;m]
        i=P.link_id
        global index_i=[index_i;i]
    end
end
    
eff_index=[(index_m[n]-1)*NL+(index_i[n]-1)*T+t for n=1:NL]
h_star=value.(h_t)
hbar_star=value.(hbar_t)
hs_star=value.(hs_t)
tl_star=value.(tl_t)
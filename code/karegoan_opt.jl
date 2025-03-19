using JuMP, GLPK ,HiGHS
using ExcelReaders
using Base.Iterators
using DataStructures
using LinearAlgebra

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

function save_network(title_text,filename,demand_nodes,ref_node,tank_loc,tank_type,T_min,T_max,UP,pipe_seg,link_list,l_star,possible_pipes)
    file_text="[TITLE]\n$(title_text)"
    junctions="\n[JUNCTIONS]\n;ID\tElev\tDemand\tPattern\n"
    for de_n in demand_nodes
        if de_n.node_id!=ref_node.node_id
        junctions*="$(de_n.node_id)\t$(de_n.elevation)\t$(de_n.water_demand)\t1\t;\n"
        end
    end
    file_text*=junctions

    resevior="\n[RESERVOIRS]\n;ID\t	Head\tPattern\n"
    resevior*="$(ref_node.node_id)\t$(ref_node.pressure)\t1\t;\n"
    file_text*=resevior

    tanks="\n[TANKS]\n;ID\tElevation\tInitLevel\tMinLevel\tMaxLevel\tDiameter\tMinVol\tVolCurve\n"
    for (k,t) in enumerate(tank_loc)
        tanks*="tank$(k)\t$(demand_nodes[k].elevation)\t$(T_min)\t$(T_max)\t$(T_max)\t$(UP[k]/T_max)\t$(0)\t\t;\n"
    end
    file_text*=tanks

    pipes="\n[PIPES]\n;ID\tNode1\tNode2\tLength\tDiameter\tRoughness\tMinorLoss\tStatus;\n"
    for (i,k) in pipe_seg
        n=link_list[i].start_node
        if n in tank_loc
            tank_k=findfirst(tloc->tloc==n,tank_loc)
            pipes*="$(i)\t$(link_list[i].start_node)\ttank$(tank_k)\t$(l_star[i][k])\t$(possible_pipes[k].roughness)\t$(0)\tOpen\t;\n"
            pipes*="$(NL+i)\ttank$(tank_k)\t$(link_list[i].end_node)\t$(l_star[i][k])\t$(0)\t$(0)\tOpen\t;\n"
        else
            pipes*="$(i)\t$(link_list[i].start_node)\t$(link_list[i].end_node)\t$(l_star[i][k])\t$(possible_pipes[k].roughness)\t$(0)\tOpen\t;\n"
        end
    end
    file_text*=pipes

    pumps="\n[PUMPS]\n;ID\tNode1\tNode2\tParameters\n"
    for i in findall(>(0),pe_star)
        n=link_list[i].start_node
        if n in tank_loc
            tank_k=findfirst(tloc->tloc==n,tank_loc)
            pumps*="$(i)\t$(link_list[i].start_node)\ttank$(tank_k)\tPump$(i)\t;\n"
            #pumps*="$(i)\t$(link_list[i].start_node)\ttank$(tank_k)\tPump$(i)"
        else
             pumps*="$(i)\t$(link_list[i].start_node)\t$(link_list[i].end_node)\tPump$(i)\t;\n"
        end
    end
    file_text*=pumps

    patterns="\n[PATTERNS]\n;ID\tMultipliers\n;General Default Demand Pattern\n"
    patterns*="1\t1.34\t1.94\t1.46\t1.44\t.76\t.92\n"
    file_text*=patterns

    curves="\n[CURVES]\n;ID\tX-Value\tY-Value;\n"
    file_text*=curves

    controls="\n[CONTROLS]\n"
    file_text*=controls

    energy="\n[ENERGY]\n"
    energy*="Global Efficiency\t$(eff*100)\nGlobal Price\t$(EP)\n"
    file_text*=energy

    times="\n[TIMES]\n"
    times*="Duration\t24:00\nHydraulic Timestep\t3:00\nQuality Timestep\t3:00\nPattern Timestep\t3:00\nPattern Start\t0:00\nReport Timestep\t3:00\nReport Start\t0:00\n"
    file_text*=times

    report="\n[REPORT]\n"
    report*="Status\tYes\nSummary\tNo\nPage\t0\n"
    file_text*=report

    coordinates="\n[COORDINATES]\n;Node\tX-Coord\tY-Coord\n"
    for de_n in demand_nodes
        coordinates*="$(de_n.node_id)\t$(ref_node.elevation*(N-de_n.node_id)/N)\t$(de_n.elevation*(1-rand(1:de_n.node_id)/(3*de_n.node_id)))\n"
    end
    for (k,t) in enumerate(tank_loc)
        coordinates*="tank$(k)\t$(ref_node.elevation*(N-demand_nodes[k].node_id)/N - 5*rand(0:1))\t$(demand_nodes[k].elevation*(1-rand(1:demand_nodes[k].node_id)/(3*demand_nodes[k].node_id))-5*rand(0:1))\n"
    end

    
    file_text*=coordinates
    dir="output"
    fpath=joinpath(dir,"$(filename).inp")
    open(fpath, "w") do file
        write(file, file_text)
    end
end

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

function tank_location(es_star)
    global tank_loc=[]
    for n=1:N
        i=findall(>(0),es_star[n])
        if length(i)>0
            tank_loc=[tank_loc;n]
        end
    end
end

function tank_d(e_star)
    global tank_type=[]
    for n in tank_loc
        k=findall(>(0),e_star[n])
        tank_type=[tank_type;k]
    end
end

function pipe_segment(l_star)
    global pipe_seg=[]
    for i=1:NL
        k_segs=findall(>(0),l_star[i])
        for k in k_segs
            pipe_seg=[pipe_seg;[(i,k)]]
        end
    end
end


function check_feasibility(model)
	optimize!(model)
	if termination_status(model)==INFEASIBLE::TerminationStatusCode
		println("model is INFEASIBLE")
		println(solution_summary(model,verbose=true))
    elseif (termination_status(model)==OPTIMAL::TerminationStatusCode)
        global es_star=value.(es)
        global es_star=[es_star[(n-1)*NL+1:n*NL] for n=1:N]
        tank_location(es_star)
        global e_star=value.(e)
        global e_star=[e_star[(n-1)*NE+1:n*NE] for n=1:N]
        tank_d(e_star)
        global l_star=value.(l)
        global l_star = [l_star[(i-1)*NP+1:i*NP] for i=1:NL]        
        pipe_segment(l_star)
        global pe_star=value.(pe)
    else    
		println("Model status : $(termination_status(model))")
	end
end



#reading source nodes 
data = readxl("input/karegon_scheme.xls", "Sheet1!A35:G35")
ref_node=reference_node(data[1],data[4],0,data[7])
water_density=1
g=9.8
println(ref_node)

#reading general data 
data = readxl("input/karegon_scheme.xls", "Sheet1!B37:B40")
#P_min::Float64=data[1]
#R_def::Float64=data[2]
PH::Int=data[4]
println(" Primary Hours : $(PH)")


# reading nodes
data = readxl("input/karegon_scheme.xls", "Sheet1!A3:F35")
#println(data)
N=33 
number_of_rows = length(data[:,1])
#println(number_of_rows)
def_link=link(0,0,0,0)
link_structure = fill(def_link,N, N)
demand_nodes=Array{demand_node}(undef,N)
for i = 1:number_of_rows 
    d1=demand_node(data[i,1],data[i,4],data[i,6],data[i,5])
    demand_nodes[d1.node_id]=d1
end 
demand_nodes[ref_node.node_id]=demand_node(ref_node.node_id,ref_node.elevation,0,0)
println(demand_nodes)


#readin pipes
data = readxl("input/karegon_scheme.xls", "Sheet1!A110:C122")
NP=13#number of pipes=number of rows
possible_pipes= Array{pipe}(undef,NP)
for i = 1:NP 
        possible_pipes[i]=pipe(i,data[i,1],data[i,2],data[i,3],0)
end
println(possible_pipes)


# reading links
data = readxl("input/karegon_scheme.xls", "Sheet1!A50:G80")
NL=31#number of links = number of rows
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
println("link structure is  : $(link_structure[33,:])")
println("Link List is : $(link_list[1])")
#println("pipe structure is : $(pipe_structure[3,7,:])")



#tank general data 
data = readxl("input/karegon_scheme.xls", "Sheet1!D100:D106")
SH::Float64=data[2]
tank_cp::Float64=data[3]
T_max::Float64=data[4]
T_min::Float64=10
#tank data 
data = readxl("input/karegon_scheme.xls", "Sheet1!A88:D97")
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

data = readxl("input/karegon_scheme.xls", "Sheet1!D136:D146")
PP_min=data[2]
eff=data[11]
CP=data[4]
EP=data[5]
Y=data[6]
DF=data[7]
INFR=data[8]

println("\nMin Pump Size : $(PP_min), Pump Efficiency : $(eff), Capital Cost : $(CP)\n
Opex : $(EP), Pump Lifetime : $(Y), Discount Rate : $(DF), Inflation Rate : $(INFR)")


# creating models 
wds_opt=Model(HiGHS.Optimizer)
set_silent(wds_opt)

#link variables
@variable(wds_opt,l[1:NP*NL].>=0)#possible lengths of pipes
@variable(wds_opt,l_p[1:NP*NL].>=0)#possible lengths of primary pipes
@variable(wds_opt,f[1:NL],Bin)#links in primary
@variable(wds_opt,hl[1:NL].>=0)#headloss across links 

#node variables
@variable(wds_opt,h[1:N].>=0) 
@variable(wds_opt,hs[1:N*NL].>=0)
@variable(wds_opt,d[1:N].>=0)

#tank variables 
@variable(wds_opt,z[1:NE*N].>=0) #demand being served by tank k at node N
@variable(wds_opt,e[1:NE*N],Bin) #if tank k is at node n 
@variable(wds_opt,s[1:N*N],Bin) #if tank at m is served by tank at N
@variable(wds_opt,t[1:N].>=0) #height of tank at node N
@variable(wds_opt,es[1:NL*N],Bin) #if tank at node n serves its immediated downstream link I

#pump variables
@variable(wds_opt,p[1:NL].>=0) #power of pump at link I
@variable(wds_opt,p_p[1:NL].>=0) #power of pump at primary networks link I
@variable(wds_opt,p_s[1:NL].>=0) #power of pump at secondary networks link I
@variable(wds_opt,pe[1:NL],Bin) #if pump is installed at link I
@variable(wds_opt,ph[1:NL].>=0) #head provided by the pump at link i

#forming the objective function
#link part
#cost vector ,roughness vector and Daimeter vectors
C=zeros(Float64,NP)
R=zeros(Float64,NP)
D=zeros(Float64,NP)

for i =1:NP
	C[i]=possible_pipes[i].cost
	R[i]=possible_pipes[i].roughness
	D[i]=possible_pipes[i].diameter/100
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

#considering tank cost 
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

L=zeros(NL)
for i=1:NL
    L[i]=link_list[i].length
end

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

L=zeros(Float64,NL)
FL_p=zeros(Float64,NL)
FL_s=zeros(Float64,NL)
for i =1:NL
	L[i]=link_list[i].length
	FL_p[i]=flow_p[Int(link_list[i].end_node)]/100
	FL_s[i]=flow_s[Int(link_list[i].end_node)]/100
end




# incorporating the previous network pipes 
ref_id=ref_node.node_id
data = readxl("input/karegon_scheme.xls", "previous network!A3:A9")
prev_N=[Int(d) for d in data]
prev_link=Set{Int}()
for prev_n in prev_N
    prev_path =S_n(prev_n,link_structure,ref_id)
    #println(prev_path)
    prev_path=get_curr_path(ref_id,prev_n,prev_path)
    global prev_link=union(prev_link,[l.link_id for l in prev_path])
end
println("\nPrevious links present : $(prev_link)")
index_i=collect(prev_link)
prev_L=length(index_i)
index_l=[(index_i[x]-1)*NP+4 for x in 1:prev_L]

@constraint(wds_opt,l[index_l].==L[index_i])
prev_cl=@expression(wds_opt,sum(C_l[index_l].*l[index_l]))

println("\nPrevious link feasibility")
check_feasibility(wds_opt)

#objective function : UN_n.*z + CT_n. *e

#objective function
@objective(wds_opt,Min,sum(C_l.*l)+sum(CP*p)+sum(UN_n.*z+CT_n.*e)+sum(PH*EP*DF*p_p)+sum(SH*EP*DF*p_s)-prev_cl)
pipe_cost=@expression(wds_opt,sum(C_l.*l))
tank_cost=@expression(wds_opt,sum(UN_n.*z+CT_n.*e))
pump_capital=@expression(wds_opt,sum(CP*p))
pump_opex=@expression(wds_opt,sum(PH*EP*DF*p_p)+sum(SH*EP*DF*p_s))
total_cost=0
#+Sum(UN_n.*z+CT_n.*e)
#link constraints 

#length constraint
#creating parameters for our problem 


println("Demand nodes : $(demand_nodes)")
println("Water flowing into each node at Primary networks : $(flow_p)\n Water Flowing into each node at Secondary Networks : $(flow_s)")


#println("Water flowing into each node : $(flow)")




path=S_n(3,link_structure,33)
#println("\nPath from 8 to 7 is $(path)\n")
#navigating the path


println("Path from 8->7 is : $(get_curr_path(33,3,path))")







#println("Length Vector : $(L)\n")
#println("Water Flow in Primary network : $(FL_p)\n")
#println("Water Flow in Secondary network : $(FL_s)\n")
for j=1:NL
	@constraint(wds_opt,sum(l[1+(j-1)*NP:(j)*NP])==L[j])
end

println("Length feasibility\n")
check_feasibility(wds_opt)

#headloss at everly link
#Headloss per unit length in each link
num_p=10.68*FL_p.^1.852
num_s=10.68*FL_s.^1.852
denom=1 ./((R.^1.852).*(D.^4.87))
HL_p=num_p*denom'
HL_s=num_s*denom'
println("Headloss for each possible assignment in the Primary network: $(HL_p[1,:]) \n its size $(size(HL_p))")
#println("Headloss for each possible assignment : in the Secondary network $(HL_s[1,:]) \n its size $(size(HL_s))")

#for i = 1:NL
#hl[i]=sum(HL_p[1+(i-1)*NP:(i)*NP]*.l_p[1+(i-1)*NP:(i)*NP]+HL_s[1+(i-1)*NP:(i)*NP]*.(l-l_p)[1+(i-1)*NP:(i)*NP])-ph[i]

for i=1:NL
	@constraint(wds_opt,hl[i]==sum(HL_p[i,:].*l_p[1+(i-1)*NP:(i)*NP]+HL_s[i,:].*(l-l_p)[1+(i-1)*NP:(i)*NP])-ph[i])
	@constraint(wds_opt,l_p[1+(i-1)*NP:(i)*NP].<=L[i]*f[i])
	@constraint(wds_opt,l[1+(i-1)*NP:(i)*NP]-l_p[1+(i-1)*NP:(i)*NP].<=L[i]*(1-f[i]))
	@constraint(wds_opt,l_p[1+(i-1)*NP:(i)*NP].<=l[1+(i-1)*NP:(i)*NP])
end

println("checking headloss check_feasibility\n")
check_feasibility(wds_opt)
#throw(DomainError)
#adding effective head due tanks constraint
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

M_p=2*maximum([L'*HL_p;L'*HL_s])

function Parent(node,link_list)
	for i=1:NL
		lk=link_list[i]
		if lk.end_node==node
			return lk
		end
	end
end

for n=1:N
	P=Parent(n,link_list)
	if P !== nothing
		m=P.start_node
		i=P.link_id
		println(" starts node $(m) ,link : $(i) ends at $(n)")
		@constraint(wds_opt,h[n]==hs[(m-1)*NL+i]-hl[i])
		@constraint(wds_opt,t[m]+E_n[m]-M_p*(1-es[(m-1)*NL+i])<=hs[(m-1)*NL+i])
		@constraint(wds_opt,hs[(m-1)*NL+i]<=M_p*(1-es[(m-1)*NL+i])+t[m]+E_n[m])
		@constraint(wds_opt,h[m]-M_p*(es[(m-1)*NL+i])<=hs[(m-1)*NL+i])
		@constraint(wds_opt,hs[(m-1)*NL+i]<=M_p*(es[(m-1)*NL+i])+h[m])
		@constraint(wds_opt,es[(m-1)*NL+i]<=s[(m-1)*N+m])
		@constraint(wds_opt,es[(m-1)*NL+i]<=1-f[i])
		@constraint(wds_opt,es[(m-1)*NL+i]>=s[(m-1)*N+m]-f[i])
	end
	@constraint(wds_opt,P_n[n]<=h[n]-(E_n[n]+t[n]))
	#@constraint(wds_opt,hs[(m)*NL+i]==(t[m]+E_n[m])*es[m*NL+i]+h[m]*(1-es[m*NL+i]))
end
ref_id=ref_node.node_id
@constraint(wds_opt,H_r==h[ref_id])
println("Tank head feasibility\n")
check_feasibility(wds_opt)

# tank capacity constraints
for n=1:N
	@constraint(wds_opt,T_min.<=t[n])
	@constraint(wds_opt,T_max.>=t[n])
end
println("Tank capacity feasibility\n")
check_feasibility(wds_opt)

DE=zeros(N)
for n = 1:N
	if isassigned(demand_nodes,n)
		DE[n]=demand_nodes[n].water_demand
	end
end
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
		O=setdiff!(union!(desc,asc),n)
		O_index=collect([(o-1)*N+o for o in O])
		for o in O_index
			#@constraint(wds_opt,s[(n-1)*N+m]+s[o]<=1)
		end
		 # each node only has 1 tank located in a path
		@constraint(wds_opt,s[(n-1)*N+n]>=s[(m-1)*N+m]) #if n can't serve itself then its descandants m can't serve itself, should be removed
		@constraint(wds_opt,s[(n-1)*N+n]>=s[(n-1)*N+m])#if n can serve itself it can serve its descandants m
	end
	desc=collect(desc)
	index_desc=[(n-1)*N+m for m in desc]
	index_asc=[(m-1)*N+n for m in asc]
	if desc !== nothing
		@constraint(wds_opt,d[n]==sum(DE[desc].*s[index_desc])) #total demand served by node n is the sum of demands of nodes that relies on it
	end
	if asc !==nothing
		@constraint(wds_opt,sum(s[index_asc])==1) #only 1 Ascending node m can serve n
	end
	if P !==nothing
		i=P.link_id
		@constraint(wds_opt,f[i]==s[(n-1)*N+n])
	end
end

println("\nTank heirarchy feasibility")
check_feasibility(wds_opt)

#setting up znk
#z[(n-1)*NE+k]=e[(n-1)*NE+k]*d[n]
M_c=6*maximum([t.UP for t in possible_tanks])
println("maximum z value is $(M_c)")

#z<=M_c*e
#z<=d
#z>=d-M_c*(1-e)
for n =1:N
	for k=1:NE
		@constraint(wds_opt,z[(n-1)*NE+k]<=M_c*e[(n-1)*NE+k])
		@constraint(wds_opt,z[(n-1)*NE+k]<=d[n])
		@constraint(wds_opt,z[(n-1)*NE+k]>=-M_c*(1-e[(n-1)*NE+k])+d[n])
	end
end

println("\nTank selection feasibility")
check_feasibility(wds_opt)

#one possible tank for a node
#e[1+(k-1)*NE:(k)*NE]==1
for n =1:N
	@constraint(wds_opt,sum(e[1+(n-1)*NE:(n)*NE])==1) # making model infeasible

end
println("\n1 Tank feasibility")
check_feasibility(wds_opt)

#d between min and max capacity 
for n =1:N
		@constraint(wds_opt,sum(LO.*e[(n-1)*NE+1:n*NE])<=d[n])
		@constraint(wds_opt,d[n]<=sum(UP.*e[(n-1)*NE+1:n*NE]))
end
println("\nTank capacity feasibility")
check_feasibility(wds_opt)
# println(wds_opt)

# pump constraints 
M_pow=2*maximum((water_density*g/eff)*M_p*FL_p)
PP_max=M_pow
@constraint(wds_opt,p_p+p_s.==p)
@constraint(wds_opt,p_p.<=6*M_p*f)
@constraint(wds_opt,p_p.<=(water_density*g/eff)*FL_p.*ph)
@constraint(wds_opt,p_p.>=((water_density*g/eff)*FL_p.*ph).-(6*M_p*(ones(NL)-f)))
@constraint(wds_opt,p_s.<=6*M_p*(ones(NL)-f))
@constraint(wds_opt,p_s.<=(water_density*g/eff)*FL_s.*ph)
@constraint(wds_opt,p_s.>=((water_density*g/eff)*FL_p.*ph).-(6*M_p*f))
@constraint(wds_opt,PP_min*pe.<=p)
@constraint(wds_opt,PP_max*pe.>=p)

println("Pump feasibility\n")
check_feasibility(wds_opt)

title_text="Optimal hydralic component allocation for Karegaon Redesign scheme"
fname="karegoan_ilp"
save_network(title_text,fname,demand_nodes,ref_node,tank_loc,tank_type,T_min,T_max,UP,pipe_seg,link_list,l_star,possible_pipes)
#throw(DomainError("Stoped by user"))
#optimize our Model
#optimize!(wds_opt)
#println(solution_summary(wds_opt))
#l_star=reshape(value.(l),NL,NP)
#println("optimal values : $(l_star)")
total_cost=value.(objective_function(wds_opt))
pipe_cost=value.(pipe_cost)
tank_cost=value.(tank_cost)
pump_capital=value.(pump_capital)
pump_opex=value.(pump_opex)
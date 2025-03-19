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
	end_note::Int
	length::Float64
end

function check_feasibility(model)
	!optimize(model)
	if termination_status(model)==INFEASIBLE::TerminationStatusCode
		println("model is INFEASIBLE")
		println(solution_summary(model,verbose=True))
	else
		println("Model status : $(termination_status(model))")
	end
end

#reading source nodes 
data = readxl("input/Sample_input_jt.xls", "General!D16:D19")
ref_node=reference_node(data[1],data[3],0,data[4])
#println(ref_node)


# reading nodes
data = readxl("input/Sample_input_jt.xls", "Node Data!A2:E10")
N=11
number_of_rows = length(data[:,1])
#println(number_of_rows)
def_link=link(0,0,0,0)
link_structure = fill(def_link,N, N)
demand_nodes=Array{demand_node}(undef,N)
for i = 1:number_of_rows 
    d1=demand_node(data[i,1],data[i,3],2*data[i,4],data[i,5])
    demand_nodes[d1.node_id]=d1
end 
demand_nodes[ref_node.node_id]=demand_node(ref_node.node_id,ref_node.elevation,0,0)
println(demand_nodes)



#readin pipes
data = readxl("input/Sample_input_jt.xls", "Pipe Data!A3:C15")
NP=13
possible_pipes= Array{pipe}(undef,NP)
for i = 1:NP 
    possible_pipes[i]=pipe(i,data[i,1],data[i,2],data[i,3],0)
end
println(possible_pipes)

# reading links
data = readxl("input/Sample_input_jt.xls", "Link Data!A2:G10")
NL=9
pipe_structure =zeros(N,N,NP)
adj_mat=Matrix{Float64}(I,N,N)
for i = 1:NL
	link1 = link(data[i,1],data[i,2],data[i,3],data[i,4])
    link_structure[Int(data[i,2]),Int(data[i,3])]=link1
	diameter_present=data[i,5]
	if diameter_present>0
		for i =1:NP
			if possible_pipes[i].diameter==diameter_present
				pipe_structure[Int(data[i,2]),Int(data[i,3]),i] = data[i,4]
			end
		end
	end
	#println(diameter_present)
end 
#println("link structure is  : $(link_structure[8,:])")
#println("pipe structure is : $(pipe_structure)")
adj_mat=map(x-> if x.link_id>0 1 else 0 end,link_structure)
# water flowing into each node
flow=zeros(N)
for i=1:N
	if isassigned(demand_nodes,i)
		flow[i]=demand_nodes[i].water_demand
	else
		flow[i]=0.0
	end
end
flow_n=flow
for i =1:NL
	global flow_n=adj_mat*flow_n+flow
end
flow=flow_n
#flow=adj_mat*(adj_mat*(adj_mat*flow+flow)+flow)+flow
#println("Water flowing into each node : $(flow)")


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

#creating parameters for our problem 
L=zeros(Float64,NL)
FL=zeros(Float64,NL)
for i =1:NL
	L[i]=data[i,4]
	FL[i]=flow[Int(data[i,3])]
end


#println("Length vector : $(L)")
#println("Flow Link vector : $(FL)")

C=zeros(Float64,NP)
R=zeros(Float64,NP)
D=zeros(Float64,NP)

for i =1:NP
	C[i]=possible_pipes[i].cost
	R[i]=possible_pipes[i].roughness
	D[i]=possible_pipes[i].diameter
end
#println("Cost vector : $(C)")
#println("Roughness vector : $(R)")
#println("Daimeter vector : $(D)")

I_c=Matrix{Float64}(I,NP,NP)
C_l=[]
for i=1:NL
    global C_l=[C_l;I_c]
end
C_l=C_l*C
#println("Cost per length : $(C_l)")
#println("Size of cost vector is $(size(C_l))")

num=10.68*FL.^1.852
denom=1 ./(R.*D.^4.87)
HL=num*denom'
#LC=reduce(vcat,C*L')
#println(LC[1:NP])
#println((C*L')[:,1])
#println("$(NP),$(NL)")

#println("Headloss for each possible assignment : $(HL[1,:]) \n its size $(size(HL))")

H_r=ref_node.pressure
println("Source node pressure : $(H_r)")
P_n=zeros(Float64,N)
E_n=zeros(Float64,N)
println
for i=1:N
	if isassigned(demand_nodes,i)
		P_n[i]=demand_nodes[i].min_pressure
		E_n[i]=demand_nodes[i].elevation
	end
end
#println("Min pressure vector : $(P_n)")
#println("Elevation vector : $(E_n)")
#setting up the model
pipe_model=Model(HiGHS.Optimizer)
set_silent(pipe_model)
#binary variable
@variable(pipe_model,l[1:NP*NL].>=0)
#println(pipe_model)
#objective
@objective(pipe_model,Min,sum(C_l.*l))

#length constraint
for j=1:NL
		@constraint(pipe_model,sum(l[1+(j-1)*NP:(j)*NP])==L[j])
end
#println(pipe_model)
for n=1:N
	path_n=S_n(n,link_structure,8)
	HL_n=[]
	index_n=[]
	for p in values(path_n)
		i=p.link_id
		HL_n=[HL_n;HL[i,:]]
		#print(HL_n)
		index_n=[index_n;collect(1+(i-1)*NP:i*NP)]		
	end
	#println("Indexes for $(n) is $(index_n)")
	@constraint(pipe_model,P_n[n]<=H_r-E_n[n]-sum([HL_n.*l[index_n];0]))
end
optimize!(pipe_model)
println(solution_summary(pipe_model))
l_star=reshape(value.(l),NL,NP)
println("optimal values : $(l_star)")
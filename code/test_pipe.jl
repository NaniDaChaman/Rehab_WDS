include("branched_pipe_lp.jl")

mutable struct demand_nodes_output
    node_id::Int
    head::Float64
    pressure::Float64
end

mutable struct pipe_segment_output
    pipe_id::Int
    link_id::Int
    length::Float64
    flow::Float64
    speed::Float64
    headloss::Float64
    cost::Float64
end

#reading general output 
data = readxl("output/Sample_output.xls", "General!D23:D25")
total_cost::Float64= data[3]
pipe_length::Float64=data[2]
network_length::Float64=data[1]

# reading the node output
#A35:G44
data = readxl("output/Sample_output.xls", "General!A36:G44")
number_of_rows = length(data[:,1])
#println(number_of_rows)
def_link=link(0,0,0,0)
link_structure_output = fill(def_link,N, N)
demand_nodes_outputs=Array{demand_nodes_output}(undef,N)
for i=1:number_of_rows
    d1=demand_nodes_output(data[i,1],data[i,5],data[i,6])
    demand_nodes_outputs[Int(d1.node_id)]=d1
end

#reading the pipe segment results
#A50:K62
data = readxl("output/Sample_output.xls", "General!A50:K62")
number_of_rows = length(data[:,1])
pipe_segment_outputs=Array{pipe_segment_output}(undef,number_of_rows)
for i =1:number_of_rows
    n=Int(data[i,2])
    m=Int(data[i,3])
    println("Starting node : $(n), Ending node : $(m)")
    ps1=pipe_segment_output(data[i,1],link_structure[n,m].link_id,data[i,4],data[i,5],data[i,6],data[i,9],data[i,11])
    pipe_segment_outputs[i]=ps1
end


index_star=findall(x-> x>0 ,l_star)   
O_N=length(index_star)
pipe_segments_model=Array{pipe_segment_output}(undef,O_N) 
for i=1:O_N
    ps2=pipe_segment_output(index_star[i][1],index_star[i][2],l_star[index_star[i]],FL[index_star[i][1]],0,0,0)
    pipe_segments_model[i]=ps2
end


println(pipe_segments_model)
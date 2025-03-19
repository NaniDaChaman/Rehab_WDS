using JuMP, GLPK ,HiGHS

mutable struct Nanobot #needed if you wanna specify that a struct can be changed
    x::Int
    y::Int
    z::Int
    radius::Int
end


bot1=Nanobot(0,0,0,0) #need to declare all the values for the bot 
bot1.x=5
println("Bot1 value : ",bot1.x)

function processinput(filename)
    puzin = [[parse(Int, m.match) for m in eachmatch(r"(\-?\d+)", line)]
             for line in readlines(filename)]
    bots = [Nanobot(l...) for l in puzin]
    return bots
end

bots=processinput("input.txt")

model=Model(HiGHS.Optimizer)
set_silent(model)
@variable(model,x,Int)
@variable(model,y,Int)
@variable(model,z,Int)
@variable(model,botsinrange[1:length(bots)],Bin)

slack = 2*4*maximum(b -> max(abs(b.x), abs(b.y), abs(b.z), b.radius), bots)
for (i,b) in enumerate(bots)
    for absmult in Iterators.product((-1,1),(-1,1),(-1,1))
        @constraint(model,absmult[1]*(b.x-x)+absmult[2]*(b.y-y)+absmult[3]*(b.z-z)<=b.radius+(1-botsinrange[i])*slack)
    end
end

@objective(model, Max, sum(botsinrange))
optimize!(model)
println("""x=$(round(Int, value(x))), y=$(round(Int, value(y))), z=$(round(Int, value(z)))
ans=$(round(Int, value(x)) + round(Int, value(y)) + round(Int, value(z)))""")

#print(model)

# building my model 
model2=Model(HiGHS.Optimizer)
set_silent(model2)
R=[bot.radius for bot in bots]
X=[bot.x for bot in bots]
Y=[bot.y for bot in bots]
Z=[bot.z for bot in bots]
#println(R)
# desc var i, this desc var is a one hot vector which specifies the nanobot with the highest range
@variable(model2,i[1:length(bots)],Bin)

# desc var j, this variable tells us which bots are in the range of the bot with the maximum range
@variable(model2,j[1:length(bots)],Bin)

# one hot vector constraint
@constraint(model2,sum(i)==1)
# selection constraints 
slack = 2*4*maximum(b -> max(abs(b.x), abs(b.y), abs(b.z), b.radius), bots)
for (k,b) in enumerate(bots)
    for absmult in Iterators.product((-1,1),(-1,1),(-1,1))
        @constraint(model2,absmult[1]*(b.x-sum(X.*i)) + absmult[2]*(b.y-sum(Y.*i))+absmult[3]*(b.z-sum(Z.*i))<=sum(R.*i)+(1-j[k])*slack)
        
    end
end

@objective(model2,Max,sum(R.*i)+sum(j))
optimize!(model2)
println(value.(i))
println(sum(value.(j)))

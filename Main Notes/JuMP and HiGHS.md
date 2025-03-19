2025-01-2710:39

Tags : ongoing

Status : #opt 

## Advent of Code : Experimental Emergency Teleportation
- [ ] What is the problem ? 
- [ ] How can it be re written as an MILP ? 

# Check the status of our solutions in our Model 
- Local or Global solution available 
```julia
is_solved_and_feasible(model)
```
- Global solution available
```julia
is_solved_and_feasible(model,allow_local=false)
#true/false
```
- Solutions Summary 
```julia
solution_sumamary(model)
```
- Why did the solver stop 
- Primal and Dual Solutions
- Workflow recommended
## JuMP and HiGHS Questions
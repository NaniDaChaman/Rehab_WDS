2024-12-2014:08

Tags :

Status : 

## Todo
- [ ] Simplex function has some problems
- [ ] Fixing the pivot operation
- [ ] updating basis and non basis after every pivot
- [ ] Displaying the solution

## Preview of the simplex method : 

- [ ] What is the optimality condition ? 
- All the non basic variables in a canonical form have non positive coeffs in the objetive function
- Optimal value : basic solution
- [ ] What is the un bounded criterion ? 
- Any of the non basic variables in the canonical form has a positive coeff in objective function and non pos coeffs in all the constraints 
- [ ] How can we improve a non optimal solution ? 
- Non optimality criterion : one of the non basic var has pos coeffs in obj but pos coeefs in atleast one of the constraints.
- 
- [ ] What is pivoting across a non basic solution ? 
- [ ] How can we reduce a general equation to its canonical form ? 

## Implementing canonical form algorithm 
- We are stuck at getting the optimal solution of progam that is at an optimal criterion
- the basic feasible solution is the optimal solution 
- I think we just : 
```cpp
float opt_val=c(-1); \\get the last element
```
- Plan is to write some test
- Note down 3 lps which are in canonical form of varying size and length
- make sure that your program compiles and there aren't any syntactical errors 
- make sure you can import and run all your functions 

## Test for LPs
```cpp
VectorXf c{{0,0,-3,-1,20}};
MatriXf A {{1,0,,-3,3}
		  {0,1,-8,4}};
VectorXf b {{6,4}};
```
## Linear program Questions
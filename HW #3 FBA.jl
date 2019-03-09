using JuMP
using GLPKMathProgInterface
using DelimitedFiles

include("HW #3 with Balance.Jl")
S_matrix

#gets length of v vector
n=length(exprxnv)

m=Int64(length(S_matrix)/n)
#C transposed
#Maximize Urea production
cT=zeros(1,n)

#Maximize Urea production
for i in 1:length(exprxnv)#length(exprxnv) #column i of matrix
    bigsplit=split(exprxnv[i],r" -> | <=> ")#splits rxn into left and right half
    right_rxn=bigsplit[2]
    if occursin("C00086",right_rxn)
        cT[i]=1
    end
end
cT
#S_matrix
S_matrix

#check up on information vectors
rxns
ec
exprxnv

#V_boundaries_vector

#buidling V_boundtries
include("V Building.jl")


V_boundaries_vector #built by hand see V Bulding.jl and read me explaining how they were made

#FBA constraint
zero_vector_FBA_constraint=zeros(m,1)

#Model Construction

P_Model=Model(solver=GLPKSolverLP())

#Variables

@variable(P_Model, V_boundaries_vector[i,1] <= V_vector[i in 1:n] <= V_boundaries_vector[i,2])

#Objective

@objective(P_Model, Max, (cT*V_vector)[1])

#Constraints

@constraint(P_Model, S_matrix*V_vector .== zero_vector_FBA_constraint)

#Output
Status=solve(P_Model)

print(P_Model)

v=zeros(n,1)

#Building solution vector
for i in 1:n
    v[i]=getvalue(V_vector[i])
end

#here is the flux's for this system
v


#Obtaining optimized value, here it is the maximum rate of urea production
opt_value=(cT*v)[1]

using LibCURL
using DelimitedFiles

#User input (need list of e.c's) put as column vector only the form "x.x.x.x" x's can
#be more than one digit but the x's represent a number
#results will not contrain exchange reactions those must be manuelly inputed
ec=["3.5.3.1";"2.1.3.3";"4.3.2.1";"6.3.4.5";"1.14.13.39"]

#rxns will contain a list of the chemical reactions
#creactants contains a list of chemical reactants in terms of
#kegg nomincalature from which matrix EB short for elemental balence is
#constructed and with the rows correctioning to the elements (can be expanded to entire periodic table if wanted)
#but here it is [C;H;N;O;P;S], as for the columns they are ordered by the metabolites with row equal to column number

reb=["C";"H";"N";"O";"P";"S"]#rows of elemental balance matrix

rxns=Set{String}() # used to prevent duplicates
creactants=Set{String}() #no used here but could be used for larger reaction networks
exprxns=Set{String}() #used to preven duplicates

#h's are temporary vectors
for i in 1:length(ec)
    global rxns
    global stp
    h=ec[i]
    h1=read(`curl -X GET http://rest.kegg.jp/link/rn/$h/`, String)

        s=split(h1,"\t")
        st=split.(s,"\n")
        stp=st[2:length(st)]
        for i in 1:length(stp)
            global stp
            global rxns
            rxns=union(rxns,[split.(stp[i][1],"rn:")[2]])
        end
end

#change set to vector for easy access
rxnv=collect(rxns)

for j in 1:length(rxnv)
    global rxnv
    global exprxns

    h=rxnv[j]
    v=read(`curl -X GET http://rest.kegg.jp/get/reaction:$h/`, String)
    t=split(v,"\nEQUATION    ")[2]
    d=split(t,"\nRCLASS")[1]
    c=split(d,"\nCOMMENT")[1]
    exprxns=union(exprxns,[c])

end

#vectorizing
exprxnv=collect(exprxns)

#adding on ex rxn givens
exprxnv=push!(exprxnv," -> C00169")
exprxnv=push!(exprxnv," -> C00049")
exprxnv=push!(exprxnv,"C00122 -> ")
exprxnv=push!(exprxnv,"C00086 -> ")

#adding on ex rxn to balance elements, Added after running the first tome
#needed an ATP exchange rxn and an ADP exchange rxn to have a full elemental balance
#assume this is an exchange required and a H20 export rxn. the none exchange rxns were already balenced
exprxnv=push!(exprxnv,"C00001 -> ")
exprxnv=push!(exprxnv,"C00002 -> ")
exprxnv=push!(exprxnv," -> C00008")

chemicals_raw=[]
for i in 1:length(exprxnv)#length(rxns) #column i of matrix
    global exprxnv
    global chemicals_raw
    bigsplit=split(exprxnv[i],r" -> | <=> ")#splits rxn into left and right half
    left_rxn=bigsplit[1]
    left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
    left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
    for j in 1:length(left_chem_and_co)#collects all unique chemicals on left side
        if length(left_chem_and_co[j])==2 && (!(in(left_chem_and_co[j][2],chemicals_raw)) && left_chem_and_co[j][2]!="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][2])
        elseif length(left_chem_and_co[j])==1 && (!(in(left_chem_and_co[j][1],chemicals_raw)) && left_chem_and_co[j][1]!="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,left_chem_and_co[j][1])
        end
    end
    right_rxn=bigsplit[2]
    right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
    right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
    for j in 1:length(right_chem_and_co)#collects all unique chemicals on right side
        if length(right_chem_and_co[j])==2 && !(in(right_chem_and_co[j][2],chemicals_raw)) && !(right_chem_and_co[j][2]=="")#if chemical has coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][2])
        elseif length(right_chem_and_co[j])==1 && !(in(right_chem_and_co[j][1],chemicals_raw)) && !(right_chem_and_co[j][1]=="")#if chemical has no coeffiecent + uniqueness constraints + protection against "" inculsion
            chemicals_raw=push!(chemicals_raw,right_chem_and_co[j][1])
        end
    end
end
chemicals_raw

#building S_matrix

S_matrix=zeros(length(chemicals_raw),length(exprxnv))

#processing of rxns into stoichometric matrix
for j in 1:length(exprxnv)#length(rxns) #column j of matrix
    global S_matrix
    global exprxnv
    global chemicals_raw
    bigsplit=split(exprxnv[j],r" -> | <=> ")#splits rxn into left and right half
    left_rxn=bigsplit[1]
    left_chem_raw=split(left_rxn," + ")#subdivides the left half into chemicals and their coeffiecent
    left_chem_and_co=split.(left_chem_raw," ")#splits up chemicals and their coeffiecent
    for k in 1:length(left_chem_and_co)#Places coeffiecents in s_matrix and since its the left side the entries recieve an negative sign
        if length(left_chem_and_co[k])==2 && left_chem_and_co[k][2]!=""#grabs coeffiecent of chemicals if present
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(left_chem_and_co[k][2],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=-1*parse(Int,left_chem_and_co[k][1])
                end
            end
        elseif length(left_chem_and_co[k])==1 && left_chem_and_co[k][1]!=""#grabs present chemicals and assigns -1 to their S_matrix possition
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(left_chem_and_co[k][1],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=-1
                end
            end
        end
    end
    right_rxn=bigsplit[2]
    right_chem_raw=split(right_rxn," + ")#subdivides the right half into chemicals and their coeffiecent
    right_chem_and_co=split.(right_chem_raw)#splits up chemicals and their coeffiecent
    for k in 1:length(right_chem_and_co)
        if length(right_chem_and_co[k])==2 && !(right_chem_and_co[k][2]=="")#grabs coeffiecent of chemicals if present
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(right_chem_and_co[k][2],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=parse(Int,right_chem_and_co[k][1])
                end
            end
        elseif length(right_chem_and_co[k])==1 && right_chem_and_co[k][1]!=""#grabs present chemicals and assigns 1 to their S_matrix position
            for i in 1:length(chemicals_raw) #possible row number
                if occursin(right_chem_and_co[k][1],chemicals_raw[i])#finds correct row number
                    S_matrix[i,j]=1
                end
            end
        end
    end
end

S_matrix

chemicals_raw

chemicals_raw_formula=[]
#building EB
for k in 1:length(chemicals_raw)
    global chemicals_raw
    global chemicals_raw_formula

    cname=chemicals_raw[k]
    h=read(`curl -X GET http://rest.kegg.jp/find/compound/$cname/`, String)
    Cname=split(h,"\t")[1]
    h1=read(`curl -X GET http://rest.kegg.jp/get/$Cname/`, String)
    h2=split(split(h1,"\nFORMULA")[2],"\nEXACT_MASS")[1]
    h3=strip(h2)
    chemicals_raw_formula=append!(chemicals_raw_formula,[h3])
end
#post building check point
chemicals_raw_formula

#prebuilding EB matrix
EB=zeros(length(reb),length(chemicals_raw_formula))

#EB reference vector
EBR="CHNOPS"

#generalizable if have more chemicals by assembling list of chemicals
#and adding in addition for loop inexing by reb instead of "Element"
for l in 1:length(chemicals_raw_formula)
    global chemicals_raw_formula
    global EB
    global EBR
    for p in 1:length(chemicals_raw_formula[l])
        global EB
        global chemicals_raw_formula
        global EBR
            k=chemicals_raw_formula[l][p]
            knp='X'
            z=1
            while (p+z)<=length(chemicals_raw_formula[l]) && !occursin(chemicals_raw_formula[l][p+z],EBR)
                if z==1
                    knp=convert(Int64,chemicals_raw_formula[l][p+z])-48#-48 to convert correcting to 0 at '0'
                    z=z+1
                elseif z==2 #(can be expanded even further but its a fair assumption that we are dealing with macromolecules)
                    knp=((convert(Int64,chemicals_raw_formula[l][p+1])-48)*10)+(convert(Int64,chemicals_raw_formula[l][p+z])-48)#-48 to convert correcting to 0 at '0' terms
                    z=z+1
                elseif z==3 #(can be expanded even further but its a fair assumption that we are dealing with macromolecules)
                    knp=((convert(Int64,chemicals_raw_formula[l][p+2])-48)*100)+((convert(Int64,chemicals_raw_formula[l][p+2])-48)*10)+(convert(Int64,chemicals_raw_formula[l][p+z])-48)#-48 to convert correcting to 0 at '0' terms
                    z=z+1
                end
            end
            #finding which row and column the multiplier goes into
            #then determining if it is just 1 or if it has a multiplier
            if isequal('C',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[1,l]=knp
                else
                    EB[1,l]=1
                end
            elseif isequal('H',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[2,l]=knp
                else
                    EB[2,l]=1
                end
            elseif isequal('N',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[3,l]=knp
                else
                    EB[3,l]=1
                end
            elseif isequal('O',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[4,l]=knp
                else
                    EB[4,l]=1
                end
            elseif isequal('P',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[5,l]=knp
                else
                    EB[5,l]=1
                end
            elseif isequal('S',chemicals_raw_formula[l][p])
                if Int64==typeof(knp)
                    EB[6,l]=knp
                else
                    EB[6,l]=1
                end
            end
    end
end
EB

#checking if elementarily balanced
ebs=EB*S_matrix

#confirming the elentary balance as the rows add sum to zero this comfirms the
#elementary balance
sum(ebs[1,:])
sum(ebs[2,:])
sum(ebs[3,:])
sum(ebs[4,:])
sum(ebs[5,:])
sum(ebs[6,:])

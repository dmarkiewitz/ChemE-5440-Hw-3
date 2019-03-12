#builds the v_boundaries_vector
#E is in units of mmol/gDW
#kcat is in units of 1/hr
#sat is the saturation terms found in excel file found in get file and the
#column containing the information needed is highlighted yellow
#column 1 is lower bound
#column 2 is upper bound
V_boundaries_vector=zeros(n,2)

#kcats
kcat3531=249*3600
kcat2133=88.1*3600
kcat4321=34.5*3600
kcat6345=203*3600
kcat1141339=13.7*3600

#Enzyme levels in units of mmol/gDW
E=10^-5

sat3531=1.42/10 #used homo
sat2133=1.1/100 #used e. coli not great assumption but only value found in paper
sat4321=8.38/100 #used homo
sat6345=9.23/10 #used Mus musculus
sat1141339=9.86/10 #used Mus musculus

#constructing V_boundaries_vectors
V_boundaries_vector[1,1]=-sat2133*E*kcat2133
V_boundaries_vector[1,2]=sat2133*E*kcat2133
V_boundaries_vector[2,1]=-sat6345*E*kcat6345
V_boundaries_vector[2,2]=sat6345*E*kcat6345
V_boundaries_vector[3,1]=-sat1141339*E*kcat1141339
V_boundaries_vector[3,2]=sat1141339*E*kcat1141339
V_boundaries_vector[4,1]=-sat1141339*E*kcat1141339
V_boundaries_vector[4,2]=sat1141339*E*kcat1141339
V_boundaries_vector[5,1]=-sat4321*E*kcat4321
V_boundaries_vector[5,2]=sat4321*E*kcat4321
V_boundaries_vector[6,1]=-sat3531*E*kcat3531
V_boundaries_vector[6,2]=sat3531*E*kcat3531
V_boundaries_vector[7,1]=-sat1141339*E*kcat1141339
V_boundaries_vector[7,2]=sat1141339*E*kcat1141339
V_boundaries_vector[8,1]=0
V_boundaries_vector[8,2]=10
V_boundaries_vector[9,1]=0
V_boundaries_vector[9,2]=10
V_boundaries_vector[10,1]=0
V_boundaries_vector[10,2]=10
V_boundaries_vector[11,1]=0
V_boundaries_vector[11,2]=10
V_boundaries_vector[12,1]=0
V_boundaries_vector[12,2]=10
V_boundaries_vector[13,1]=0
V_boundaries_vector[13,2]=10
V_boundaries_vector[14,1]=0
V_boundaries_vector[14,2]=10

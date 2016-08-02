
#setwd("G://work_2013/Papers/Paper1/GHG_kong/2009")
#source("G://work_2013/Papers/Paper1/ghg_function1.r")
# We load aggregation function to obtain GHG fit for the model.
#source("G://work_2013/Papers/Paper1/agg.r")
# loading data

#setwd("G://work_2013/Papers/Paper1/GHG_kong/2009")
source("ghg_function1.r")
# We load aggregation function to obtain GHG fit for the model.
source("agg.r")



#Energy Balance in Gross Calorific Value
EB_G=read.csv(file="EB_G_2009.csv",header=T)

#Conversion Factor from Gross Calorific Value to Net Calorific Value
CF=read.csv(file="CF_2009.csv",header=T)

#Row names of IO table
RowName=read.csv(file="EB_row_2009.csv",header=F)

# IO tables in basic price (YR 2009)
IO_Domestic=read.csv(file="IO_Domestic_2009.csv", header=T)
IO_Import=read.csv(file="IO_Import_2009.csv",header=T)
IO_whole=read.csv(file="IO_whole_2009.csv",header=T)
EC=read.csv(file="EC1996.csv",header=T)

# Coefficients : buring rate of Crude Oil/Change in Oil Product stock/Oil product domestic production
Oil.Data=read.csv(file="OilData_2009.csv",header=T)
YR=2009
#YR_1=YR-1
ECO2_2009=Gen_CO2.1(EB_G,CF,RowName,IO_Domestic,IO_Import,IO_whole,EC,Oil.Data,YR)

ghg_2009_m=ECO2_2009[[1]]
ghg_2009_T=ECO2_2009[[2]]

## ghg fit for model
# loading index file to aggregate CO2 generation : feeding files for mocel
#(1) GTAP equivalent set
# row index is index for fuels and electricities. column index is the same with regular column index , except we don't have import and import taxes

gindr=read.csv(file="gdix_r2.csv",header=F)
#gindc=read.csv(file="G://work_2013/Papers/Paper1/kdix_c2.csv",header=F)
gindc=read.csv(file="kdix_c2.csv",header=F)
# we cut down import taxes
gindc=gindc[(1:65),]
# we cut industry index column and the total row
gdata2009=ghg_2009_m[-22,-1]
# change into matrix
gdata2009.m=as.matrix(gdata2009)
# obtain gtap equivalence in column, model equivalance in row
ghg2009_prelim=agg2(gindr,gindc,gdata2009.m)

# (2) grown up ghg
# row is already aggregated to fit for the model, column needs some more aggregation
g22_r_g=as.matrix(1:7)
#g22_c_g=read.csv(file="G://work_2013/Papers/Paper1/g22_c.csv",header=F)
g22_c_g=read.csv(file="g22_c.csv",header=F)
# Again we cut down column indices for import and import taxes
g22_c_g=g22_c_g[(1:31),]
ghg2009_prelim.m=as.matrix(ghg2009_prelim)
ghg2009_grown=agg2(g22_r_g,g22_c_g,ghg2009_prelim.m)
# We cut down final demand, only keep CO2 generated via intermediate demand consumption 
ghg2009_grown=ghg2009_grown[,(1:22)]

# (3) baby ghg
br_g=g22_r_g
#bc_g=read.csv(file="G://work_2013/Papers/Paper1/bc.csv",header=F)
bc_g=read.csv(file="bc.csv",header=F)
bc_g=bc_g[(1:21),]
ghg2009_baby=agg2(br_g,bc_g,ghg2009_prelim.m)
ghg2009_baby=ghg2009_baby[,(1:12)]

#(4) add name of industries
rownames(ghg2009_grown)=c("Elec-c","Oil-c","Roil-c","Coal-c","Gas-c","Nuclear-C","Hydro-C")
colnames(ghg2009_grown)=c("AFF-a","Paper-a","Chemi-a","I_S-a","Nonfm-a","Mineral-a","Food-a","Othr_m-a","Vehicle-a","Water-a","Cons-a","Wood-a","Trans-a","Serv-a","Elec-a","Oil-a","ROil-a","Coal-a","Gas-a","Nuclear-a","Hydro-a","Dwelling-a")

rownames(ghg2009_baby)=c("Elec-c","Oil-c","Roil-c","Coal-c","Gas-c","Nuclear-C","Hydro-C")
colnames(ghg2009_baby)= c("agri-a","eint-a","othr-a","serv-a","trans-a","elec-a","oil-a","roil-a","coal-a","Gas-a","Nuclear-a","Hydro-a")

#(5) take out zero cells. Oil/Elec/Hydro/Nuclear are out. They don't generate CO2
ghg2009_grown_core=ghg2009_grown[c("Roil-c","Coal-c","Gas-c"),]/100 #convert units in 100,000 ton 
ghg2009_grown_core["Coal-c","Coal-a"]=0 # CO2 generated in coal mining is ignored. It is rather small (300,000 ton)
ghg2009_grown_core=ghg2009_grown_core*(ghg2009_grown_core>=1) #Ignore emission less than 100,000 ton

#(6) add process CO2
ghg2009_grown_core=rbind(ghg2009_grown_core, rep(0,(dim(ghg2009_grown)[2])))
ghg2009_grown_core[(dim(ghg2009_grown_core)[1]),c("Mineral-a","Chemi-a","I_S-a")]=c(287.6894,1.4893,1.6919)
rownames(ghg2009_grown_core)[(dim(ghg2009_grown_core)[1])]="Process"
ghg2009_grown_core=rbind(ghg2009_grown_core,colSums(ghg2009_grown_core))
rownames(ghg2009_grown_core)[(dim(ghg2009_grown_core)[1])]="Total"

#(7) take out zero cells. add process CO2 (baby)
ghg2009_baby_core=ghg2009_baby[c("Roil-c","Coal-c","Gas-c"),]/100
ghg2009_baby_core["Coal-c","coal-a"]=0
ghg2009_baby_core=ghg2009_baby_core*(ghg2009_baby_core>=1)
ghg2009_baby_core=rbind(ghg2009_baby_core, rep(0,(dim(ghg2009_baby_core)[2])))
ghg2009_baby_core[(dim(ghg2009_baby_core)[1]),"eint-a"]=287.6894
rownames(ghg2009_baby_core)[(dim(ghg2009_baby_core)[1])]="Process"
ghg2009_baby_core=rbind(ghg2009_baby_core,colSums(ghg2009_baby_core))
rownames(ghg2009_baby_core)[(dim(ghg2009_baby_core)[1])]="Total"

#write.csv(ghg2009_grown_core,file="G://work_2013/Papers/Paper1/ghg2009_gw_core.csv")
#write.csv(ghg2009_baby_core,file="G://work_2013/Papers/Paper1/ghg2009_baby_core.csv")





Gen_CO2.1=function(EB_G,CF,RowName,IO_Domestic,IO_Import,IO_whole,EC,Oil.Data,YR){
# step 1 : Gross Calorific Value E balance => Net Calorific Value Energy E balance
## Multiply Gross Calorific Value E balance (EB_G) with  conversion factor = [Net Calorific Value Energy Coversion Factor/Gross Calorific Value Energy Conversion Factor] (CF) for each fule/Electricity
## Energy Balance is obtained from "Yearbook of Energy Statistics" (KEEI)
## Conversion Factor is obatiend from 'Greenhouse Gas Inventrory and Research Center (2011). "National Greenhouse Gas Inventory Report of Korea"'
## For Solvent/Asphalt/Parafiin-Wax, conversion factor is zero/For 'Other products' Conversion Factor is 0.95935 for no good reason.
## For JA-1, JA-2, AVI-G, conversion factor for Jet is used
## For Hydro, Nuclear, Electricity,Heat, Renewable, conversion factor is set as 1

# Setting time
YR_1=YR-1

### step 1-1 : Loading data, conversion factor
### step 1-2: convert into matrix and post muliply with diag CF (scale up columnes with CF)
EB_N=(as.matrix(EB_G))%*%(diag(as.numeric(CF)))
### step 1-3: add subsum and total columnes 
EB_N=data.frame(EB_N)
colnames(EB_N)=colnames(EB_G)

Anthra_N=EB_N$Anthra_D+EB_N$Anthra_Im
Bit_N=EB_N$Bit_C+EB_N$Bit_S                        
Coal_N=Anthra_N+Bit_N
E_USE_N=rowSums(EB_N[,5:13])
LPG_N=rowSums(EB_N[,14:15])
NE_USE_N=rowSums(EB_N[,16:22])
Petro_N=E_USE_N+LPG_N+NE_USE_N
Total_N= Coal_N+Petro_N+rowSums(EB_N[,23:29])
EB_N_total=cbind(Coal_N,Anthra_N,EB_N[,1:2],Bit_N,EB_N[,3:4],Petro_N,E_USE_N,EB_N[,5:13],LPG_N,EB_N[,14:15],NE_USE_N,EB_N[,16:29],Total_N)
 
Anthra_G=EB_G$Anthra_D+EB_G$Anthra_Im
Bit_G=EB_G$Bit_C+EB_G$Bit_S
Coal_G=Anthra_G+Bit_G
E_USE_G=rowSums(EB_G[,5:13])
LPG_G=rowSums(EB_G[,14:15])
NE_USE_G=rowSums(EB_G[,16:22])
Petro_G=E_USE_G+LPG_G+NE_USE_G
Total_G= Coal_G+Petro_G+rowSums(EB_G[,23:29])
EB_G_total=cbind(Coal_G,Anthra_G,EB_G[,1:2],Bit_G,EB_G[,3:4],Petro_G,E_USE_G,EB_G[,5:13],LPG_G,EB_G[,14:15],NE_USE_G,EB_G[,16:29],Total_G)

### setp 1-4: add row names

row.names(EB_N_total)=c("	Domestic","Import","PetroPro","PetroIm","Export","IntBunk","Stock_d","f_Stock	","e_Stock","Diff","P_con","Transform","E_Gene","Heating","Gas_Manu","OwnUse","F_con","Ind","AFF","Mine","Manu","Food","Textile","Wood","PPP","PetChem","NonMetal","IS","Nonferr","FabMetal","Other_Manu","Other_E","	Const","Transport","Rail_t","Land_t","Water_t","Air_t","Residential","Commercial","Public")	

row.names(EB_G_total)=c("	Domestic","Import","PetroPro","PetroIm","Export","IntBunk","Stock_d","f_Stock	","e_Stock","Diff","P_con","Transform","E_Gene","Heating","Gas_Manu","OwnUse","F_con","Ind","AFF","Mine","Manu","Food","Textile","Wood","PPP","PetChem","NonMetal","IS","Nonferr","FabMetal","Other_Manu","Other_E","	Const","Transport","Rail_t","Land_t","Water_t","Air_t","Residential","Commercial","Public")	

row.names(EB_N)=row.names(EB_N_total)
row.names(EB_G)=row.names(EB_G_total)

# step 2 : Net Calorific Value E balance => Total E demand by type of fuel/Electricity generation

##step 2-1: Construction E_Demand by pick EB_N["P_con",]

E_Demand=EB_N["P_con",]

## step 2-2: Modification
### (i) TownGas Primanry consumption is set as zero b/c Towngas is generated using LNG. 
###     But, burning Towngas is a primary source of GHG emission. 
###      So, we assign Primary Consumption of LNG to Towngas 

E_Demand$TownGas=E_Demand$LNG

### (ii) Electricity generation : E demand for Hydro and nuclear is discounted with (860/2150) Why(?)

E_Demand$Hydro=(860/2150)*E_Demand$Hydro
E_Demand$Nuclear=(860/2150)*E_Demand$Nuclear

### (iii) Electricity generation: E_demand for Electricity (Fossile fule) 
### a. P_con for Electricity doesn't exists, b/c Elec is generated from natural resources
### b. F_con for Electricity contains E_Demand for Hydro and Nuclear. We want to get rid of them
### c. F_con usually contains E_Demand for internatinal use. We take Export and IntBunk out to eleminate it
### d. To get F_con, we sum P_con and Transformation. The Transformation has Electricity generation, District Heating, GAS manufacturing.
###    If any of this item being negative, that means some of the Electricity generated has lost.
###    But lost Electricity still generates GHG. We want to put those lost Elecricity back to Elec Demand
###    If any of this item being positive, that means there is some byproduct with same amount of fuels are used. We ignore them
### e. Transfomration also includes OwnUse and Loss, usually negative. They are not demanded but used and generate GHG. We want them back.
### The formula for Electricity Demand become as follows
### Elect Demand=F_con-Export-IntBunk-Stock_d-E_Gene*I(E_Gene<0)-Heating*I(Heating<0)-Gas_manu*I(Gas_manu<0)-Ownuse&Loss-E_Demand$Hydro-E_Demand$Nuclear
 
 

#### step a. construct selection vector for relevent elements N. E_Demand
N.E_Demand=c("F_con","Export","IntBunk","Stock_d","E_Gene","Heating","Gas_Manu","OwnUse")
#### step b. obtain related values needed to construct formula
Elect_D=EB_N[N.E_Demand,]$Elect
#### step c. obtain indicator function 
condElect=c(rep(1,4),as.numeric(Elect_D[5:7]<0),1)
#### step d. obtain operation vector
funElect=c(1, rep(-1,7))
#### step e. obtain Elect Deamnd before hydro and nuclear adjustment
S.Elect_D=sum(Elect_D*condElect*funElect)
#### step f. substract Hydro and Nuclear
S.Elect_D=S.Elect_D-E_Demand$Hydro-E_Demand$Nuclear
#### step g. replace E_Demand$Elect
E_Demand$Elect=S.Elect_D

### (iv) Heating = P_con of Renewable (?)

E_Demand$Heat=E_Demand$Renewable

### (v) Renewable is set 0 (No GHG comes from Renewable)

E_Demand$Renewable=0 


##step 2-3: Add subsum

Anthra_E=E_Demand$Anthra_D+E_Demand$Anthra_Im
Bit_E=E_Demand$Bit_C+E_Demand$Bit_S
Coal_E=Anthra_E+Bit_E
E_USE_E=rowSums(E_Demand[,5:13])
LPG_E=rowSums(E_Demand[,14:15])
NE_USE_E=rowSums(E_Demand[,16:22])
Petro_E=E_USE_E+LPG_E+NE_USE_E
Total_E= Coal_E+Petro_E+rowSums(E_Demand[,23:29])
E.Demand_total=cbind(Coal_E,Anthra_E,E_Demand[,1:2],Bit_E,E_Demand[,3:4],Petro_E,E_USE_E,E_Demand[,5:13],LPG_E,E_Demand[,14:15],NE_USE_E,E_Demand[,16:29],Total_E)

## step 2-4: Odd Items. Yuntan/Crude Oil/Other Coal Product are in I-O, not in Energy Balance. 
## We assgin 'temporary energy demand' for these three items. 
## (Personally, I have doubts about this approach. 
##  a. Crude oil is mostly used in Oil production. Most of Energies are used  in Oil Product form. 
##  What we are missing here is the crude oil energy used in refineries. That wouldn't be much.
##  b. Yuntan, is a form of Anthracite coal. Already the domestic use of Anthracite coal is included in Energy balance. Add yuntan could be double counting
##  c. Other coal products. Very samll.  Can be ingonred.
##  d. Overall. The error in 'guessing' the proper energy formula for these three item could be larger than the error in ignoring these three item.
## Anyway, according to Kong, we assign energy demand for these three item as follows.

### (i). Crude Oil

#### a. Obtain the adjusted sum of [PetroProd + Diff*[PetroProd/Import]]*(1-I(Petroprod=0)) for coal products 
##### a.1. select appropriate statistics
Petroitems=c("Gasoline","Kerosene","Diesel",	"B.A",	"B.B",	"B.C",	"JA.1",	"JP.4",	"AVI.G",	"Propane",	"Butane",	"Naphta",	"Solvent",	"Asphalt",	"Lubricant",	"ParaWax",	"PetroCoke",	"OtherProd")
CrudeEitems=c("Import","PetroPro","Diff")
CrudeEdata=EB_N[CrudeEitems,Petroitems]

#####  a.2. Excluded data with zero Import 
##### Since Crude Oil is 100% imported and converted into Petroleum Products, the energy demand in the form of crude oil is included in imports of each petroleum product. 
##### Itmes with zero Petroproduct are actually not used at all.
S.CrudeEdata=CrudeEdata[,(CrudeEdata["PetroPro",]!=0)]

##### a.3. Obtain Crude Oil E demand included in each petroleum product.        
##### P_Con= Domestic Product+Import+Export+IntBunk+Stock_d+Diff. We exclude Export and IntBunk (not in Korea), and Stock_d is added afterwords. 
##### The rest is Domestic Product, Import, Diff. But we don't have Domeistic Product.
##### Import consists of Petro Product and Petro Import. 
##### We assume Petro Product E demand is E demand of Crude oil used to produce each Petro Product.
##### And we assume that Stat.difference mainly conmes from Petro product and Petro Import. We also assume that the size of each fraction in Stat.difference is proportional to the E_demand  
S.CrudeEdata=S.CrudeEdata["PetroPro",]+S.CrudeEdata["Diff",]*(S.CrudeEdata["PetroPro",]/S.CrudeEdata["Import",])

##### a.4. take sum of the Crude Oil E demand included in each petroleum product.
E.Demand_Crude_petroprod=sum(S.CrudeEdata)

##### a.5. adjust for efficiency in burning Crude oil. Only 99% of Crude Oil input are actually burn. The rest is wasted. 
######  We adjsut E.Demand_Crude by diving it with 0.99. More Crude Oil is burned than the energy. (Emit CO2)

E.Demand_Crude_petroprod=E.Demand_Crude_petroprod/0.99

#### b. Add Stock change using YES 2010. Lower oil stock implies more oil than production is used to fulfil the demand. Higher oil stock implies some oil producec is not used. 
#### We adjuste E.Deamnd_Crude with Stock change.
#### The 'stock chnage' here is the negative value of (2009 stock of Oil and Petro Product[7244 (1,000bbl)] -2008 stock of Oil and Petro Product[10723(1,000bbl)]  (2010 Yearbook of Energy Statistics, p.70)  
#### The data of Oil and Petro Product stock is recorded in  Annual Energy Statistics Yearbook. 
#### Since the Product stock is recorded in barel (1000 bbl), we should devided E.Demand_Crude with 7.33 to convert bbl into toe. 1ton=1/7.33 bbl not 1toe. But Crude Oil 1 ton generates 1toe.       

E.Demand_Crude_Stock_d=-1*(Oil.Data$Stocks[which(Oil.Data$year==YR)]-Oil.Data$Stocks[which(Oil.Data$year==YR_1)])/7.33

#### c. Add domestic production
#### We assume Domestic Production of crude oil = 321 (1000 bbl) (2010 Yearbook of Energy Statistics, p.70) 
#### Convert bbl into teo by dividing with 7.33

E.Demand_Crude_domestic=Oil.Data$Domestic[which(Oil.Data$year==YR)]/7.33
                                                                      

#### d. Obtain Import.
#### We exploit this identity P_con=Domestic Prod+Import+Stock_c+Diff.
#### Wa assume that P_con=E.Demand_Crude_petroprod, Diff=0
#### Then Import= P_con-Domestic Prod-Stock_c

E.Demand_Crude_Im= E.Demand_Crude_petroprod-E.Demand_Crude_domestic-E.Demand_Crude_Stock_d

#### e. obtain E.Demand for Crude Oil = Domestic Production + Import
E.Demand_Crude=E.Demand_Crude_Im+E.Demand_Crude_domestic
#### f. obtain E.Supply for Crude Oil = E.Demand-(stock change)
E.Supply_Crude=E.Demand_Crude_petroprod-E.Demand_Crude_Stock_d

### (ii). Coal Briquette (Yuntan) 
#### a. Residential Demand for Anthracite is used for the Demand for Yuntan.
E.Demand_Coalbriquette = sum(EB_N["Residential",1:2])
#### b. Inflate  Domestic Demand with [1+(Total Final Demand of Yuntan|Import I-O)/(Total Final Demand of Yuntan|Domestic I-O)] to take Coal briguette import into account.
####    The multiplication factor shows how much coal briqutte is imported when one unit of domesitc coal briqutte is used. 
##### load I-O information
##### Find appropriate value to construct inflating number
E.Demand_Coalbriquette =E.Demand_Coalbriquette*(1+(IO_Import[131,"Total.demand"]/IO_Domestic[131,"Total.demand"]))
#E.Demand_Coalbriquette =E.Demand_Coalbriquette*(1+(13458/396625))

### (iii) Other Coal Product (132 in 403)
#### a. E.Demand.Other_c_dom is E.Demand of 'domestically produced' Other coal products. 
#### All imported Bit coal used in the form of other coal products are assumed to be used in domestic other coal product production. (Same as Yuntan)
#### E.Demand.Other_c_dom is part of E.Demand of Bituminous Coal (Bit_N) not used in E_Gene
#### E.Demand.Other_c_dom =(E.Demand Bit Coal not used in E_Gene)*(proportion of E.Demand of Bit Coal used in Other Coal Product form among E.Demand of Bit Coal unsed in E_Gene)
#### (E.Demand of Bit Coal Other coal product/E.Demand of Bit Coal not used in E_Gene)
#### =(Money spent in Bit Coal in Othe Coal Product production|(I-O)/[Money Spent in Bit Coal-Money Spend in Bit Coal used in E_Gene](I-O))|(Domestic+Import)
E.Demand.Bit_residual=Bit_E-(-1)*min(c(EB_N_total["E_Gene","Bit_N"],0))
E.Demand.Other_c_dom=E.Demand.Bit_residual*(IO_whole[31,"Coke.and.other.coal.products"]/(IO_whole[31,"Total.demand"]-IO_whole[31,"Fire.power.generation"]))
#E.Demand.Other_c_dom=E.Demand.Bit_residual*(3650132/6032043)

#### b. Inflate E.Demand.Other_c_dom with [1+(Total Final Demand of Other Coal Product|Import I-O)/(Total Final Demand of Other Coal Product|Domestic I-O)) to take Other Coal Products into account.
E.Demand.Other_c=E.Demand.Other_c_dom*(1+(IO_Import[132,"Total.demand"]/IO_Domestic[132,"Total.demand"])) 


## step 2-5: Arrange to obtain Total E.Demand_total matching I-O (E.Demand_IO)
## For Jet.Oil (135), We only use JA.1
## For Heavy.Oil(138), We use the sum of B.A, B.B, B.C
## For Misc. petroleum.refinery products(141), we use the sum of PetoCoke and OtherProd
## For Fire.power.generation, we use Elect*((Total Demand of Fire Power generation(299)|IO_whole/(Total Demand of Fire Power Generation(299) +Total Demand of Other generation(301)]
## For Other.generation, we use Elect*((Total Demand of Other generation(299)|IO_whole/(Total Demand of Fire Power Generation(299) +Total Demand of Other generation(301)]
## (Since we don't have Other generation in EB, and we don't have Renewable in IO, 
##  we assume that Electricity demand excluding hydro and nuclear is the sum of Fire generataion and Other generation)

### pick indexes of related items saved in E.Demand_total (30~33,131~141,298~303)

index1.c=c("Anthra_E","Bit_E")#Anthracite(30)/Bituminous.coal(31)
index1=match(index1.c,colnames(E.Demand_total))

index2.c=c("Naphta","Gasoline","JA.1","Kerosene","Diesel") #Naphta(133)/Gasoline(134)/Jet.Oil(135)/Kerosene(136)/Light.oil(137)
index2=match(index2.c,colnames(E.Demand_total))
  
index3.c=c("B.A","B.B","B.C")#B.A+B.B+B.C = Heavy.Oil (138)
index3=match(index3.c,colnames(E.Demand_total))

index4.c=c("LPG_E","Lubricant") #Liqufied.petroleum.gas(139)/Lubricants(140)
index4=match(index4.c,colnames(E.Demand_total))

index5.c=c("PetroCoke","OtherProd")#PetroCoke+OtherProd = Misc.petroleum.refinery. products(141)
index5=match(index5.c,colnames(E.Demand_total))

index6.c=c("TownGas","Heat") #Manufactured.gas.supply(302)/Steam.and.hot.water.supply(303)
index6=match(index6.c,colnames(E.Demand_total))

index_IO=c(30:33,131:141,298:303)

### calculate share of Fire/Other generation to split E.Demand_total$Elect
Elect.Fire=IO_whole[299,"Total.demand"]/(IO_whole[299,"Total.demand"]+IO_whole[301,"Total.demand"])
Elect.Other=IO_whole[301,"Total.demand"]/(IO_whole[299,"Total.demand"]+IO_whole[301,"Total.demand"])

### Construdt E.Demand_IO
E.Demand_IO=c(E.Demand_total[index1],E.Demand_Crude,E.Demand_total$LNG,E.Demand_Coalbriquette,E.Demand.Other_c,E.Demand_total[index2],sum(E.Demand_total[index3]),E.Demand_total[index4],sum(E.Demand_total[index5]),E.Demand_total$Hydro,(E.Demand_total$Elect*Elect.Fire),E.Demand_total$Nuclear,(E.Demand_total$Elect*Elect.Other),E.Demand_total[index6]) 
E.Demand_IO=data.frame(E.Demand_IO)

### Provide col names
IO_name=colnames(IO_whole)
colnames(E.Demand_IO)=IO_name[index_IO]

### Save as k*1 matrix
E.Demand_IO.m=as.matrix((as.numeric(E.Demand_IO)),nrow=length(E.Demand_IO))
rownames(E.Demand_IO.m)=IO_name[index_IO]
colnames(E.Demand_IO.m)=c("Total.demand")


# step 3:  Allocated intermediate demand/final demand of Total E generated by each type of fuel/Electricity generation                                                                                                    
## In this step, we redistribute "Total.demand" according to I-O intermediate demand /final demand. 
## [Residuals]For five fuels, we put some important components of E demand into specified I-O cells, and redistribute the residual
##  acoording to the rest of I-O. These Five Fuels are Anthracite, Bituminous.coal, Naphta, Heavy.oil, Towngas. 
## [Overalls] For the rest, we redistribute "Total.Demand" itself according to I-O

## step 3-1: Obtain relevent data,  add row names and numbers
Energy_IO=IO_whole[index_IO,]
Energy_IO=cbind(data.frame(index_IO),Energy_IO)
rownames(Energy_IO)=rownames(E.Demand_IO.m)
## step 3-2: Separate Residual and Overall
index_residual=c("Anthracite", "Bituminous.coal","Naphtha","Heavy.oil","Manufactured.gas.supply")
index.residual.m=match(index_residual,rownames(Energy_IO))
#index.residual.m=c(1,2,7,12,20)

Energy_IO_residual=Energy_IO[index_residual,]
Index_IO_residual=Energy_IO_residual$index_IO
Energy_IO_overall=Energy_IO[-index.residual.m,]
Index_IO_overall=Energy_IO_overall$index_IO

E.Demand_IO.residual=E.Demand_IO.m[index.residual.m,]
E.Demand_IO.overall=E.Demand_IO.m[-index.residual.m,]

## step 3-3. overall. 
### (i) select "Total.supply" column from Energy_IO_residual (TS_IO_overall)
TS_IO_overall=Energy_IO_overall[,"Total.supply"]
###  (ii) devide E.Demand_IO_overall by "Total. supply" (optain E.demand per moneytary unit. E.Demand_overall_money)
E.Demand_IO.overall_money=E.Demand_IO.overall/TS_IO_overall
###   (iii) multiply Energy_IO_overall with E.demand per monetary unit
Energy_BE_overall=Energy_IO_overall*E.Demand_IO.overall_money
####   (iv) obtain sum of intermediate demand and sum of final demand. check matching
Energy_BE_overall[,"Total.intermediate.demand"]=rowSums(Energy_BE_overall[,2:404])
Energy_BE_overall[,"Total.final.demand"]=rowSums(Energy_BE_overall[,406:411])

## step 3-4 residual.
### (i) prepare monetary residual
#### a. identify designated cells for 'residual members'
index.anth=c("Coal.briquettes","Fire.power.generation")
index.bit="Fire.power.generation"
index.Naphtha=143:172
index.heavyoil="Fire.power.generation"
index.towngas=c("Fire.power.generation","Steam.and.hot.water.supply")
#### b. take out the money for designated cells from total money spent 
money.out=c(sum(Energy_IO_residual["Anthracite",index.anth]),Energy_IO_residual["Bituminous.coal",index.bit],sum(Energy_IO_residual["Naphtha",index.Naphtha]),Energy_IO_residual["Heavy.oil", index.heavyoil], sum(Energy_IO_residual["Manufactured.gas.supply",index.towngas]))
TS_IO_residual=Energy_IO_residual[,"Total.supply"]-money.out
#### c.indentify energy for designated cells
anth.desig.bit=EB_N_total["Residential","Anthra_N"]
anth.desig.E_Gene=-1*min(c(EB_N_total["E_Gene","Anthra_N"],0))
anth.desig=c(anth.desig.bit,anth.desig.E_Gene)
bit.desig=-1*min(c(EB_N_total["E_Gene","Bit_N"],0))
Naphtha.desig=EB_N_total["PetChem","Naphta"]
Heavoil.desig=-1*min(c(EB_N_total["E_Gene","B.A"],0))-1*min(c(EB_N_total["E_Gene","B.B"],0))-1*min(c(EB_N_total["E_Gene","B.C"],0))
Towngas.desig.E_Gene=-1*min(c(EB_N_total["E_Gene","LNG"],0))-1*min(c(EB_N_total["E_Gene","TownGas"],0))
Towngas.desig.Heating=-1*min(c(EB_N_total["Heating","LNG"],0))-1*min(c(EB_N_total["Heating","TownGas"],0))
Towngas.desig=c(Towngas.desig.E_Gene,Towngas.desig.Heating)
#### d.take designated energy from Total E. Demand 
E.Demand_IO.residual_d=c(sum(anth.desig),bit.desig,Naphtha.desig,Heavoil.desig,sum(Towngas.desig))
E.Demand_IO.residual_r=E.Demand_IO.residual-E.Demand_IO.residual_d
### (ii) devide E.Demand_IO_residual_r by TS_IO_residual (obtain E. demand per monetary unit for residual. E.Demand_residual_money)
E.Demand_IO.residual_money=E.Demand_IO.residual_r/TS_IO_residual
### (iii) multiply Energy_IO_residual wit E.Demand_IO.residual_money (distribute residual energy)
Energy_BE_residual=Energy_IO_residual*E.Demand_IO.residual_money
### (iv) Replace designated cell with designated Energy demand
Energy_BE_residual["Anthracite",index.anth]=anth.desig
Energy_BE_residual["Bituminous.coal",index.bit]=bit.desig
Energy_BE_residual["Heavy.oil", index.heavyoil]=Heavoil.desig
Energy_BE_residual["Manufactured.gas.supply",index.towngas]=Towngas.desig
Energy_BE_residual["Naphtha",index.Naphtha]=Naphtha.desig*(Energy_IO_residual["Naphtha",index.Naphtha]/(sum(Energy_IO_residual["Naphtha",index.Naphtha]))) 
### (v) Replace subsum columns with actual subsums
Energy_BE_residual[,"Total.intermediate.demand"]=rowSums(Energy_BE_residual[,2:404])
Energy_BE_residual[,"Total.final.demand"]=rowSums(Energy_BE_residual[,406:411])
Energy_BE_residual[,"Total.demand"]=Energy_BE_residual[,"Total.intermediate.demand"]+Energy_BE_residual[,"Total.final.demand"]
Energy_BE_residual[,"Total.supply"]=Energy_BE_residual[,"Total.demand"] 

## step 3-5  merge and Replace "Gross.domestic.output.Basic price", "Import.Basic.price" with 
###  Gross.domestic.output.Basic.price = sum of intermediate demand for domestic energy + sum of final demand for domestic energy
###  * intermediate demand for domestic energy : 0 if intermediate demand in total IO is 0, 
###  otherwise intermediate demand for energy *(intermediate domesitc demand in domestic IO/intermediate domestic demand in total IO(Energy_IO))
###  * final demand for domestic energy: 0 if final demand in total IO is 0,
###  otherwise final demand for energy*(final demand in doemstic IO/final demand in total IO)

#### (i) merge Engery_BE_residual and Energy_BE_overall
##### a. put back indexes
Energy_BE_overall$index_IO=Index_IO_overall        
Energy_BE_residual$index_IO=Index_IO_residual
##### b.merge and sort
Energy_BE=rbind(Energy_BE_overall,Energy_BE_residual)
Energy_BE=Energy_BE[order(Energy_BE$index_IO),]

##### (ii). calculate intermediate demand_dometic/intermeidate demand_whole
###### a. Obtain Energy_IO_domestic
Energy_IO_Domestic=IO_Domestic[index_IO,]
Energy_IO_Domestic=cbind(data.frame(index_IO),Energy_IO_Domestic)
rownames(Energy_IO_Domestic)=rownames(E.Demand_IO.m)
###### b.  Obtain Domestic/Whole share
Energy_IO.m=as.matrix(Energy_IO)
Energy_IO_Domestic.m=as.matrix(Energy_IO_Domestic)
Energy_IO_Domestic_share=Energy_IO.m
Energy_IO_Domestic_share=ifelse(Energy_IO.m==0,0,(Energy_IO_Domestic.m/Energy_IO.m)) # set Energy_IO_Domestic_share as zero if Energy_IO is zero
Energy_IO_Domestic_share=data.frame(Energy_IO_Domestic_share)

#### (iii) Obtain E demand for domestic production/Imports

###### a. Multiply Energy demand with (Domestic/Total) factor
Energy_BE_Domestic=Energy_BE*Energy_IO_Domestic_share
###### b.  Obtain subsums
Energy_BE_Domestic$Total.intermediate.demand=rowSums(Energy_BE_Domestic[,2:404])
Energy_BE_Domestic$Total.final.demand=rowSums(Energy_BE_Domestic[,406:411])

Energy_BE_Domestic$Total.demand=Energy_BE_Domestic$Total.intermediate.demand+Energy_BE_Domestic$Total.final.demand
##### c. Replace GDP and Imports in Energy_BE with calculated subsums
Energy_BE$Gross.domestic.output.Basic.price.=Energy_BE_Domestic$Total.demand
Energy_BE$Imports..Basic.price.=Energy_BE$Total.demand-Energy_BE$Gross.domestic.output.Basic.price.
                                                                 
# step 4: Convert intermediate demand of E/final demand of E/ E demand assigned to domestic production(import) into CO2 emission
## In this step, we multiply Energy demand with emission factor (unit =toe)
## For Naphtha and Lubricant, some of carbon contained in these materials would sink to the final output. The rest would be released.
## For generation (Electricity, Heat), carbon is released in the generation process only. 
## What is emitted after generation is assumed as the hidden emission in generation. Using electricity or heat wouldn't emit additional GHG.
## So, we set zero emisson factor for (power.generation, Steam.and.hot.water.supply)
## After multiplication, we set zero from emission in inrease.in.stock/Exports. Those are not emitted within our territory.
## Finally, adjuste Total.intermediate.demand/Total.final.demand/Total.demand. These subsomes should be recalculated.

## Step 4-1. prepare emission factor 

### (ii) convert carbon emission factor(ton/toe) to CO2 emission factor
EC$CO2.emission=(44/12)*EC$C.emission
### (iii)adjust for sinking in Naphtha and Lubricants
EC$EF=(1-EC$Sinking)*EC$CO2.emission

##Step 4-2. Obtaion emission factor*Energy
### Multiply Energy demand with emission factor
ECO2=Energy_BE*EC$EF
### setting increase.in.stock/Exports. as zero
ECO2$Increase.in.stocks=0
ECO2$Exports=0
### Adjusting subsums and obtain totals.
ECO2$Total.intermediate.demand=rowSums(ECO2[,2:404])
ECO2$Total.final.demand=rowSums(ECO2[,406:411])
ECO2$Total.demand=ECO2$Total.intermediate.demand+ECO2$Total.final.demand

## Step 4-3. Excluding double counting
### indexing double counting cells.
index.Anth.Exc=132:133 # Coal.briquettes and Coke.and.other.coal.products  (CO2 is emitted when 2nd products are used)
index.Bit.Exc=132:133  # Coal.briquettes and Coke.and.other.coal.products  (CO2 is emitted when 2nd products are used)
index.Crude.Exc=134:142 # Petro products  (CO2 is emitted when 2nd products are used)
index.Ngas.Exc=303 #Manufactured.gas.supply (CO2 is emitted when Manufactured gas is used)
index.LPG.Exc=303  #Manufactured.gas.supply (CO2 is emitted when Manufactured gas is used)

### setting indexed cells as zero
ECO2["Anthracite",index.Anth.Exc]=0                        
ECO2["Bituminous.coal",index.Bit.Exc]=0                  
ECO2["Crude.petroleum",index.Crude.Exc]=0                   
ECO2["Natural.gas",index.Ngas.Exc]=0                      
ECO2["Liquefied.petroleum.gas",index.LPG.Exc]=0    

### Adjusting subsums and obtain totals.
ECO2$Total.intermediate.demand=rowSums(ECO2[,2:404])
ECO2$Total.final.demand=rowSums(ECO2[,406:411])
ECO2$Total.demand=ECO2$Total.intermediate.demand+ECO2$Total.final.demand

## Step 4-4. Obtain CO2 emission by industry
ECO2.slim=cbind(Energy_BE$index_IO,ECO2[,2:413]) # cut off final demand section after Total.demand/shed industry index
ECO2.ind=colSums(ECO2.slim[-1])
ECO2.ind.5dg=as.numeric(formatC(ECO2.ind,format="f",digits=5))
names(ECO2.ind.5dg)=names(ECO2.ind)
ECO2.slim=rbind(ECO2.slim,c(0,ECO2.ind.5dg))
ECO2.result=list(ECO2.slim,ECO2.ind.5dg) 
return(ECO2.result)
} 

# step 5: add CO2 in process                                   
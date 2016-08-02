sam_br=read.csv(file="b_sam_br_ng_1202.csv",header=T,as.is=T)
ghg=read.csv(file="GHG_BR_p_1202.csv",header=T,as.is=T)
samsize=c(7,1,2,1,1,1,4,1,1,1)
filename="set_br_20151211.txt"

source("HC.r")
sink(file=filename)
# (1) preperation


## set column and low names with variable names
### take the first column containing variable names and remove empty spaces
samnames=sam_br[,1]
samnames=sub("[[:space:]]+$","",samnames)
ghgnames=ghg[,1]
ghgnames=sub("[[:space:]]+$","",ghgnames)
#ghgnames=c(paste(ghgnames[1:(length(ghgnames)-2)], "-c",sep=""),ghgnames[(length(ghgnames)-1):length(ghgnames)])
#ghgnames=c("GAS-c",  "ROIL-c", "CoalPro-c", "Coal-c",   "Process", "Total")
### remove the first column from the data seet
sam_br=sam_br[,(2:(dim(sam_br)[2]))]
ghg=ghg[,(2:(dim(ghg)[2]))]
### replace column names with variable names
colnames(sam_br)=samnames
### replace row names with variable names
rownames(sam_br)=samnames
rownames(ghg)=ghgnames
sam_br.raw=sam_br
### indicating size
ind=7
green=1
fac=2
#hh=1
hh=1
gov=1
Nres=1
tax=4
S_I=1
ROW=1
Total=1

ind=samsize[1]
green=samsize[2]
fac=samsize[3]
hh=samsize[4]
gov=samsize[5]
Nres=samsize[6]
tax=samsize[7]
S_I=samsize[8]
ROW=samsize[9]
Total=samsize[10]

#ACTP: all entry
ACTP=unique(c(samnames,ghgnames))
setwords=function(setname,domain,x){paste(setname,"(",domain,") /", HC.com(x)," /",sep='')}
ACTP.words=paste("ACTP"," /",HC.com(ACTP)," /",sep='')
cat(ACTP.words, sep='\n')


#ACT=ACT-"process"
#ACT=setdiff(ACTP,setdiff(ghgnames,samnames))
ACT=ACTP[ACTP!="Process"]
ACT.words=setwords("ACT", "ACTP",ACT)
#ACNT: null set
ACNT.words="ACNT(ACT)"
cat(ACT.words, sep='\n')
cat(ACNT.words, sep='\n')

#ghg dataset entry
#GOT=ghgnames
#GO=ghgnames[ghgnames!="Total"]
#GOT.words=setwords("GOT","ACTP",GOT)
#GO.words=setwords("GO","GOT",GO)
#cat(GOT.words, sep='\n')
#cat(GO.words, sep='\n')


#modification of whole entry
AC=ACT[ACT!="Total"]
ACNGT=ACT[(ACT!="CO2-a")&(ACT!="CO2-c")]
ACNGA=ACNGT[ACNGT!="Total"]
AC.words=setwords("AC","ACT",AC)
ACNGT.words=setwords("ACNGT","ACT",ACNGT)
ACNGA.words=setwords("ACNGA","AC",ACNGA)
cat(AC.words, sep='\n')
cat(ACNGT.words, sep='\n')
cat(ACNGA.words, sep='\n')


#Activity and Commodity
Activity=AC[1:ind]
Commodity=AC[(ind+green+1):(2*(ind)+green)]
Activity.words=setwords("A","AC",Activity)
Commodity.words=setwords("C","AC",Commodity)
cat(Activity.words,sep='\n')
cat(Commodity.words,sep='\n')

# Commodities with positive exports
#CE=Commodity[(sam_br[Commodity,"ROW"]>0)]
#CEN=setdiff(Commodity,CE)
#CE.words=setwords("CE","C",CE)
#CEN.words=setwords("CEN","C",CEN)
#cat(CE.words,sep='\n')
#cat(CEN.words,sep='\n')

# Commodites with positive imports
#CM=Commodity[t(sam_br["ROW",Commodity]>0)]
#CMN=setdiff(Commodity,CM)
#CM.words=setwords("CM","C",CM)
#CMN.words=setwords("CMN","C",CMN)
#cat(CM.words,sep='\n')
#cat(CMN.words,sep='\n')


#CD: Commodities with postivie domestic production(Activity)
#IO AV section
AV=sam_br[,Activity]
#getting rid of total row
#AV=AV[1:(dim(AV)[1]-1),]
#take out activitties with total production zero
#CD=Commodity[colSums(AV)>0]
#CDN=setdiff(Commodity,CD)
#CD.words=setwords("CD","C",CD)
#CDN.words=setwords("CDN","C",CDN)
#cat(CD.words, sep='\n')
#cat(CDN.words, sep='\n')

#CGov: Commodity with positive government consumption
#CInv: Commodity with positive investment demand
#CPC: Commodity with positive household consumption
#Final consumption could be negative (ex. Inventory)
#CGov=Commodity[(sam_br[Commodity,"GoV"]!=0)]
#CInv=Commodity[(sam_br[Commodity,"S-I"]!=0)]
#CPC=Commodity[(sam_br[Commodity,"Household"]!=0)]
#CGov.words=setwords("CGov","C",CGov)
#CInv.words=setwords("CInv","C",CInv)
#CPC.words=setwords("CPC","C",CPC)
#cat(CGov.words, sep='\n')
#cat(CInv.words, sep='\n')
#cat(CPC.words, sep='\n')

#ENC : Energy commodities
ENC=c("GASHeat-c", "OIL-c","COAL-c" )
ENCN=setdiff(Commodity,ENC)
ENC.words=setwords("ENC","C",ENC)
ENCN.words=setwords("ENCN","C",ENCN)
cat(ENC.words, sep='\n')
cat(ENCN.words, sep='\n')

#Sfuel : Solid Fuel commodities
Sfuel=c("COAL-c")
Sfuel.words=setwords("Sfuel","C",Sfuel)
cat(Sfuel.words,sep='\n')

#Lfuel : Liqoud Fuel commodities
Lfuel=setdiff(ENC,Sfuel)
Lfuel.words=setwords("Lfuel","C",Lfuel)
cat(Lfuel.words,sep='\n')


#ELECC: Electricity, ENERGY=ENC and ELECC, Material (M)=commodity - ENergy

ELECC=c("ELEC-c")
ENERGY=c(ELECC,ENC)
M=setdiff(Commodity,ENERGY)
M.words=setwords("M","C",M)
cat(M.words, sep='\n')

#For each commodity, define uique set
#Comset_name_1=paste(substring(setdiff(Commodity,ENC),1,6),"C(C) /",sep="")
Comset_name_1=paste(substring(setdiff(Commodity,ENC),1,(nchar(setdiff(Commodity,ENC))-2)),"C(C) /",sep="")
Comset_name_2=paste(substring(ENC,1,(nchar(ENC)-2)),"ac(C) /",sep="")
Comset_name=c(Comset_name_1,Comset_name_2)
Comset_set=c(setdiff(Commodity,ENC),ENC)
Comset.words=paste(paste(Comset_name,Comset_set,sep=""),"/",sep="")
cat(Comset.words, sep='\n')

                            
#Contingency set of fuel mix for each Activity

QINT_FUEL=sam_br[ENC,Activity]

fuelpattern=mapply(function(x){rownames(QINT_FUEL)[x>0]},QINT_FUEL)

Fuel.A={}
for (i in (1:ind)){
                            #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
                            Fuel.A_i=paste(paste(names(fuelpattern)[i],HC.com(fuelpattern[[i]]),sep=".("),")")
                            Fuel.A=rbind(Fuel.A,Fuel.A_i)
}
cat("FuelA(A,C)",sep="\n")
cat("/",sep="\n")
cat(Fuel.A,sep="\n")
cat("/",sep="\n")



#Greenhouse gas.Commodity
GC="CO2-c"
GC.words=setwords("GC","AC",GC)
cat(GC.words, sep='\n')

# contingency set of sfuel mix for Activity

QINT_SFUEL=sam_br[Sfuel,Activity]
sfuelpattern=lapply(QINT_SFUEL,function(x){rownames(QINT_SFUEL)[x>0]})
sfuelpattern=sfuelpattern[mapply(FUN=length,sfuelpattern)>0]
SFuel.A={}
for (i in (1:length(sfuelpattern))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  SFuel.A_i=paste(paste(names(sfuelpattern)[i],HC.com(sfuelpattern[[i]]),sep=".("),")")
  SFuel.A=rbind(SFuel.A,SFuel.A_i)
}
cat("SFuelA(A,C)",sep="\n")
cat("/",sep="\n")
cat(SFuel.A,sep="\n")
cat("/",sep="\n")

# contingency set of Lfuel mix for Activity

QINT_LFUEL=sam_br[Lfuel,Activity]
Lfuelpattern=lapply(QINT_LFUEL,function(x){rownames(QINT_LFUEL)[x>0]})
Lfuelpattern=Lfuelpattern[mapply(FUN=length,Lfuelpattern)>0]
LFuel.A={}
for (i in (1:length(Lfuelpattern))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  LFuel.A_i=paste(paste(names(Lfuelpattern)[i],HC.com(Lfuelpattern[[i]]),sep=".("),")")
  LFuel.A=rbind(LFuel.A,LFuel.A_i)
}
cat("LFuelA(A,C)",sep="\n")
cat("/",sep="\n")
cat(LFuel.A,sep="\n")
cat("/",sep="\n")

# commodity-activity mapping
QINT_AC=sam_br[Activity,Commodity]
ACpattern=lapply(QINT_AC,function(x){rownames(QINT_AC)[x>0]})
ACpattern=ACpattern[mapply(FUN=length,ACpattern)>0]
AC.A={}
for (i in (1:length(ACpattern))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  AC.A_i=paste(paste(names(ACpattern)[i],HC.com(ACpattern[[i]]),sep=".("),")")
  AC.A=rbind(AC.A,AC.A_i)
}
cat("XPXC(C,A)",sep="\n")
cat("/",sep="\n")
cat(AC.A,sep="\n")
cat("/",sep="\n")


# Activities holing each M and Electricity as intermediate input

QINT=sam_br[c(ELECC,M),Activity]
T_QINT=data.frame(t(QINT))
colnames(T_QINT)=rownames(QINT)
Positive.int.demand=lapply(T_QINT,function(x){rownames(T_QINT)[x>0]})
Positive.int.demand=Positive.int.demand[mapply(FUN=length,Positive.int.demand)>0]
CA.A={}
for (i in (1:length(Positive.int.demand))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  CA.A_i=paste(paste(names(Positive.int.demand)[i],HC.com(Positive.int.demand[[i]]),sep=".("),")")
  CA.A=rbind(CA.A,CA.A_i)
}
cat("XAPA(C,A)",sep="\n")
cat("/",sep="\n")
cat(CA.A,sep="\n")
cat("/",sep="\n")

# Activities holing each fuel as intermediate input
QCE=sam_br[ENC,Activity]
T_QCE=data.frame(t(QCE))
colnames(T_QCE)=rownames(QCE)
Positive.enc.demand=lapply(T_QCE,function(x){rownames(T_QCE)[x>0]})
Positive.enc.demand=Positive.enc.demand[mapply(FUN=length,Positive.enc.demand)>0]
ENCA.A={}
for (i in (1:length(Positive.enc.demand))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  ENCA.A_i=paste(paste(names(Positive.enc.demand)[i],HC.com(Positive.enc.demand[[i]]),sep=".("),")")
  ENCA.A=rbind(ENCA.A,ENCA.A_i)
}
cat("XEPA(C,A)",sep="\n")
cat("/",sep="\n")
cat(ENCA.A,sep="\n")
cat("/",sep="\n")

#Non household (Domestic) Final demand: FD

FD=c("GoV","S-I")
cat(paste(paste("FD(ACT) ",HC.com(FD),sep="/"),"/"),sep="\n")

#Final demand mix for non household institutions   

XFA=sam_br[Commodity,match(FD,colnames(sam_br))]
T_XFA=data.frame(t(XFA))
colnames(T_XFA)=rownames(XFA)
Positive.fin.demand=lapply(T_XFA,function(x){rownames(T_XFA)[x!=0]})
Positive.fin.demand=Positive.fin.demand[mapply(FUN=length,Positive.fin.demand)>0]
XFA.A={}
for (i in (1:length(Positive.fin.demand))){
  #Fuel.A_i=HC.com(names(fuelpattern.v)[fuelpattern.v==fuelpattern_u.v[i]])
  XFA.A_i=paste(paste(names(Positive.fin.demand)[i],HC.com(Positive.fin.demand[[i]]),sep=".("),")")
  XFA.A=rbind(XFA.A,XFA.A_i)
}
cat("FD_C(C,FD)",sep="\n")
cat("/",sep="\n")
cat(XFA.A,sep="\n")
cat("/",sep="\n")


#Roil singlton as a subset of activity
#ROILA="ROil-a"
#ROILA.word=paste(sub("-a","a(A)",ROILA),"/",ROILA,"/",sep="")
#cat(ROILA.word,sep='\n')

#Factors as a subset of AC
FACTOR.endow=ACT[(2*(ind+green)+1):(2*(ind+green)+fac)]
Labor=FACTOR.endow[1]
Capital=FACTOR.endow[2]
FACTOR.words=paste("F(AC) /",HC.com(FACTOR.endow),"/",sep='')
Labor.words=paste("L(F) /",Labor,"/",sep='')
Capital.words=paste("K(F) /",Capital,"/",sep='')

FACTOR.words=setwords("F","AC",FACTOR.endow)
Labor.words=setwords("L","F",Labor)
Capital.words=setwords("K","F",Capital)

cat(FACTOR.words, sep='\n')
cat(Labor.words, sep='\n')
cat(Capital.words, sep='\n')




#Institution as a subset of AC

INS=ACT[(2*(ind+green)+fac+hh):(2*(ind+green)+fac+hh+gov+Nres+tax+S_I+ROW)]
INSD=setdiff(INS,INS[length(INS)])
INSDN=INS[length(INS)]
H=ACT[(2*(ind+green)+fac+1):(2*(ind+green)+fac+hh)]
INS.words=setwords("INS","AC",INS)
INSD.words=setwords("INSD","INS",INSD)
INSDN.words=setwords("INSDN","INS",INSDN)
H.words=setwords("H","INSD",H)

cat(INS.words,sep='\n')
cat(INSD.words,sep='\n')
cat(INSDN.words,sep='\n')
cat(H.words,sep='\n')

# Comodities with positive household demand XACH(H,C)

HC=data.frame(sam_br[Commodity,match(H,colnames(sam_br))])
colnames(HC)=H
rownames(HC)=Commodity
cat("XACH(H,C)",sep="\n")
cat("/",sep="\n")
cat(paste("Household",rownames(HC)[HC>0],sep="."),sep="\n")
cat("/",sep="\n")



#Alias
Alias1=c("ACT","INS","AC","A","C","F","H","GC","ENC")
Alias2=c("ACTPP","INSP","ACP","AP","CP","FP","HP","GCP","ENCPP")
Alias=cbind(Alias1,Alias2)
#Alias
Alias.words=paste("Alias(", Alias1, ",", Alias2, ");",sep="")
cat(Alias.words,sep='\n')

#set modification
Modification.words=c("ACNT(ACT)=yes;","ACNT('Total')=no;")
cat(Modification.words,sep='\n')
#}
sink()



$TITLE Hybrid model top down module trial version
$ONTEXT
- Oct 31. 2015
This model is for testing hybrid alogrithm only.
We have only 7 industries.
Production function and Utility function don't have any nesting structure.
Data : 2010 IO table. Basic price data.

- Sep. 22. 2015
Non nested version test-run is complete

-Sep. 23.24. 2015
Changing notation. Simplify equation using conditional set fuelA,XAPA,XEPA,FD_C

-Sep. 24. 2015 .Depreciation cost is subtracted from household taxable income.

Dep(H)=shr(K,H)*delta*K(S) (deprecation cost is proportional to household capital endowment)
Tax is imposed on factor income + residue sales income - deprecitaion cost
YD becomes (1-Tinsr)*(LY+KY+ResinC-Dep(H))
XAC is still mpc*YD
SH becomes mps*YD+Dep(H)
YH=sum(mpc.mpc*(1-TINSR)*(LY+KY+ResinC-Dep(H)))+mps*(1-TINSR)*(LY+KY+ResinC-Dep(H))+TINSR(LY+KY+ResinC-Dep(H))+Dep(H)

delta: deprecation ratio =Return rate *[[depreciation cost]/[depreciation cost+operation surplus]]=Return rate*0.39 (2010 IO basic price table)

- Sep 25. 2015. nesting begins.
step 1. Leontief input (Materials)
step 2. Labor

- Sep 26. 2015. nesting goes on
Step 3. Kapital
Setp 4. Electricity

-Setp.27. 2015. Final nest
Setp 5. Coal

-Oct. 7. Impose BR(2009) nesting structure
-Oct. 8 (Lafayett, IN USA). calibration /initiation are completed.



$OFFTEXT
OPTION SYSOUT=ON;
*OPTION limrow=50;
*SETS=================================
parameter
operation;
operation=0;
*operation=1;
SET
* All entry
ACTP /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a,CO2-a,ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c,CO2-c,Labor,Capital,Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW,process,Total /
* All entry -process
ACT(ACTP) /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a,CO2-a,ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c,CO2-c,Labor,Capital,Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW,process,Total /
* Null set
ACNT(ACT)
* ACT-Total
AC(ACT) /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a,CO2-a,ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c,CO2-c,Labor,Capital,Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW,process /
*ACT-{CO2-a,CO2-c}
ACNGT(ACT) /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a,ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c,Labor,Capital,Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW,process,Total /
*ACNGT-{Total}
ACNGA(AC) /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a,ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c,Labor,Capital,Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW,process /

A(AC) /ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a /
C(AC) /ELEC-c,GASHeat-c,OIL-c,COAL-c,ENIT-c,NEINT-c,AGRI-c /
*Fuel commodities
ENC(C) /GASHeat-c,OIL-c,COAL-c /
*Non Fuel commodities
ENCN(C) /ELEC-c,ENIT-c,NEINT-c,AGRI-c /
*Material commodities
M(C) /ENIT-c,NEINT-c,AGRI-c /
*Unique set for each commodity
ELECC(C) /ELEC-c/
ENITC(C) /ENIT-c/
NEINTC(C) /NEINT-c/
AGRIC(C) /AGRI-c/
GASHeatac(C) /GASHeat-c/
OILac(C) /OIL-c/
COALac(C) /COAL-c/
*Unique pattern of fuel mix and Activities using each mix
fuelA(A,C)
/
ELEC-a.(GASHeat-c,OIL-c,COAL-c)
GASHeat-a.(GASHeat-c,OIL-c,COAL-c)
OIL-a.(GASHeat-c,OIL-c)
COAL-a.(GASHeat-c,OIL-c,COAL-c)
ENIT-a.(GASHeat-c,OIL-c,COAL-c)
NEINT-a.(GASHeat-c,OIL-c,COAL-c)
AGRI-a.(GASHeat-c,OIL-c,COAL-c)
/
Sfuel(C)
/coal-c/
Lfuel(C)
/Oil-c,GASHeat-C/

SfuelA(A,C)
/
ELEC-a.(COAL-c)
GASHeat-a.(COAL-c)
COAL-a.(COAL-c)
ENIT-a.(COAL-c)
NEINT-a.(COAL-c)
AGRI-a.(COAL-c)
/
*Fuel0(C)
*/GASHeat-c,OIL-c,COAL-c/
*Fuel1(C)
*/GASHeat-c,OIL-c/
LfuelA(A,C)
/
ELEC-a.(GASHeat-c,OIL-c)
GASHeat-a.(GASHeat-c,OIL-c)
OIL-a.(GASHeat-c,OIL-c)
COAL-a.(GASHeat-c,OIL-c)
ENIT-a.(GASHeat-c,OIL-c)
NEINT-a.(GASHeat-c,OIL-c)
AGRI-a.(GASHeat-c,OIL-c)
/

*Greenhouse gas commodity
GC(AC) /CO2-c /
*Commodity-Activity mapping
XPXC(C,A)
/ELEC-c.ELEC-a,
GASHeat-c.GASHeat-a,
OIL-c.OIL-a,
COAL-c.COAL-a,
ENIT-c.ENIT-a,
NEINT-c.NEINT-a,
AGRI-c.AGRI-a/
* intermediate demand of C in A production
XAPA(C,A)
/
Elec-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
ENIT-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
NEINT-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
AGRI-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
/
* fuel demand of C in A produciton
XEPA(C,A)
/
GASHeat-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
Oil-c.(ELEC-a,GASHeat-a,OIL-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
COAL-c.(ELEC-a,GASHeat-a,COAL-a,ENIT-a,NEINT-a,AGRI-a)
/

FD(ACT) /Gov,S-I/
*Final demand mix for non household institutions
FD_C(C,FD)
/
GASHeat-c.S-I
OIL-c.S-I
COAL-c.S-I
ENIT-c.S-I
NEINT-c.(Gov,S-I)
AGRI-c.S-I
/
*Factor
F(AC) /Labor,Capital /
L(F) /Labor /
K(F) /Capital /
*Institutes
INS(AC) /Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I,ROW /
*Domestic institutes
INSD(INS) /Household,GoV,NRES,Ptaxin,Ptaxetc,Tarrif,YTAX,S-I /
*Foreign institute
INSDN(INS) /ROW /
*Household
H(INSD) /Household /
* comsumption mix for each household
XACH(H,C)
/Household.ELEC-c
Household.GASHeat-c
Household.OIL-c
Household.COAL-c
Household.ENIT-c
Household.NEINT-c
Household.AGRI-c /
*t/0*1/
t/0*21/
Alias(ACT,ACTPP);
Alias(INS,INSP);
Alias(AC,ACP);
Alias(A,AP);
Alias(C,CP);
Alias(F,FP);
Alias(L,LP);
Alias(K,KP);
Alias(H,HP);
Alias(GC,GCP);
Alias(ENC,ENCPP);
ACNT(ACT)=yes;
ACNT('Total')=no;


PARAMETERS
alpha_nres(A) net residue to output ratio
ta_in(A)        net producer's tax rate in production a (PTAXin)
ta_ex(A)        etc producer's tax rate in production a (PTAXex)
ica(C,A)      Material(C) intermediate demand coefficient in production A:  XAP(C_A) over XC(A)

thetaP(GC,A)  GHG emission per unit of Activity a output

alphaq(C)       Armington CES function shifting parameter
deltaq(C)       Armington CES function share parameter
rhoq(C)                Armington CES function exponent
sigmaq(C)       1 over (1+rhoq(CM))

alphat(C)       CET function shifting parameter
deltat(C)       CET function share parameter
rhot(C)                CET function exponent
sigmat(C)       1 over (1-rhot(CE))

alphaaVAE(A)    VAE CES function shifting parameter
deltaXEP(A)     VAE CES function share parameter for Energy composite XEP
deltaVA(A)      VAE CES function share parameter for Value Added composite VA
deltaf(F,A)     VA CES function share parameter for factor
deltac(C,A)     Leontief perameter for non energy intermediate input and CES share parameter for intermediate energy input (ENC and ENCN)
rhoaVAE(A)      VAE CES function exponent
sigmaaVAE(A)    1 over (1+rhoaVAE(A))
rhoaXFL(A)      XFL CES function exponent (fuel)
sigmaaXFL(A)    1 over (1+rhoaXLF(A))

thetaE(GC,C,A)   GHG emmision per unit of Energy good C input in A: valid if C in ENC

tm(C)           Import tax rate
tm_in           Net producers' tax levied on import
te(C)           Export tax rate

pwm(C)          CIF price of good C in foreign currency
pwe(C)          FOB price of good C in foreign currency
fsav0           base year current account balance
FSAD            Foreign saving adjustment factor exogenous


tr0_per(H)           base year transfer payment to Oldpopulation ratio of Household H
shr(F,H)         Household H share of endownment F
thetaRES_c       share of net residue payment to consumption
* ResinC (Residue payment to household) =thetaRES_c*(sum(A, alphaNRES(A)*XC(A)))

mu(C,H)           marginal propensity to consume Commodity C
mus(H)             marginal propensity to save
epsilon_L(L)     uncompensated elasticity of labor to real wage
Lw0(L)           Labor supply adjusting constant Ls=Lw0*(realwage)^epsilon_L(L)
tc_in            Net Producers' tax in Household consumption to Total Household consumption


*Government
qg0(C)           Base year government consumption of Commodity c
qgr0(C)          ratio of Base year gov consumption of commodity c to Absortion
sg0              Base year government savings.
tg_in            Net Producers' tax revenue in government consumption to Gov consumption
*tax policy paremeters
TINS0            Base year income tax rate
YD0(H)          Base Year Household Disposable Income

*carbon tax policy paremeters
gtax(GC)        carbon tax rate
gtax_policy(GC,t) carbon tax rate evolution
*dgtax(t)      Optional GHG tax rate evolution when all tax increse would same
CrevH         Revenue Recycling share of Household transfer
CrevH_share(H)   Household share of Revenue Recycling share of Household transfer
CrevC         Revenue Recycling share of Government Consumption
CrevI         Revenue Recycling share of Total Industry subsidy
CrevIw(A)     Revenue Recycling industry share of Industry subsidy
Crevtax       Revenue Recycling share of income tax cut

*initialization parameter replicating bau. Set as zero
gtaxrate_o
crevI_o

cwrt(C)         Price index weight
theta(A,C)      Yield of commodity C form one unit of Activity A
qinv_o(C)        Initial investment level
tiv_in          Net Producers' tax revenue in government consumption to Investment demand
thetaRes_iv     share of net residue payment to investment

cpi_o             Consumer price index
lambdat         overall labor productivity
lambda(A)       Activity specific labor productivity
AEEI(C,A)       Automatic energy efficiency imporvment index.
lambdak         overal capital productivity
lambdaka(A)     Activity specific capital productivity
Oldpop          Population age 65 and older at base year
Oldpopg(t)      Population age 65 and older growth rate
TBg(t)          Foreing saving growth rate
Scenario        Policy scenarios BAU TR GE GS LCUT NCUT
delta           Depreciation rate
SFuelmix(A) Size of Solid fuel input set
/ELEC-a   0
GASHeat-a 0
OIL-a     0
COAL-a    0
ENIT-a    0
NEINT-a   0
AGRI-a    0/
LFuelmix(A) Size of Liquid fuel input set and Solid fuel input set
/ELEC-a   1
GASHeat-a 1
OIL-a     1
COAL-a    1
ENIT-a    1
NEINT-a   1
AGRI-a    1/;
*VARIABLES========================================================================
positive variables
PC(A)           Activity price
XC(A)           Activity supply
XMT(C)           Import composite of good C
PMT(C)           Import composite price of good C
XD(C)           Domestic supply of domestic production C
PD(C)           Price of domestic supply of domestic production C

XA(C)           Supply of C
PA(C)           Market price of C
ES(C)           Export of good C
PET(C)           Export price of C
XP(C)           Export domestic supply composite of domestic production
PP(C)           Price of export domestic supply composite of domestic production
QNEG(C,A)  Demand of commodity c GHG composite for non electricity composite
PNEG(C,A)       Price of commodity c GHG composite


Variables

PVAE(A)                Value added energy composite price
QVAE(A)                Demand for Value added energy composite
PVA(A)                 Value Added composite price
VA(A)                  Demand for Value Added composite
PEP (A)                Energy composite price
XEP (A)                Energy Composite demand
PFL (A)                Fuel composite price
XFL (A)                Fuel composite demand
PLFL (A)               Liquid Fuel composite price
XLFL (A)               Liquid Fuel composite demand

XAP(C,A)       Intermediate demand of commodity C
QINTG(GC,A)   Intermediate demand of GHG G

R(K)           rental price of Capital
W(L)           wage
LD(L,A)            labor demand for production A
KD(K,A)            capial demand for production A
Ks(K)           Capital supply
Ls(L)           Labor supply

QCE(C,A)      Demand of commodity c for Electricity composite in A production
QGE(GC,C,A)      Demand of GHG g for commodity c-GHG composite in A production

EXR             The ratio of domestic currency to foeign currency
YG              Government revenue

TR(H)           Transfer payment to Household H
XAF(FD,C)       Final Demand of good C government and Investment
SG              Government savings



MPS(H)          Household Marginal Propensity to save
LY(H)           Household labor income
KY(H)           Household capital income
YH(H)           Household income
YD(H)            Household disposable income
XAC(C,H)                Consumption C of household H
SH(H)            Household Savings

TINSR            income tax rate


FSAV            Net import or Foreign Saving
IVAD            Investment adjustment factor
Warlas          S-I eq dummy


CREV(A)       Carbon tax revenue collected from Activity A
TCREV         Total Carbon tax revenue
ResinC        Residue income from consumption
ResinI        Residue income from investment
CPI           Consumer price index
*** Loading Data

table sam(ACT,ACTP) data in CSV format
$Ondelim
$include b_sam_br_g.csv
$Offdelim

table samng(ACT,ACTP) data in CSV format
$Ondelim
$include b_sam_br_ng.csv
$Offdelim

table ghg(ACT,ACTP) data in CSV format
$Ondelim
$include GHG_BR_p.csv
$Offdelim


Equations

*Price Block
ImPr(C) Import Price
ExPr(C) Export price
AspPr(C) Absorption Price PA f of PD and PMT
AspPrni(C) Absorption price PA without import PA=PD
AspPrnd(C) Absorption price PA without domestic production PA=PMT
ProdPr(C) Supply Price PP f of PET and PD
ProdPrne(C) Supply Price PP without export PP=PD
ProdPrnd(C) Supply Price PP without domestic supply PP=PET


ComPr(C) Activity Price and Commodity Price
ActR(A) Activity revenue and costs with C02 not in process
ActRp(A) Activity Revenue and costs with C02 in process

VAEPr0(A) KLEM composite price with VA only (without EXP)
VAEPr1(A) KLEM composite price with VA and XEP
VAEPr2(A) KLEM composite price with XEP only (without VA)

VAPr0(A)  VA composite price with labor only
VAPr1(A)  VA composite price
VAPr2(A)  VA composite price with capital only
XEPr0 (A) Energy composite price (XEP without fuel)
XEPr1 (A) Energy composite price (XEP with fuel)
XEPr2 (A) Energy composite price (XEP without electricity)
XFLPr0(A) Fuel composite price type 0 (XFL without solid fuel)
XFLPr1(A) Fuel composite price type 1 (XFL with solid fuel)
XFLPr2(A) Fuel composite price type 2 (XFL solid fuel only with multiple solid fuel)
XFLPr3(A) Feul composite price type 3 (XFL solid fuel only with single solid fuel)
XLFLPr0 (A) Liquid Fuel composite price (single Liquid fuel)
XLFLPr1 (A) Liquid Fuel composite price (multiple Liquid fuel)
NEGPr(C,A) Non Electricity GHG composite price with CO2 PNEG f of PA of fuel and PG

Norm Definition of Consummer Price Index
CPIfix Fixing CPI
Labsup(L) Labor supply

*Production and Trade block
QVAED(A) value added - Energy composite demand QVAE f of XC
XVAD0(A) value addede demand VA f of QVAE without XEP
XVAD1(A) valude added demand VA f of QVAE
XEPD1(A)  Energy composite demand XEP function of QVAE
XEPD2(A)  Energy composite demand XEP function of QVAE without VA
XFLD1(A)  Fuel composite demand XFL function of XEP (XEP with electricity and fuel)
XFLD2(A)  Fuel composite demand XFL function of XEP (XEP without electricity)
XLFLD0(A) liquid fuel composite demand XFLD function of XFL(Single Liquid fuel)
XLFLD1(A) liquid fuel composite demand XFLD function of XFL(multiple Liqudi fuel)

INTDM(C,A) Material intermeidate demand for XC(A): Leontief of XC
INTDE1(C,A) Electricity intermediate demand for XC(A): CES of VAE
INTGD(GC,A) intermediate demand of CO2 CO2 emission in process

LDA1(L,A) labor demand LD f of QVAE
KDA1(K,A) Capital demand  KD f of QVAE

NEGDL0(C,A) Fuel-CO2 composite Demand: Liquid Singleton. f of XLFL
NEGDL1(C,A) Fuel-CO2 composite Demand: Liquid. Multi input . f of XLFL

NEGDS1(C,A) Fuel-CO2 composite Demand: Solid . f of XFL
NEGDS2(C,A) Fuel-CO2 composite Demand: Solid . f of XFL (XFL without Liquid fuel demand. Multiple solid fuel)
NEGDS3(C,A) Fuel-CO2 composite Demand: Solid . f of XFL (XFL without Liquid fuel demand. single solid fuel)

NELQCED(ENC,A) fuel demand QCE f of QNEG

GD(GC,ENC,A) CO2 Demand(emission) in fuel use QGE f of QCE


ActDC(A) Activity demand by commodity

XDD(C) Domestic Product Demand f of XA
XDDni(C) Domestic Product Demand without import
XMTD(C) Import Demand  f of XA
XMTDnd(C) Import Demand without domestic product

XDS(C) Domestic Produciton Supply
XDSne(C) Domestic Production Supply without export
ESS(C) Export Supply
ESSnd(C) Export without Domestic Production.

*Institution
InvD(C) Investment Demand
HouseLY(H) Household labor Income
HouseKY(H) Household capital income
HouseY(H) Household Income
HouseYD(H) Household Disposable Income

HouseD(C,H) Household Demand
Hsave(H) Household Savings
Saver(H) Marginal propensity to save


GovE(C) Government Spending
Tras(H) transfer payment
Ytax Income tax rate

ForS Foreign Saving
GovI Government Income

*Market clearing
LabM(L) Labor market clearing
CapM(K) Capital market clearing


ComMENCN(C) Market clearing(ENCN)
ComMENC(C)  Market clearing (ENC)


CREVE(A)  carbon tax revenue without process emission (A)
CREVP(A)  carbon tax revenue with process emission(A)
 TCREVsum   Total carbon tax revenue

InvM Savings and Investment clearing
CAB Current Accoutn Balance
GovB Governmetn Balance

*Residue income for household and investment
ResC Residue income due to consumption
ResI Residuw income due to investment;

*Price Block
ImPr(C)$(SAM('ROW',C) ne 0)..pwm(C)*(1+tm(C))*EXR=g=PMT(C);
ImPr.m(C)=1;
ExPr(C)$(SAM(C,'ROW') ne 0)..PET(C)=g=pwe(C)*(1-te(C))*EXR;
ExPr.m(C)=1;
AspPr(C)$((SAM('ROW',C) ne 0) and (sum(A,SAM(A,C))>0))..(1/alphaq(C))*(((deltaq(C))**(sigmaq(C))*(PD(C))**(1-sigmaq(c))+ (1-deltaq(C))**sigmaq(C)*(PMT(C))**(1-sigmaq(c)))**(1/(1-sigmaq(C))))=g=PA(C);
AspPr.m(C)=1;
AspPrni(C)$(SAM('ROW',C)=0 and (sum(A,SAM(A,C))>0))..PA(C)=e=PD(C);
AspPrni.m(C)=1;
AspPrnd(C)$(sum(A,SAM(C,A))+sum(H,SAM(C,H))+sum(FD,SAM(C,FD))-SAM('ROW',C)=0)..PA(C)=e=PMT(C);
AspPrnd.m(C)=1;


ProdPr(C)$((SAM(C,'ROW') ne 0) and (sum(A,SAM(A,C))>0))..PP(C)=g=(1/alphat(C))*(((deltat(C))**(sigmat(C))*(PET(C))**(1-sigmat(c))+ (1-deltat(C))**(sigmat(C))*(PD(C))**(1-sigmat(c)))**(1/(1-sigmat(C))));
ProdPr.m(C)=1;
ProdPrne(C)$(SAM(C,'ROW')=0 and (sum(A,SAM(A,C))>0))..PP(C)=e=PD(C);
ProdPrne.m(C)=1;
ProdPrnd(C)$(sum(A,SAM(A,C))-SAM(C,'ROW')=0)..PP(C)=e=PET(C);
ProdPrnd.m(C)=1;


ComPr(C)$(sum(A,SAM(A,C))>0)..sum(A$XPXC(C,A),theta(A,C)*PC(A))=g=PP(C);
ComPr.m(C)=1;

*ActR(A)$(sum(C,SAM(A,C))>0 and not ghg('process',A))..PC(A)*alpha_nres(A)*XC(A)+sum(C$M(C),PA(C)*ica(C,A)*XC(A))+PVAE(A)*QVAE(A)=g=PC(A)*(1-ta_in(A)-ta_ex(A))*XC(A)+crevI*crevIw(A)*CREV(A);
*ActR(A)$(sum(C,SAM(A,C))>0 and not ghg('process',A))..PC(A)*alpha_nres(A)*XC(A)+sum(C$M(C),PA(C)*ica(C,A)*XC(A))+PVAE(A)*QVAE(A)=g=PC(A)*(1-ta_in(A)-ta_ex(A))*XC(A);
ActR(A)$(sum(C,SAM(A,C))>0 and not ghg('process',A))..PC(A)*alpha_nres(A)+sum(C$M(C),PA(C)*ica(C,A))+PVAE(A)=g=PC(A)*(1-ta_in(A)-ta_ex(A));
ActR.m(A)=1;
*ActRp(A)$(ghg('process',A) ne 0)..PC(A)*alpha_nres(A)*XC(A)+sum(C$M(C),PA(C)*ica(C,A)*XC(A))+sum(GC,thetaP(GC,A)*XC(A)*gtax(GC)*cpi)+PVAE(A)*QVAE(A)=g=PC(A)*(1-ta_in(A)-ta_ex(A))*XC(A)+crevI*crevIw(A)*CREV(A);
ActRp(A)$(ghg('process',A) ne 0)..PC(A)*alpha_nres(A)+sum(C$M(C),PA(C)*ica(C,A))+sum(GC,thetaP(GC,A)*gtax(GC)*cpi)+PVAE(A)=g=PC(A)*(1-ta_in(A)-ta_ex(A));
ActRP.m(A)=1;


VAEPr0(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))=0 and sum(F,SAM(F,A))>0)..PVAE(A)=e=(1/alphaaVAE(A))*PVA(A);
VAEPr0.m(A)=1;


VAEPr1(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))>0 and sum(F,SAM(F,A))>0)..PVAE(A)=e=(1/alphaaVAE(A))*
(
deltaXEP(A)**sigmaaVAE(A)*PEP(A)**(1-sigmaaVAE(A))
+ deltaVA(A)**sigmaaVAE(A)*PVA(A)**(1-sigmaaVAE(A))
)**(1/(1-sigmaaVAE(A)));
VAEPr1.m(A)=1;

VAEPr2(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))>0 and sum(F,SAM(F,A))=0)..PVAE(A)=e=(1/alphaaVAE(A))*PEP(A);
VAEPr2.m(A)=1;


VAPr0(A)$(sum(L,SAM(L,A))>0 and sum(K,SAM(K,A))=0)..PVA(A)=e=
(
     prod(L,(W(L)/(lambdat*lambda(A)*deltaf(L,A)))**(deltaf(L,A)))
);
VAPr0.m(A)=1;


VAPr1(A)$(sum(L,SAM(L,A))>0 and sum(K,SAM(K,A))>0)..PVA(A)=e=
(
     prod(K,(R(K)/(lambdak*lambdaka(A)*deltaf(K,A)))**(deltaf(K,A)))*
     prod(L,(W(L)/(lambdat*lambda(A)*deltaf(L,A)))**(deltaf(L,A)))
);
VAPr1.m(A)=1;

VAPr2(A)$(sum(L,SAM(L,A))=0 and sum(K,SAM(K,A))>0)..PVA(A)=e=
(
     prod(K,(R(K)/(lambdak*lambdaka(A)*deltaf(K,A)))**(deltaf(K,A)))
);
VAPr2.m(A)=1;


XEPr0(A)$(sum(ELECC,SAM(ELECC,A))>0 and sum(ENC,SAM(ENC,A))=0 )..PEP(A)=e=prod(C$ELECC(C),(PA(C)/deltaC(C,A))**(deltaC(C,A)));
XEPr0.m(A)=1;



XEPr1(A)$(sum(ELECC,SAM(ELECC,A))>0 and sum(ENC,SAM(ENC,A))>0)..PEP(A)=e=
  prod(C$ELECC(C),(PA(C)/deltaC(C,A))**(deltaC(C,A)))*
(
PFL(A)/
(1-sum(C$ELECC(C),deltaC(C,A)))
)**
(
1-sum(C$ELECC(C),deltaC(C,A))
)
;
XEPr1.m(A)=1;

XEPr2(A)$(sum(ELECC,SAM(ELECC,A))=0 and sum(ENC,SAM(ENC,A))>0)..PEP(A)=e=PFL(A);
XEPr2.m(A)=1;


XFLPr0(A)$(sum(C$sfuel(C),SAM(C,A)=0) and sum(C$Lfuel(C), SAM(C,A))>0)..PFL(A)=e=PLFL(A);
XFLPr0.m(A)=1;

XFLPr1(A)$(sum(C$sfuel(C),SAM(C,A) >0) and sum(C$Lfuel(C), SAM(C,A))>0)..PFL(A)=e=
(
sum(C$SfuelA(A,C),deltaC(C,A)**sigmaaXFL(A)*(PNEG(C,A)/AEEI(C,A))**(1-sigmaaXFL(A)))
+
(1-sum(C$SfuelA(A,C),deltaC(C,A)))**(sigmaaXFL(A))*PLFL(A)**(1-sigmaaXFL(A))
)**(1/(1-sigmaaXFL(A)));
XFLPr1.m(A)=1;

XFLPr2(A)$(sum(C$sfuel(C),SAM(C,A) >0) and sum(C$Lfuel(C), SAM(C,A))=0 and Sfuelmix(A) eq 1)..PFL(A)=e=
(
sum(C$SfuelA(A,C),deltaC(C,A)**sigmaaXFL(A)*(PNEG(C,A)/AEEI(C,A))**(1-sigmaaXFL(A)))
)**(1/(1-sigmaaXFL(A)));
XFLPr2.m(A)=1;

XFLPr3(A)$(sum(C$sfuel(C),SAM(C,A) >0) and sum(C$Lfuel(C), SAM(C,A))=0 and Sfuelmix(A) eq 0)..PFL(A)=e=sum(C$(Sfuel(C) and SfuelA(A,C)),PNEG(C,A))/sum(C$(Sfuel(C) and SfuelA(A,C)),AEEI(C,A));
XFLPr3.m(A)=1;


XLFLPr0(A)$(Lfuelmix(A) eq 0)..PLFL(A)=e=sum(C$(ENC(C) and LfuelA(A,C)),PNEG(C,A))/sum(C$(ENC(C) and LfuelA(A,C)),AEEI(C,A));
XLFLPr0.m(A)=1;

XLFLPr1(A)$(Lfuelmix(A) eq 1)..PLFL(A)=e=
prod(C$LfuelA(A,C),(PNEG(C,A)/(deltaC(C,A)*AEEI(C,A)))**deltaC(C,A))
;
XLFLPr1.m(A)=1;



NEGPr(C,A)$(ghg(C,A)>0)..PA(C)+sum(GC,thetaE(GC,C,A)*gtax(GC)*cpi)=e=PNEG(C,A);
NEGPr.m(C,A)=1;

*Numeria
Norm..cpi=e=sum(C,PA(C)*cwrt(C));
Norm.m=1;
CPIfix..cpi=e=cpi_o;
CPIfix.m=1;

*Production and Trade block
QVAED(A)$(sum(C,SAM(A,C))>0)..QVAE(A)=e=XC(A);
QVAED.m(A)=1;

XVAD0(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))=0 and sum(F,SAM(F,A))>0)..VA(A)=e=(1/alphaaVAE(A))*QVAE(A);
XVAD0.m(A)=1;

XVAD1(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))>0 and sum(F,SAM(F,A))>0)..VA(A)=e=(deltaVA(A)**sigmaaVAE(A))*((PVAE(A)/PVA(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE(A);
XVAD1.m(A)=1;

XEPD1(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))>0 and sum(F,SAM(F,A))>0)..XEP(A)=e=(deltaXEP(A)**sigmaaVAE(A))*((PVAE(A)/PEP(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE(A);
XEPD1.m(A)=1;

XEPD2(A)$(sum(ELECC,SAM(ELECC,A))+sum(ENC,SAM(ENC,A))>0 and sum(F,SAM(F,A))=0)..XEP(A)=e=(1/alphaaVAE(A))*QVAE(A);
XEPD2.m(A)=1;



LDA1(L,A)$(sum(LP,SAM(LP,A))>0)..LD(L,A)=e=deltaf(L,A)*(PVA(A)*VA(A))/W(L);
LDA1.m(L,A)=1;

KDA1(K,A)$(sum(KP,SAM(KP,A))>0)..KD(K,A)=e=deltaf(K,A)*(PVA(A)*VA(A))/R(K);
KDA1.m(K,A)=1;

XFLD1(A)$(sum(ELECC,SAM(ELECC,A))>0 and sum(ENC,SAM(ENC,A))>0)..XFL(A)=e=(1-sum(C$ELECC(C),deltaC(C,A)))*(PEP(A)/PFL(A))*XEP(A);
XFLD1.m(A)=1;
XFLD2(A)$(sum(ELECC,SAM(ELECC,A))=0 and sum(ENC,SAM(ENC,A))>0)..XFL(A)=e=XEP(A);
XFLD2.m(A)=1;


INTDE1(C,A)$(ELECC(C) and SAM(C,A)$ELECC(C)>0)..XAP(C,A)=e= deltaC(C,A)*(PEP(A)/PA(C))*XEP(A);
INTDE1.m(C,A)=1;



XLFLD0(A)$(sum(C$sfuel(C),SAM(C,A)=0) and (sum(CP$Lfuel(CP),SAM(CP,A)) >0))..XLFL(A)=e=XFL(A);
XLFLD0.m(A)=1;

NEGDS1(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0) and (sum(CP$Lfuel(CP),SAM(CP,A)) >0))..QNEG(C,A)=e=                (deltaC(C,A)**sigmaaXFL(A))*((PFL(A)/PNEG(C,A))**sigmaaXFL(A))*((AEEI(C,A))**(sigmaaXFL(A)-1))*XFL(A);
NEGDS1.m(C,A)=1;

NEGDS2(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0) and (sum(CP$Lfuel(CP),SAM(CP,A))=0) and Sfuelmix(A) eq 1)..QNEG(C,A)=e=                (deltaC(C,A)**sigmaaXFL(A))*((PFL(A)/PNEG(C,A))**sigmaaXFL(A))*((AEEI(C,A))**(sigmaaXFL(A)-1))*XFL(A);
NEGDS2.m(C,A)=1;

NEGDS3(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0) and (sum(CP$Lfuel(CP),SAM(CP,A))=0) and Sfuelmix(A) eq 0)..QNEG(C,A)=e=(1/AEEI(C,A))*XFL(A);
NEGDS3.m(C,A)=1;


XLFLD1(A)$(sum(C$sfuel(C),SAM(C,A) >0) and (sum(CP$Lfuel(CP),SAM(CP,A)) >0))..XLFL(A)=e=(1-sum(C$sfuel(C),deltaC(C,A)))**sigmaaXFL(A)*(PFL(A)/PLFL(A))**(sigmaaXFL(A))*XFL(A);
XLFLD1.m(A)=1;

NEGDL0(C,A)$(Lfuel(C) and Lfuelmix(A) eq 0)..QNEG(C,A)=e=(1/AEEI(C,A))*XLFL(A);
NEGDL0.m(C,A)=1;

NEGDL1(C,A)$(Lfuel(C) and Lfuelmix(A) eq 1)..QNEG(C,A)=e=deltaC(C,A)*(PLFL(A)/PNEG(C,A))*XLFL(A);
NEGDL1.m(C,A)=1;







INTDM(C,A)$(M(C) and SAM(C,A)$M(C)>0)..XAP(C,A)=e=ica(C,A)*XC(A);
INTDM.m(C,A)=1;
INTGD(GC,A)$(ghg('process',A)>0)..QINTG(GC,A)=e=thetaP(GC,A)*XC(A);
INTGD.m(GC,A)=1;




NELQCED(ENC,A)$(SAM(ENC,A)>0)..QCE(ENC,A)=e=QNEG(ENC,A);
NELQCED.m(ENC,A)=1;

GD(GC,ENC,A)$(ghg(ENC,A)>0)..QGE(GC,ENC,A)=e=thetaE(GC,ENC,A)*QNEG(ENC,A);
GD.m(GC,ENC,A)=1;

ActDC(A)$(sum(C,SAM(A,C))>0)..XC(A)=g=sum(C$XPXC(C,A),theta(A,C)*XP(C));
ActDC.m(A)=1;

XDD(C)$((SAM('ROW',C) ne 0) and (sum(A,SAM(A,C))>0))..XD(C)=g=(deltaq(C)**sigmaq(C))*((PA(C)/(PD(C)))**sigmaq(C))*(alphaq(C)**(sigmaq(C)-1))*XA(C);
XDD.m(C)=1;
XDDni(C)$(SAM('ROW',C)=0 and (sum(A,SAM(A,C))>0)).. XD(C)=e=XA(C);
XDDni.m(C)=1;

XMTDnd(C)$(sum(A,SAM(C,A))+sum(H,SAM(C,H))+sum(FD,SAM(C,FD))-SAM('ROW',C)=0)..XMT(C)=e=XA(C);
XMTDnd.m(C)=1;
XMTD(C)$((SAM('ROW',C) ne 0) and (sum(A,SAM(A,C))>0))..XMT(C)=g=((1-deltaq(C))**sigmaq(C))*((PA(C)/(PMT(C)))**sigmaq(C))*(alphaq(C)**(sigmaq(C)-1))*XA(C);
XMTD.m(C)=1;

ESS(C)$(SAM(C,'ROW') ne 0 and (sum(A,SAM(A,C))>0))..(deltat(C)**sigmat(C))*((PP(C)/PET(C))**sigmat(C))*(alphat(C)**(sigmat(C)-1))*XP(C)=g=ES(C);
ESS.m(C)=1;
ESSnd(C)$(sum(A,SAM(A,C))-SAM(C,'ROW')=0)..ES(C)=e=XP(C);
Essnd.m(C)=1;

XDS(C)$(SAM(C,'ROW') ne 0 and (sum(A,SAM(A,C))>0))..((1-deltat(C))**sigmat(C))*((PP(C)/PD(C))**sigmat(C))*(alphat(C)**(sigmat(C)-1))*XP(C)=g=XD(C);
XDS.m(C)=1;
XDSne(C)$(SAM(C,'ROW')=0 and (sum(A,SAM(A,C))>0)).. XD(C)=e=XP(C);
XDSne.m(C)=1;


Labsup(L)..Ls(L)=e=Lw0(L)*(((1-TINSR)*(W(L)/cpi))**epsilon_L(L));
Labsup.m(L)=1;
*Institution
InvD(C)$(SAM(C,'S-I') ne 0)..XAF('S-I',C)=e=qinv_o(C)*IVAD;
InvD.m(C)=1;
HouseLY(H)..LY(H)=e=sum(L,shr(L,H)*Ls(L)*W(L));
HouseLY.m(H)=1;
HouseKY(H)..KY(H)=e=sum(K,shr(K,H)*Ks(K)*R(K));
HouseKY.m(H)=1;
HouseYD(H)..YD(H)=e=(1-TINSR)*(LY(H)+KY(H)+TR(H)+ResinC-delta*sum(K,shr(K,H)*Ks(K)));
HouseYD.m(H)=1;
HouseY(H)..YH(H)=e=LY(H)+KY(H)+TR(H)+ResinC;
HouseY.m(H)=1;
Saver(H)..MPS(H)=e=mus(H);
Saver.m(H)=1;
HouseD(C,H)$(SAM(C,'Household') ne 0)..XAC(C,H)=e=mu(C,H)*(YD(H)/(PA(C)*(1+tc_in)));
HouseD.m(C,H)=1;
Hsave(H)..SH(H)=e=MPS(H)*YD(H)+delta*sum(K,shr(K,H)*Ks(K));
Hsave.m(H)=1;

GovE(C)$(SAM(C,'Gov') ne 0)..XAF('Gov',C)=e=qgr0(C)*sum(CP,PA(CP)*XA(CP))*(1/(PA(C)*(1+tg_in)))+(qg0(C)/sum(CP,qg0(C)))*crevc*TCREV/(PA(C)*(1+tg_in));
GovE.m(C)=1;
Tras(H)..Tr(H)=e=tr0_per(H)*cpi*Oldpop+crevh*crevh_share(H)*TCREV;
Tras.m(H)=1;

Ytax..TINSR=e=(sum(H,(YH(H)-delta*sum(K,shr(K,H)*Ks(K))))*TINS0-crevtax*TCREV)/sum(H,(YH(H)-delta*sum(K,shr(K,H)*Ks(K))));
Ytax.m=1;

ForS..FSAV=e=fsav0*FSAD;
ForS.m=1;
GovI..YG=e=
sum(A$(sum(C,SAM(A,C))>0),(ta_in(A)+ta_ex(A))*PC(A)*XC(A))
+sum(C$(SAM('ROW',C) ne 0),tm(C)*pwm(C)*XMT(C)*EXR)
+sum(C$(SAM(C,'ROW') ne 0),te(C)*pwe(C)*ES(C)*EXR)
+(TINSR)*sum(H,(YH(H)-delta*sum(K,shr(K,H)*Ks(K))))
+sum((GC,A)$(ghg('process',A)>0),gtax(GC)*cpi*QINTG(GC,A))
+sum((GC,C,A)$(ghg(C,A)>0),gtax(GC)*cpi*QGE(GC,C,A))
+(tc_in)*sum((C,H)$(SAM(C,'Household') ne 0),PA(C)*XAC(C,H))
+tg_in*sum(C$(SAM(C,'Gov') ne 0),PA(C)*XAF('Gov',C))
+tiv_in*sum(C$(SAM(C,'S-I') ne 0),PA(C)*XAF('S-I',C))
+tm_in*sum(C$(SAM('ROW',C) ne 0),PMT(C)*XMT(C));
GovI.m=1;
******Market clearing
LabM(L)..Ls(L)=g=sum(A$(SAM(L,A)>0),LD(L,A));
LabM.m(L)=1;
CapM(K)..Ks(K)=g=sum(A$(SAM(K,A)>0),KD(K,A));
CapM.m(K)=1;


ComMENCN(C)$(ENCN(C))..XA(C)=g=sum(A$XAPA(C,A),XAP(C,A))      +sum(H$XACH(H,C),XAC(C,H))+sum(FD$FD_C(C,FD),XAF(FD,C));
ComMENCN.m(C)=1;
ComMENC(C)$(ENC(C))..XA(C)=g=sum(A$XEPA(C,A) ,QCE(C,A ))+sum(H$XACH(H,C),XAC(C,H))+sum(FD$FD_C(C,FD),XAF(FD,C));
ComMENC.m(C)=1;





CREVE(A)$(ghg('process',A) eq 0)..CREV(A)=e=sum(GC,gtax(GC)*cpi*sum(C$fuelA(A,C),QGE(GC,C,A)));
CREVE.m(A)=1;
CREVP(A)$(ghg('process',A) ne 0)..CREV(A)=e=sum(GC,gtax(GC)*cpi*sum(C$fuelA(A,C),QGE(GC,C,A)))+sum(GC,gtax(GC)*cpi*QINTG(GC,A));
CREVP.m(A)=1;
TCREVsum..TCREV=e=sum(A,CREV(A));
TCREVsum.m=1;

InvM..sum(C$(SAM(C,'S-I') ne 0),PA(C)*(1+tiv_in)*XAF('S-I',C))=e=Warlas+sum(H,SH(H))+SG+FSAV*EXR+ResinI;
InvM.m=1;
CAB..sum(C$(SAM('ROW',C) ne 0),pwm(C)*XMT(C))=e=sum(C$(SAM(C,'ROW') ne 0),pwe(C)*ES(C))+FSAV+(tm_in)*(sum(C$(SAM('ROW',C) ne 0),PMT(C)*XMT(C))/EXR);
CAB.m=1;
GovB..SG=e=YG-sum(C$(SAM(C,'Gov') ne 0),PA(C)*(1+tg_in)*XAF('Gov',C))-sum(H,TR(H))-crevI*TCREV;
GovB.m=1;

ResC..ResinC=e=thetaRes_c*sum(A$(sum(C,SAM(A,C))>0),PC(A)*alpha_nres(A)*XC(A));
ResC.m=1;
ResI..ResinI=e=thetaRes_iv*sum(A$(sum(C,SAM(A,C))>0),PC(A)*alpha_nres(A)*XC(A));
ResI.m=1;

**** Initialization and calibration

**** Declare initial values
Parameters
PC0(A)
XC0(A)
XMT0(C)
PMT0(C)
XD0(C)
PD0(C)

XA0(C)
PA0(C)
ES0(C)
PET0(C)
XP0(C)
PP0(C)
QNEG0(C,A)
PNEG0(C,A)

PVAE0(A)
QVAE0(A)
PVA0(A)
VA0(A)
*HKTE0(A)
*PHKTE0(A)
XEP0(A)
PEP0(A)
XFL0(A)
PFL0(A)
XLFL0(A)
PLFL0(A)
XAP0(C,A)
QINTG0(GC,A)

R0(K)
W0(L)
LD0(L,A)
KD0(K,A)
Ks0(K)
Ls0(L)

QCE0(C,A)
QGE0(GC,C,A)

EXR0
YG0

TR0(H)
XAF0(FD,C)
SG0

MPS0(H)
LY0(H)
KY0(H)
YH0(H)
YD0(H)
XAC0(C,H)
SH0(H)

TINSR0


FSAV0
IVAD0
Warlas0


CREV0(A)
TCREV0
ResinC0
ResinI0;

**initial point creating paremeters set as 0


CrevI_o=0;
gtaxrate_o=0.0;

CREV0(A)=gtaxrate_o*sam('CO2-c',A)+gtaxrate_o*ghg('process',A);
TCREV0=sum(A,CREV0(A));
*display CREV0;

*** Policy Scenario ***
* Scenario=0 : BAU
* Scenario=1 : TR
* Scenario=2 : GE
* Scenario=3 : GS
* Scenario=4 : LCUT
* Scenario=5 : NCUT
** Dynamic adjustments**

**carbon tax policy parameters
gtax(GC)=0;
CrevH_share(H)=sam(H,'Gov')/sum(HP,sam(HP,'Gov'));
parameter
gtaxp_BAU(GC,t)
*/CO2-c.0 0/
*/CO2-c.0*1 0/
/CO2-c.0*21 0/
gtaxp_ctax(GC,t)
*/CO2-c.0 0/
*/CO2-c.0 0.18
* CO2-c.1 0.18/
/CO2-c.0*10 0.18
 CO2-c.(11,21) 0
 CO2-c.(12,16,20) 0.1
 CO2-c.13*15 0
 CO2-c.17*19 0/
;

*** policy simulation setting****
Scenario=0;
if (scenario=0,
gtax_policy(GC,t)=gtaxp_BAU(GC,t);
CrevH=0;
CrevC=0;
CrevI=0;
CrevIw(A)=0;
Crevtax=0;
elseif (scenario=1),
gtax_policy(GC,t)=gtaxp_ctax(GC,t);
CrevH=1;
CrevC=0;
CrevI=0;
CrevIw(A)=0;
Crevtax=0;

elseif (scenario=2),
gtax_policy(GC,t)=gtaxp_ctax(GC,t);
CrevH=0;
CrevC=1;
CrevI=0;
CrevIw(A)=0;
Crevtax=0;

elseif (scenario=3),
gtax_policy(GC,t)=gtaxp_ctax(GC,t);
CrevH=0;
CrevC=0;
CrevI=0;
CrevIw(A)=0;
Crevtax=0;

else
gtax_policy(GC,t)=gtaxp_ctax(GC,t);
CrevH=0;
CrevC=0;
CrevI=0;
CrevIw(A)=0;
Crevtax=1;
);

*display
*CrevI_o,gtaxrate_o,CREV0,TCREV0,gtax,CrevH_share,gtaxp_BAU,gtaxp_ctax,scenario,gtax_policy,CrevH,CrevC,CrevI,CrevIW,Crevtax;

sigmaaVAE(A)=0.5;
sigmaaXFL(A)=0.5;
sigmat(C)=-3;
sigmaq('ELEC-c')=0.3;
sigmaq('GASHeat-c')=0.3;
sigmaq('OIL-c')=2.1;
sigmaq('COAL-c')=3.05;
sigmaq('ENIT-c')=3;
sigmaq('NEINT-c')=1.9;
sigmaq('AGRI-c')=3.29;

rhoaVAE(A)=(1/sigmaaVAE(A))-1;
rhoaXFL(A)=(1/sigmaaXFL(A))-1;
rhot(C)$((SAM(C,'ROW') ne 0) and (sum(A,SAM(A,C))>0))=1-(1/sigmat(C));
rhoq(C)$((SAM('ROW',C) ne 0) and (sum(A,SAM(A,C))>0))=(1/sigmaq(C))-1;

*display rhot,rhoq;

*labor productivity
lambdat=1;
lambda(A)=1;
epsilon_L(L)=0.4;

*capital productivity
lambdak=1;
lambdaka(A)=1;

*AEEI
AEEI(C,A)=1;

*base year capital rent : IO capital payment + IO depreciation/ KPI DB 2010 Total capital stock evaluated in 2010 ppi (Unit: 1,000,000,000 won)

R0(K)=sum(A,SAM(K,A))/3403090.255;
KD0(K,A)=SAM(K,A)/R0(K);
*base year wage : IO payroll / 2010 employment in Thousand
W0(L)=sum(A,SAM(L,A))/23890;
LD0(L,A)=SAM(L,A)/W0(L);
*depreciation
** delta=0.046;
delta=R0('Capital')*0.390061;

**Price
PET0(C)$(SAM(C,'ROW') ne 0)=1;
PP0(C)$(sum(A,SAM(A,C))>0)=1;
PC0(A)=1;
PMT0(C)$(SAM('ROW',C) ne 0)=1;
PD0(C)$(sum(A,SAM(A,C))>0)=1;
PA0(C)=1;
EXR0=1;

XC0(A)$(sum(ACTpp,SAMng(A,ACTpp))>0)=sum(ACT,SAMng(ACT,A))/PC0(A);
QCE0(ENC,A)$(sum(ACTpp,SAMng(A,ACTpp))>0)=SAM(ENC,A)/PA0(ENC);
XAP0(C,A)$(not ENC(C))=SAM(C,A)/PA0(C);

ica(C,A)$(M(C))=XAP0(C,A)/XC0(A);


QVAE0(A)$(sum(ACT,SAM(A,ACT))>0)=XC0(A);
PVAE0(A)$(sum(ACT,SAM(A,ACT))>0)=(sum(F,SAMng(F,A))+sum(C$ELECC(C),SAMng(C,A))+sum(C$ENC(C),SAMng(C,A))+sum(ENC,gtaxrate_o*ghg(ENC,A)))/QVAE0(A);



QGE0(GC,C,A)=ghg(C,A);
thetaE(GC,C,A)$(QCE0(C,A)>0)=QGE0(GC,C,A)/QCE0(C,A);

QINTG0(GC,A)$(ghg('process',A)>0)=ghg('process',A);
thetaP(GC,A)$(XC0(A)>0 and ghg('process',A)>0)=QINTG0(GC,A)/XC0(A);

PNEG0(C,A)$(ENC(C) and (sum(ACT,SAM(A,ACT))>0))=PA0(C)+sum(GC,thetaE(GC,C,A)*gtaxrate_o);
QNEG0(C,A)$(ENC(C) and (sum(ACT,SAM(A,ACT))>0))=(QCE0(C,A)*PA0(C)+sum(GC,gtaxrate_o*QGE0(GC,C,A)))/PNEG0(C,A);

deltaF(L,A)=SAM(L,A)/sum(F,SAM(F,A));
deltaF(K,A)=SAM(K,A)/sum(F,SAM(F,A));


VA0(A)$(sum(ACT,SAM(ACT,A)>0))=prod(L,(lambdat*lambda(A)*LD0(L,A))**deltaF(L,A))*prod(K,(lambdak*lambdaka(A)*KD0(K,A))**deltaF(K,A));
PVA0(A)$(sum(ACT,SAM(ACT,A)>0))=sum(F,SAM(F,A))/VA0(A);


parameter
LD1(L,A)
KD1(K,A)
PVA1(A)
;

PVA1(A)
=(
     prod(K,(R0(K)/(lambdak*lambdaka(A)*deltaf(K,A)))**(deltaf(K,A)))*
     prod(L,(W0(L)/(lambdat*lambda(A)*deltaf(L,A)))**(deltaf(L,A)))
);


LD1(L,A)=deltaf(L,A)*(PVA0(A)*VA0(A))/W0(L);

KD1(K,A)=deltaf(K,A)*(PVA0(A)*VA0(A))/R0(K);



deltaC(C,A)$(Lfuel(C) and (SAM(C,A)>0))=SAM(C,A)/sum(CP$Lfuel(CP),sam(CP,A));

parameter

XLFL0(A)
PLFL0(A);

XLFL0(A)=prod(C$Lfuel(C),(AEEI(C,A)*QNEG0(C,A))**deltac(C,A));

PLFL0(A)=(sum(C$Lfuel(C),SAMng(C,A))+sum(C$Lfuel(C),gtaxrate_o*ghg(C,A)))/XLFL0(A);



parameter

XLFL1(A)
PLFL1(A);

PLFL1(A)=prod(C$LfuelA(A,C),(PNEG0(C,A)/(deltaC(C,A)*AEEI(C,A)))**deltaC(C,A))
;



parameter

QNEG1(C,A);

QNEG1(C,A)=deltaC(C,A)*(PLFL0(A)/PNEG0(C,A))*XLFL0(A);




parameter

*KC_XLFL(C,A)
deltaXLFL(A)
KC_XFL(C,A)
KET(C)
KDM(C);

KC_XFL(C,A)$(Sfuel(C) and (sam(C,A)>0))=
((PLFL0(A)/PNEG0(C,A))**(1-sigmaaXFL(A))*((SAM(C,A)+gtaxrate_o*ghg(C,A))/sum(CP$LFuel(CP),SAM(CP,A)+gtaxrate_o*ghg(CP,A)))*(1/AEEI(C,A))**(1-sigmaaXFL(A)))**(1/sigmaaXFL(A));
deltaXLFL(A)=1/(1+sum(C$(sfuel(C)and(sam(C,A)>0)),KC_XFL(C,A)));
deltaC(C,A)$(Sfuel(C) and (sam(C,A)>0))=KC_XFL(C,A)*deltaXLFL(A);


XFL0(A)$(sum (C$sfuel(C), sam(C,A))>0)=
(sum(C$(sfuel(C) and (sam(C,A)>0)),deltaC(C,A)*(QNEG0(C,A)*AEEI(C,A))**(-rhoaXFL(A)))
+(1-sum(C$(sfuel(C) and (sam(C,A)>0)),deltaC(C,A)))*XLFL0(A)**(-rhoaXFL(A)))**(-1/rhoaXFL(A));

XFL0(A)$(sum (C$sfuel(C), sam(C,A)=0))=XLFL0(A);

PFL0(A)$(sum (C$sfuel(C), sam(C,A))>0)=
(sum(C$ENC(C),SAMng(C,A))+sum(ENC,gtaxrate_o*ghg(ENC,A)))/XFL0(A);

PFL0(A)$(sum(C$sfuel(C),SAM(C,A)=0))=PLFL0(A);


parameter
PFL1(A);

PFL1(A)$(sum(C$sfuel(C),SAM(C,A)=0))=PLFL0(A);


PFL1(A)$(sum(C$sfuel(C),SAM(C,A) >0))=
(
sum(C$SfuelA(A,C),deltac(C,A)**sigmaaXFL(A)*(PNEG0(C,A)/AEEI(C,A))**(1-sigmaaXFL(A)))
+
(1-sum(C$SfuelA(A,C),deltac(C,A)))**(sigmaaXFL(A))*PLFL0(A)**(1-sigmaaXFL(A))
)**(1/(1-sigmaaXFL(A)));




XLFL1(A)$(sum(C$sfuel(C),SAM(C,A)=0))=XFL0(A);

QNEG1(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0))=(deltaC(C,A)**sigmaaXFL(A))*((PFL0(A)/PNEG0(C,A))**sigmaaXFL(A))*((AEEI(C,A))**(sigmaaXFL(A)-1))*XFL0(A);
XLFL1(A)$(sum(C$sfuel(C),SAM(C,A) >0))=(1-sum(C$sfuel(C),deltaC(C,A)))**sigmaaXFL(A)*(PFL0(A)/PLFL0(A))**(sigmaaXFL(A))*XFL0(A);



deltaC(C,A)$(ELECC(C) and (SAM(C,A)>0))=SAM(C,A)/(sum(CP$ENC(CP),SAMng(CP,A))+sum(CP$ELECC(CP),SAMng(CP,A))+sum(ENC,gtaxrate_o*ghg(ENC,A)));

XEP0(A)=prod(C$ELECC(C),XAP0(C,A)**deltaC(C,A))*(XFL0(A)**(1-sum(C$ELECC(C),deltaC(C,A))));

PEP0(A)=(sum(C$ELECC(C),SAMng(C,A))+sum(C$ENC(C),SAMng(C,A))+sum(ENC,gtaxrate_o*ghg(ENC,A)))/XEP0(A);


parameter
PEP1(A);

PEP1(A)=
  prod(C$ELECC(C),(PA0(C)/deltaC(C,A))**(deltaC(C,A)))*
(
PFL0(A)/
(1-sum(C$ELECC(C),deltaC(C,A)))
)**
(
1-sum(C$ELECC(C),deltaC(C,A))
)
;



parameter
XFL1(A)
XAP1(C,A);

XFL1(A)$(sum(C,SAM(A,C))>0)=(1-sum(C$ELECC(C),deltaC(C,A)))*(PEP0(A)/PFL0(A))*XEP0(A);
XAP1(C,A)$(ELECC(C) and SAM(C,A)$ELECC(C)>0)= deltaC(C,A)*(PEP0(A)/PA0(C))*XEP0(A);

XAP1(C,A)$(M(C) and SAM(C,A)$M(C)>0)=ica(C,A)*XC0(A);


*Resume here
* Oct. 08. 15:34 (Korean Time)
parameter
KEP(A)
;

KEP(A)=
(
 (PEP0(A)/PVA0(A))**(sigmaaVAE(A)-1)*
 (
  (
   sum(C$ELECC(C),SAMng(C,A))+sum(C$ENC(C),SAMng(C,A))+sum(ENC,gtaxrate_o*ghg(ENC,A))
   )
   /
    sum(F,SAMng(F,A))
  )
)**(1/sigmaaVAE(A));
deltaXEP(A)=KEP(A)/(1+KEP(A));
deltaVA(A)=1/(1+KEP(A));

alphaaVAE(A)$(sum(ACT,SAM(A,ACT))>0)=QVAE0(A)/(deltaXEP(A)*XEP0(A)**(-rhoaVAE(A))+deltaVA(A)*VA0(A)**(-rhoaVAE(A)))**(-1/rhoaVAE(A));

parameter
PVAE1(A);

PVAE1(A)$(sum(C,SAM(A,C))>0)=(1/alphaaVAE(A))*
(
deltaXEP(A)**sigmaaVAE(A)*PEP0(A)**(1-sigmaaVAE(A))
+ deltaVA(A)**sigmaaVAE(A)*PVA0(A)**(1-sigmaaVAE(A))
)**(1/(1-sigmaaVAE(A)));



parameter
XEP1(A)
VA1(A);

VA1(A)$(sum(C,SAM(A,C))>0)=(deltaVA(A)**sigmaaVAE(A))*((PVAE0(A)/PVA0(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE0(A);
XEP1(A)$(sum(C,SAM(A,C))>0)=(deltaXEP(A)**sigmaaVAE(A))*((PVAE0(A)/PEP0(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE0(A);

*display PVAE0,PVAE1;
*display XEP0,XEP1,VA0,VA1;
*display PEP0,PEP1,PVA0,PVA1;
*display XFL0,XFL1,XAP0,XAP1;
*display PFL0,PFL1;
*display XLFL0,XLFL1,QNEG0,QNEG1;
*display PLFL0,PLFL1;
*display KD1,KD0,LD1,LD0;



theta(A,C)$((sum(ACT,SAM(A,ACT))>0) and (sum(AP,SAM(AP,C))>0))=(SAM(A,C)/PP0(C))/XC0(A);


** tax parameters

ta_in(A)$((sum(ACT,SAM(A,ACT))>0))=SAMng('Ptaxin',A)/sum(ACT,samng(ACT,A));
ta_ex(A)$((sum(ACT,SAM(A,ACT))>0))=SAMng('Ptaxetc',A)/sum(ACT,samng(ACT,A));

te(C)$(SAM(C,'ROW') ne 0)=0;
tm(C)$(SAM('ROW',C) ne 0)=SAM('Tarrif',C)/SAM('ROW',C);

** Trade related parameters
*(i) world price
pwm(C)$(SAM('ROW',C) ne 0)=PMT0(C)/((1+tm(C))*EXR0);
pwe(C)$(SAM(C,'ROW') ne 0)=PET0(C)/((1-te(C))*EXR0);

*(ii) CET coefficient**

ES0(C)$(SAM(C,'ROW') ne 0)=SAM(C,'ROW')/PET0(C);
XD0(C)$(sum(A,SAM(A,C))>0)=(sum(A,SAM(A,C))-SAM(C,'ROW'))/PD0(C);

KET(C)$(SAM(C,'ROW') ne 0 and (sum(A,SAM(A,C))>0))=((ES0(C)/XD0(C))*((PET0(C)/PD0(C))**sigmat(C)))**(1/sigmat(C));
deltat(C)$(SAM(C,'ROW') ne 0 and (sum(A,SAM(A,C))>0))=KET(C)/(1+KET(C));

XP0(C)$((sum(A,SAM(A,C))>0))=(PET0(C)*ES0(C)+XD0(C)*PD0(C))/PP0(C);

alphat(C)$(SAM(C,'ROW') ne 0)=XP0(C)/((deltat(C)*(ES0(C)**(rhot(C)))+(1-deltat(C))*(XD0(C))**(rhot(C)))**(1/rhot(C)));


**(iii) Armington Coefficient**
XMT0(C)$(SAM('ROW',C) ne 0)=(SAM('ROW',C)+SAM('Tarrif',C))/PMT0(C);
XA0(C)=(sum(ACT,SAM(ACT,C))-SAM(C,'ROW'))/PA0(C);
KDM(C)$((SAM('ROW',C) ne 0) and (sum(A,SAM(A,C))>0))=((XD0(C)/XMT0(C))*((PD0(C)/PMT0(C))**sigmaq(C)))**(1/sigmaq(C));
deltaq(C)$((SAM('ROW',C) ne 0)and (sum(A,SAM(A,C))>0))=KDM(C)/(1+KDM(C));

alphaq(C)$((SAM('ROW',C) ne 0)and (sum(A,SAM(A,C))>0))=XA0(C)/((deltaq(C)*(XD0(C)**(-rhoq(C)))+(1-deltaq(C))*(XMT0(C))**(-rhoq(C)))**(-1/rhoq(C)));




cwrt(C)=sum(H,SAM(C,H)/PA0(C))/sum((CP,H),SAM(CP,H));
cpi_o=sum(C,cwrt(C)*PA0(C));

*display deltat, deltaq, alphat, alphaq, cwrt, cpi_o;

tm_in=SAM('Ptaxin','ROW')/(sum(C,SAM('ROW',C))+sum(C,SAM('Tarrif',C)));
FSAV0=SAM('S-I','ROW');
* FSAD is exogenous. It will be fit to forcasted trade balance growth rate
FSAD=1;

** Household **

TR0(H)=SAM(H,'Gov');
* Oldpop: 2010 65+ popluation (Census, not projection) in thousand 5424.667 (Kosis)
Oldpop=5424.667;
tr0_per(H)=TR0(H)/Oldpop;

shr(F,H)=SAM(H,F)/sum(HP,SAM(HP,F));
Ls0(L)=sum(H,SAM(H,L))/W0(L);
Ks0(K)=sum(H,SAM(H,K))/R0(K);
LY0(H)=sum(L,SAM(H,L));
KY0(H)=sum(K,SAM(H,K));
ResinC0=sum(H,SAM(H,'NRES'));
YH0(H)=TR0(H)+LY0(H)+KY0(H)+ResinC0;
TINS0=sum(H,SAM('Ytax',H))/(sum(H,YH0(H))-delta*(sum(K,Ks0(K))));
TINSR0=(sum(H,YH0(H)*TINS0)-TINS0*delta*(sum(K,Ks0(K)))-crevtax*TCREV0)/(sum(H,YH0(H))-delta*(sum(K,Ks0(K))));
YD0(H)=(1-TINSR0)*(YH0(H)-sum(K,shr(K,H)*delta*Ks0(K)));
tc_in=sum(H,SAM('Ptaxin',H))/sum((C,H),SAM(C,H));
XAC0(C,H)=SAM(C,H)/PA0(C);
SH0(H)=SAM('S-I',H);


mu(C,H)=XAC0(C,H)*PA0(C)*(1+tc_in)/YD0(H);
mus(H)=(SH0(H)-sum(K,shr(K,H)*delta*Ks0(K)))/YD0(H);
*mu(C)=SAM22(C,'Household')/(SAM22('Household','Total')-SAM22('inc_tax','Household'));
*mus=SAM22('S-I','Household')/(SAM22('Household','Total')-SAM22('inc_tax','Household'));
MPS0(H)=mus(H);

Lw0(L)=Ls0(L)/(((1-TINSR0)*W0(L)/cpi_o)**epsilon_L(L));


**Investment
XAF0('S-I',C)=SAM(C, 'S-I')/PA0(C);
qinv_o(C)=SAM(C, 'S-I')/PA0(C);
IVAD0=sum(C,XAF0('S-I',C))/sum(CP,qinv_o(CP));
tiv_in=SAM('Ptaxin','S-I')/sum(C,SAM(C, 'S-I'));
ResinI0=SAM('S-I','NRES');
*display QINV0,qinv_o,IVAD0,tiv_in,ResinI0;


** Government
qg0(C)=SAM(C,'Gov')/PA0(C);
sg0=SAM('S-I','Gov');
XAF0('Gov',C)=qg0(C)*cpi_o;
SG0=sg0*cpi_o;
tg_in=SAM('Ptaxin','Gov')/sum(C,PA0(C)*XAF0('Gov',C));
qgr0(C)=PA0(C)*(1+tg_in)*qg0(C)/sum(CP,PA0(CP)*XA0(CP));
YG0=SAM('Gov','Ptaxin')+SAM('Gov','Ptaxetc')+SAM('Gov','Tarrif')+SAM('Gov','Ytax');

*display qg0,sg0,YG0,tg_in;

**Warlas dummy initial value

alpha_nres(A)=SAMng('NRES',A)/sum(ACT,samng(ACT,A));
thetaRes_c=sum(H,SAMng(H,'NRES'))/sum(A,SAMng('NRES',A));
thetaRes_iv=SAMng('S-I','NRES')/sum(A,SAMng('NRES',A));

Warlas0=sum(C$(sam(C,'S-I') ne 0),PA0(C)*(1+tiv_in)*XAF0('S-I',C))-sum(H,SH0(H))-SG0-FSAV0*EXR0-ResinI0;

*parameter
*VAEPr1_LHS(A)
*VAEPr1_RHS(A)
*;

*VAEPr1_LHS(A)$(sum(C,SAM(A,C))>0)=PVAE0(A);

*VAEPr1_RHS(A)$(sum(C,SAM(A,C))>0)=(1/alphaaVAE(A))*
*(
*deltaXEP(A)**sigmaaVAE(A)*PEP0(A)**(1-sigmaaVAE(A))
*+ deltaVA(A)**sigmaaVAE(A)*PVA0(A)**(1-sigmaaVAE(A))
*)**(1/(1-sigmaaVAE(A)));

*display VAEPr1_LHS, VAEPr1_RHS;

*parameter
*VAPr1_LHS(A)
*VAPr1_RHS(A)
*;
*VAPr1_LHS(A)$(sum(C,SAM(A,C))>0)=PVA0(A);
*VAPr1_RHS(A)$(sum(C,SAM(A,C))>0)=
*(
*     prod(K,(R0(K)/(lambdak*lambdaka(A)*deltaf(K,A)))**(deltaf(K,A)))*
*     prod(L,(W0(L)/(lambdat*lambda(A)*deltaf(L,A)))**(deltaf(L,A)))
*);

*display VAPr1_LHS,VAPr1_RHS;
*parameter
*XEPr1_LHS(A)
*XEPr1_RHS(A);

*XEPr1_LHS(A)$(sum(C,SAM(A,C))>0)=PEP0(A);
*XEPr1_RHS(A)$(sum(C,SAM(A,C))>0)=
*  prod(C$ELECC(C),(PA0(C)/deltaC(C,A))**(deltaC(C,A)))*
*(
*PFL0(A)/
*(1-sum(C$ELECC(C),deltaC(C,A)))
*)**
*(
*1-sum(C$ELECC(C),deltaC(C,A))
*)
*;

*display XEPr1_LHS, XEPr1_RHS;

*parameter
*XELPr0_LHS(A)
*XELPr0_RHS(A);

*XELPr0_LHS(A)$(sum(C$sfuel(C),SAM(C,A)=0))=PFL0(A);
*XELPr0_RHS(A)$(sum(C$sfuel(C),SAM(C,A)=0))=PLFL0(A);


*display XELPr0_LHS,XELPr0_RHS;
*parameters
*XFLPr1_LHS(A)
*XFLPr1_RHS(A);

*XFLPr1_LHS(A)$(sum(C$sfuel(C),SAM(C,A) >0))=PFL0(A);
*XFLPr1_RHS(A)$(sum(C$sfuel(C),SAM(C,A) >0))=
*(
*sum(C$SfuelA(A,C),deltaC(C,A)**sigmaaXFL(A)*(PNEG0(C,A)/AEEI(C,A))**(1-sigmaaXFL(A)))
*+
*(1-sum(C$SfuelA(A,C),deltaC(C,A)))**(sigmaaXFL(A))*PLFL0(A)**(1-sigmaaXFL(A))
*)**(1/(1-sigmaaXFL(A)));

*display XFLPr1_LHS, XFLPr1_RHS;

*parameter
*XLFLPr0_LHS(A)
*XLFLPr0_RHS(A);

*XLFLPr0_LHS(A)$(Lfuelmix(A) eq 0)=PLFL0(A);
*XLFLPr0_RHS(A)$(Lfuelmix(A) eq 0)=sum(C$(ENC(C) and LfuelA(A,C)),PNEG0(C,A))/sum(C$(ENC(C) and LfuelA(A,C)),AEEI(C,A));

*display
*XLFLPr0_LHS, XLFLPr0_RHS;

*parameter
*XLFLPr1_LHS(A)
*XLFLPr1_RHS(A);

*XLFLPr1_LHS(A)$(Lfuelmix(A) eq 1)=PLFL0(A);
*XLFLPr1_RHS(A)$(Lfuelmix(A) eq 1)=
*prod(C$LfuelA(A,C),(PNEG0(C,A)/(deltaC(C,A)*AEEI(C,A)))**deltaC(C,A))
*;

*display XLFLPr1_LHS,XLFLPr1_RHS;

*parameter
*XVAD1_LHS(A)
*XVAD1_RHS(A);

*XVAD1_LHS(A)$(sum(C,SAM(A,C))>0)=VA0(A);
*XVAD1_RHS(A)$(sum(C,SAM(A,C))>0)=(deltaVA(A)**sigmaaVAE(A))*((PVAE0(A)/PVA0(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE0(A);

*display XVAD1_LHS,XVAD1_RHS;
*parameter
*XEPD1_LHS(A)
*XEPD1_RHS(A)
*;
*XEPD1_LHS(A)$(sum(C,SAM(A,C))>0)=XEP0(A);
*XEPD1_RHS(A)$(sum(C,SAM(A,C))>0)=(deltaXEP(A)**sigmaaVAE(A))*((PVAE0(A)/PEP0(A))**sigmaaVAE(A))*( alphaaVAE(A)**(sigmaaVAE(A)-1))*QVAE0(A);

*display XEPD1_LHS,XEPD1_RHS;

*parameter
*LDA1_LHS(L,A)
*LDA1_RHS(L,A);

*LDA1_LHS(L,A)$(sum(C,SAM(A,C))>0)=LD0(L,A);
*LDA1_RHS(L,A)$(sum(C,SAM(A,C))>0)=deltaf(L,A)*(PVA0(A)*VA0(A))/W0(L);

*display LDA1_LHS,LDA1_RHS;
*parameter
*KDA1_LHS(K,A)
*KDA1_RHS(K,A);

*KDA1_LHS(K,A)$(sum(C,SAM(A,C))>0)=KD0(K,A);
*KDA1_RHS(K,A)$(sum(C,SAM(A,C))>0)=deltaf(K,A)*(PVA0(A)*VA0(A))/R0(K);

*display KDA1_LHS,KDA1_RHS;
*parameter
*XFLD1_LHS(A)
*XFLD1_RHS(A);

*XFLD1_LHS(A)$(sum(C,SAM(A,C))>0)=XFL0(A);
*XFLD1_RHS(A)$(sum(C,SAM(A,C))>0)=(1-sum(C$ELECC(C),deltaC(C,A)))*(PEP0(A)/PFL0(A))*XEP0(A);

*display XFLD1_LHS,XFLD1_RHS;

*parameter
*INTDE1_LHS(C,A)
*INTDE1_RHS(C,A)
*;
*INTDE1_LHS(C,A)$(ELECC(C) and SAM(C,A)$ELECC(C)>0)=XAP0(C,A);
*INTDE1_RHS(C,A)$(ELECC(C) and SAM(C,A)$ELECC(C)>0)= deltaC(C,A)*(PEP0(A)/PA0(C))*XEP0(A);

*display INTDE1_LHS,INTDE1_RHS;


*parameter
*XLFLD0_LHS(A)
*XLFLD0_RHS(A);

*XLFLD0_LHS(A)$(sum(C$sfuel(C),SAM(C,A)=0))=XLFL0(A);
*XLFLD0_RHS(A)$(sum(C$sfuel(C),SAM(C,A)=0))=XFL0(A);

*display XLFLD0_LHS,XLFLD0_RHS;
*parameter
*NEGDS1_LHS(C,A)
*NEGDS1_RHS(C,A);

*NEGDS1_LHS(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0))=QNEG0(C,A);
*NEGDS1_RHS(C,A)$(Sfuel(C) and (sum(CP$sfuel(CP),SAM(CP,A)) >0))=(deltaC(C,A)**sigmaaXFL(A))*((PFL0(A)/PNEG0(C,A))**sigmaaXFL(A))*((AEEI(C,A))**(sigmaaXFL(A)-1))*XFL0(A);

*display NEGDS1_LHS,NEGDS1_RHS;

*parameter
*XLFLD1_LHS(A)
*XLFLD1_RHS(A);

*XLFLD1_LHS(A)$(sum(C$sfuel(C),SAM(C,A) >0))=XLFL0(A);
*XLFLD1_RHS(A)$(sum(C$sfuel(C),SAM(C,A) >0))=(1-sum(C$sfuel(C),deltaC(C,A)))**sigmaaXFL(A)*(PFL0(A)/PLFL0(A))**(sigmaaXFL(A))*XFL0(A);

*display XLFLD1_LHS,XLFLD1_RHS;

*parameter
*NEGDL0_LHS(C,A)
*NEGDL0_RHS(C,A);

*NEGDL0_LHS(C,A)$(Lfuel(C) and Lfuelmix(A) eq 0)=QNEG0(C,A);
*NEGDL0_RHS(C,A)$(Lfuel(C) and Lfuelmix(A) eq 0)=XLFL0(A);

*display NEGDL0_LHS,NEGDL0_RHS;
*parameter
*NEGDL1_LHS(C,A)
*NEGDL1_RHS(C,A);


*NEGDL1_LHS(C,A)$(Lfuel(C) and Lfuelmix(A) eq 1)=QNEG0(C,A);
*NEGDL1_RHS(C,A)$(Lfuel(C) and Lfuelmix(A) eq 1)=deltaC(C,A)*(PLFL0(A)/PNEG0(C,A))*XLFL0(A);

*display NEGDL1_LHS,NEGDL1_RHS;
*parameter
Model
BR 7 ind model
/ImPr.XMT
ExPr.ES
AspPr.XA
AspPrni.XA
AspPrnd.XA
ProdPr.XD
ProdPrne.XD
ProdPrnd.XD
ComPr.XP
ActR.XC
ActRp.XC
VAEPr0
VAEPr1
VAEPr2
*HKTEPr
VAPr0
VAPr1
VApr2
XEPr0
XEPr1
XEPr2
XFLPr0
XFLPr1
XFLPr2
XFLPr3
XLFLPr0
XLFLPr1
NEGPr.QNEG
Norm
CPIfix
QVAED
*HKTED
XVAD0
XVAD1
XEPD1
XEPD2
XFLD1
XFLD2
XLFLD0
XLFLD1
INTDM
INTDE1
INTGD
LDA1
KDA1
NEGDL0.PNEG
NEGDL1.PNEG
NEGDS1.PNEG
NEGDS2.PNEG
NEGDS3.PNEG
NELQCED
GD
ActDC.PC
XDD.PD
XDDni.PD

XMTD.PMT
XMTDnd.PMT

ESS.PET
ESSnd.PET

XDS.PP
XDSne.PP

Labsup

InvD
HouseLY
HouseKY
HouseYD
HouseY
Saver
HouseD
Hsave
GovE
Tras
Ytax
ForS
GovI

LabM.W
CapM.R
ComMENCN.PA
ComMENC.PA
CREVE
CREVP
TCREVsum

InvM
CAB
GovB
ResC
ResI/;


*** setting up initial values *****

PC.L(A)        =        PC0(A)        ;
XC.L(A)        =        XC0(A)        ;
XMT.L(C)        =        XMT0(C)        ;
PMT.L(C)        =        PMT0(C)        ;
XD.L(C)        =        XD0(C)        ;
PD.L(C)        =        PD0(C)        ;
XA.L(C)        =        XA0(C)        ;
PA.L(C)        =        PA0(C)        ;
ES.L(C)        =        ES0(C)        ;
PET.L(C)        =        PET0(C)        ;
XP.L(C)        =        XP0(C)        ;
PP.L(C)        =        PP0(C)        ;
QNEG.L(C,A)        =        QNEG0(C,A)        ;
PNEG.L(C,A)        =        PNEG0(C,A)        ;
PVAE.L(A)        =        PVAE0(A)        ;
QVAE.L(A)        =        QVAE0(A)        ;
PVA.L(A)=PVA0(A);
VA.L(A)=VA0(A);
XEP.L(A)=XEP0(A);
PEP.L(A)=PEP0(A);
XFL.L(A)=XFL0(A);
PFL.L(A)=PFL0(A);
XLFL.L(A)=XLFL0(A);
PLFL.L(A)=PLFL0(A);
XAP.L(C,A)        =        XAP0(C,A)        ;
QINTG.L(GC,A)        =        QINTG0(GC,A)        ;
R.L(K)        =        R0(K)        ;
W.L(L)        =        W0(L)        ;
LD.L(L,A)        =        LD0(L,A)        ;
KD.L(K,A)        =        KD0(K,A)        ;
Ks.Fx(K)        =        Ks0(K)        ;
Ls.L(L)        =        Ls0(L)        ;
QCE.L(C,A)        =        QCE0(C,A)        ;
QGE.L(GC,C,A)        =        QGE0(GC,C,A)        ;
EXR.L        =        EXR0        ;
YG.L        =        YG0        ;
TR.L(H)        =        TR0(H)        ;
XAF.L(FD,C)        =        XAF0(FD,C)        ;
SG.L        =        SG0        ;
MPS.L(H)        =        MPS0(H)        ;
LY.L(H)        =        LY0(H)        ;
KY.L(H)        =        KY0(H)        ;
YH.L(H)        =        YH0(H)        ;
YD.L(H)        =        YD0(H)        ;
XAC.L(C,H)        =        XAC0(C,H)        ;
SH.L(H)        =        SH0(H)        ;
TINSR.L        =        TINSR0        ;
FSAV.L        =        FSAV0        ;
IVAD.L        =        IVAD0        ;
Warlas.L        =        Warlas0        ;
CREV.L(A)        =        CREV0(A)        ;
TCREV.L        =        TCREV0        ;
ResinC.L        =        ResinC0        ;
ResinI.L        =        ResinI0        ;
CPI.L = cpi_o;
BR.HOLDFIXED=1;

SET ACGDP GDP Items
/GDPMP1  GDP Aggregate Demand_market price
CON_P private consumption
CON_G government consumption
INV   Investment
EXP   EXPORT
IMP   IMPORT
NITAX Indirect tax
NRES  Net Residue
GDPFP GDP Factor Price and Residue
GDPMP2 GDP Aggregate income_market price
GAP GDP_AD-GDP_Y/
ACGDP1(ACGDP)
/CON_P private consumption
CON_G government consumption
INV   Investment
EXP   EXPORT
IMP   IMPORT/

parameters
PCREP(A,t)
XCREP(A,t)
XMTREP(C,t)
PMTREP(C,t)
XDREP(C,t)
PDREP(C,t)
XAREP(C,t)
PAREP(C,t)
ESREP(C,t)
PETREP(C,t)
XPREP(C,t)
PPREP(C,t)
QNEGREP(C,A,t)
PNEGREP(C,A,t)
PVAEREP(A,t)
QVAEREP(A,t)
PVAREP(A,t)
VAREP(A,t)
PEPREP(A,t)
XEPREP(A,t)
PFLREP(A,t)
XFLREP(A,t)
PLFLREP(A,t)
XLFLREP(A,t)
XAPREP(C,A,t)
QINTGREP(GC,A,t)
RREP(K,t)
WREP(L,t)
LDREP(L,A,t)
KDREP(K,A,t)
KsREP(K,t)
LsREP(L,t)
QCEREP(C,A,t)
QGEREP(GC,C,A,t)
EXRREP(t)
YGREP(t)
TRREP(H,t)
XAFREP(FD,C,t)
SGREP(t)
*QINVREP(C,t)
MPSREP(H,t)
LYREP(H,t)
KYREP(H,t)
YHREP(H,t)
YDREP(H,t)
XACREP(C,H,t)
SHREP(H,t)
TINSRREP(t)
FSAVREP(t)
IVADREP(t)
WarlasREP(t)
CREVREP(A,t)
TCREVREP(t)
ResinCREP(t)
ResinIREP(t)
BALCHK (t,ACT)
GDPREP(*,t)
SAMREP(t,ACT,ACTPP)
GCREP(t,GC,ACT,ACTP)
DGCREP(t,GC,A)

TGC(t,GC)

SAMDIFF(*,ACT,ACTPP)
DiFFT
gtaxrep(GC,t)
CPIREP(t)
;

**====================Growth parameters=======================================**

Parameter
lgrow(t)
*/0       0.017416
*1        0.01
*/


/0       0.017416
1        0.01
2        0.01
3        0.01
4        0.01
5        0.01
6        0.006
7        0.006
8        0.006
9        0.006
10        0.006
11        0.002
12        0.002
13        0.002
14        0.002
15        0.002
16       -0.002
17       -0.002
18       -0.002
19       -0.002
20       -0.002
21        0
/
lpgrow(t)

*/0       0.0170
*1        0.01155
*/
/0       0.0170
1        0.01155
2        0.01155
3        0.01155
4        0.01155
5        0.01155
6        0.0138
7        0.0138
8        0.0138
9        0.0138
10        0.0138
11        0.0122
12        0.0122
13        0.0122
14        0.0122
15        0.0122
16        0.0124
17        0.0124
18        0.0124
19        0.0124
20        0.0124
21        0
/

Oldpopg(t)
*/
*0        0.0373
*1        0.0374
*/

/
0        0.0373
1        0.0374
2        0.0414
3        0.0421
4        0.0404
5        0.0373
6        0.0362
7        0.0372
8        0.0389
9        0.0433
10        0.0477
11        0.0496
12        0.0502
13        0.0504
14        0.0506
15        0.0505
16        0.0493
17        0.0461
18        0.0417
19        0.0379
20        0.0351
21        0
/

TBg(t)
*/0        0.073
*1        0.073
*/

/0        0.073
1        0.073
2        0.056
3        -0.101
4        -0.009
5        0.097
6        0.081
7        0.049
8        0.066
9        0.060
10        0.057
11        0.028
12        0.041
13        0.023
14        0.042
15        0.020
16        0.021
17        0.037
18        0.015
19        0.012
20        0.012
21        0
/
;

**============================================================================**
set TO(t) /0/;
parameter GDP1(t),ygrow(t);

Loop (t,
SOLVE BR Using MCP;
*)
PCREP(A,t)        =        PC.L(A)        ;
XCREP(A,t)        =        XC.L(A)        ;
XMTREP(C,t)        =        XMT.L(C)        ;
PMTREP(C,t)        =        PMT.L(C)        ;
XDREP(C,t)        =        XD.L(C)        ;
PDREP(C,t)        =        PD.L(C)        ;
XAREP(C,t)        =        XA.L(C)        ;
PAREP(C,t)        =        PA.L(C)        ;
ESREP(C,t)        =        ES.L(C)        ;
PETREP(C,t)        =        PET.L(C)        ;
XPREP(C,t)        =        XP.L(C)        ;
PPREP(C,t)        =        PP.L(C)        ;
QNEGREP(C,A,t)        =        QNEG.L(C,A)        ;
PNEGREP(C,A,t)        =        PNEG.L(C,A)        ;
PVAEREP(A,t)        =        PVAE.L(A)        ;
QVAEREP(A,t)        =        QVAE.L(A)        ;
PVAREP(A,t)=PVA.L(A);
VAREP(A,t)=VA.L(A);
PEPREP(A,t)=PEP.L(A);
XEPREP(A,t)=XEP.L(A);
PFLREP(A,t)=PFL.L(A);
XFLREP(A,t)=XFL.L(A);
PLFLREP(A,t)=PLFL.L(A);
XLFLREP(A,t)=XLFL.L(A);
XAPREP(C,A,t)        =        XAP.L(C,A)        ;
QINTGREP(GC,A,t)        =        QINTG.L(GC,A)        ;
RREP(K,t)        =        R.L(K)        ;
WREP(L,t)        =        W.L(L)        ;
LDREP(L,A,t)        =        LD.L(L,A)        ;
KDREP(K,A,t)        =        KD.L(K,A)        ;
KsREP(K,t)        =        Ks.L(K)        ;
LsREP(L,t)        =        Ls.L(L)        ;
QCEREP(C,A,t)        =        QCE.L(C,A)        ;
QGEREP(GC,C,A,t)        =        QGE.L(GC,C,A)        ;
EXRREP(t)        =        EXR.L        ;
YGREP(t)        =        YG.L        ;
TRREP(H,t)        =        TR.L(H)        ;
XAFREP(FD,C,t)        =        XAF.L(FD,C)        ;
SGREP(t)        =        SG.L        ;
MPSREP(H,t)        =        MPS.L(H)        ;
LYREP(H,t)        =        LY.L(H)        ;
KYREP(H,t)        =        KY.L(H);
YHREP(H,t)        =        YH.L(H)        ;
YDREP(H,t)        =        YD.L(H)        ;
XACREP(C,H,t)        =        XAC.L(C,H)        ;
SHREP(H,t)        =        SH.L(H)        ;
TINSRREP(t)        =        TINSR.L        ;
FSAVREP(t)        =        FSAV.L        ;
IVADREP(t)        =        IVAD.L        ;
WarlasREP(t)        =        Warlas.L        ;
CREVREP(A,t)        =        CREV.L(A)        ;
TCREVREP(t)        =        TCREV.L        ;
ResinCREP(t)        =        ResinC.L        ;
ResinIREP(t)        =        ResinI.L        ;
CPIREP(t)           = CPI.L;
gtaxrep(GC,t)=gtax(GC);

SAMREP(t,ENCN,A)=XAP.L(ENCN,A)*PA.L(ENCN);
SAMREP(t,ENC,A)=QCE.L(ENC,A)*PA.L(ENC);
SAMREP(t,K,A)=KD.L(K,A)*R.L(K);
SAMREP(t,L,A)=LD.L(L,A)*W.L(L);
SAMREP(t,'NRES',A)=alpha_nres(A)*XC.L(A)*PC.L(A);
SAMREP(t,'Ptaxin',A)=ta_in(A)*XC.L(A)*PC.L(A);
SAMREP(t,'Ptaxetc',A)=ta_ex(A)*XC.L(A)*PC.L(A);


SAMREP(t,A,C)=theta(A,C)*XC.L(A)*PC.L(A);
SAMREP(t,'Tarrif',C)=tm(C)*pwm(C)*EXR.L*XMT.L(C);
SAMREP(t,'ROW',C)=pwm(C)*EXR.L*XMT.L(C);



SAMREP(t,H,L)=W.L(L)*shr(L,H)*Ls.L(L);
SAMREP(t,H,K)=R.L(K)*shr(K,H)*Ks.L(K);



SAMREP(t,C,H)=PA.L(C)*XAC.L(C,H);
SAMREP(t,'Ptaxin',H)=tc_in*sum(C,PA.L(C)*XAC.L(C,H));
SAMREP(t,'Ytax',H)=TINSR.L*(YH.L(H)-sum(K,delta*shr(K,H)*KS.L(K)));
SAMREP(t,'S-I',H)=SH.L(H);

SAMREP(t,C,'Gov')=PA.L(C)*XAF.L('Gov',C);
SAMREP(t,H,'Gov')=TR.L(H);
SAMREP(t,'Ptaxin','Gov')=sum(C,tg_in*PA.L(C)*XAF.L('Gov',C));
SAMREP(t,'S-I','Gov')=SG.L;
*
SAMREP(t,'Household','NRES')=ResinC.L;
SAMREP(t,'S-I','NRES')=ResinI.L;

SAMREP(t,'Gov','Ptaxin')=sum(A,ta_in(A)*PC.L(A)*XC.L(A))+(tc_in)*sum((C,H)$(sam(C,H) ne 0),PA.L(C)*XAC.L(C,H))+tg_in*sum(C$(sam(C,'Gov')ne 0),PA.L(C)*XAF.L('Gov',C))+tiv_in*sum(C$(sam(C,'S-I') ne 0),PA.L(C)*XAF.L('S-I',C))+tm_in*sum(C$(sam('ROW',C) ne 0),PMT.L(C)*XMT.L(C));
SAMREP(t,'Gov','Ptaxetc')=sum(A,ta_ex(A)*PC.L(A)*XC.L(A));
SAMREP(t,'Gov','Tarrif')=sum(C$(sam('ROW',c)>0),tm(C)*EXR.L*pwm(C)*XMT.L(C));
SAMREP(t,'GOV','Ytax')=sum(H, TINSR.L*(YH.L(H)-sum(K,delta*shr(K,H)*KS.L(K))));



SAMREP(t,C,'S-I')=PA.L(C)*XAF.L('S-I',C);
SAMREP(t,'Ptaxin','S-I')=sum(C,tiv_in*PA.L(C)*XAF.L('S-I',C));
SAMREP(t,C,'ROW')=PET.L(C)*ES.L(C);
SAMREP(t,'Ptaxin','ROW')=tm_in*sum(C$(sam('ROW',c)>0),PMT.L(C)*XMT.L(C));
SAMREP(t,'S-I','ROW')=FSAV.L*EXR.L;

GDPREP('CON_P',t)=sum((C,H), PA.L(C)*XAC.L(C,H));
GDPREP('CON_G',t)=sum(C,PA.L(C)*XAF.L('Gov',C));
GDPREP('INV',t)=sum(C,PA.L(C)*XAF.L('S-I',C));
GDPREP('EXP',t)=sum(C$(sam(C,'ROW')>0),EXR.L*pwe(C)*ES.L(C));
GDPREP('IMP',t)=-1*sum(C$(sam('ROW',c)>0),EXR.L*pwm(C)*XMT.L(C));
GDPREP('GDPMP1',t)=sum(ACGDP1,GDPREP(ACGDP1,t));
GDPREP('GDPFP',t)=sum((K,A),R.L(K)*KD.L(K,A))+sum((L,A),W.L(L)*LD.L(L,A));
GDPREP('NITAX',t)=
sum(A$(sum(ACT,SAM(A,ACT))>0),(ta_in(A)+ta_ex(A))*PC.L(A)*XC.L(A))
+sum(C$(sam('ROW',c)>0),tm(C)*pwm(C)*XMT.L(C)*EXR.L)
+sum(C$(sam(C,'ROW')>0),te(C)*pwe(C)*ES.L(C)*EXR.L)
+TCREV.L
-sum(A,crevI*crevIw(A)*CREV.L(A));
GDPREP('NRES',t)=sum(A$(sum(ACT,SAM(A,ACT))>0),alpha_nres(A)*PC.L(A)*XC.L(A));
GDPREP('GDPMP2',t)=GDPREP('GDPFP',t)+GDPREP('NITAX',t)+GDPREP('NRES',t);
GDPREP('GAP',t)=GDPREP('GDPMP1',t)-GDPREP('GDPMP2',t);

GCREP(t,GC,C,A)=QGE.L(GC,C,A);
GCREP(t,GC,'process',A)=QINTG.L(GC,A);
GCREP(t,GC,'Total',A)=sum(ACT,GCREP(t,GC,ACT,A));

DGCREP(t,GC,A)=GCREP(t,GC,'Total',A);
TGC(t,GC)=sum(A,GCREP(t,GC,'Total',A));
BALCHK(t,ACT)=sum(ACTPP,SAMREP(t,ACT,ACTPP))-sum(ACTPP,SAMREP(t,ACTPP,ACT));
**==== Law of Motion =====*
Ks.Fx('Capital')=Ks.L('Capital')*(1-delta)+sum(C,XAF.L('S-I',C));
Lw0(L)=Lw0(L)*(1+lgrow(t));
Oldpop=Oldpop*(1+Oldpopg(t));
FSAD=FSAD*(1+TBg(t));
lambdat=lambdat+lpgrow(t);
gtax(GC)=gtax(GC)+gtax_policy(GC,t);
* =========================== *


);
SAMDIFF('0',ACT,ACTPP)=samng(ACT,ACTPP)-SAMREP('0',ACT,ACTPP);
DiFFT=sum((ACT,ACTPP),SAMDIFF('0',ACT,ACTPP));

GDP1(t)=GDPREP('GDPMP1',t);
loop (t$(not TO(t)),ygrow(t)=(GDP1(t)-GDP1(t-1))/GDP1(t-1));

display scenario;
display SAMDIFF,DiFFT;
display BALCHK;
*display CREVREP,TCREVREP;
display warlasrep;
*display PVAEREP,QVAEREP;
*display PHKTEREP,HKTEREP;
*display SAMREP;
display GDPrep,GDP1;
*display ygrow;
display TGC;
display gtaxrep;
display CPIrep;
*display R0,delta;


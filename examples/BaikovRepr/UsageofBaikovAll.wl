(* ::Package:: *)

(*we list all four examples in one .wl file for the convenience of demonstration.*) 
(*However, we recommand each example in one file.*)


(*using same way of loading as AMflow and Blade*)
current = If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName;
Get[current<>"../../"<>"BaikovAll.wl"];


(* ::Section:: *)
(*Outer-Massive double box*)


dlist={
	PD[k1,0,-m2],
	PD[k1,-p1,-m2],
	PD[k1,-p1-p2,-m2],
	PD[k1-k2,0,0],
	PD[k2,-p1-p2,-m2],
	PD[k2,-p1-p2-p3,-m2],
	PD[k2,0,-m2],
	PD[k2,-p1,-m2],
	PD[k1,-p1-p2-p3,-m2]
};(*input of propagators and independent scalar products*)
conservation={p4->-p1-p2-p3};(*if dependent variables does not appear, this can be an empty set*)
kinematics={
	SProd[p1]->0,
	SProd[p2]->0,
	SProd[p3]->0,
	SProd[p4]->0,
	SProd[p1+p2]->s,
	SProd[p2+p3]->t
};(*define independent scalar products, SProd[a]=SProd[a,a]*)


kinerep=GetIndExtSProd[kinematics,conservation];(*represent all scalar products of external momenta with independent ones*)
krep=Join[BaikovTrans[dlist,kinerep][[2]],kinerep];(*a replacement rule for all scalar products involving loop momenta*)


stdrep=BaikovRep[dlist,{k1,k2},{p1,p2,p3},Abst->True];(*the standard representation with an abstract form*)
resultmat=AllSectorBaikovMat[stdrep,krep,Exc->{}];(*derive all representations from the standard one*)


Export[current<>"db_rep.wxf",resultmat,PerformanceGoal->"Size"];


(*we can now explore some function in this family*)


(*1. Get the loop by loop representation from the whole representations*)
(*now we need to specify the sector we want, {1,2,3,4,5,6,7} specifies which propagators exist*)
(*two independent loop-by-loop representation will be extracted*)
toprep=GetBaikovMatRep[resultmat,{1,2,3,4,5,6,7}]


(*2. get the explicit expressions for the representation*)
(*we can also get the dimension of the top topology under maximal cut*)
(*dim means the number of master integrals under maximal cut without considering symmetry*)
toprepexp = Times@@Power@@@toprep[[1,2,1]]//Gram2Poly[#,krep]&//Simplify
dim = GetDimension[toprepexp/.{Subscript[x, i_/;(i<8)]:>0},{Subscript[x, 8]}]


(*3. considering the symmetry of top sector by analyzing the polynomial in standard representation*)
poly=stdrep[[1,2,1]]//Gram2Poly[#,krep]&;
symrule=PolySym[poly,Join[Subscript[x,#]&/@Range[$BvNum],{s,t,m2}],$BvNum]


(*check these are indeed the symmetry of this polynomial*)
PolySymCheck[poly,symrule[[2]]]


(*4. get all zero sectors derived from Baikov representation*)
zerolist=GetMatZeroSector[resultmat,$BvNum,{8,9}](*the last list is where the isps locate in dlist*)
PadRight[Reverse[IntegerDigits[#,2]],$BvNum]&/@zerolist (*look up these sector numbers as binary lists*)


(* ::Section:: *)
(*Linear propagator involved*)


(*now we present above example with linear propagator involved*)
(*just replace PD with LPD: LPD[k,p,-m2]=k.p-m2, PD[k,p,-m2]=(k+p)^2-m2*)
dlist={
	PD[k1,0,-m2],
	LPD[k1,-p1,-m2],
	LPD[k1,-p1-p2,-m2],
	LPD[k1,k2,0],
	PD[k2,-p1-p2,-m2],
	PD[k2,-p1-p2-p3,-m2],
	PD[k2,0,-m2],
	LPD[k2,-p1,-m2],
	PD[k1,-p1-p2-p3,-m2]
};(*input of propagators and independent scalar products*)
conservation={p4->-p1-p2-p3};(*if dependent variables does not appear, this can be an empty set*)
kinematics={
	SProd[p1]->0,
	SProd[p2]->0,
	SProd[p3]->0,
	SProd[p4]->0,
	SProd[p1+p2]->s,
	SProd[p2+p3]->t
};(*define independent scalar products, SProd[a]=SProd[a,a]*)


kinerep=GetIndExtSProd[kinematics,conservation];(*represent all scalar products of external momenta with independent ones*)
nkrep=Join[BaikovTrans[dlist,kinerep][[2]],kinerep];(*a replacement rule for all scalar products involving loop momenta*)


(*one can compare the result with the former example*)
nkrep(*new one with linear propagator*)
krep(*former one with all normal propagators*)


(*derive the representations in the same way*)
stdrep=BaikovRep[dlist,{k1,k2},{p1,p2,p3},Abst->True];(*the standard representation with an abstract form*)
nresultmat=AllSectorBaikovMat[stdrep,nkrep,Exc->{}];(*derive all representations from the standard one*)


(* ::Section:: *)
(*additional parameters in definition of propagators*)


(*in some cases, there may be parameters between the combinations of momenta*)(*We thank Xing Wang for providing this example*)
(*we need to declare these variables as numbers*)
DeclareVarAsNum[{r,z}]
dlist={
	PD[l1, -(1 - z)*k1, 0],
	PD[l1, z*k1, 0],
	PD[l2, 0, 0],
	PD[l3, 0, 0],
	PD[l1 - l2, -(1 - z)*k1, 0],
	PD[l1 + l3, z*k1, 0],
	PD[l1 + l3, z*k1 + k2, 0],
	PD[l1 - l2 + l3, z*k1 + k2, 0],
	PD[l1 - l2, z*k1 + k2, 0],
	PD[l2, k1, 0],
	PD[l3, k1, 0],
	PD[l1, k2, 0]
};(*input of propagators and independent scalar products*)
conservation={};(*if dependent variables does not appear, this can be an empty set*)
kinematics={
	SProd[k1]->0,
	SProd[k2]->0,
	SProd[k1,k2]->s/2
};(*define independent scalar products, SProd[a]=SProd[a,a]*)


kinerep=GetIndExtSProd[kinematics,conservation];(*represent all scalar products of external momenta with independent ones*)
krep=Join[BaikovTrans[dlist,kinerep][[2]],kinerep];(*a replacement rule for all scalar products involving loop momenta*)


krep


stdrep=BaikovRep[dlist,{l1,l2,l3},{k1,k2},Abst->True];(*the standard representation with an abstract form*)
resultmat=AllSectorBaikovMat[stdrep,krep];(*derive all representations from the standard one*)


Export[current<>"para_rep.wxf",resultmat,PerformanceGoal->"Size"];


Options[GetBaikovMatRep]


GetBaikovMatRep[resultmat,{1,2,3,4,5,6,7,8,9}]//Simplify


(* ::Section:: *)
(*five-loop banana: Derive representation starting from somewhere not being the standard representation*)


dlist={
	PD[k1,0,-m2],
	PD[k2-k1,0,0],
	PD[k3-k2,0,0],
	PD[k4-k3,0,0],
	PD[k5-k4,0,0],
	PD[k5,p,-m2],
	PD[k2,0,-m2],
	PD[k3,0,-m2],
	PD[k4,0,-m2],
	PD[k5,0,-m2],
	PD[k1,p,-m2],
	PD[k2,p,-m2],
	PD[k3,p,-m2],
	PD[k4,p,-m2],
	PD[k5-k1,0,0],
	PD[k5-k2,0,0],
	PD[k5-k3,0,0],
	PD[k4-k1,0,0],
	PD[k4-k2,0,0],
	PD[k3-k1,0,0]
};(*input of propagators and independent scalar products*)
conservation={};(*if dependent variables does not appear, this can be an empty set*)
kinematics={
	SProd[p]->s
};(*define independent scalar products, SProd[a]=SProd[a,a]*)


kinerep=GetIndExtSProd[kinematics,conservation];(*represent all scalar products of external momenta with independent ones*)
krep=Join[BaikovTrans[dlist,kinerep][[2]],kinerep];(*a replacement rule for all scalar products involving loop momenta*)


(*we use Exc->{...} to first specify which variables are not integrated out.*)
(*Then we get all possible loop-by-loop representations for the top sector*)
(*at last, we start from these loop-by-loop representations and derive all the relevant representations for this family*)
AbsoluteTiming[
stdrep=BaikovRep[dlist,{k1,k2,k3,k4,k5},{p},Abst->True];
resultmat=AllSectorBaikovMat[stdrep,krep,Exc->{Subscript[x, 1],Subscript[x, 2],Subscript[x, 3],Subscript[x, 4],Subscript[x, 5],Subscript[x, 6]}];
resultc=AllSectorBaikovMatC[resultmat,Subscript[x,#]&/@Range[20],krep];
Export[current<>"banana5l_rep.wxf",resultc,PerformanceGoal->"Size"];
]


(*there are totally 42 representations for this single top topology*)
(*but they are symmetric as we will see*)
toprep=GetBaikovMatRep[resultc,{1,2,3,4,5,6}];
toprep//Length


(*half in integer power, half in half-integer power*)
toprep[[1]]//Simplify
toprep[[-1]]//Simplify


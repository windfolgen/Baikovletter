(* ::Package:: *)

(*
Package Name: BaikovLetter
Version: 1.1.0
Description: Generate all rational letters and algebraic letters for families which belong to MPL functions. 
*)

(*
If you find any bug or suggest for this program, please contact: phamapku14@gmail.com
*)


(* ::Subsection:: *)
(*Begin*)


$BaikovPath = DirectoryName[$InputFileName];
If[FreeQ[$Path,$BaikovPath],PrependTo[$Path,$BaikovPath]];(*Load Baikov` package first*)
Get["Baikov`"];


WriteString["stdout","BaikovLetter: a package for analysing singularity for all kinds of families and \n generating candidate rational letters and algebraic letters for families which belong to MPL functions.\n"];
BeginPackage["BaikovLetter`",{"Baikov`"}];


dlogForm::usage="dlogForm[a,b] gives Log[(a+Sqrt[b])/(a-Sqrt[b])].";
Lambda::usage="Lambda[x,y,z] is the Kallen function x^2+y^2+z^2-2xy-2yz-2zx.";


AllRuleQ::usage="AllRuleQ[list] returns True if all elements in the list are rules.";


AbsGram::usage="Remove the minus sign before a Gram mat G[]";
EquivalentGramQ::usage="EquivalentGramQ[G1,G2,krep] finds whether two Gram determinants are equal to each other.";


Sector2Digits::usage="Sector2Digits[sector] transforms a sector list to a number. For example, Sector2Digits[{1,2,4}]=11.";
Digits2Sector::usage="Digits2Sector[digits] transforms a digital number to a sector. For example, Digits2Sector[11]={1,2,4}";


MatrixPosition::usage="MaxtriPosition[mat,fun] gives positions of element which satisfy fun[element]===True in matrix 'mat'.";
ManifestFactorized::usage="ManifestFactorized[mat,gram] return the factorized result if a matrix 'mat' is manifestly factorized. 'gram' is the corresponding Gram representation of this matrix.";

ExistRelationQ::usage="ExistRelationQ[list] find whether there are nontrivial linear relations with integer coefficients among elements in 'list'.";
CongruenceTrans::usage="CongruenceTrans[mat,gram] tries to find a congruence transformation that makes the matrix 'mat' manifestly factorized. 'gram' is the corresponding Gram representation of this matrix. It returns another gram which is congruent to the original one. The congruence transformation lies in integer domain.";

IsReducible::usage="IsReducible[gram,cut,krep] checks whether 'gram' is reducible under the cut condition 'cut','krep' is the replacement rule for Baikov variables and scalar products. If 'gram' is reducible, it will return the reduced result which should be a product (perfect square) of Grams. If 'gram' is not reducible, it will return original gram.\n The output is in the form {Ture/False,results}";


IsReducible::warning="The gram matrix `1` is 0 under the cut condition `2`.";


LeadingSingularities::usage="LeadingSingularities[rep,cut,krep] gives all the leading singularities related to a representation 'rep'. 'cut' is a set of conditions which specify propagators being cut, 'krep' is the replacement rule for scalar products. ";


Begin["`Private`"]


(* ::Subsection:: *)
(*General Functions*)


(*some primary functions*)
dlogForm[g1_,g2_]:=Log[(g1+Sqrt[g2])/(g1-Sqrt[g2])];
Lambda[x_,y_,z_]:=x^2+y^2+z^2-2x*y-2y*z-2x*z;


AllRuleQ[list_]:=If[list==={},True,AllTrue[list,MatchQ[#,Rule[_,_]]&]];


AbsGram[exp_]:=If[MatchQ[exp,-G[__]],Return[-exp],Return[exp]];

EquivalentGramQ[G1_,G2_,krep_]:=Catch@Module[{tem,tem1,var,var1,numrep},
If[Length[G1[[1]]]-Length[G2[[1]]]!=0,Throw[False]];
tem=G1/.{G[a_,b_]:>GramMat[a,b,krep]};
var=Variables[tem]//Sort;
tem1=G2/.{G[a_,b_]:>GramMat[a,b,krep]};
var1=Variables[tem1]//Sort;
If[var=!=var1,Throw[False]];
If[var==={},Throw[False]];(*if there are no variables at all, we keep them*)
Do[
numrep=Thread@Rule[var,RandomPrime[{100,1000},Length[var]]];
If[(Abs@Det[tem/.numrep]-Abs@Det[tem1/.numrep])=!=0,Throw[False]]
,{i,1,3}];
Throw[True];
];


Sector2Digits[sector_]:=Sum[Power[2,sector[[i]]-1],{i,1,Length[sector]}];
Digits2Sector[digits_]:=Position[IntegerDigits[digits,2]//Reverse,1]//Flatten;


(* ::Subsection:: *)
(*Core functions*)


MatrixPosition[mat_,fun_]:=Position[Map[fun,mat,{2}],True];


Options[ManifestFactorized]={"debug"->False};
ManifestFactorized[omat_,gram_G,OptionsPattern[]]:=Catch@Module[{mat,pos,len,int,row,max,temc,temr},
	mat=omat//Factor;
	len=Length[mat];(*the dimension of a square matrix*)
	pos=MatrixPosition[mat,#===0&]//GatherBy[#,First]&;(*gather by rows*)
	If[OptionValue["debug"],Print["pos: ",pos]];
	Do[
		If[Length[pos]<i,Throw[gram]];(*if the number of rows containing 0 are less than i, return original results*)
		row=Subsets[#[[All,2]]&/@pos,{i}];(*extract indices for rows and generate i-pairs of them*)
		int=Length/@(Intersection@@@row);(*check length of intersection for each pair*)
		max=Max[int];
		If[max<len-i,
			Continue[],(*if number of 0's is less than the remaining dimension, then the matrix can not be factorized*)
			If[max>len-i,Throw[0]](*if number of 0's is larger than the remaining dimension, then the matrix is 0*)
		];
		temc=Extract[row,Position[int,max][[1]]];
		temr=Select[pos,MemberQ[temc,#[[All,2]]]&,i];(*we only need to find i such rows that can give the max intersection number*)
		temr={temr[[All,1,1]],Complement[Range[len],temr[[All,1,1]]]};(*{rows containing 0, rows not containing 0}*)
		temc={Intersection@@temc,Complement[Range[len],Intersection@@temc]};(*{columns containing 0, columns not containing 0}*)
		If[OptionValue["debug"],Print["row info: ",temr];Print["column info: ",temc];];
		Throw[Power[-1,Total[temc[[2]]]+Total[temr[[1]]]]*ManifestFactorized[mat[[temr[[1]],temc[[2]]]],G[gram[[1,temr[[1]]]],gram[[2,temc[[2]]]]]]*ManifestFactorized[mat[[temc[[1]],temr[[2]]]],G[gram[[2,temc[[1]]]],gram[[1,temr[[2]]]]]]](*calculate the factorization recursively*)
	,{i,1,len}];
	Throw[gram];
];


ExistRelationQ[list_]:=Catch@Module[{len,sys,c,cv,var,numrep,sol},
	(*find if there are elements that can be combined to 0 in the list*)
	len=Length[list];
	var=Variables[list];
	numrep=Thread@Rule[var,RandomPrime[{10000,100000},Length[var]]];(*generate a set of random values for variables*)
	cv=(c/@Range[len]);
	sys={(cv) . (list/.numrep)==0};
	sys=Join[sys,Thread@LessEqual[Abs/@(cv),1],Thread@Equal[Complement[cv,Cases[sys,c[_],Infinity]],0],{cv . cv>0}];
	sol=FindInstance[And@@sys,cv,Integers];(*find a linear combination that gives 0*)
	If[sol==={},Throw[False],sol=sol[[1]]];
	Throw[cv/.sol];(*return the combination*)
];


Options[CongruenceTrans]={"debug"->False};
CongruenceTrans[mat_,gram_G,OptionsPattern[]]:=Module[{trans,cong,tem,tem1},
	cong=IdentityMatrix[Length[mat]];
	trans=mat//Transpose;
	tem1=gram;
	Do[
			tem=ExistRelationQ[trans[[i]]];
			If[OptionValue["debug"],Print["ExistRelationQ: ",tem]];
			If[tem===False,Continue[]];
			cong[[FirstPosition[tem,_?(#!=0&)][[1]]]]=tem;
			tem1=ManifestFactorized[cong . mat . Transpose[cong],G[cong . gram[[1]],cong . gram[[2]]]];
			If[tem1=!=gram,Break[],Continue[]];
	,{i,1,Length[trans]}];
	Return[{cong,tem1}];
];


Options[IsReducible]={"debug"->False};
IsReducible[gram_G,cut_,krep_,OptionsPattern[]]:=Module[{mat,cutsys,cutsol,tem},
	mat=gram//Gram2Mat[#,krep]&;(*transform the gram to matrix form*)
	cutsys=Thread@Equal[cut,0]//Gram2Mat[#,krep]&;(*solve the cut system for Baikov variables*)
	cutsol=Solve[cutsys,Cases[cutsys,Subscript[x,_],Infinity]//DeleteDuplicates][[1]];
	mat=mat/.cutsol//Factor;
	If[OptionValue["debug"],Print["mat: ",MatrixForm[mat]]];
	
	(*check whether the gram is manifestly factorizable*)
	tem=ManifestFactorized[mat,gram];
	If[tem===0,Message[IsReducible::warning,gram,cut]];
	If[tem=!=gram,Return[{True,tem}]];
	If[OptionValue["debug"],Print["tem: ",tem]];
	
	(*if it is not manifestly factorized, it can still possibly be factorzied after congruence transformation*)
	tem=CongruenceTrans[mat,gram];
	If[tem[[2]]=!=gram,
		(*try to reduce it recursively so that its simplest form achieved*)
		If[Head[tem[[2]]]===G,Return[{True,IsReducible[tem[[2]],cut,krep][[2]]}],Return[{True,tem[[2]]}]],
		Return[{False,gram}]
	];
];


LeadingSingularities[rep_,icut_,krep_,OptionsPattern[]]:=Module[{cut,tem,newrep,var,isp,temrep,temkrep,ls},
	(*'rep' is a basic element from the output of AllSectorBaikovMat[]. Its form is {{variables already integrated out},{glist,const}}*)
	(*for the first step, we check whether new representation can be generated*)
	var=rep[[2,1]]//Gram2Mat[#,krep]&//Cases[#,Subscript[x,_],Infinity]&//DeleteDuplicates;(*all Baikov variables involved*)
	If[Not@FreeQ[icut,x],cut=icut/.{Subscript[x,a_]:>a,x[a_]:>a},cut=icut];(*the input can be either a list of Subscript[x, i] or just the numbers*)
	isp=Complement[var,Subscript[x,#]&/@cut];(*all isps for this cut*)
	tem=NewReducibleRep[rep[[2,1]],isp,krep,"sector"->cut];
	If[tem[[1]],
		Print["    New representations found in rep id: ",rep[[1]]," with cut: ",cut];
		(*in this case, we need also consider the leading singularities of these new representations*)
		newrep=Table[{Append[rep[[1]],newrep[[i,2]]],{newrep[[i,1,1]],newrep[[i,1,2]]*rep[[2,2]]//FullSimplify},newrep[[i,3]]},{i,1,Length[newrep]}];
		newrep=Prepend[newrep,rep],
		newrep={{rep[[1]],rep[[2]],{}}};(*if no new representation found, then we just analysis the input representation*)
	];
	
	(*in the second step, we analyse every element of newrep. Note that the krep should be replaced if new representations have been obtained*)
	ls=Reap[
		Do[
			temrep=newrep[[i]][[{1,2}]];(*its form will be {{variables already integrated out},{glist,const}}*)
			temkrep=krep/.newrep[[i,3]];(*new krep we should use*)
			Sow[sLeadingSingularities[temrep,cut,temkrep]];
		,{i,1,Length[newrep]}]
	][[2,1]]//Flatten;(*collect all leading singularities*)
];


(* ::Subsection:: *)
(*End*)


End[];


EndPackage[];

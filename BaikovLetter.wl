(* ::Package:: *)

(*
Package Name: BaikovLetter
Version: 2.0.0
Description: Generate all rational letters and algebraic letters for planar families which belong to MPL functions. 
*)

(*
If you find any bug or suggest for this program, please contact: xhjiang@itp.ac.cn
*)


BeginPackage["BaikovLetter`",{"Baikov`"}]


dlogForm::usage="dlogForm[a,b] gives Log[(a+Sqrt[b])/(a-Sqrt[b])].";
Lambda::usage="Lambda[x,y,z] is the Kallen function x^2+y^2+z^2-2xy-2yz-2zx.";


Shrink::usage="Shrink[list] reduce a set of lists with first two same and third term different to one list with the third element to be a list the third item of all lists.";
ReformRep::usage="ReformRep[result] collects all Baikov Gram determinants in the representation and keep records of their power and path. Its output will be in the form: {{Gram,power,path},...}";
EquivalentGramQ::usage="EquivalentGramQ[G1,G2,krep] finds whether two Gram determinants are equal to each other.";


CompatibleGramQ::usage="CompatibleGramQ[G1,G2] detects whether G1 and G2 are permitable combination in the construction and return its type. 0 is not compatible, 1 is the first type of compatible combiantion GXiXi and GXjXj, 2 is the second type of compatible combination GXX and GXijXij.";
CompatibleGramQ::err="Input must be in the form G[__] or -G[__]";
CompatiblePathQ::usage="CompatiblePathQ[pathl,path] detects whether there is a path in path list pathl compatible with the second argument path. There is an option TopOnly, when it is true, only those paths which are parents of the second path will be detected. For example, {} is a parent path of {x9}. {x9,x8} is compatible with {x9} but not a parent path.";


AllSubRep::usage="AllSubRep[reform,rep] gets all the grams in a representation which belong to the top-sector representation rep.";
ExtractRelevantGramSim::usage="ExtractRelevantGramSim[rep,reform] extracts all compatible grams from a specific representation rep. rep is the output of GetBaikovMatRep[]. reform is the output of ReformRep[].";
ExtractRelevantGramSingle::usage="ExtractRelevantGramSingle[g,reform] extracts all compatible grams for a single gram list g. g is in the form {Gram,power,path}. reform is the output of ReformRep[].";


Sector2Digits::usage="Sector2Digits[sector] transforms a sector list to a number. For example, Sector2Digits[{1,2,4}]=11.";
Digits2Sector::usage="Digits2Sector[digits] transforms a digital number to a sector. For example, Digits2Sector[11]={1,2,4}";
PerfectSquareSplit::usage="PerfectSquareSplit[exp] split the perfect square in an expression exp. For internal use only.";
PerfectSquareSplit::warning="Power odd appear under the square root: `1`.";

TakeSquareRoot::usage="TakeSquareRoot[exp] take the square roots of an expression.";
EnhancedFactor::usage="EnhancedFactor[list,z] factor the expression with respect to z.";

ResolveSingleVariable::usage="ResolveSingleVariable[intset,hintset,z] resolves singularities of z from the integer power set: intset and half-integer power set: hintset.";
ResolveSingleVariable::warning="Odd power under square: `1` !";
ResolveSingleVariable::err="In ResolveSingleVariable[intset,hintset,z], intset and hintset can't both be {}";
ResolveSingleVariable::elliptic="Elliptic case encountered!";
ResolveSingleVariable::additional="This construction may be viable by proper variable transformation but we won't include it: {poly, poly list under square root} `1`.";

ResolveSingularities::usage="ResolveSingularities[intset,hintset] resolves singularities for all the variables in intset and hintset. It is the multivariate version of ResolveSingleVariable[].";
ResolveSingularities::err="There is something wrong with the arguments. They are all irrelevant to Baikov variables.";
ResolveSingularities::warning="In this case it is equivalent to integrating this variable out, so we won't consider its pole here but will include relevant Landau varieties. `1`.";
ResolveSingularities::oddpower="There is an odd power higher than 1 of variable under square root or it's a linear expression but with no other poles in the denominator! The infinity is a branch cut. Its expression is `1` and the input is `2`.";
ResolveSingularities::elliptic="Elliptic case encountered!";

Trans2Inf::usage="Trans2Inf[intset,hintset,xl] transforms the polynomial to infinity plane which corresponds to second Landau singularities.";


CheckValidity::usage="CheckValidity[brep,subset] checks whether a reprentation is valid under the maximal cut, if there is no valid representation, it will return {}.";


ExistRelationQ::usage="ExistRelationQ[poly,sqpolylist] detects whether poly has some relation with a list of other polynomials under square roots sqpolylist.";


FindCover::usage="FindCover[list] find a cover of the set of all variables present from a set of sets of these variables.";


PolesAnalyze::usage="PolesAnalyze[result,topsector,krep,n] gives the possible poles, their paths and leading singularities related. topsector is specified by a number list like {1,2,3,4} where 1-st, 2-nd, 3-rd and 4-th propagators are present in this sector. n is the total number of Baikov variables in a family. Actually, this serves as a limited case of dlog construction. The output will be a list of forms like {{{pole,poly,path},...},polylist1,polylist2,sector}. poles are classified by their sectors.";
PolesAnalyze::misp = "We need consider at least next-to-next-to-maximal cut for this sector: `1` and only `2` of the propagators can be set to 0.";


RIntersection::usage="RIntersection[list1,list2] calculates the intersection of two lists where terms with different signs are taken as the same.";
AllRuleQ::usage="AllRuleQ[list] returns True if all elements in the list are rules.";


ExtractPoleInfo::usage="ExtractPoleInfo[exp,n] extracts the information of pole from the output of PolesAnalyze[] function. n is the total number of Baikov variables in an integral family. There is an option OutputLevel, it is set to 1 by default. When it is 2, the output will include some second type of poles.";

CompatiblePoleQ::usage="CompatiblePoleQ[p1,p2] checks whether two poles p1 and p2 are compatible.";
MaxPole::usage="MaxPole[p1,p2] compares two poles. If they don't belong to each other then it will return 0. If p1 contains p2, it will return 1 and otherwise it will return 2.";
RemoveSPole::usage="RemoveSPole[poles] removes those poles that is already contained in other poles.";
RemoveSPoleBeta::usage="Enhanced version of RemoveSPole.";
CompatiblePoleGraph::usage="CompatiblePoleGraph[poles] constructs a graph for compatible poles and find all cliques in this graph. It will return {graph,cliques}";
MergePoles::usage="MergePoles[poleslist] merges a list of compatible poles, poleslist.";

ExtractSquareRoots::usage="ExtractSquareRoots[exp] extracts the leading singularities which are square roots. These are the possible square roots appearing in UT basis. exp is the output of PolesAnalyze[]";


FindPoleMaps::uage="FindPoleMaps[poles,var,sector] finds all possible pole maps that can map var to kinematic variables in a given sector. poles is the output of ExtractPoleInfo[]. var is a list of Baikov variables involved in a construction of algebraic letter. The sector information is used to set the constraints that all propagators in this sector should be cut.";


ConstructLetter::usage="ConstructLetter[G1,G2,type] gives the form of letter of two compatible Gram. When type=1, it is 1/Sqrt[G1G2], when type=2, it is 1/(G1G2).";
ConstructLetter::warning="This two grams are same with each other! `1`";
ConstructLetter::err="The following case has failed: `1`.";
ConstructLetter::fail="Fail to find the relations between these two list: `1`";

ConstructFromGram::usage="ConstructFromGram[gramlist] constructs all possible letters from the compatible Gram list.";

FindGram::usage="FindGram[{gram},glist,krep] finds the Gram in the glist, note that two grams can be equivalent to each other if they evaluate to the same polynomial.";

LetterPathRefined::usage="LetterPathRefined[result,letters,krep] gives a refined path for every letter, so that only few pole with the same path can be substituted into the letter.";

ExtractAlgLetter::usage="ExtractAlgLetter[result,topsector] constructs all possible letters from the result which is the output of AllSectorBaikovMat[]. topsector is the top sector specified. All subsectors of it will be analyzed.";


AllRationalLetters::usage="AllRationalLetters[polestructure] gives all rational letters for an integrall family where polestructure is the output of PolesAnalyze[]. It has an option OutputLevel which is set to 1 by default. When it is set to another value, the output will include some spurious expressions.";

SpecialSimplify::usage="SpecialSimplify[{gr,ga}] simplifies expressions like (gr+Sqrt[ga])/(gr-Sqrt[ga]) and returns simplified {gr,ga}";
ApplyPoleToLetter::usage="ApplyPoleToLetter[pole,letter,krep] apply one pole to the letter constructed";
ApplyPoleToGram::usage="ApplyPoleToGram[pole,gram,krep] apply one pole to the gram matrix";

PerfectSquareQ::usage="PerfectSquareQ[monomial] find whether monomial given is a perfect square of algebraic variables.";

AdmissiblePoleQ::usage="AdmissiblePoleQ[glist,pole,permsq,krep] decides whether a pole can be substituted into glist. permsq is a list of polynomials the square root of which are leading singularities appearing in the construction of UT basis.";
PathsDistance::usage="PathsDistance[p1,p2] calculate the minimum distance of two path list.";


SelectAlgLetter::usage="SelectAlgLetter[alg,rl] selects algebraic letters in the list alg according to whether this algebraic letter is related to some rational letter in list rl.";


GramPCQ::err="There must be two grams in the arguments";
GramPCQ::usage="GramPCQ[glist] finds whether two grams in glist are parent Gram and child Gram. For example, GramPCQ[{G[{1,3,4}],G[{3,4}]}] will be True";


ApplyPolesToAlgLetter1::usage="ApplyPolesToAlgLetter1[poles,letters,reform,krep] apply all possible poles to the letters constructed. reform is the output of ReformRep[result]. This will give all first type letters. The output consists of three part: nletters (which is the normal letter), kletters (which consists of only pure kinematics grams) and sletters (which are not permissible because the distance between the selection rule we put on)";


GetQuadMatrix::usage="GetQuadMatrix[exp] get the quadratic matrix from a one-loop baikov polynomial.";
ApplyPoleToQM::usage="ApplyPoleToQM[poles,G,krep] apply two poles to one Q matrix. G is the Gram G[]. The ouput will be {{Log[...],G,{pole1,pole2}},...}.";
ApplyPolesToAlgLetter2::usage="ApplyPolesToAlgLetter2[poles,gram,reform,krep] apply poles to one Gram to construct the second type of letter. When gram={}, reform will give all possible gram. The output will be in the form {{{Log[],Gram1,{pole1,pole2}},{Log[],Gram1,{pole3,pole4}},...},{{Log[],Gram2,{pole1,pole2}},{Log[],Gram2,{pole3,pole4}},...},...}";


AllAlgLetters::usage="AllAlgLetters[poles,algletter,reform,krep] get all possible algbraic letter for a family. algletter is the output of ExtractAlgLetter[]. reform is the output of ReformRep[result]";
AllAlgLettersPL::usage="The parallel version of AllAlgLetters[poles,algletter,reform,krep]";


SameAlgLetterQ::usage="SameAlgLetterQ[a1,a2] determines whether a1 and a2 are two equivalent algebraic letters.";
GetCoord::usage="GetCoord[exp] gives {sq,r} from an algebraic letter (b+Sqrt[sq])/(b-Sqrt[sq]). r= b^2-sq.";
DeleteSameAlgLetter::usage="DeleteSameAlgLetter[alglist] delete the same algletter in a list.";


SearchIndepLetterNum::usage="SearchIndepLetterNum[letter,repnum] finds the independent letter with all its variables have been replaced by random num.";
GetAllAlgIndepLetter::usage="GetAllAlgIndepLetter[alglist] gets all independent algebraic letters in alglist. The result is classified by the square roots appearing in the letters.";


LetterInfo::usage="LetterInfo[letter,algresult,polestructure] show information about all this letter is constructed and where it comes from.";


FindLetterLinearRelation::usage="FindLetterLinearRelation[basis,target] expand the target letter on a set of letter basis.";


FindGramFromPoly::usage="FindGramFromPoly[poly,poles,reform,krep] finds the poles and Grams corresponding to the poly.";


RegularizeSquareRoots::usage="RegularizeSquareRoots[list] makes the square roots in the algebraic letter appearing in a canonical form which is defined by mathematica function Factor[]. This avoids the same square root being taken as different ones due to their different forms caused by the simplification process.";


CheckPosition::usage="CheckPosition[pos] checks the positions pos of one variable appearing in the Gram matrix are whether in one row and column or not."


Begin["`Private`"]


(*some primary functions*)
dlogForm[g1_,g2_]:=Log[(g1+Sqrt[g2])/(g1-Sqrt[g2])];
Lambda[x_,y_,z_]:=x^2+y^2+z^2-2x*y-2y*z-2x*z;


AllRuleQ[list_]:=If[list==={},True,AllTrue[list,MatchQ[#,Rule[_,_]]&]];


AbsGram[exp_]:=If[MatchQ[exp,-G[__]],Return[-exp],Return[exp]];

Shrink[list_]:=Module[{un},
If[Length[list]==1,Return[{{list[[1,1]],list[[1,2]],{list[[1,3]]}}}]];
un=GatherBy[list,IntegerQ[(#[[2]]/.{\[Epsilon]->0})]&];
Return[Table[{AbsGram[un[[i,1,1]]],un[[i,1,2]],un[[i,All,3]]},{i,1,Length[un]}]];
];

Options[ReformRep]={ZeroSector->False};
ReformRep[result_,OptionsPattern[]]:=Module[{rep,glist},
If[OptionValue[ZeroSector],
	rep=result//Flatten[#,1]&,
	rep=result//Flatten[#,1]&//DeleteCases[#,{_,{_,0}}]&
];
glist=Table[Append[#,rep[[i,1]]]&/@rep[[i,2,1]],{i,1,Length[rep]}]//Flatten[#,1]&//GatherBy[#,AbsGram[First[#]]&]&;
Return[MapAt[AbsGram,(Shrink/@glist)//Flatten[#,1]&,{All,1}]];
];

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
If[(Det[tem/.numrep]-Det[tem1/.numrep])=!=0,Throw[False]]
,{i,1,3}];
Throw[True];
];


Options[CompatibleGramQ]={deBug->False};
CompatibleGramQ[Gram1_,Gram2_,OptionsPattern[]]:=Catch@Module[{G1,G2,l1,l2,sol,var1,var2,mat,mat1,t,t1,sys,s1,s2,tem,ntem},
If[Head[Gram1]===G,G1=Gram1,If[MatchQ[Gram1,-G[__]],G1=-Gram1,Message[CompatibleGramQ::err];Throw[$Failed]]];
If[Head[Gram2]===G,G2=Gram2,If[MatchQ[Gram2,-G[__]],G2=-Gram2,Message[CompatibleGramQ::err];Throw[$Failed]]];
(*when two Grams are the same, we just skip this case*)
If[(G1[[1]]//Sort)===(G2[[1]]//Sort),Throw[0]];
l1=Length[G1[[1]]];
l2=Length[G2[[1]]];
(*only two kind of Grams are compatible: G(q,q1) and G(q,q2); G(q) and G(q,q1,q2); q is a list of momenta*)
If[Abs[l1-l2]!=0&&Abs[l1-l2]!=2,Throw[0]];
If[Abs[l1-l2]==0,
	var1=Union[Variables[G1[[1]]],Variables[G2[[1]]]];
	s1=G1[[1]];s2=G2[[1]];
	mat=Table[t[i,j],{i,1,Length[s2]},{j,1,Length[s1]}];
	sys=Thread@Equal[(Coefficient[mat . s1-s2,#]&/@var1)//Flatten,0];
	sol=FindInstance[sys,mat//Flatten,Integers];
	If[OptionValue[deBug],Print["sol full length: ",sol]];
	If[sol=!={},Throw[0]];(*In this case, these two grams are the same with each other*)
	(*We require there are at least l-1 common entries, so will solve the remaining entries to see whether there are so many relations*)
	tem=Intersection[G1[[1]],G2[[1]]];
	If[Length[tem]==(l1-1),Throw[1]];
	s1=Complement[G1[[1]],tem];
	s2=Complement[G2[[1]],tem];
	If[OptionValue[deBug],Print["s1,s2: ",{s1,s2}]];
	mat=Table[t[j],(*{i,1,l1-1},*){j,1,Length[s1]}];
	mat1=Table[t1[j],(*{i,1,l1-1},*){j,1,Length[s2]}];
	sys=(Coefficient[(Coefficient[mat . s1-mat1 . s2,#]&/@var1)//Flatten,#]&/@Join[mat,mat1])//Transpose;
	sol=NullSpace[sys];
	If[OptionValue[deBug],Print["sys: ",sys];Print["sol: ",sol]];
	If[sol==={},Throw[0]];
	sol=Table[Thread@Rule[Join[mat,mat1],sol[[i]]],{i,1,Length[sol]}];
	ntem=mat/.sol;
	If[Length[sol]==(l1-1-Length[tem]),
		If[OptionValue[deBug],Print["ntem: ",ntem];Print["s1,s2: ",{s1,s2}]];
		(*we can constrain the possible recombination of momenta, the assignment of loop momentum is not arbitrary, but there are no reason we should do this*)
		(*If[AnyTrue[ntem,(Length[DeleteCases[#,0]]>2)&],Throw[0],Throw[1]]*)Throw[1],
		Throw[0]
	],
	(*Case when the length is not equal to each other*)
	If[l1>l2,s1=G1[[1]];s2=G2[[1]],s1=G2[[1]];s2=G1[[1]]];
	var1=Variables[s1];
	var2=Variables[s2];
	If[OptionValue[deBug],Print["case length not equal. {s1,s2} ",{s1,s2}]];
	If[ContainsAll[var1,var2],
		mat=Table[t[i,j],{i,1,Length[s2]},{j,1,Length[s1]}];
		sys=Thread@Equal[(Coefficient[mat . s1-s2,#]&/@var1)//Flatten,0];
		sol=FindInstance[sys,mat//Flatten,Integers];
		If[sol==={},
			(*Message[CompatibleGramQ::warning,{s1,s2}];*)Throw[0],
			tem=mat/.sol[[1]];
			If[OptionValue[deBug],Print["tem: ",tem];Print["s1,s2: ",{s1,s2}]];
			(*Do[If[Length[DeleteCases[tem[[j]],0]]>2,Throw[0]],{j,1,Length[tem]}];*)(*we can constrain the possible recombination of momenta, the assignment of loop momentum is not arbitrary, but there are no reason we should do this*)
			If[MatrixRank[tem]<Min[l2,l1],Throw[0],Throw[2]]
		],
		Throw[0]
	];
];
];


Options[PathsDistance]={PathRelax->False};
PathsDistance[p1i_,p2i_,OptionsPattern[]]:=Module[{p1,p2,dis,min,pos,un},
If[p1i==={},p1={{}},If[Head[p1i[[1]]]=!=List,p1={p1i},p1=p1i]];
If[p2i==={},p2={{}},If[Head[p2i[[1]]]=!=List,p2={p2i},p2=p2i]];
(*we require the input in the form: {{path 1},{path 2},...}*)
dis=Reap[
	Do[
		Sow[(Length[Union[p1[[i]],#]]-Length[Intersection[p1[[i]],#]])&/@p2]
	,{i,1,Length[p1]}]
	];
min=Min[Union@@(dis[[2,1]])];
If[!OptionValue[PathRelax],Return[min]];
If[min>2,Return[min]];
If[min==2,
	pos=Position[dis[[2,1]],2];
	Do[
		un=Union[p1[[pos[[j,1]]]],p2[[pos[[j,2]]]]];
		If[Complement[un,p1[[pos[[j,1]]]]]==={}||Complement[un,p2[[pos[[j,2]]]]]==={},min=min-1;Break[]]
	,{j,1,Length[pos]}]
];
Return[min];
];


Options[CompatiblePathQ]={TopOnly->False,SubOnly->False};
CompatiblePathQ[pathl_,path_,OptionsPattern[]]:=Catch@Module[{},
If[path=!={}&&Head[path[[1]]]===List,
	(*in case path is a list of paths like {{1,2},{1,2,4}}*)
	Do[
		If[OptionValue[TopOnly],
			(*pathl can be a parent of path*)
			If[AnyTrue[path,ContainsAll[#,pathl[[i]]]&],Throw[True],Continue[]],
			If[OptionValue[SubOnly],
				(*pathl can be a child of path*)
				If[AnyTrue[path,ContainsAll[pathl[[i]],#]&],Throw[True],Continue[]],
				(*either case is allowed*)
				If[AnyTrue[path,ContainsAll[#,pathl[[i]]]&]||AnyTrue[path,ContainsAll[pathl[[i]],#]&],Throw[True],Continue[]]
			]
		]
	,{i,1,Length[pathl]}],
	(*in case path is a single path like {1,2}*)
	Do[
		If[OptionValue[TopOnly],
			(*pathl can be a parent of path*)
			If[ContainsAll[path,pathl[[i]]],Throw[True],Continue[]],
			If[OptionValue[SubOnly],
				(*pathl can be a child of path*)
				If[ContainsAll[pathl[[i]],path],Throw[True],Continue[]],
				(*either case is ok*)
				If[ContainsAll[path,pathl[[i]]]||ContainsAll[pathl[[i]],path],Throw[True],Continue[]]
			]
		]
	,{i,1,Length[pathl]}]];
	Throw[False];
];


AllSubRep[reform_,rep_]:=Module[{toppath,l,glist={}},
toppath=rep[[All,1]];
l=Length[toppath];
Do[
If[AnyTrue[reform[[i,3]],(Or@@Table[ContainsAll[#,toppath[[j]]],{j,1,l}])&],AppendTo[glist,reform[[i]]]]
,{i,1,Length[reform]}];
Return[glist];
];


Options[ExtractRelevantGram]={TopOnly->True,deBug->False};
ExtractRelevantGram[result_,sector_,n_,OptionsPattern[]]:=Module[{reform,flag,rep,glist,tem,cmGram,cutGram,eq,xl,x1,sol,sol1,pow,disc,fac},
rep=GetBaikovMatRep[result,sector,n]//Simplify;
reform=ReformRep[result];(*Note that the structure of reform is {{g,c,{paths}},{g,c,{paths}},...}*)
If[OptionValue[TopOnly],
glist=Table[Append[#,rep[[i,1]]]&/@rep[[i,2,1]],{i,1,Length[rep]}]//Flatten[#,1]&(*//GatherBy[#,IntegerQ[(#[[2]]/.{\[Epsilon]->0}//Factor)]&]&*),
glist=AllSubRep[reform,rep]
];(*classify all the Gram determinant.Note that the structure of glist is {{g,c,path},{g,c,path},...}*)
(*Now for each Gram, we find their compatible Gram in the whole representation*)
If[OptionValue[deBug],Print["glist: ",glist]];
cmGram=Table[
tem={glist[[i]]};
Do[
If[IntegerQ[(reform[[j,2]]-glist[[i,2]])/.{\[Epsilon]->0}],
(*If[OptionValue[deBug],Print["reform: ",reform[[j]]]];*)
If[CompatiblePathQ[reform[[j,3]],glist[[i,3]],TopOnly->True],
If[OptionValue[deBug],Print["Gram: ",reform[[j]]]];
If[CompatibleGramQ[reform[[j,1]],glist[[i,1]]]>0,
AppendTo[tem,reform[[j]]],
Continue[]]
,Continue[]]
,Continue[]]
,{j,1,Length[reform]}];
tem
,{i,1,Length[glist]}];(*get the compatible Gram for every Gram appearing in the list*)
Return[cmGram];
];

Options[ExtractRelevantGramSim]={TopOnly->True,deBug->False,KineVar->{},KineRestrict->False,LoopRestrict->False,OneLoop->False,LEN->2};
ExtractRelevantGramSim[rep_,reform_,n_,OptionsPattern[]]:=Module[{flag,glist,tem,cmGram,cutGram,eq,xl,x1,sol,sol1,pow,disc,fac},
If[OptionValue[TopOnly],
	glist=Table[Append[#,rep[[i,1]]]&/@rep[[i,2,1]],{i,1,Length[rep]}]//Flatten[#,1]&(*//GatherBy[#,IntegerQ[(#[[2]]/.{\[Epsilon]->0}//Factor)]&]&*),
	glist=AllSubRep[reform,rep]
];(*classify all the Gram determinant.Note that the structure of glist is {{g,c,path},{g,c,path},...}*)
(*Now for each Gram, we find their compatible Gram in the whole representation*)
If[OptionValue[deBug],Print["glist: ",glist]];
cmGram=Table[
		tem={glist[[i]]};
		Do[
			If[IntegerQ[(reform[[j,2]]-glist[[i,2]])/.{\[Epsilon]->0}],(*two grams are relevant only when their powers are all integer or half integer*)
				(*If[OptionValue[deBug],Print["reform: ",reform[[j]]]];*)
				If[CompatiblePathQ[reform[[j,3]],glist[[i,3]],TopOnly->True],
					(*note here we put on the constraint that the Grams relevant should not be too far away from the current Gram*)
					(*Print["var: ",Variables[Join[Abs[glist[[i,1]]],reform[[j,1,1]]]]];*)
					If[ContainsAll[OptionValue[KineVar],Variables[Join[AbsGram[glist[[i,1]]][[1]],reform[[j,1,1]]]]],
						(*the case when two Grams only involve external kinematics*)
						(*if we put constraints on the paths of kinematic Grams, the we remove those too far away from each other*)
						If[OptionValue[KineRestrict],If[PathsDistance[reform[[j,3]],glist[[i,3]]]>OptionValue[LEN],Continue[]]];
						If[OptionValue[OneLoop],If[Abs[Length[AbsGram[glist[[i,1]]][[1]]]-Length[reform[[j,1,1]]]]==2&&IntegerQ[glist[[i,2]]/.{\[Epsilon]->0}],(*Print["glist, Gram: ",{glist[[i,1]],reform[[j]]}];*)Continue[]]];
						If[CompatibleGramQ[reform[[j,1]],glist[[i,1]]]>0,AppendTo[tem,reform[[j]]],Continue[]]
						,
						If[OptionValue[LoopRestrict],If[PathsDistance[reform[[j,3]],glist[[i,3]]]>OptionValue[LEN],Continue[]]];
						If[CompatibleGramQ[reform[[j,1]],glist[[i,1]]]>0,AppendTo[tem,reform[[j]]],Continue[]]
					],
					Continue[]
				],
				Continue[]
			]
		,{j,1,Length[reform]}];
		(*now we add the special case where the other Gram is G[{},{}]=1. This is always possible as long as glist[[i]] is a length-2 Gram*)
		If[Length[glist[[i,1,1]]]==2,AppendTo[tem,{G[{},{}],glist[[i,2]],{glist[[i,3]]}}]];
		tem
		,{i,1,Length[glist]}];
(*get the compatible Gram for every Gram appearing in the list*)
Return[cmGram];
];

Options[ExtractRelevantGramSingle]={deBug->False};
ExtractRelevantGramSingle[g_,reform_,OptionsPattern[]]:=Module[{tem={g}},
(*If[OptionValue[deBug],Print["glist: ",glist]];*)
Do[
If[IntegerQ[(reform[[j,2]]-g[[2]])/.{\[Epsilon]->0}],
(*If[OptionValue[deBug],Print["reform: ",reform[[j]]]];*)
If[CompatiblePathQ[reform[[j,3]],g[[3]]],
If[OptionValue[deBug],Print["Gram: ",reform[[j]]]];
If[CompatibleGramQ[reform[[j,1]],g[[1]]]>0,
AppendTo[tem,reform[[j]]],
Continue[]]
,Continue[]]
,Continue[]]
,{j,1,Length[reform]}];
(*get the compatible Gram for every Gram appearing in the list*)
Return[tem];
];


FindGram[Gram_,glist_,krep_]:=Module[{gmat,bench,var,temmat,temvar,check,temc,rep,pos},
gmat=Gram/.{G[a_,a_]:>GramMat[a,a,krep]};
pos=Table[{},{i,1,Length[gmat]}];
var=gmat//Variables//Sort;
rep=Table[Thread@Rule[var,RandomPrime[{100,1000},Length[var]]],{i,1,3}];
bench=Table[Abs[Det[#]]&/@(gmat/.rep[[i]]),{i,1,3}]//Transpose;
temmat=glist/.{G[a_,a_]:>GramMat[a,a,krep]};
Do[
	temvar=temmat[[i]]//Variables;
	(*If[FreeQ[temvar,x],Continue[]];*)
	temc=Table[Abs[Det[temmat[[i]]/.rep[[j]]]],{j,1,3}];
	Do[
		check=bench[[j]]-temc;
		If[check==={0,0,0},AppendTo[pos[[j]],i]]
	,{j,1,Length[bench]}]
,{i,1,Length[temmat]}];
(*If[!FreeQ[pos,{}],Message[FindGram::warning,Gram]];*)
Return[pos];
];


Sector2Digits[sector_]:=Sum[Power[2,sector[[i]]-1],{i,1,Length[sector]}];
Digits2Sector[digits_]:=Position[IntegerDigits[digits,2]//Reverse,1]//Flatten;


PerfectSquareQ[mon_]:=Module[{fl,tem,ntem},
fl=FactorList[mon]//DeleteCases[#,_?(NumericQ[First[#]]&)]&;
If[fl==={},Return[True]];
If[AllTrue[fl[[All,2]],EvenQ],Return[True],Return[False]];
];


PerfectSquareSplit[expr_]:=Module[{exp,list,sl={},nsl={}},
exp=expr//Factor;
If[Head[exp]===Power,If[MatchQ[exp,Power[_,_?EvenQ]],Return[{{exp/.{Power[z_,n_]:>z}},{}}],Return[{{exp/.{Power[z_,n_]:>Power[z,Quotient[n,2]]}},{exp/.{Power[z_,n_]:>z}}}]]];
If[Head[exp]===Times,list=List@@exp,Return[{{},{exp}}]];
Do[
	If[MatchQ[list[[i]],Power[_,_?EvenQ]],AppendTo[sl,list[[i]]/.{Power[z_,n_]:>z}],AppendTo[nsl,list[[i]]]]
,{i,1,Length[list]}];
Return[{sl,nsl}];(*sl is the perfect square after removing power, nsl is non-perfect square and it will remains in the square root*)
];


ExistRelationQ[p1_,sqpl_]:=Catch@Module[{var,pow,pos,sol,tem},
If[sqpl==={},Throw[False]];
var=Intersection[Variables[p1],Variables[sqpl]];
pow=Exponent[p1,#]&/@var;
If[FreeQ[pow,1],Throw[False]];(*find whether there is a variable where it is linear in p1*)
(*We haven't exclude all the possibility that there may exist relations between the polynomials since there may exists a variable transformation to make it possible, but in most cases, this will work*)
pos=FirstPosition[pow,1][[1]];
sol=Solve[p1==0,var[[pos]]][[1]];
tem=sqpl/.sol//Factor//Numerator;
If[AnyTrue[tem,PerfectSquareQ],Throw[True],Throw[False]];
];


RemoveCoeff[exp_]:=If[exp===0,0,Times@@Power@@@(FactorList[exp]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)];


TakeSquareRoot[exp_]:=Module[{list,a,s,id,pos},
	list=FactorList[exp];
	a=Product[Power[list[[i,1]],Quotient[list[[i,2]],2]],{i,1,Length[list]}];
	s=Product[Power[list[[i,1]],Mod[list[[i,2]],2]],{i,1,Length[list]}];
	If[s=!=1,
		If[!FreeQ[s,x],Return[0]];
		If[!MemberQ[Rlist,s],
		Unprotect[R,Rlist];
		AppendTo[Rlist,s];
		id=Length[Rlist];
		R/:R[id]*R[id]:=Rlist[[-1]];
		R/:Power[R[id],n_/;(EvenQ[n])]:=Power[Rlist[[-1]],n/2];
		R/:Power[R[id],n_/;(OddQ[n])]:=Power[Rlist[[-1]],Quotient[n,2]]*R[id];
		Protect[R,Rlist];
		,
		id=Position[Rlist,s][[1,1]]
		];
		s=R[id];
	];
	Return[a*s];
];


EnhancedFactor[intset_,z_]:=Module[{tem,result,pow,sol,sq,rep},
	result=Reap[
	Do[
		If[FreeQ[intset[[i]],z],
			Sow[intset[[i]]];Continue[],
			pow=Exponent[intset[[i]],z];
			If[pow==1||pow>2,
				Sow[intset[[i]]];Continue[],
				sol=Solve[intset[[i]]==0,z];
				sq=Cases[sol,Power[_,1/2],Infinity]//DeleteDuplicates;
				If[sq==={},Sow[(Power@@@(FactorList[intset[[i]]]))//DeleteCases[#,_?NumericQ]&];Continue[]];
				tem=TakeSquareRoot/@(sq/.{Power[a_,1/2]:>a});
				If[tem==={0},Sow[intset[[i]]];Continue[]];(*in this case, there are still Baikov variables under square roots*)
				rep=Thread@Rule[sq,tem];
				tem=z-(z/.(sol/.rep))//Factor//Numerator;
				Do[Sow[tem[[j]]],{j,1,Length[tem]}];
			]
		]
	,{i,1,Length[intset]}]
	][[2]];
	If[result=!={},Return[result[[1]]//Flatten],Return[result]];
];


Options[ResolveSingleVariable]={deBug->False,LastVar->False,HigherPower->False,DLog->False,AdInfo->True,SelectQ->True};
ResolveSingleVariable[iintset_,ihintset_,z_,OptionsPattern[]]:=Module[{intset,hintset,p,tem,tem1,tem2,coe,irintset={},rintset={},irhintset={},rhintset={},result={},sresult={},sol,hpintset={},dlog={},a1,b1},
intset=iintset;hintset=ihintset;
If[intset==={}&&hintset==={},Message[ResolveSingleVariable::err];Return[$Failed]];

(*some repeated terms in hintset should be moved to intset*)
tem=Times@@hintset;
tem1=PerfectSquareSplit[tem];
intset=Join[intset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates;
hintset=tem1[[2]]//DeleteCases[#,_?NumericQ]&;
intset=(FactorList[#][[All,1]]&/@intset)//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates;(*keep intset to be irreducible polynomials*)
(*Print["intset: ",intset];Print["hintset: ",hintset];*)

(*split the set into relevant set and irrelevant set w.r.t variable z*);
If[intset=!={},
	tem=intset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,z]&]&;
	If[!FreeQ[tem[[1]],z],rintset=tem[[1]];irintset=Drop[tem,1]//Flatten[#,1]&,irintset=tem[[1]];rintset=Drop[tem,1]//Flatten[#,1]&],
	(*if there cannot be rational polynomials in the denominator*)
	rintset={};irintset={}
];
If[hintset=!={},
	tem=hintset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,z]&]&;
	If[!FreeQ[tem[[1]],z],rhintset=tem[[1]];irhintset=Drop[tem,1]//Flatten[#,1]&,irhintset=tem[[1]];rhintset=Drop[tem,1]//Flatten[#,1]&],
	(*if there are no square roots in the denominator*)
	rhintset={};irhintset={}
];
(*irintset and irhintset will be polynomial free of variable z, they will be kept to the resolution of next variable*)
If[OptionValue[deBug],Print["rintset,rhintset,irintset,irhintset: ",{rintset,rhintset,irintset,irhintset}]];

(*we will include the 'intersection' from two polynomials under square root in elliptic case*)
If[OptionValue[LastVar]&&!OptionValue[DLog],If[Length[rhintset]>1&&Exponent[Times@@rhintset,z]>2,Table[irhintset=Join[irhintset,{Resultant[rhintset[[i]],rhintset[[j]],z]//Factor}];,{i,1,Length[rhintset]},{j,i+1,Length[rhintset]}]]];

(*resolve singularities relevant to square roots*)
If[rhintset=!={},
	tem=Times@@rhintset;(*expression under square root*)
	tem1=Exponent[tem,z];
	If[OddQ[tem1],(*if polynomial under square root is of odd power*)
		If[tem1>1,
			Message[ResolveSingleVariable::warning,tem1];(*polynomial under square root of odd power higher than 2*)
			If[OptionValue[LastVar],Message[ResolveSingularities::elliptic];Print["tem: ",tem]];(*when there is only one integration variable, it is elliptic or hyperelliptic*)
			coe=0,(*we don't consider this case since it is hardly related to MPL*)
			coe=Prime[2024](*when this is a linear polynomial it can be related to dlog like Sqrt[b*a-c]/(z-a)/Sqrt[b*z-c], but we need to keep in mind that Sqrt[b]/Sqrt[b*z-c] is not a dlog*)
			(*here we use coe=Prime[2024] which is chosen randomly to distinguish it from the following case Sqrt[a]/Sqrt[az^2+bz+c] which is a dlog form*)
		],
		(*even power case*)
		If[tem1>2,
			If[OptionValue[LastVar],Message[ResolveSingularities::elliptic];Print["tem: ",tem]];(*when there is only one integration variable, it is elliptic or hyperelliptic*)
			coe=0,
			coe=Coefficient[tem,z,tem1](*in this case, there is a quadratic polynomial under square root*)
		]
	],
	(*if there are no square roots in the integrand.*)
	coe=Prime[2024](*here we use coe=Prime[2024] which is chosen randomly, it should be understood as an identification code*)
];

(*one type of dlog constructed with only one square root in the denominator Sqrt[a]/Sqrt[a*z^2+bz+c]*)
If[coe=!=0,
	If[coe===Prime[2024],
		coe=1(*redefine coe to be 1*),
		(*the following is the case Sqrt[a]/Sqrt[a*z^2+bz+c]*)
		tem1=PerfectSquareSplit[coe];
		If[OptionValue[SelectQ],
			tem1={{},tem1[[2]]};(*--------------------(^_^)--------------------*)
			(*here we remove those terms from perfect square roots under square root which we believe to be spurious letters*)
		];
		AppendTo[sresult,{{Sqrt[Times@@rhintset]},Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}];
		AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}];
	];
];

If[OptionValue[deBug],Print["coe: ",coe];Print["result: ",result]];
(*If[coe=!=0&&coe=!=1,
	tem1=PerfectSquareSplit[coe];
	(*AppendTo[result,{{z->Infinity},Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}];*)
	If[Length[rintset]>0,
		tem=DeleteCases[rintset,_?(Exponent[#,z]>1&)];
		Do[
			sol=(Solve[tem[[k]]==0,z]//DeleteDuplicates);
			AppendTo[result,{sol[[1]],Join[irintset,((FactorList[#][[All,1]]&/@(tem/.sol[[1]]//Factor//Numerator))//Flatten),tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}]
		,{k,1,Length[tem]}];
	];
];*)

(*resolve singularities from rational part*)
rintset=EnhancedFactor[rintset,z];(*Enhanced factor can factor some quadratic polynomials and abbreviate some square roots of kinematic variables*)
p=Exponent[#,z]&/@rintset;(*power of rational polynomials w.r.t variable z*)
sol=Table[If[p[[i]]<3,(Solve[rintset[[i]]==0,z]//DeleteDuplicates),Null],{i,1,Length[rintset]}];(*only solve linear polynomials and quadratic polynomials*)
tem2={};
Do[
	(*when there are polynomials with power of z higher than 2, we don't consider them since they may generate singularities with power other than 1/2*)
	If[p[[i]]>2,AppendTo[hpintset,rintset[[i]]];Continue[]];
	
	If[Length[sol[[i]]]==1,(*if this pole is a linear rational pole, we will keep this pole*)
		If[rhintset==={},
			(*if there are no square roots*)
			(*take pole at this solution*)
			If[!OptionValue[DLog],
				AppendTo[result,{sol[[i,1]],Join[irintset,(*{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)},*)(RemoveCoeff/@(NumeratorDenominator[(rintset/.sol[[i,1]]//Factor)]//Flatten))]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,irhintset}],
				(*dlog type construction*)
				AppendTo[dlog,{Join[irintset,{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)}],irhintset}];(*a/(a*z-b) type dlog*)
				Table[If[k!=i,AppendTo[dlog,{Join[irintset,{RemoveCoeff[Times@@NumeratorDenominator[rintset[[k]]/.sol[[i,1]]//Factor]]}],irhintset}]],{k,1,Length[rintset]}](*(bc-ad)/(az-b)/(cz-d) type dlog*)
			]
			,
			(*if there are square roots*)
			If[OptionValue[LastVar]&&Exponent[Times@@rhintset,z]>2,Continue[]];(*when it is the last variable and the power is higher than 2 in square root, then this is elliptic or hyperelliptic function*)
			tem=Times@@@NumeratorDenominator[rhintset/.sol[[i,1]]//Factor];(*move denominator to numerator*)
			(*when the polynomials under square root vanish, then this is a branch cut, not a pole*)
			If[!FreeQ[tem,0],Continue[]];
			If[OptionValue[SelectQ],
				tem=tem//DeleteCases[#,_?PerfectSquareQ]&;(*------------------(^_^)-------------------*)
				(*remove perfect squares under square roots, they are conjectured to be spurious*)
			];
			tem1=PerfectSquareSplit[Times@@tem];
			AppendTo[tem2,{i,Times@@@tem1}];(*for later construction*)
			If[!OptionValue[DLog],
				(*AppendTo[result,{sol[[i,1]],Join[irintset,tem1[[1]]]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}]*)
				AppendTo[result,{sol[[i,1]],Join[irintset,tem1[[1]](*,FactorList[Coefficient[rintset[[i]],z,p[[i]]]][[All,1]]//Flatten*),((((rintset(*//DeleteCases[#,_?(Exponent[#,z]>1&)]&*))/.sol[[i,1]]//Factor)//Numerator))//Flatten]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}],
				(*dlog of type (Sqrt[ad^2+bd+c]/(z-d)/(Sqrt[az^2+bz+c]))*)
				(*Note that when two polynomials are rational, we can always perform the variable transformation, so we can set some poly in rintset at this pole in the same time*)
				AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}]
			]
		],
		(*---------------------------when there is an irreducible quadratic polynomial in the denominator----------------------*)
		If[rhintset=!={}&&!ExistRelationQ[rintset[[i]],rhintset],If[OptionValue[AdInfo],Message[ResolveSingleVariable::additional,{rintset[[i]],rhintset}];Print[Short[{rintset[[i]],rhintset},20]]];Continue[]];
		tem=Discriminant[rintset[[i]],z]*coe//Factor;
		If[tem===0,Continue[]];(*if tem is 0, this construction will stop here*)
		tem1=PerfectSquareSplit[tem];
		If[OptionValue[SelectQ],
			tem1={{},tem1[[2]]};(*---------------------------(^_^)-------------------------------*)
			(*here we remove those terms from perfect squares under square roots which we believe to be spurious letters*)
		];
		If[!OptionValue[DLog],
			AppendTo[sresult,{{rintset[[i]]},Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}],
			(*a special kind of dlog, we slightly enlarge its possible realization.*)
			AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}]
		]
	]
,{i,1,Length[rintset]}];

(*a sepcial kind of dlog, they are actually combinations of above dlogs already constructed, but they may give new varieties*)
If[Length[tem2]>1,
	(*we search for the possible combinations of dlogs, they may give new varieties.*)
	tem=Gather[tem2,(((Last[#1][[2]]-Last[#2][[2]])//Expand)===0||((Last[#1][[2]]+Last[#2][[2]])//Expand)===0)&]//DeleteCases[#,_?(Length[#]<2&)]&;(*we also consider the cases where the square roots are the same except for a minus sign *)
	Do[
		tem1=tem[[k]];(*tem1 will be like {{1,{..,..}},{3,{..,..}}}*)
		Table[
		a1=(rintset[[tem1[[i,1]]]]*tem1[[j,2,1]]-rintset[[tem1[[j,1]]]]*tem1[[i,2,1]])//Factor//Numerator;
		b1=(rintset[[tem1[[i,1]]]]*tem1[[j,2,1]]+rintset[[tem1[[j,1]]]]*tem1[[i,2,1]])//Factor//Numerator;
		If[FreeQ[a1,z],
			AppendTo[result,{sol[[i,1]],Join[irintset,{a1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,{tem1[[i,2,2]]}//Flatten]//DeleteCases[#,_?NumericQ]&}];
			AppendTo[dlog,{Join[irintset,{a1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,{tem1[[i,2,2]]}//Flatten]//DeleteCases[#,_?NumericQ]&}];
		];
		If[FreeQ[b1,z],
			AppendTo[result,{sol[[i,1]],Join[irintset,{b1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,{tem1[[i,2,2]]}//Flatten]//DeleteCases[#,_?NumericQ]&}];
			AppendTo[dlog,{Join[irintset,{b1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,{tem1[[i,2,2]]}//Flatten]//DeleteCases[#,_?NumericQ]&}];
		];
		,{i,1,Length[tem1]},{j,i+1,Length[tem1]}];
	,{k,1,Length[tem]}];
];

(*higher power polynomials may involve higher power singularity, we don't consider them by default*)
If[OptionValue[HigherPower]&&!OptionValue[DLog],
	If[Length[hpintset]>1,Table[AppendTo[sresult,{{hpintset[[i]]},irintset,irhintset,Resultant[hpintset[[i]],hpintset[[j]],z]}];,{i,1,Length[hpintset]},{j,i+1,Length[hpintset]}]];
	Table[AppendTo[sresult,{{hpintset[[i]]},irintset,irhintset,Resultant[hpintset[[i]],D[hpintset[[i]],z],z]}],{i,1,Length[hpintset]}]
];
If[!OptionValue[DLog],
	Return[{result(*normal solution*),sresult(*special solution*)}],
	Return[dlog](*dlog type construction*)
];
];


RemoveSquareRootsSol[sol_]:=Module[{keys,values,den,tem,temq,result},
	result=Reap[
		Do[
			keys=Keys[sol[[i]]];
			values=Values[sol[[i]]];
			den=Quiet[Denominator[values]//.sol[[i]]//Factor];
			If[!FreeQ[den,0],Continue[]];
			tem=keys//.sol[[i]]//Factor;
			temq={Plus@@tem,Times@@tem}//Factor;
			If[FreeQ[temq,R],
				Sow[Thread@Rule[keys,tem]],
				Continue[]
			]
		,{i,1,Length[sol]}]
	][[2]];
	If[result=!={},Return[result[[1]]],Return[{}]];
];


Options[ResolveSingularities]={deBug->False,SortQ->True,AdInfo->True,SelectQ->True,RemoveSquareRoots->True};
ResolveSingularities[iintset_,ihintset_,OptionsPattern[]]:=Module[{intset,hintset,xl,nl,lv={},sq={},sol={},tem,p,tem1,tem2,sol1={},ssol={},result,a1,b1},
intset=iintset;hintset=ihintset;
xl=Cases[{intset,hintset},Subscript[x,_],Infinity]//DeleteDuplicates;(*number of remaining ISPs*)
If[xl==={},Message[ResolveSingularities::err];Return[$Failed]];

(*sort the variable by their appearance in polynomials*)
If[OptionValue[SortQ],
	nl=Table[{Count[FreeQ[#,xl[[i]]]&/@Join[(*intset,*)hintset],False],(Exponent[#,xl[[i]]]&/@intset)//DeleteCases[#,0]&//Min,Count[FreeQ[#,xl[[i]]]&/@intset,True],xl[[i]]},{i,1,Length[xl]}]//SortBy[#,{First,#[[2]]&,#[[3]]&}]&;
	xl=nl[[All,-1]];
	(*sort the variables by their appearing frequences*)
	(*here we set the rule for the order. In principle, we should consider every order. First, we choose those that appear less in the square root with most priority. Then we choose those with less power outside the square root. At last, we choose those appear more outside the square root *)
];
If[OptionValue[deBug],Print["xl: ",xl]];

(*we consider the single ISP case separately*)
If[Length[xl]==1,(*if there is only one ISP remaining*)
	Do[(*this will give all possible pole given by linear and quadratic equations*)
		If[Exponent[intset[[i]],xl[[1]]]<2,
			(*linear pole*)
			AppendTo[sol,Solve[intset[[i]]==0,xl][[1]]];
			AppendTo[lv,Coefficient[intset[[i]],xl[[1]],1]]
			If[OptionValue[deBug],Print["lv: ",lv," intset:",intset[[i]]]],
			(*quadratic pole*)
			AppendTo[ssol,{intset[[i]]}];
			If[hintset=!={}&&!ExistRelationQ[intset[[i]],hintset],
				If[OptionValue[AdInfo],Message[ResolveSingleVariable::additional,{intset[[i]],hintset}];Print[{intset[[i]],hintset}]],
				If[OptionValue[SelectQ],
					tem=PerfectSquareSplit[Discriminant[intset[[i]],xl[[1]]]];
					tem={{},tem[[2]]};(*---------------------(^_^)-------------------------*)
					(*here we remove those terms which arise from the perfect squares under square roots which we conjecture to be spurious*)
					AppendTo[lv,Times@@tem[[2]]];
					AppendTo[sq,Times@@tem[[2]]]
					,
					tem=PerfectSquareSplit[Discriminant[intset[[i]],xl[[1]]]];
					lv=Join[lv,{Times@@tem[[1]],Times@@tem[[2]]}];
					AppendTo[sq,Times@@tem[[2]]]
				];
			];
		]
	,{i,1,Length[intset]}];
	If[hintset==={},
		(*if there are no square roots*)
		lv=Join[lv,((*FactorList[#][[All,1]]&/@*)(Table[RemoveCoeff/@(NumeratorDenominator[intset/.sol[[i]]//Factor]//Flatten//DeleteDuplicates),{i,1,Length[sol]}]//Flatten))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates];
		(*here, it need to stress that (a*b)/(z*(z+a)(z-b)) can not be a dlog, only two linear term can be combined to a dlog form a/(z(z-a)) or ((a-b)/((z-a)(z-b))). but familiar construction can work for multivariate case because we can perform a change of variables*)
		If[OptionValue[deBug],Print["lv: ",lv," intset:",intset]];
		AppendTo[sq,1];(*this corresponds to the case 1/x, that is one single pole of sol[[i]], it will multiply a factor from coefficients*)
		sol1=sol;
		(*AppendTo[sol1,{xl[[1]]->Infinity,xl[[1]]->0}](*in this case Infinity pole is present, not a branch cut*)*),
		
		(*when there are square roots in the denominator*)
		(*here we only concern square roots in dlog, so we won't consider the combinations of dlog any more since this won't give us new square root*)
		tem=Times@@hintset;
		p=Exponent[tem,xl[[1]]];
		If[p>2,
			Message[ResolveSingularities::elliptic];Print["only one ISP! tem: ",tem];
			If[Length[hintset]>1,Table[lv=Join[lv,{Resultant[hintset[[i]],hintset[[j]],xl[[1]]]//Factor}],{i,1,Length[hintset]},{j,i+1,Length[hintset]}]]
		];(*elliptic case*)
		AppendTo[sq,1];
		If[Length[intset]>0,lv=Join[lv,Table[intset/.sol[[i]],{i,1,Length[sol]}]//NumeratorDenominator//Flatten//Factor//DeleteCases[#,_?NumericQ]&//DeleteDuplicates]];
		If[EvenQ[p],
			lv=Join[lv,{Coefficient[tem,xl[[1]],p]}]//DeleteCases[#,0]&//DeleteDuplicates;
			If[OptionValue[deBug],Print["lv: ",lv," tem:",tem]];
			sq=sq*Coefficient[tem,xl[[1]],p]//DeleteCases[#,0]&//DeleteDuplicates;
			(*pay attention that in this case, this square root coefficients should be multiplied to the former ones, because we are not considering a different construction here, it is in one construction process. Then in the last case, we are considering a different construction*)
			sol1={};
			tem2={};
			Do[
				tem1=Times@@@NumeratorDenominator[hintset/.sol[[i]]//Factor];
				(*when the polynomials square root vanish under the outside pole, then it is a branch cut actually*)
				If[FreeQ[tem1,0],
					If[OptionValue[SelectQ],
					tem1=tem1//DeleteCases[#,_?PerfectSquareQ]&(*-----------------(^_^)-------------------*)
					(*here we remove those terms which arising from perfect square under square root, they are conjectured to be spurious letters*)
					];
					tem1=PerfectSquareSplit[Times@@tem1];
					lv=Join[lv,{Times@@tem1[[1]],Times@@tem1[[2]]}];
					If[OptionValue[deBug],Print["lv: ",lv," tem1:",tem1]];
					sq=Join[sq,{Times@@tem1[[2]]}];
					AppendTo[sol1,sol[[i]]];
					AppendTo[tem2,{Numerator[(xl[[1]]-(xl[[1]]/.sol[[i]]))//Factor],Times@@@tem1}];(*for later construction*)
					,
					Continue[]
				]
			,{i,1,Length[sol]}];			
			(*AppendTo[sol1,{xl[[1]]->Infinity,xl[[1]]->0}](*in this case Infinity pole is present, not a branch cut*)*)
			,
			If[p==1&&Length[sol]!=0,(*when there is a linear expression in the square root, Sqrt[ac-b]/(z-c)/Sqrt[az-b] is still a dlog*)
			tem2={};
			Do[
				tem1=Times@@@NumeratorDenominator[hintset/.sol[[i]]//Factor];
				If[FreeQ[tem1,0],
					If[OptionValue[SelectQ],
						tem1=tem1//DeleteCases[#,_?PerfectSquareQ]&(*-----------------(^_^)-------------------*)
						(*here we remove those terms which arising from perfect square under square root, they are conjectured to be spurious letters*)
					];
					tem1=PerfectSquareSplit[Times@@tem1];
					lv=Join[lv,{Times@@tem1[[1]],Times@@tem1[[2]]}];
					If[OptionValue[deBug],Print["lv: ",lv," tem1:",tem1]];
					sq=Join[sq,{Times@@tem1[[2]]}];
					AppendTo[sol1,sol[[i]]];
					AppendTo[tem2,{Numerator[(xl[[1]]-(xl[[1]]/.sol[[i]]))//Factor],Times@@@tem1}];(*for later construction*)
					,
					Continue[]
				]
			,{i,1,Length[sol]}],
			Message[ResolveSingularities::oddpower,tem,{iintset,ihintset}];(*the remaining cases should be elliptic*)
			];
		];
		
		(*a sepcial kind of dlog, they are actually combinations of above dlogs already constructed, but they may give new varieties*)
		If[Length[tem2]>1,
			(*we search for the possible combinations of dlogs, they may give new varieties.*)
			tem=Gather[tem2,(((Last[#1][[2]]-Last[#2][[2]])//Expand)===0||((Last[#1][[2]]+Last[#2][[2]])//Expand)===0)&]//DeleteCases[#,_?(Length[#]<2&)]&;(*we also consider the cases where the square roots are the same except for a minus sign *)
			Do[
				tem1=tem[[k]];(*tem1 will be like {{a*z+b,{..,..}},{c*z+d,{..,..}}}*)
				Table[
					a1=(tem1[[i,1]]*tem1[[j,2,1]]-tem1[[j,1]]*tem1[[i,2,1]])//Factor//Numerator;
					b1=(tem1[[i,1]]*tem1[[j,2,1]]+tem1[[j,1]]*tem1[[i,2,1]])//Factor//Numerator;
					If[FreeQ[a1,xl[[1]]],
						lv=Join[lv,{a1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates;
						If[OptionValue[deBug],Print["lv: ",lv," a1:",a1]];
					];
					If[FreeQ[b1,xl[[1]]],
						lv=Join[lv,{b1}//Flatten//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates;
						If[OptionValue[deBug],Print["lv: ",lv," b1:",b1]];
					];
				,{i,1,Length[tem1]},{j,i+1,Length[tem1]}];
			,{k,1,Length[tem]}];
		];
	],
	
	(*if there are more than one ISPs*)
	tem={{{{},intset,hintset}},{}};
	(*Print["RSS tem: ",tem];*)
	Do[
		(*if one element of tem doesn't contain the variable to be studied at all, then this is not what we want*)
		tem[[1]]=DeleteCases[tem[[1]],_?(FreeQ[#[[{2,3}]],xl[[i]]]&)];
		tem[[2]]=DeleteCases[tem[[2]],_?(FreeQ[#[[{2,3}]],xl[[i]]]&)];
		If[OptionValue[deBug],Print["tem: ",tem]];
		tem1={{},{}};(*structure: {{first type},{second type}}, element:{{},intset,hintset}*)
		Do[
			tem2=ResolveSingleVariable[tem[[1,j,2]],tem[[1,j,3]],xl[[i]],LastVar->If[i==Length[xl],True,False],SelectQ->OptionValue[SelectQ]];
			If[tem2[[1]]=!={},tem1[[1]]=Join[tem1[[1]],MapAt[Join[tem[[1,j,1]],#]&,tem2[[1]],{All,1}]]];
			If[tem2[[2]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[1,j,1]],#]&,tem2[[2]],{All,1}]]]
		,{j,1,Length[tem[[1]]]}];(*Analyze the first part of tem*)
		Do[
			tem2=ResolveSingleVariable[tem[[2,j,2]],tem[[2,j,3]],xl[[i]],LastVar->If[i==Length[xl],True,False],SelectQ->OptionValue[SelectQ]];
			If[tem2[[1]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[2,j,1]],#]&,tem2[[1]],{All,1}]]];
			If[tem2[[2]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[2,j,1]],#]&,tem2[[2]],{All,1}]]]
		,{j,1,Length[tem[[2]]]}];(*Analyze the second part of tem, these will all be classified as second type*)
		tem=tem1;(*renew the original list*)
	,{i,1,Length[xl]}];
	If[tem[[1]]=!={},
		sol1=tem[[1,All,1]];(*form like {{a->1,b->2},{a->2,b->1}}*)
		lv=Join[lv,tem[[1,All,{2,3}]]//Flatten];
		tem2={};
		Do[(*move the abbreviate R[i] from rational letter list to square root list*)
			tem1=Times@@Cases[(FactorList[#][[All,1]]&/@tem[[1,s,2]])//Flatten,R[_]]/.{R[a_]:>Rlist[[a]]};
			AppendTo[tem2,tem1];
		,{s,1,Length[tem[[1]]]}];
		sq=Join[sq,Table[Times@@(Join[tem[[1,s,{3}]],{tem2[[s]]}]//Flatten),{s,1,Length[tem[[1]]]}]];
		If[OptionValue[RemoveSquareRoots],(*remove those solution which will result in square roots in the end*)
			sol1=RemoveSquareRootsSol[sol1];
			lv=DeleteCases[lv,_?(!FreeQ[#,R]&)];
			sq=DeleteCases[sq,_?(!FreeQ[#,R]&)];
		];
	];
	If[tem[[2]]=!={},
		ssol=tem[[2,All,1]];
		lv=Join[lv,tem[[2,All,{2,3}]]//Flatten];
		tem2={};
		Do[(*move the abbreviate R[i] from rational letter list to square root list*)
			tem1=Times@@Cases[(FactorList[#][[All,1]]&/@tem[[2,s,2]])//Flatten,R[_]]/.{R[a_]:>Rlist[[a]]};
			AppendTo[tem2,tem1];
		,{s,1,Length[tem[[2]]]}];
		sq=Join[sq,Table[Times@@(Join[tem[[2,s,{3}]],{tem2[[s]]}]//Flatten),{s,1,Length[tem[[2]]]}]];
		If[OptionValue[RemoveSquareRoots],(*remove those solution which will result in square roots in the end*)
			lv=DeleteCases[lv,_?(!FreeQ[#,R]&)];
			sq=DeleteCases[sq,_?(!FreeQ[#,R]&)];
		];
	];
	(*(*then we analyze the possible square roots in dlog construction*)
	tem={{intset,hintset}};
	Do[
		(*if one element of tem doesn't contain the variable to be studied at all, then this is not what we want*)
		tem=DeleteCases[tem,_?(FreeQ[#,xl[[i]]]&)];
		If[OptionValue[deBug],Print["tem: ",tem]];
		tem1={};(*structure: {dlogs}, element:{intset,hintset}*)
		Do[
			tem2=ResolveSingleVariable[tem[[j,1]],tem[[j,2]],xl[[i]],LastVar->If[i==Length[xl],True,False],DLog->True];
			If[tem2=!={},tem1=Join[tem1,tem2]]
		,{j,1,Length[tem]}];
		tem=tem1;(*renew the original list*)
	,{i,1,Length[xl]}];
	sq=Join[sq,Table[Times@@(tem[[i,{2}]]//Flatten),{i,1,Length[tem]}]];*)
];
Return[{{sol1,ssol},lv//DeleteCases[#,_?NumericQ]&//Factor//DeleteDuplicates,sq//Factor//DeleteDuplicates}];
(*we need to note that here we cannot remove number in sq, since it has a meaning that there are cases where no square roots are present*)
];


Options[Trans2Inf]={deBug->False};
Trans2Inf[intset_,hintset_,xl_,OptionsPattern[]]:=Module[{tem,tem1,nintset={},nhintset={},LVlocal={},sqlocal={}},
tem=(Numerator/@(intset/.Thread@Rule[xl,xl/Subscript[x, 0]]//Factor))/.{Subscript[x, 0]->0}//Factor;
nintset=(FactorList[#][[All,1]]&/@tem)//Flatten//DeleteCases[#,_?NumericQ]&;
If[OptionValue[deBug],Print["nintset: ",nintset]];
tem=(Numerator/@(hintset/.Thread@Rule[xl,xl/Subscript[x, 0]]//Factor))/.{Subscript[x, 0]->0}//Factor;
nhintset={};
Do[
	tem1=PerfectSquareSplit[tem[[k]]];
	nintset=Join[nintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&];
	nhintset=Join[nhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&];
,{k,1,Length[tem]}];
If[OptionValue[deBug],Print["nintset: ",nintset];Print["nhintset: ",nhintset]];
nintset=nintset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,x]&]&;
nhintset=nhintset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,x]&]&;
If[nintset=!={},
	If[!FreeQ[nintset[[1]],x],
		LVlocal=Join[LVlocal,Drop[nintset,1]//Flatten[#,1]&];nintset=nintset[[1]],
		LVlocal=Join[LVlocal,nintset[[1]]];nintset=Drop[nintset,1]//Flatten[#,1]&
	]
];
If[nhintset=!={},
	If[!FreeQ[nhintset[[1]],x],
		LVlocal=Join[LVlocal,Drop[nhintset,1]//Flatten[#,1]&];sqlocal=Join[sqlocal,Drop[nhintset,1]//Flatten[#,1]&];nhintset=nhintset[[1]],
		LVlocal=Join[LVlocal,nhintset[[1]]];sqlocal=Join[sqlocal,nhintset[[1]]];nhintset=Drop[nhintset,1]//Flatten[#,1]&
	]
];
(*This step further splits possible varieties independent of Baikov variable x*)
Return[{nintset//DeleteDuplicates,nhintset,LVlocal//DeleteDuplicates,sqlocal}];
];

CheckValidity[brep_,subset_,krep_]:=Catch@Module[{cut,tem},
cut=Thread@Rule[Subscript[x,#]&/@subset,0];
Do[
	tem=(brep[[j,2,1,All,1]]//Gram2Poly[#,krep]&)/.cut//Factor;
	If[!FreeQ[tem,0],Throw[{}]](*if one in the brep is 0 under cut, then we will return {}. This is designed to keep all the representations equal*)
,{j,1,Length[brep]}];
Throw[brep];
];


Options[FindCover]={deBug->False};
FindCover[list_,OptionsPattern[]]:=Module[{var,subset,flag=0,L=5,tem,temset,pos,result={},k=1},
	var=list//Flatten//Union;
	If[OptionValue[deBug],Print["var: ",var]];
	If[Length[var]<=6,Return[{{var,{}}}]];
	While[flag==0&&k<10,(*we will find a set of subsets which can cover the original set but with length smaller than the total set*)
		subset=Subsets[var,{L}];
		If[OptionValue[deBug],Print["subset: ",subset]];
		tem=Table[Intersection[subset[[i]],#]&/@list,{i,1,Length[subset]}];
		pos=Position[tem,_?(!FreeQ[#,{}]&),1]//DeleteCases[#,{0}]&;
		subset=Delete[subset,pos];(*keep those subsets that has non-zero intersection with each set in the list*)
		If[OptionValue[deBug],Print["subset: ",subset]];
		If[Length[Union[subset//Flatten]]<Length[var],L=L+1;If[L>Length[var],Break[],Continue[]]];(*if remaining sets can not form a cover of original set, then we increase the length of these subsets*)
		
		(*then we pick a collection of sets that can cover the original set*)
		AppendTo[result,{subset[[1]],Thread@Rule[#,Table[Prime[20+i*2],{i,1,Length[#]}]]&[Complement[var,subset[[1]]]]}];
		temset=subset[[1]];
		While[flag==0&&k<10,
			tem=Length[Complement[#,temset]]&/@subset;
			If[Max[tem]==0,flag=1;Break[]];
			pos=PositionLargest[tem][[1]];
			temset=Join[temset,subset[[pos]]]//Union;
			AppendTo[result,{subset[[pos]],Thread@Rule[#,Table[Prime[20+i*2],{i,1,Length[#]}]]&[Complement[var,subset[[pos]]]]}];
			k=k+1;
		];
		k=k+1
	];
	Return[result];
];


Options[PolesAnalyze]={deBug->False,AugAna->True,SelectQ->False,SelectAllQ->False,AdInfo->False,RemoveSquareRoots->True};
PolesAnalyze[result_,topsector_,krep_,n_,OptionsPattern[]]:=Module[{subset,zerolist,brep,supersec,xl,cutr,cut,adcut,singular={},sol={},LVlocal,sqlocal,sqt,intset,hintset,tem,tem1,tem2,tem3,tsingular={},cc,len,var,groebner,pos,pow,flag,brepr,bench,u,tsol,vrep,cover},
subset=Subsets[topsector]//ReverseSortBy[#,Length]&//DeleteCases[#,{}]&;
zerolist=GetMatZeroSector[result,n,Complement[Range[n],topsector]];(*all zero sectors*)
subset=DeleteCases[subset,_?(MemberQ[zerolist,Sector2Digits[#]]&)];(*remove all zero sectors*)
Print["totally ",Length[subset]," sectors need to be analyzed!"];

(*subset={{1,2,3,4,5,6,7}};*)
Monitor[Do[(*analyze sector by sector*)
	(*Print["subset: ",subset[[i]]];(*//////////////////////////////////*)*)
	singular={};
	brep=GetBaikovMatRep[result,subset[[i]],n];
	(*If[OptionValue[deBug],Print["path: ",brep]];*)
	If[OptionValue[AugAna],
		If[CheckValidity[brep,{},krep]==={},Continue[]];(*if some Grams in the representation equal to 0 before cut, then this sector is actually reducible already. We don't consider it.*)
		tem={CheckValidity[brep,subset[[i]],krep],Thread@Rule[Subscript[x,#]&/@subset[[i]],0]};(*check validity of the representations*)
		cc=Length[subset[[i]]];
		While[tem[[1]]==={},
			(*when maximal cut is 0, we consider the next-to-maximal cut and iterate this procedure*)
			(*we assume that there is at least one next-to-maximal cut is not 0, when it is not the case, a warning will be generated*)
			cc=cc-1;
			If[cc<=Length[subset[[i]]]-2,If[OptionValue[AdInfo],Message[PolesAnalyze::misp,subset[[i]],cc]]];
			cutr=Subsets[subset[[i]],{cc}];
			Do[
				tem={CheckValidity[brep,cutr[[j]],krep],Thread@Rule[Subscript[x,#]&/@cutr[[j]],0]};
				If[tem[[1]]=!={},Break[]]
			,{j,1,Length[cutr]}];
		];
		brep=tem[[1]];
		cutr=tem[[2]]
		,
		brep=GetBaikovMatRep[result,subset[[i]],n];
		cutr=Thread@Rule[Subscript[x,#]&/@subset[[i]],0];
	];
	If[OptionValue[SelectQ],If[OptionValue[SelectAllQ],len=2,len=Length[brep](*the number of minimal Baikov representations*)],len=1];
	Do[(*iteration for different representations of one sector*)
		intset={};
		hintset={};
		sol={};
		LVlocal={};
		sqlocal={};
		adcut=Complement[Subscript[x,#]&/@subset[[i]],Keys[cutr]];
		If[OptionValue[deBug],Print["adcut: ",adcut]];
		If[adcut=!={},adcut=Thread@Rule[adcut,0]];
		(*when this is a next-to-maximal cut, we need to remember that the propagators should be always set to 0, this is done by adcut.*)
		cut=Join[cutr,Thread@Rule[brep[[j,1]],0]];(*these variables brep[[j,1]] has been integrated out from the representation, we set it to 0 here*)
		Do[(*this iteration performs the classification of polynomials*)
			If[IntegerQ[brep[[j,2,1,k,2]]/.{\[Epsilon]->0}],
				tem=(brep[[j,2,1,k,1]]//Gram2Poly[#,krep]&)/.cut//Factor;
				If[FreeQ[tem,x],
					AppendTo[LVlocal,tem],
					intset=Join[intset,(FactorList[tem][[All,1]])/.adcut//DeleteCases[#,_?NumericQ]&]
				],
				tem=(brep[[j,2,1,k,1]]//Gram2Poly[#,krep]&)/.cut//Factor;
				If[FreeQ[tem,x],
					AppendTo[LVlocal,tem];AppendTo[sqlocal,tem],
					tem1=PerfectSquareSplit[tem];
					intset=Join[intset,(tem1[[1]])/.adcut//DeleteCases[#,_?NumericQ]&];
					tem1=PerfectSquareSplit[(Times@@(tem1[[2]])/.adcut)];
					intset=Join[intset,(tem1[[1]])//DeleteCases[#,_?NumericQ]&];
					hintset=Join[hintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&];
				]
			]
		,{k,1,Length[brep[[j,2,1]]]}];
		(*The following step further splits possible varieties independent of Baikov variable x*)
		(*some repeated terms in hintset should be moved to intset*)
		tem=Times@@hintset;
		tem1=PerfectSquareSplit[tem];
		intset=Join[intset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates;
		hintset=tem1[[2]]//DeleteCases[#,_?NumericQ]&;
		intset=(FactorList[#][[All,1]]&/@intset)//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates;(*keep intset to be irreducible polynomials*)
		(*Print["intset: ",intset];Print["hintset: ",hintset];*)
		intset=intset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,x]&]&;
		hintset=hintset//DeleteCases[#,_?NumericQ]&//GatherBy[#,FreeQ[#,x]&]&;
		(*If[OptionValue[deBug],Print["intset: ",intset];Print["hintset: ",hintset]];*)
		If[intset=!={},
			If[!FreeQ[intset[[1]],x],
				LVlocal=Join[LVlocal,Drop[intset,1]//Flatten[#,1]&];intset=intset[[1]],
				LVlocal=Join[LVlocal,intset[[1]]];intset=Drop[intset,1]//Flatten
			]
		];
		If[hintset=!={},
			If[!FreeQ[hintset[[1]],x],
				LVlocal=Join[LVlocal,Drop[hintset,1]//Flatten[#,1]&];sqlocal=Join[sqlocal,Drop[hintset,1]//Flatten[#,1]&];hintset=hintset[[1]],
				LVlocal=Join[LVlocal,hintset[[1]]];sqlocal=Join[sqlocal,hintset[[1]]];hintset=Drop[hintset,1]//Flatten
			]
		];
		If[intset==={}&&hintset==={},(*if there are no Baikov variable remaining*)
			AppendTo[singular,{{},LVlocal//DeleteDuplicates,(Times@@sqlocal),brep[[j,1]]}];Continue[]
		];
		intset=intset//DeleteDuplicates;
		sqt=(Times@@sqlocal);
		If[OptionValue[deBug],Print["intset,hintset: ",{intset,hintset}]];
		tem=ResolveSingularities[intset,hintset,SortQ->True,SelectQ->If[len>1,True,False],AdInfo->OptionValue[AdInfo],RemoveSquareRoots->OptionValue[RemoveSquareRoots]];(*when len>1, it means that we should always apply our selection rule*)
		If[OptionValue[deBug],Print["tem: ",tem];Print["sqlocal: ",sqlocal]];
		tem1=tem[[1]];
		sol=Join[sol,tem1//Flatten[#,1]&];
		LVlocal=Join[LVlocal,tem[[2]]];
		sqlocal=Join[{},sqt*tem[[3]]];
		(*the square roots from the coefficients should be multiplied to those constructed from integrands*)

		(*now we consider sending all the isps to infinity plane which corresponds to the second Landau singularity*)
		xl=Cases[{intset,hintset},Subscript[x,_],Infinity]//DeleteDuplicates;
		If[OptionValue[deBug],If[brep[[j,1]]==={Subscript[x, 9],Subscript[x, 8],Subscript[x, 7],Subscript[x, 5],Subscript[x, 2]},Print["subset: ",subset[[i]]];Print["intset,hintset: ",{intset,hintset}]]];
		If[Length[xl]>0,(*Length[xl]=1 case has been include before*)
			tem=Trans2Inf[intset,hintset,xl];
			intset=tem[[1]];
			hintset=tem[[2]];
			LVlocal=Join[LVlocal,tem[[3]]];
			sqt=sqt*Times@@(tem[[4]]);
			tem1=ResolveSingularities[intset,hintset,SortQ->True,SelectQ->If[len>1,True,False],AdInfo->OptionValue[AdInfo],RemoveSquareRoots->OptionValue[RemoveSquareRoots]];
			(*If[OptionValue[deBug],Print["tem1: ",tem1]];*)
			tem=Join[Thread@Rule[xl,Infinity],#]&/@((tem1[[1]]//Flatten[#,1]&));
			sol=Join[sol,tem(*//Flatten[#,1]&*)];
			LVlocal=Join[LVlocal,tem1[[2]]];
			sqlocal=Join[sqlocal,sqt*tem1[[3]]]
		];
		AppendTo[singular,{sol//DeleteDuplicates,LVlocal//DeleteDuplicates,sqlocal//DeleteDuplicates,brep[[j,1]]}]
	,{j,1,Length[brep]}];
	tem=singular[[All,2]];
	If[!OptionValue[SelectQ],
		(*now we try to remove all spurious letters in a more rigorous way, we identify them by two steps, the first step will rule out some and the second step rule all of them out in principle*)
		If[Length[tem]>1,
			tem3=Reap[Do[
				tem2=Delete[tem,k];
				var=Variables/@(tem[[k]]);
				cover=FindCover[var];(*find a cover of variables set. We will only take some subset of all variables as symbols and all others will be set to some numeric values*)
				(*when there are too many variables, the calculation of Groebner basis is difficult, here we adopt a semi-numerical method to speed up the calculation*)
				tem1={};
				Do[
					groebner=GroebnerBasis[tem[[k]]/.cover[[j,2]],cover[[j,1]],Modulus->Prime[2024]];(*calculate the Groebner basis from one representation*)
					tem1=Append[tem1,Table[Position[(PolynomialReduce[#/.cover[[j,2]],groebner,cover[[j,1]],Modulus->Prime[2024]]&/@tem2[[l]])[[All,2]],_?(#=!=0&),1]//DeleteCases[#,{0}]&,{l,1,Length[tem2]}]];(*Find those polynomials that can not be reduced by the groebner basis calculated before*)
				,{j,1,Length[cover]}];
				Sow[Join[{tem[[k]]},Table[Delete[tem2[[l]],tem1[[All,l]]//Flatten[#,1]&//DeleteDuplicates],{l,1,Length[tem2]}]]//Flatten]
				,{k,1,Length[tem]}];
			][[2]];
			If[tem3=!={},tem3=tem3[[1]]];
			tem3=Intersection[Sequence@@tem3],
			tem3=tem[[1]]
		];
		tem1=(FactorList[#][[All,1]]&/@tem3)//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates;
		(*in the last step, we check again the letters which can be spurious by the dimension calculation*)
		tem=Table[(FactorList[#][[All,1]]&/@tem[[k]])//Flatten//DeleteCases[#,_?NumericQ]&,{k,1,Length[tem]}];
		tem2=If[tem=!={},Intersection[Sequence@@(tem)],{}];
		tem3=Complement[tem1,tem2];
		(*Print["subset: ",subset[[i]]];*)
		(*If[subset[[i]]==={1,2,3,4,5,6,8},Print[" tem3: ",tem3];Print["cutr: ",cutr]];*)
		Do[
			flag=0;
			pos=Position[singular[[All,2]],_?(MemberQ[#,tem3[[k]]]&),1];
			brepr=(Delete[brep,pos]);
			var=Variables[tem3[[k]]];
			pow=Exponent[tem3[[k]],#]&/@var;
			pos=FirstPosition[pow,1];
			If[pos===Missing["NotFound"],
				If[subset[[i]]==={2,3,4,5,6,8},Print["tem3,var: ",{tem3,var}]];(*//////////////////////////////////*)
				sol=TimeConstrained[FindInstance[Join[{tem3[[k]]==0,var[[1]]>0},Thread@Unequal[var,0]],var,Rationals,2],60,0];(*find two numeric solutions to check*)
				If[sol===0,Print["haven't find solutions in given time. relax the condition and add this singularity: ",tem3[[k]]];Continue[]];
				If[sol==={},Print["this polynomial has no non-zero solution: ",tem3[[k]]," sec: ",i];sol=FindInstance[{tem3[[k]]==0},var]],
				tsol=Solve[tem3[[k]]==0,var[[pos]]][[1]];
				sol=Table[vrep=Thread@Rule[Delete[var,pos],Table[Prime[2024+l*m],{m,1,Length[var]-1}]];Join[Thread@Rule[var[[pos]],var[[pos]]/.tsol/.vrep],vrep],{l,1,2}]//DeleteDuplicates
			];
			(*Print["dimension calculation"];(*//////////////////////////////////*)*)
			Do[
				u=(((((brepr[[l,2,1]])//Gram2Poly[#,krep]&)/.cutr//Factor//DeleteCases[#,_?(FreeQ[First[#],x]&)]&)));
				xl=Cases[u,Subscript[x,_],Infinity]//DeleteDuplicates;
				If[xl==={},Continue[]];
				u=Times@@Power@@@u;
				bench=GetDimension[u,xl];
				If[bench===$Failed,Continue[]];(*overtime, we will pass this case*)
				(*If[subset[[i]]==={1,3,4,5,6,7,8},Print["brep: ",brepr[[l,2,1]]," cut: ",cut," u: ",u,"sol: ",sol," tem3[[k]]: ",tem3[[k]]]];*)
				If[Max[Table[GetDimension[u/.sol[[m]],xl],{m,1,Length[sol]}]]>=bench,flag=1;Break[]]
			,{l,1,Length[brepr]}];
			(*Print["dimension calculation finished!"];(*//////////////////////////////////*)*)
			If[flag==1,tem1=Complement[tem1,{tem3[[k]]}]]
		,{k,1,Length[tem3]}],
		tem=Table[(FactorList[#][[All,1]]&/@tem[[k]])//Flatten//DeleteCases[#,_?NumericQ]&,{k,1,Length[tem]}];
		tem2=If[tem=!={},Intersection[Sequence@@(tem)],{}];
		tem1=If[tem=!={},Union[Sequence@@(tem)],{}](*our old way to remove spurious letters which is experimental*)
	];
	tsingular=AppendTo[tsingular,{singular[[All,{1,2,3,4}]],tem2,tem1,singular[[All,3]]//Flatten//DeleteDuplicates,subset[[i]]}]
,{i,1,Length[subset]}],ProgressIndicator[i,{1,Length[subset]}]];
Return[tsingular];
];


RIntersection[l1_,l2_]:=Module[{int1={},int2={}},
Do[
If[MemberQ[l1,l2[[i]]|(-l2[[i]])],AppendTo[int2,l2[[i]]]]
,{i,1,Length[l2]}];
Do[
If[MemberQ[l2,l1[[i]]|(-l1[[i]])],AppendTo[int1,l1[[i]]]]
,{i,1,Length[l1]}];
Return[{int1,int2}];
];


(*Options[ExtractPoleInfo]={deBug->False,OutputLevel->1};
ExtractPoleInfo[exp_,n_,OptionsPattern[]]:=Module[{tem,xl,isp={},result={},rep,spresult={}},
xl=Subscript[x,#]&/@Range[n];
tem=(exp[[All,1]]//Flatten[#,1]&)[[All,{1,4}]];
Do[
isp=Cases[tem[[i,1]],Subscript[x,_],Infinity]//DeleteDuplicates;
rep=Thread@Rule[Complement[xl,isp],0];
If[OptionValue[deBug],Print["isp: ",isp]];
If[tem[[i,1]]==={},AppendTo[result,{{rep},tem[[i,2]]}]];
Do[
If[AllRuleQ[tem[[i,1,j]]],
AppendTo[result,{{rep,tem[[i,1,j]]},tem[[i,2]]}],
(*-------------This is the quadratic polynomial list--------------*)
AppendTo[spresult,{{rep,tem[[i,1,j]]},tem[[i,2]]}]
]
,{j,1,Length[tem[[i,1]]]}];
,{i,1,Length[tem]}];
tem=result//GatherBy[#,First]&;
result=Table[{tem[[i,1,1]],tem[[i,All,2]]},{i,1,Length[tem]}];
tem=spresult//GatherBy[#,First]&;
spresult=Table[{tem[[i,1,1]],tem[[i,All,2]]},{i,1,Length[tem]}];
If[OptionValue[OutputLevel]==1,Return[result]];
If[OptionValue[OutputLevel]==2,Return[{result,spresult}]];
];*)


Options[CompatiblePoleQ]={deBug->False,CutConstraint->False};
CompatiblePoleQ[p1_,p2_,OptionsPattern[]]:=Module[{pt1,pt2,tem,xint},
	If[p1==={}||p2==={},Return[True]];
	If[OptionValue[CutConstraint],(*if we set constraint on the propagators that are cut*)
		xint={Cases[p1[[1]],Subscript[x,_],Infinity]//DeleteDuplicates,Cases[p2[[1]],Subscript[x,_],Infinity]//DeleteDuplicates}//SortBy[#,Length]&;
		(* x variables which are cut variables*)
		If[!ContainsAll[xint[[2]],xint[[1]]],Return[False]];
	];
	xint=Intersection[Cases[p1,Subscript[x,_],Infinity]//DeleteDuplicates,Cases[p2,Subscript[x,_],Infinity]//DeleteDuplicates];
	If[FreeQ[p1,Rule[_,Infinity]],
		pt1=p1//Flatten,
		tem=p1[[2]]//DeleteCases[#,Rule[_,Infinity]]&;
		If[tem=!={},
			pt1=Join[p1[[1]],Thread@Rule[Keys[tem],g/@Values[tem]]],
			pt1=p1//Flatten
		]
	];(*distinguish replacement containing infinity*)
	If[FreeQ[p2,Rule[_,Infinity]],
		pt2=p2//Flatten,
		tem=p2[[2]]//DeleteCases[#,Rule[_,Infinity]]&;
		If[tem=!={},
			pt2=Join[p2[[1]],Thread@Rule[Keys[tem],g/@Values[tem]]],
			pt2=p2//Flatten
		]
	];(*distinguish replacement containing infinity*)
	If[SameQ[xint/.pt1,xint/.pt2],Return[True],Return[False]];
];


Options[MaxPole]={deBug->False};
MaxPole[p1_,p2_]:=Module[{a,b},
	If[!CompatiblePoleQ[p1,p2],Return[0]];
	a=Cases[p1,Subscript[x,_],Infinity];
	b=Cases[p2,Subscript[x,_],Infinity];
	If[ContainsAll[a,b],Return[1]];
	If[ContainsAll[b,a],Return[2]];
	Return[0];
];


Options[RemoveSPole]={deBug->False};
RemoveSPole[poles_,OptionsPattern[]]:=Module[{tem,tem1,l1=1,l2=2},
	tem=poles;
	While[l1<l2,
		tem1={};
		Do[
			If[MaxPole[tem[[l1]],tem[[i]]]==0,Continue[],
				If[MaxPole[tem[[l1]],tem[[i]]]==1,AppendTo[tem1,i](*keep track of those poles that are contained in this pole*),
					tem1={1};
					If[OptionValue[deBug],Print[tem[[l1]]]];
					Break[](*this pole is contained in another pole*)
				]
			]
		,{i,l1+1,Length[tem]}];
		If[tem1==={},l1=l1+1;l2=Length[tem],(*if there are no poles contain or be contained w.r.t this pole*)
			If[tem1==={1},tem=Delete[tem,l1];l2=Length[tem],(*if this pole is contained in another pole*)
				l1=l1+1;tem=Delete[tem,List/@tem1];l2=Length[tem](*if this pole contains other poles, remove them*)
			]
		];
		If[OptionValue[deBug],Print[tem1];Print[tem[[l1-1]]]];
	];
	Return[tem];
];


(*enhanced version of RemoveSPole*)
RemoveSPoleBeta[poles_]:=Module[{tem,tem1,den},
	tem=poles;
	(*now if the isp part of poles do not involve infinity, we directly substitute its values*)
	tem1=Reap[
		Do[
			If[!FreeQ[tem[[i]],Infinity],Sow[tem[[i]]];Continue[]];
			den=Denominator[Values[tem[[i,2]]]]//.tem[[i,2]]//Factor;
			If[!FreeQ[den,0],Sow[tem[[i]]];Continue[]];(*avoid the situation where some variable becomes infinity after substituting values*)
			Sow[{tem[[i,1]],Thread@Rule[Keys[tem[[i,2]]],Values[tem[[i,2]]]//.tem[[i,2]]//Factor]}]
		,{i,1,Length[tem]}];
	][[2]];
	If[tem1=!={},tem1=RemoveSPole[tem1[[1]]]];
	Return[tem1];
];


Options[CompatiblePoleGraph]={deBug->False,OutputLevel->2};
CompatiblePoleGraph[poles_,OptionsPattern[]]:=Module[{remain,graph,tem},
	remain=poles;
	(*Print["original total ",Length[remain]," poles"];*)
	graph=Reap[
			Do[
				Do[
					If[CompatiblePoleQ[remain[[i]],remain[[j]]],(*if this two poles are compatible*)
						Sow[UndirectedEdge[i,j]](*then we add an undirected edge to this graph*)
					]
				,{j,i+1,Length[remain]}]
			,{i,1,Length[remain]}];
	][[2,1]];
	If[OptionValue[OutputLevel]==1,Return[{graph}]];
	(*edge=EdgeList[graph]/.{UndirectedEdge->List};*)
	(*If[OptionValue[deBug],Print[FindClique[graph,{3},All]]];*)
	tem=FindClique[graph,Infinity,All]//Reverse;(*we find complete subgraphs formed by compatible poles*)
	If[tem=!={},Print["maximal clique: ",Length[tem[[-1]]],"; # of cliques: ",Length[tem]]];
	If[OptionValue[OutputLevel]==2,Return[{graph,tem}]];
];


Options[ExtractPoleInfo]={deBug->False,OutputLevel->1,Simp->True,Enhanced->False};
ExtractPoleInfo[exp_,n_,OptionsPattern[]]:=Module[{tem,tem1,xl,isp={},result={},poles},
result=Reap[Do[
	tem=exp[[i]];(*for every sector*)
	Sow[{Thread@Rule[Subscript[x,#]&/@tem[[-1]],0],#}&/@(tem[[1]][[All,1]]//Flatten[#,1]&//Cases[#,_?(AllRuleQ)]&),a];
	Sow[{Thread@Rule[Subscript[x,#]&/@tem[[-1]],0],#}&/@(tem[[1]][[All,1]]//Flatten[#,1]&//DeleteCases[#,_?(AllRuleQ)]&),b];
,{i,1,Length[exp]}];
];(*All poles appearing in the Baikov representations. For every pole, there may be some x variable missing (integrated out)*)
xl=Subscript[x,#]&/@Range[n];
If[OptionValue[OutputLevel]==1,
	tem1=result[[2,1]]//Flatten[#,1]&,(*poles only involve simple pole*)
	If[OptionValue[OutputLevel]==2,Return[result[[2,2]]//Flatten[#,1]&]];(*poles that involve quadratic polynomial*)
];
Print["original total ",Length[tem1]," poles"];
If[!OptionValue[Simp],Return[tem1]];
poles=RemoveSPole[tem1];(*remove poles that are already contained in other poles*)
If[OptionValue[Enhanced],poles=RemoveSPoleBeta[poles]];(*in enhanced version, we further remove some poles which can be identified to each other by substitute the values of the variables. For example, {x1->s-x2,x2->s} is equivalent to {x1->0,x2->s}*)
(*tem=CompatiblePoleMerge[poles];(*merge compatible poles*)
poles=RemoveSPole[tem];*)
Print["final total ",Length[poles]," poles"];
Return[{poles,Sequence@@CompatiblePoleGraph[poles]}];
];


ExtractSquareRoots[exp_]:=((Times@@Power@@@(MapAt[Mod[#,2]&,DeleteCases[FactorList[#],{_?NumericQ,_}],{All,2}]))&/@(exp[[All,4]]//Flatten//DeleteCases[#,_?NumericQ]&))//DeleteCases[#,1]&//DeleteDuplicates;


Options[ConstructLetter]={deBug->False};
ConstructLetter[Gl1_,Gl2_,type_,OptionsPattern[]]:=Module[{G1,G2,l1,l2,od,mm,MM,var,var1,var2,s1,s2,tem,tem1,sys,sol,mat,mat1,t,t1,pos,r},
If[Head[Gl1[[1]]]===G,G1=Gl1[[1]],If[MatchQ[Gl1[[1]],-G[__]],G1=-Gl1[[1]],Return[$Failed]]];
If[Head[Gl2[[1]]]===G,G2=Gl2[[1]],If[MatchQ[Gl2[[1]],-G[__]],G2=-Gl2[[1]],Return[$Failed]]];
l1=Length[G1[[1]]];
l2=Length[G2[[1]]];
(*if this two Grams are the same with each other, we give a warning*)
If[(G1[[1]]//Sort)===(G2[[1]]//Sort),Message[ConstructLetter::warning,{Gl1,Gl2}]];
If[Abs[l1-l2]==0,
	var=Union[Variables[G1[[1]]],Variables[G2[[1]]]];
	tem=Intersection[G1[[1]],G2[[1]]];
	If[Length[tem]==(l1-1),
		(*if all but one are common*)
		mat=Table[t[j],(*{i,1,l1-1-Length[tem]},*){j,1,Length[G1[[1]]]}];
		sys=Thread@Equal[(Coefficient[mat . G1[[1]]- (Complement[G2[[1]],tem][[1]]),#]&/@var)//Flatten,0];
		sol=FindInstance[sys,mat,Integers];
		If[sol==={},
			mm=tem;od={G1[[1]],G2[[1]]};MM=Union[G1[[1]],G2[[1]]];Goto[FINAL],
			Message[ConstructLetter::warning,{Gl1,Gl2}]
		],
		
		(*if common term is less than l-1, then we solve their relations to find more common entries*)
		s1=Complement[G1[[1]],tem];
		s2=Complement[G2[[1]],tem];
		r=l1-1-Length[tem];
		(*try to solve the relations existing in the remaining entries*)
		mat=Table[t[j],(*{i,1,l1-1-Length[tem]},*){j,1,Length[s1]}];
		mat1=Table[t1[j],(*{i,1,l1-1-Length[tem]},*){j,1,Length[s2]}];
		sys=(Coefficient[(Coefficient[mat . s1-mat1 . s2,#]&/@var)//Flatten,#]&/@Join[mat,mat1])//Transpose;
		sol=NullSpace[sys];
		If[OptionValue[deBug],Print["sys: ",sys];Print["sol: ",sol]];
		sol=Table[Thread@Rule[Join[mat,mat1],sol[[i]]],{i,1,Length[sol]}];
		mat=(mat/.sol)//RowReduce;(*note that there are sevel set solutions, so this is a matrix*)
		mat1=(mat1/.sol)//RowReduce;
		If[Length[sol]<r,Message[ConstructLetter::fail,{G1[[1]],G2[[1]]}]];
		(*If[OptionValue[deBug],Print["s1,s2: ",{s1,s2}];Print["mat: ",mat]];*)
		tem1={};
		Do[
			pos=FirstPosition[mat[[i]],_?(#!=0&)];
			AppendTo[tem1,pos[[1]]];
		,{i,1,r}];
		If[OptionValue[deBug]&&!FreeQ[pos,"NotFound"],Print["tem1: ",tem1];Print["mat: ",mat];Print["sol: ",sol];Print["G1 G2: ",{Gl1,Gl2}]];
		tem1=Delete[s1,Partition[tem1,1]];(*there exists linear combinations that these entries in s1 can be transformed to equivalent form with s2*)
		mm=Join[mat . s1,tem];(*this is the common entries after the transformation*)
		od={Join[mm,tem1]};
		tem1={};
		Do[
			pos=FirstPosition[mat1[[i]],_?(#!=0&)];
			AppendTo[tem1,pos[[1]]];
		,{i,1,r}];
		tem1=Delete[s2,Partition[tem1,1]];(*delete those entries that can be brought to be common with s1*)
		od=Append[od,Join[mm,tem1]];
		MM=Union[od//Flatten];
		(*usually this kind of letter is in a more complicated form, we add a minus sign before type to distinguish them*)
		Goto[FINAL]
	]
	,
	
	If[l1>l2,s1=G1[[1]];s2=G2[[1]],s1=G2[[1]];s2=G1[[1]]];
	If[s2==={},mm=s2;MM=s1;od=Partition[s1,1];Goto[FINAL]];(*this is the special case where one of the Gram is G[{},{}]=1*)
	mm=s2;
	var1=Variables[s1];
	var2=Variables[s2];
	mat=Table[t[i,j],{i,1,Length[s2]},{j,1,Length[s1]}];
	sys=Thread@Equal[(Coefficient[mat . s1-s2,#]&/@var1)//Flatten,0];
	sol=FindInstance[sys,mat//Flatten,Integers];
	(*If[OptionValue[deBug],Print["mm: ",mm];Print["sol: ",sol]];*)
	mat=(mat/.sol[[1]])//RowReduce;
	tem={};
	Do[
		pos=FirstPosition[mat[[i]],_?(#!=0&)];
		AppendTo[tem,pos[[1]]];
		If[Length[tem]==Length[s2],Break[]]
	,{i,1,Length[mat]}];
	tem=Delete[s1,Partition[tem,1]];
	If[Length[tem]!=2,Message[ConstructLetter::err,{s1,s2}];Return[$Failed]];
	od={Join[s2,{tem[[1]]}],Join[s2,{tem[[2]]}]};
	MM=Union[od//Flatten];
];
Label[FINAL];
If[type==1,
	If[Abs[l1-l2]==0,
		Return[{dlogForm[G[Sequence@@od],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]],-G[mm,mm]*G[MM,MM],{(*Gl1[[3]],Gl2[[3]]*)},1}],
		Return[{dlogForm[G[Sequence@@od],-G[mm,mm]*G[MM,MM]],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]],{(*Gl1[[3]],Gl2[[3]]*)},1}]
	],
	If[Abs[l1-l2]==0,
		Return[{dlogForm[G[Sequence@@od],-G[mm,mm]*G[MM,MM]],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]],{(*Gl1[[3]],Gl2[[3]]*)},2}],
		Return[{dlogForm[G[Sequence@@od],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]],-G[mm,mm]*G[MM,MM],{(*Gl1[[3]],Gl2[[3]]*)},2}]
	]
];
];

ConstructFromGram[gramlist_]:=Table[ConstructLetter[gramlist[[i,1]],gramlist[[i,j]],If[IntegerQ[gramlist[[i,1,2]]/.{\[Epsilon]->0}],2,1]],{i,1,Length[gramlist]},{j,2,Length[gramlist[[i]]]}]//Flatten[#,1]&//DeleteDuplicatesBy[#,({First[#],Last[#]}&)]&;

LetterPathRefined[reform_,letters_,krep_]:=Module[{glist,tem,path,pos,letterref},
(*Note that the structure of reform is {{g,c,{paths}},{g,c,{paths}},...}*)
glist=reform[[All,1]];
letterref=Reap[
Do[
tem=Cases[Cases[letters[[i,1]],Power[_,1/2],Infinity]//DeleteDuplicates,_G,Infinity]//DeleteDuplicates;
tem=DeleteCases[tem,_?(FreeQ[GramMat[#[[1]],#[[2]],krep],x]&)];
If[tem==={},Sow[Append[letters[[i]],Infinity]];Continue[]];
pos=FindGram[tem,glist,krep];
path=Intersection@@((Union@@(reform[[#,3]]))&/@pos);
Sow[Append[letters[[i]],path]]
,{i,1,Length[letters]}]
];
Return[letterref[[2,1]]];
];


Options[ExtractAlgLetter]={deBug->False,KineVar->{},KineRestrict->False,LoopRestrict->False,OneLoop->False};
ExtractAlgLetter[result_,topsector_,n_,OptionsPattern[]]:=Module[{subset,zerolist,reform,flag,rep,glist,tem,cmGram,letter,letterref,summary={},a},
reform=ReformRep[result];(*Note that the structure of reform is {{g,c,{paths}},{g,c,{paths}},...}*)
subset=Subsets[topsector]//ReverseSortBy[#,Length]&;(*all zero sectors*)
zerolist=GetMatZeroSector[result,n,Complement[Range[n],topsector]];
subset=DeleteCases[subset,_?(MemberQ[zerolist,Sector2Digits[#]]&)];(*remove all zero sectors*)
Print["totally ",Length[subset]," sectors need to be analyzed!"];

letterref={};(*keep track of algbraic letters so that repeated ones will not included in the result*)
Monitor[Do[
	If[OptionValue[deBug],Print["subset: ",subset[[a]]]];
	rep=GetBaikovMatRep[result,subset[[a]],n]//Simplify;
	cmGram=ExtractRelevantGramSim[rep,reform,n,KineVar->OptionValue[KineVar],KineRestrict->OptionValue[KineRestrict],LoopRestrict->OptionValue[LoopRestrict],OneLoop->OptionValue[OneLoop](*,deBug->If[subset[[a]]==={2,3,4,6},True,False]*)];
	(*If[OptionValue[deBug],Print["cmGram: ",cmGram]];*)
	letter=ConstructFromGram[cmGram];(*//DeleteCases[#,_?(!AdmissibleLetterQ[#]&)]&*);(*element: {Log[_],_,1|2}*)
	letterref=Join[letterref,letter]//DeleteDuplicatesBy[#,({First[#],Last[#]}&)]&;(*Delete duplicates by the first and last term.*)
	letter=Intersection[letter,letterref];(*take only unique terms that haven't appeared so far*)
	(*letterref=LetterPathRefined[reform,letter,krep];*)
	AppendTo[summary,{letter,subset[[a]]}];
	(*If[a==3,Break[]]*)
	,{a,1,Length[subset]}],
ProgressIndicator[a,{1,Length[subset]}]];
Return[summary];
];


Options[AllRationalLetters]={OutputLevel->2};
AllRationalLetters[polestructure_,OptionsPattern[]]:=
If[OptionValue[OutputLevel]==1,
Return[(FactorList[#][[All,1]]&/@(polestructure[[All,2]]//Flatten//DeleteDuplicates))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates],
Return[(FactorList[#][[All,1]]&/@(polestructure[[All,3]]//Flatten//DeleteDuplicates))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates]
];


Options[SpecialSimplify]={AdRep->{},deBug->False};
SpecialSimplify[{gr_,ga_},OptionsPattern[]]:=Catch@Module[{rep,rp,tem,sq,numsq,densq},
sq=ga//Factor;
If[sq===0,Throw[{0,0}]];
rp=gr//Factor;
If[rp===0,Throw[{0,0}]];(*If there is one zero, then this letter will be degenerate*)
rep=OptionValue[AdRep];
(*we apply the replacement to the expression one-by-one*)
If[rep==={},
	If[Head[sq]===Times,sq=List@@sq,sq={sq}];
	tem=Table[If[MatchQ[sq[[i]],Power[_,_?EvenQ]],Replace[sq[[i]],Power[z_,a_]:>Power[z,a/2]],If[MatchQ[sq[[i]],Power[_,_?OddQ]],Replace[sq[[i]],Power[z_,a_]:>Power[R*z,1/2]*Power[z,Quotient[a,2]]],Power[sq[[i]],1/2]]],{i,1,Length[sq]}];(*take the square root, R is introduced to avoid things like Power[a,3/2]*)
	tem=((Times@@tem)/(rp))//Factor;(*cancel common parts between gr and Sqrt[ga]*)
	If[tem===0,Throw[{0,0}]];
	numsq=Cases[{Numerator[tem]},Power[z_,1/2]->z,Infinity];(*treat the square root in numerator and denomiantor separately*)
	densq=Cases[{Denominator[tem]},Power[z_,1/2]->z,Infinity];
	sq=Join[numsq,densq];
	densq=Times@@densq;
	If[OptionValue[deBug],Print["tem ",tem];Print["sq ",sq];Print["densq: ",densq]];
	tem=tem/.{Power[_,1/2|-1/2]->1}//Factor;(*remove square root part*)
	Throw[{densq/(tem)/.{R->1}//Factor,(Times@@(sq(*/.{Power[z_,1/2]:>z}*)/.{R->1}))//Factor}],

	Do[
		sq=sq/.rep[[i]]//Factor;
		If[sq===0,Throw[{0,0}]];
		If[(Denominator[rp]/.rep[[i]]//Factor)===0,rp=0;sq=0;Break[]];(*when the rational part becomes infinity after replacement*)
		rp=rp/.rep[[i]]//Factor;
		If[rp===0,Throw[{0,0}]];(*when the rational part or square root part become 0 after replacement*)
		If[Head[sq]===Times,sq=List@@sq,sq={sq}];
		tem=Table[If[MatchQ[sq[[i]],Power[_,_?EvenQ]],Replace[sq[[i]],Power[z_,a_]:>Power[z,a/2]],If[MatchQ[sq[[i]],Power[_,_?OddQ]],Replace[sq[[i]],Power[z_,a_]:>Power[R*z,1/2]*Power[z,Quotient[a,2]]],Power[R*sq[[i]],1/2]]],{i,1,Length[sq]}];
		If[OptionValue[deBug],Print["sq: ",sq];Print["rp: ",rp];Print["tem: ",tem]];
		tem=((Times@@tem)/(rp))//Factor;
		If[tem===0,rp=0;sq=0;Break[]];
		numsq=Cases[{Numerator[tem]},Power[z_,1/2]->z,Infinity];(*treat the square root in numerator and denomiantor separately*)
		densq=Cases[{Denominator[tem]},Power[z_,1/2]->z,Infinity];
		sq=Join[numsq,densq];
		densq=Times@@densq;
		If[OptionValue[deBug],Print["tem ",tem];Print["sq ",sq];Print["densq: ",densq]];
		tem=tem/.{Power[_,1/2|-1/2]->1}//Factor;
		rp=densq/(tem)//Factor;
		sq=(Times@@(sq(*/.{Power[z_,1/2]:>z}*)))//Factor;
		If[Denominator[sq]=!=1,rp=rp*Denominator[sq];sq=Numerator[sq]*Denominator[sq]];(*when expression under square root has denominator, then multiply a common factor to remove it*)
		rp=rp/.{R->1};sq=sq/.{R->1};
		If[OptionValue[deBug],Print["{rp,sq}: ",{rp,sq}]]
	,{i,1,Length[rep]}];
	Throw[{rp//Factor,sq//Factor}]
];
];


Options[ApplyPoleToLetter]={deBug->False,AdRep->{}};
ApplyPoleToLetter[pole_,letter_,krep_,OptionsPattern[]]:=Module[{tem,ntem,rep,rep1,gr,ga,isp,pos},
(*the input of pole should be in the form like {{Subscript[x, 1]->0,Subscript[x, 2]->0,...},{Subscript[x, 5]->a,Subscript[x, 6]->b,Subscript[x, 7]->Infinity,...,Subscript[x, 8]->Infinity,Subscript[x, 7]->Subscript[x, 8],Subscript[x, 8]->0}}. After the replacement of Infinity, the curve should be projective.*)
ga=((Cases[letter,Power[z_,1/2]->z,Infinity]//DeleteDuplicates)[[1]])/.{G[{},{}]->1};(*algebraic part of letter*)
gr=(letter[[1]]//Numerator)/.{Power[_,1/2]->0};(*rational part of letter*)
If[pole=!={},
	tem=SpecialSimplify[({gr,ga}//Gram2Poly[#,krep]&)/.pole[[1]]];
	If[OptionValue[deBug],Print["gr: ",gr];Print["ga: ",ga];Print["tem: ",tem]];
	Do[(*transform the poles of infinity to the procedure of taking x0 to 0*)
		If[!FreeQ[pole[[k]],Infinity],
			rep=Split[pole[[k]],(Last[#1]===Last[#2])&];
			ntem=Table[If[FreeQ[rep[[i]],Infinity],rep[[i]],{rep[[i]],Subscript[x, 0]->0}/.{Rule[Subscript[x,a_],Infinity]:>Rule[Subscript[x,a],Subscript[x,a]/Subscript[x,0]]}],{i,1,Length[rep]}]//Flatten[#,1]&;
			isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;(*keep track of variables taken to infinity*)
			(*Do[
				pos=Position[Keys[ntem],isp[[i]]];
				If[Length[pos]>1,Continue[],If[ntem[[pos[[1,1]]+1,1]]=!=Subscript[x, 0],ntem=Insert[ntem,Subscript[x, 0]->0,pos[[1,1]]+1]]];
			,{i,1,Length[isp]}];(*it is different whether taking variables to infinity one by one or taking them to infinity at one time. *)*)
			If[OptionValue[deBug],Print["ntem: ",ntem];Print["poles: ",pole[[k]]];Print["isp: ",isp]];
			(*rep=Thread@Rule[isp,isp/Subscript[x, 0]];
			ntem=DeleteCases[pole[[k]],Rule[_,Infinity]];*)
			rep={};
			rep1=If[ntem=!={},If[MemberQ[isp,ntem[[-1,1]]],Drop[ntem,-1],Print["not projective curve!: ",ntem];ntem],{}](*when projecting to infinity space, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*)
			,
			rep=pole[[k]];rep1={}
		];
		If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1];Print["tem: ",tem]];
		tem=SpecialSimplify[tem,AdRep->Join[rep,rep1,{OptionValue[AdRep]}]];
	,{k,2,Length[pole]}],
	(*if pole==={}, then this only involves pure kinematics grams*)
	tem=SpecialSimplify[({gr,ga}//Gram2Poly[#,krep]&)]
];
If[((tem[[1]]^2-tem[[2]]//Factor)===0),Return[0]];
Return[Log[(tem[[1]]+Sqrt[tem[[2]]])/(tem[[1]]-Sqrt[tem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&]]
];

Options[ApplyPoleToGram]={deBug->False,AdRep->{}};
ApplyPoleToGram[pole_,gram_,krep_,OptionsPattern[]]:=Module[{tem,ntem,fac,rep,rep1,isp,pos},
	tem=Times@@(Det/@((gram/.{G[a_,b_]:>GramMat[a,b,krep]})/.pole[[1]]))/.OptionValue[AdRep]//Factor;
	If[Length[pole]==1,tem=Times@@(PerfectSquareSplit[tem][[2]])];
	Do[
		If[!FreeQ[pole[[k]],Infinity],
			isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
			If[OptionValue[deBug],Print["isp: ",isp]];
			rep=Thread@Rule[isp,isp/Subscript[x, 0]];
			ntem=DeleteCases[pole[[k]],Rule[_,Infinity]];
			rep1=If[ntem=!={},If[MemberQ[isp,ntem[[-1,1]]],Drop[ntem,-1],ntem],{}]
			(*when projecting to infinity plane, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*)
			,
			pos=Position[pole[[k]],{__},1];
			If[pos==={},rep={};rep1=pole[[k]],rep=Take[pole[[k]],pos[[1,1]]]//Flatten;rep1=Drop[pole[[k]],pos[[1,1]]]];
		];
		(*tem=SpecialSimplify[{1,Times@@tem},AdRep->Join[rep,rep1,{OptionValue[AdRep]}]];*)
		ntem=tem;
		If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
		If[OptionValue[deBug],Print["tem: ",tem]];
		If[!FreeQ[rep,Subscript[x, 0]],
			ntem=Numerator[ntem/.rep//Factor]/.{Subscript[x, 0]->0};(*here we take the variables to infinity first*)
			Do[
				ntem=Times@@NumeratorDenominator[ntem/.rep1[[i]]//Factor];
				ntem=Times@@(PerfectSquareSplit[ntem][[2]])(*remove perfect square*)
			,{i,1,Length[rep1]}],
			Do[
				ntem=Times@@NumeratorDenominator[ntem/.rep[[i]]//Factor];
				ntem=Times@@(PerfectSquareSplit[ntem][[2]])(*remove perfect square*)
			,{i,1,Length[rep]}];
		];
		(*If[rep1=!={},rep1=Thread@Rule[Keys[rep1],Values[rep1]//.Drop[rep1,-1]]//Factor];
		If[rep1=!={},rep1=Thread@Rule[Keys[rep1],Values[rep1]//.rep1]//Factor];*)
		tem=ntem;
	,{k,2,Length[pole]}];
Return[tem];
];

(*Options[ApplyPoleToGram]={deBug->False,AdRep->{}};
ApplyPoleToGram[pole_,gram_,krep_,OptionsPattern[]]:=Module[{tem,ntem,rep,rep1,isp,den,pos},
tem=(Det/@((gram/.{G[a_,b_]:>GramMat[a,b,krep]})/.pole[[1]]))/.OptionValue[AdRep]//Factor;
ntem=tem;
Do[
If[!FreeQ[pole[[k]],Infinity],
isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
If[OptionValue[deBug],Print["isp: ",isp]];
rep=Thread@Rule[isp,isp/Subscript[x, 0]];
ntem=DeleteCases[pole[[k]],Rule[_,Infinity]];
rep1=If[ntem=!={},If[MemberQ[isp,ntem[[-1,1]]],Drop[ntem,-1],ntem],{}]
(*when projecting to infinity plane, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*)
,
rep=pole[[k]];rep1={}];
ntem=tem;
If[rep=!={}&&FreeQ[rep,Subscript[x, 0]],rep=Thread@Rule[Keys[rep],Values[rep]//.Drop[rep,-1]]//Factor]//Quiet;
If[rep=!={}&&FreeQ[rep,Subscript[x, 0]],rep=Thread@Rule[Keys[rep],Values[rep]//.rep]//Factor]//Quiet;
If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
If[OptionValue[deBug],Print["tem: ",tem]];
Do[
If[FreeQ[rep[[i]],ComplexInfinity],ntem=ntem/.rep[[i]]//Factor,ntem=(Series[ntem,Keys[rep[[i]]]->Infinity]//Normal)/.{Keys[rep[[i]]]->1}//Factor;(*Print["rep: ",rep];Print["ntem: ",ntem]*)];
,{i,1,Length[rep]}];
If[rep1=!={},rep1=Thread@Rule[Keys[rep1],Values[rep1]//.Drop[rep1,-1]]//Factor];
If[rep1=!={},rep1=Thread@Rule[Keys[rep1],Values[rep1]//.rep1]//Factor];
den=(Denominator/@(ntem))/.{Subscript[x, 0]->1}//Factor;(*the denominator*)
tem=(Numerator/@(ntem))/.{Subscript[x, 0]->0}//Factor;
ntem=tem/den//Factor;
Do[
den=(Denominator/@(ntem))/.rep1[[i]]//Factor;(*the denominator*)
tem=(Numerator/@(ntem))/.rep1[[i]]//Factor;
ntem=tem/den//Factor
,{i,1,Length[rep1]}];
,{k,2,Length[pole]}];
Return[ntem];
];*)


Options[AdmissiblePoleQ]={deBug->False};
AdmissiblePoleQ[glist_,pole_,permsq_,krep_,OptionsPattern[]]:=Module[{tem,ntem},
If[pole==={},
	(*cases where glist is list of Grams of external kinematics*)
	tem=(Times@@glist)//Gram2Poly[#,krep]&//Factor,
	tem=ApplyPoleToGram[pole,glist,krep]/.{Subscript[x,_]->1};
];

If[OptionValue[deBug],Print["tem: ",tem];];
If[NumericQ[tem],Return[False]];(*when the polynomials under square root are all perfect square, it is actually a rational letter, we don't count them here*)
If[AnyTrue[permsq,NumericQ[tem/#//Cancel]&],Return[True],Return[False]];
];


(*Options[AdmissiblePoleQ]={deBug->False};
AdmissiblePoleQ[glist_,pole_,permsq_,krep_,OptionsPattern[]]:=Module[{tem,ntem},
If[pole==={},
	(*cases where glist is list of Grams of external kinematics*)
	tem=glist//Gram2Poly[#,krep]&//Factor,
	tem=ApplyPoleToGram[pole,glist,krep]/.{Subscript[x,_]->1};
];

ntem=DeleteCases[tem,_?PerfectSquareQ];
If[OptionValue[deBug],Print["tem: ",tem];Print["ntem: ",ntem]];
If[ntem==={},Return[False]];(*when the polynomials under square root are all perfect square, it is actually a rational letter, we don't count them here*)
If[And@@Table[AnyTrue[permsq,NumericQ[ntem[[i]]/#//Cancel]&],{i,1,Length[ntem]}],Return[True],Return[False]];
];*)


(*Options[ApplyPolesToAlgLetter]={Sector->{},AdRep->{},deBug->False,PermSq->{},PathDis->True};
ApplyPolesToAlgLetter[poles_,letters_,reform_,krep_,OptionsPattern[]]:=Module[{start,polepath,letterpath,kletters={},tem,ntem,glist,Glist,ll,nletters={},sletters={},var,rep,rep1,gr,ga,isp,pos,path},
start=SessionTime[];
Monitor[
Do[
letterpath=letters[[i,3]];
ga=(Cases[letters[[i,1]],Power[z_,1/2]->z,Infinity]//DeleteDuplicates)[[1]];(*algebraic part of letter*)
gr=(letters[[i,1]][[1]]//Numerator)/.{Power[_,1/2]->0};(*rational part of letter*)
tem=({gr,ga}//Cases[#,_G,Infinity]&)/.{G[a_,b_]:>GramMat[a,b,krep]};
If[FreeQ[tem,x],
ntem=SpecialSimplify[({gr,ga}//Gram2Poly[#,krep]&)//Factor];
If[((ntem[[1]]-ntem[[2]]//Factor)===0)||((ntem[[1]]+ntem[[2]]//Factor)===0),Continue[]];
glist=Cases[{ga},_G,Infinity]//DeleteDuplicates;
ll=(Length[#[[1]]]&/@glist);
Glist=Select[reform,MemberQ[ll,Length[#[[1,1]]]]&];
pos=FindGram[glist,Glist[[All,1]],krep];
If[!FreeQ[pos,{}],Continue[](*Message[FindGram::warning,Gram]*)];
AppendTo[kletters,{Log[(ntem[[1]]+Sqrt[ntem[[2]]])/(ntem[[1]]-Sqrt[ntem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&],letters[[i,1]],letters[[i,4]]}];
(*If[OptionValue[deBug],Print["time consuming: ",SessionTime[]-start]];*)
Continue[],
glist=Cases[{ga},_G,Infinity]//DeleteDuplicates;
(*If[OptionValue[deBug],Print["glist: ",glist]];*)
ll=(Length[#[[1]]]&/@glist);
Glist=Select[reform,MemberQ[ll,Length[#[[1,1]]]]&];
pos=FindGram[glist,Glist[[All,1]],krep];
If[!FreeQ[pos,{}],Continue[](*Message[FindGram::warning,Gram]*)];
path=(Union@@(Glist[[#,3]]))&/@pos;(*paths for grams in the square root*)
If[OptionValue[PathDis]&&PathsDistance[path[[1]],path[[2]]]>1,AppendTo[sletters,letters[[i]]];Continue[]];
ntem=GramMat[#[[1]],#[[2]],krep]&/@glist;
pos=Table[If[FreeQ[ntem[[k]],x],{k}],{k,1,Length[ntem]}]//DeleteCases[#,Null]&;(*delete Gram which doesn't depend on Baikov variables*)
path=Delete[path,pos];
(*glist=Delete[glist,pos];*)
(*If[OptionValue[deBug],Print["time consuming: ",SessionTime[]-start]];*)
(*If[OptionValue[deBug],Print["path: ",path];Print["glist: ",glist]];*)
];
If[OptionValue[PermSq]==={},
Do[
polepath=poles[[j,2]]//SortBy[#,Length]&;
(*polepath=RefinedPolePath[glist,poles[[j]],poles,krep]//SortBy[#,Length]&;*)(*some poles are special and they can appear in subsectors*)
If[AllTrue[polepath,ContainsAny[#,OptionValue[Sector]]&],(*if there is any element in the path, it means this pole only appear in subsectors, so if we only want pole appearing top sector they will be discarded*)
Continue[],
If[AnyTrue[polepath,(*MemberQ[letterpath[[2]],#]&*)(MemberQ[Join[{letterpath[[1]]},letterpath[[2]]],#])&],
(*If[!AnyTrue[path,IntersectingQ[#,polepath]&],Continue[]];*)
If[!AllTrue[path,CompatiblePathQ[#,polepath,TopOnly->True]&],Continue[]];
tem=SpecialSimplify[({gr,ga}/.{G[a_,b_]:>Det[GramMat[a,b,krep]/.poles[[j,1,1]]]})//Factor];
(*If[OptionValue[deBug],Print["gr ga: ",{gr,ga}]];*)
Do[
If[!FreeQ[poles[[j,1,k]],Infinity],
isp=Cases[poles[[j,1,k]],Subscript[x,_],Infinity]//DeleteDuplicates;
(*If[OptionValue[deBug],Print["isp: ",isp]];*)
rep=Thread@Rule[isp,isp/Subscript[x, 0]];
ntem=DeleteCases[poles[[j,1,k]],Rule[_,Infinity]];
rep1=If[ntem=!={},Drop[ntem,-1],{}](*when projecting to infinity space, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*),
rep=poles[[j,1,k]];rep1={}];
(*If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
If[OptionValue[deBug],Print["tem: ",tem]];*)
tem=SpecialSimplify[tem,AdRep->Join[rep,{Subscript[x, 0]->0},rep1,{OptionValue[AdRep]}]];
,{k,2,Length[poles[[j,1]]]}];
If[((tem[[1]]-tem[[2]]//Factor)===0)||((tem[[1]]+tem[[2]]//Factor)===0),Continue[]];
AppendTo[nletters,{Log[(tem[[1]]+Sqrt[tem[[2]]])/(tem[[1]]-Sqrt[tem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&],letters[[i,1]],poles[[j,1]],letters[[i,4]]}];
]
]
,{j,1,Length[poles]}],
(*In this case, we use additional information to constrain the form*)
Do[
polepath=poles[[j,2]]//SortBy[#,Length]&;
(*polepath=RefinedPolePath[glist,poles[[j]],poles,krep]//SortBy[#,Length]&;*)(*some poles are special and they can appear in subsectors*)
If[!AdmissiblePoleQ[glist,poles[[j,1]],OptionValue[PermSq],krep],
Continue[],
If[False(*!AllTrue[path,CompatiblePathQ[#,polepath]&]*),
Continue[],
(*If[OptionValue[deBug],Print["time consuming: ",SessionTime[]-start]];*)
tem=SpecialSimplify[({gr,ga}/.{G[a_,b_]:>Det[GramMat[a,b,krep]/.poles[[j,1,1]]]})//Factor];
(*If[OptionValue[deBug],Print["gr ga: ",{gr,ga}]];*)
Do[
If[!FreeQ[poles[[j,1,k]],Infinity],
isp=Cases[poles[[j,1,k]],Subscript[x,_],Infinity]//DeleteDuplicates;
(*If[OptionValue[deBug],Print["isp: ",isp]];*)
rep=Thread@Rule[isp,isp/Subscript[x, 0]];
ntem=DeleteCases[poles[[j,1,k]],Rule[_,Infinity]];
rep1=If[ntem=!={},Drop[ntem,-1],{}](*when projecting to infinity space, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*),
rep=poles[[j,1,k]];rep1={}];
(*If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
If[OptionValue[deBug],Print["tem: ",tem]];*)
tem=SpecialSimplify[tem,AdRep->Join[rep,{Subscript[x, 0]->0},rep1,{OptionValue[AdRep]}]];
,{k,2,Length[poles[[j,1]]]}];
If[((tem[[1]]-tem[[2]]//Factor)===0)||((tem[[1]]+tem[[2]]//Factor)===0),Continue[]];
AppendTo[nletters,{Log[(tem[[1]]+Sqrt[tem[[2]]])/(tem[[1]]-Sqrt[tem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&],letters[[i,1]],poles[[j,1]],letters[[i,4]]}];
(*If[OptionValue[deBug],Print["time consuming: ",SessionTime[]-start]];*)
]
]
,{j,1,Length[poles]}]]
,{i,1,Length[letters]}],i];
Return[{nletters//DeleteDuplicatesBy[#,First]&,kletters//DeleteDuplicatesBy[#,First]&,sletters}];
];*)


(*Options[AllSectorAlgLetter]={deBug->False,PermSq->{},PathDis->True};
AllSectorAlgLetter[poles_,algletter_,reform_,krep_,OptionsPattern[]]:=Module[{l,tem,ntem,result={},spresult={}},
l=Length[algletter];
Print["totally ",l," sectors need to be analyzed"];
Do[
If[OptionValue[deBug],Print["subset: ",algletter[[i,2]]]];
tem=ApplyPolesToAlgLetter[poles,algletter[[i,1]],reform,krep,PermSq->OptionValue[PermSq],PathDis->OptionValue[PathDis],deBug->True];
ntem=Flatten[tem[[{1,2}]],1]//DeleteDuplicatesBy[#,First]&;
AppendTo[result,{ntem,algletter[[i,2]]}];
AppendTo[spresult,{tem[[3]],algletter[[i,2]]}]
,{i,1,Length[algletter]}];
If[OptionValue[deBug],Return[{result,spresult}]];
Return[result];
];*)


RefinedPolePath[gram_,pole_,poles_,krep_]:=Module[{start,tem,ntem,int,path={}},
(*Print["start:",start=SessionTime[]];*)
(*If[FreeQ[pole[[1]],Infinity],Return[pole[[2]]]];*)
tem=ApplyPoleToGram[pole[[1]],gram,krep]/.{Subscript[x,_]->1};(*now baikov variable can be scaled to 1*)
Do[
(*Print["poles i: ",i,SessionTime[]-start];*)
If[CompatiblePathQ[poles[[i,2]],pole[[2]],SubOnly->True],(*if there is a pole in subsector which will give the same expression for one gram as this pole, then this pole can be applied to subsector.*)
(*Print["time before ",SessionTime[]-start];*)
ntem=ApplyPoleToGram[poles[[i,1]],gram,krep]/.{Subscript[x,_]->1};
(*Print["time after ",SessionTime[]-start];*)
int=Intersection[tem,ntem];
If[int=!={}&&!AllTrue[int,PerfectSquareQ],AppendTo[path,poles[[i,2]]]],
Continue[]
]
,{i,1,Length[poles]}];
Return[path//Flatten[#,1]&//DeleteDuplicates];
];


MonomialQ[exp_]:=Module[{var},
var=Variables[exp];
If[Length[var]==1&&Length[exp]==0,Return[True],Return[False]];
];

Options[SelectAlgLetter]={deBug->False};
SelectAlgLetter[alg_,rl_,OptionsPattern[]]:=Module[{tem,int,select={}},
Do[
	tem=FactorList[(Numerator[#[[1]]]*Denominator[#[[1]]]&[alg[[i]]])//Expand][[All,1]]//DeleteCases[#,_?NumericQ]&;
	int=RIntersection[tem,rl];
	If[OptionValue[deBug]==True,Print["int: ",int];Print["tem: ",tem]];
	If[Complement[tem,int[[1]]]=!={},Continue[],AppendTo[select,alg[[i]]]];
,{i,1,Length[alg]}];
Return[select];
];


GramPCQ[glist_]:=Module[{list,l1,l2},
list={glist}//Cases[#,G[__],Infinity]&//DeleteDuplicates;
If[Length[list]!=2,Message[GramPCQ::err];Return[$Failed]];
l1=list[[1,1]];
l2=list[[2,1]];
If[Complement[Union[l1,l2],l1]==={}||Complement[Union[l1,l2],l2]==={},Return[True]];
Return[False];
];


MyContainsAny[a_,list_]:=If[AnyTrue[list,ContainsAll[a,#]&],True,False];


Options[MergePoles]={deBug->False};
MergePoles[poleslist_,OptionsPattern[]]:=Module[{rep1,rep2,tem1,tem2,isp,int,intv},
	(*poleslist should be in the form {{{Subscript[x, 1]->0,Subscript[x, 2]->0,...},{Subscript[x, 6]->Infinity,...}},{...}}*)
	rep1=poleslist[[All,1]]//Flatten//DeleteDuplicates//Sort;
	tem1=Cases[poleslist[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//Flatten//DeleteDuplicates;(*Note that the second part has order so we can't sort them*)
	tem2=DeleteCases[poleslist[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//DeleteDuplicates;
	If[Length[tem2]==0,(*if there are no infinity poles*)
		rep2=tem1,
		(*if there are infinity poles*)
		tem2=Table[If[Length[tem2[[i]]]==1,Append[tem2[[i]],Rule[tem2[[i,1,1]],0]],tem2[[i]]],{i,1,Length[tem2]}];(*complete some single poles like {x->Infinity} to {x->Infinity,x->0}*)
		If[Length[tem2]==1,(*if there is only one infinity pole map*)
			isp=Cases[tem2,Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
			rep2=Join[tem1,{Thread@Rule[isp,isp/Subscript[x, 0]],Subscript[x, 0]->0},DeleteCases[tem2[[1]],Rule[_,Infinity]]//Delete[#,-1]&],
			(*if there are more than one infinity pole map*)
			isp=(Cases[#,Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates)&/@tem2;
			tem2=(DeleteCases[#,Rule[_,Infinity]]&)/@tem2//Flatten;(*the replacement rule after taking x0 to 0*)
			int=Reap[
					Do[
						Do[
							intv=Intersection[isp[[i]],isp[[j]]];
							If[intv=!={},Sow[intv]]
						,{j,i+1,Length[isp]}]
					,{i,1,Length[isp]}]
			][[2]];
			If[int=!={},int=int//Flatten//DeleteDuplicates];
			If[OptionValue[deBug],Print["int: ",int];Print["isp: ",isp]];
			isp=isp//Flatten//DeleteDuplicates;
			rep2=Join[tem1,{Thread@Rule[isp,isp/Subscript[x, 0]],Subscript[x, 0]->0},DeleteCases[tem2,_?(MemberQ[int,Keys[#]]&)],(Thread@Rule[int,int/.tem2])]//Delete[#,-1]&;
		]
	];
	Return[{rep1,rep2}];
	(*result=Table[
		tem1=remain[[edge[[i]]]];
		tem=DeleteCases[tem1[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//ReverseSortBy[#,Length]&;
		tem2=Cases[tem1[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//Flatten//DeleteDuplicates;
		If[tem==={},(*if there are no infinity poles*)
			{tem1[[All,1]]//Flatten//DeleteDuplicates,tem2//Flatten//DeleteDuplicates},
			{Join[tem1[[All,1]],{Thread@Rule[Keys[tem2],Values[tem2]//.tem2//Factor]}]//Flatten//DeleteDuplicates,If[Length[tem]==1,tem//Flatten,Riffle[tem[[1]],tem[[2]]]//DeleteDuplicates]}
		]
		,{i,1,Length[edge]}];
		(*here we define how to merge two compatible pole*)
	Return[Join[remain[[Complement[Range[Length[remain]],edge//Flatten//DeleteDuplicates]]],result]//DeleteDuplicates];
	(*we include those poles which are not compatible to any other poles separately*)
	(*the final data form will be {{{a->0,b->0},{c->1}},{...},...}*)*)
];


Options[ApplyPoleToVar]={deBug->True};
ApplyPoleToVar[poles_,var_]:=Module[{tem},
	tem=poles//Flatten//DeleteDuplicates//GatherBy[#,(Last[#]===Infinity)&]&;
	If[Length[tem]==1,Return[var/.tem]];
	Return[Table[If[MemberQ[Keys[tem[[2]]],var[[i]]],g[var[[i]]/.tem[[1]]],var[[i]]/.tem[[1]]],{i,1,Length[var]}]];
];


Options[FindPoleMaps]={deBug->False};
FindPoleMaps[poles_,var_,sector_:{},OptionsPattern[]]:=Module[{alist={},xl,tem,atem={},alisttem,l=2,max,subset,n,pos},
	(*poles consists of three parts: poles, graph, cliques of the graph.*)
	(*we find the minimal pole maps that can cover var*)
	(*minimal pole map means this map only involve the variables in var with extra variables as less as possible*)
	(*two poles can be merged only when they share the same cut*)
	If[var==={},Return[{}]];
	xl={Thread@Rule[Subscript[x,#]&/@sector,0],{}};
	tem=poles[[1]];
	alist=Reap[
			Do[
				If[ContainsAll[Keys[tem[[i]]//Flatten],var]&&CompatiblePoleQ[tem[[i]],xl],Sow[i]]
			,{i,1,Length[tem]}]
	][[2]];(*all the single pole maps that can cover var*)
	If[OptionValue[deBug],Print["alist: ",alist]];
	If[alist=!={},alist=alist[[1]]];
	tem=DeleteCases[poles[[3]],_?(ContainsAny[#,alist]&)];(*remove cliques that already contain poles in alist*)
	(*then we search in these remaining cliques for the minimal cover of var*)
	(*Apply every possible cases to the same set of variables*)
	tem=Join[List/@alist,tem];
	atem=Reap[
		Do[
			If[ContainsAll[Keys[poles[[1]][[tem[[i]]]]//Flatten],var],Sow[{ApplyPoleToVar[poles[[1]][[tem[[i]]]],var],tem[[i]]}]];
		,{i,1,Length[tem]}]
	][[2]];
	If[atem=!={},atem=atem[[1]]];
	(*then collect all different maps and remove duplicate cases*)
	atem=GatherBy[atem,First];
	Return[atem[[All,1,2]]];(*for each group, we only need the first element, the remaining map will give the same results*)
	
	(*max=poles[[3,-1]]//Length;(*maximal length of clique*)
	atem=tem;
	While[l<=max,
		(*l is number of poles we choose to form a cover. We start from the 2-pole cover*)
		alisttem={};
		pos={};
		Do[
			If[MyContainsAny[atem[[j]],alisttem],Continue[]];(*if this term is already a cover of the covers we have found*)
			If[!AllTrue[poles[[1]][[atem[[j]]]],CompatiblePoleQ[#,xl]&],AppendTo[pos,j];Continue[]];(*if in the pole map, the required constraint not be satisfied*)
			n=Length[atem[[j]]];(*length of this clique*)
			If[n==l,(*if the length of clique equals to the number n of the n-pole cover we want*)
				If[ContainsAll[Keys[poles[[1]][[atem[[j]]]]//Flatten],var],AppendTo[alisttem,atem[[j]]]],
				subset=Subsets[atem[[j]],{l,n}];(*consider all its subset with at least length l*)
				Do[
					If[ContainsAll[Keys[poles[[1]][[subset[[k]]]]//Flatten],var],
						AppendTo[alisttem,subset[[k]]];Break[](*if one already finds a subset which is a cover, then it is not necessary to search the rest of them*)
					]
				,{k,1,Length[subset]}]
			]
		,{j,1,Length[atem]}];
		atem=Delete[atem,List/@pos];(*remove candidates that are not compatible with constraints*)
		If[alisttem=!={},
			atem=DeleteCases[atem,_?(MyContainsAny[#,alisttem]&)];(*remove those terms that is already a cover of the covers we have found*)
		];
		l=l+1;
		alist=Join[alist,alisttem]//DeleteDuplicates;
		If[OptionValue[deBug],Print["length: ",l-1," alisttem: ",alisttem]];
		atem=DeleteCases[atem,_?(Length[#]<l&)];
	];*)
	(*alist is all the maps we find for the set var*)
	(*then we merge multiple pole map to single pole map*)
	(*Return[Table[If[Head[alist[[i]]]===List,{MergePoles[poles[[1]][[alist[[i]]]]],alist[[i]]},{poles[[1,alist[[i]]]],{alist[[i]]}}],{i,1,Length[alist]}]];*)
	Return[alist];
];


Options[ApplyPolesToAlgLetter1]={Sector->{},AdRep->{},deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True,NoPathInfo->True};
ApplyPolesToAlgLetter1[opoles_,letters_,reform_,krep_,OptionsPattern[]]:=Module[{start,poles,polepath,letterpath,kletters={},tem,ntem,glist,glistd,Glist,ll,nletters={},sletters={},var,assoc,rep,rep1,gr,ga,isp,pos,path,dis},
start=SessionTime[];
(*first find all possible pole maps for the algebraic letters*)
tem=Reap[
	Do[
		Sow[(letters[[i,1]]//Cases[#,_G,Infinity]&//DeleteDuplicates)/.{G[a_,b_]:>GramMat[a,b,krep]}//Cases[#,Subscript[x,_],Infinity]&//DeleteDuplicates//Sort];(*all Baikov variables involved in an algebraic letter*)		
	,{i,1,Length[letters]}]
][[2]];
If[var=!={},var=tem[[1]]//DeleteDuplicates];
assoc=Association[Table[Rule[var[[i]],FindPoleMaps[opoles,var[[i]]]],{i,1,Length[var]}]];(*the association that associate a baikov variable list to a list of pole maps*)
If[OptionValue[deBug],Print["all poles map found ",SessionTime[]-start];Print[Keys[assoc]]];

Do[
	If[OptionValue[deBug],Print[SessionTime[]-start]];
	(*letterpath=letters[[i,3]];*)
	ga=((Cases[letters[[i,1]],Power[z_,1/2]->z,Infinity]//DeleteDuplicates)[[1]])/.{G[{},{}]->1};(*algebraic part of letter*)
	gr=(letters[[i,1]][[1]]//Numerator)/.{Power[_,1/2]->0};(*rational part of letter*)
	tem=({gr,ga}//Cases[#,_G,Infinity]&)/.{G[a_,b_]:>GramMat[a,b,krep]};
	var=Cases[tem,Subscript[x,_],Infinity]//DeleteDuplicates//Sort;(*get all Baikov variables involved in this letter*)
	glist=Cases[{ga},_G,Infinity]//DeleteDuplicates;(*grams under square root*)
	glistd=Cases[{letters[[i,2]]},_G,Infinity]//DeleteDuplicates;(*grams corresponding to the multiplication of numerator and denominator of algebraic letter*)
	If[OptionValue[deBug],Print["ga: ",ga];Print["gr: ",gr];Print["var: ",var]];
	
	(*-------------------the following is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
	(*If[OptionValue[PathDis],
		var=glist/.{G->List}//Variables;(*the variables involved*)
		If[OptionValue[deBug],Print["glist: ",glist];Print["glistd: ",glistd]];
		ll=(Length[#[[1]]]&/@Join[glist,glistd]);
		Glist=Select[reform,(MemberQ[ll,Length[#[[1,1]]]]&&ContainsAll[var,Variables[#[[1,1]]]])&];(*select possible relevant Gram from the whole list*)
		pos=FindGram[Join[glist,glistd],Glist[[All,1]],krep];
		If[!FreeQ[pos,{}],Continue[]];
		(*if one gram is not in above list, then we don't use this ansatz*)
		(*------------------find the path of grams which are used to construct this letter----------------*)
		If[Abs[letters[[i,4]]]==2,
			(*second type construction: 1/(ab)*)
			path=(Union@@(Glist[[#,3]]))&/@Drop[pos,Length[glist]],
			(*first type construction: 1/Sqrt[ab]*)
			path=(Union@@(Glist[[#,3]]))&/@Take[pos,Length[glist]]
		];
		(*------------------find the path of grams which are used to construct this letter----------------*)
		If[OptionValue[deBug],Print["pos: ",pos];Print["path: ",path]];
	];*)
	(*-------------------above is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
	
	If[OptionValue[deBug],Print[SessionTime[]-start]];
	If[var==={},
		(*when this letter consists only of pure kinematics gram*)
		
		(*-------------------the following is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
		(*If[OptionValue[PathDis]&&OptionValue[KinePath],
			(*consider two different types of letters separately*)
			(*Print["letter: ",letters[[i]]];*)
			If[Abs[letters[[i,4]]]==1,
				dis=If[Length[path]==1,0,PathsDistance[path[[1]],path[[2]]]];
				If[dis>2,AppendTo[sletters,letters[[i]]];Continue[],If[dis==2&&!GramPCQ[glist],(*Print["letter: ",letters[[i]]];*)AppendTo[sletters,letters[[i]]];Continue[]]];
			];
			If[Abs[letters[[i,4]]]==2,
				dis=PathsDistance[path[[1]],path[[2]]];
				If[dis>2,AppendTo[sletters,letters[[i]]];Continue[],If[dis==2&&!GramPCQ[glist],(*Print["letter: ",letters[[i]]];*)AppendTo[sletters,letters[[i]]];Continue[]]];
			]
		];*)
		(*-------------------above is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
		
		If[!AdmissiblePoleQ[glist,{},OptionValue[PermSq],krep],Continue[],AppendTo[kletters,{letters[[i,1]],{},letters[[i,4]]}]];
		Continue[],
		
		(*when this letter consists of gram which is not free from loop moemnta*)
		(*If[OptionValue[PathDis]&&PathsDistance[path[[1]],path[[2]]]>1,AppendTo[sletters,letters[[i]]];Continue[]];*)
		
		(*-------------------the following is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
		(*If[OptionValue[PathDis],
			(*consider two different types of letters separately*)
			If[Abs[letters[[i,4]]]==1,
				(*note that there is a special case where one Gram is G[{},{}], in this case path has only one element.*)
				If[Length[path]==1,dis=0,dis=PathsDistance[path[[1]],path[[2]]]],
				If[Abs[letters[[i,4]]]==2,dis=PathsDistance[path[[1]],path[[2]]];]
			];
			If[OptionValue[LoopPath]&&dis>1,
				AppendTo[sletters,letters[[i]]];Continue[],
				If[dis>2,AppendTo[sletters,letters[[i]]];Continue[],If[dis==2&&!GramPCQ[glist],AppendTo[sletters,letters[[i]]];Continue[]]]
			];
			(*If[letters[[i,4]]==1&&PathsDistance[path[[1]],path[[2]]]>1,AppendTo[sletters,letters[[i]]];Continue[]];
			If[letters[[i,4]]==2&&PathsDistance[{letterpath[[1]]},letterpath[[2]]]>1,AppendTo[sletters,letters[[i]]];Continue[]]*)
		];*)
		(*-------------------above is a historical part, we have abandoned putting on the pathdistance constraint for the letters----------------*)
		
		(*ntem=GramMat[#[[1]],#[[2]],krep]&/@glist;
		pos=Table[If[FreeQ[ntem[[k]],x],{k}],{k,1,Length[ntem]}]//DeleteCases[#,Null]&;(*delete Gram which doesn't depend on Baikov variables*)
		path=Delete[path,pos];*)
		
		pos=assoc[var];
		poles=Table[If[Head[pos[[i]]]===List,{MergePoles[opoles[[1]][[pos[[i]]]]],pos[[i]]},{MergePoles[{opoles[[1,pos[[i]]]]}],{pos[[i]]}}],{i,1,Length[pos]}];(*the pole maps needed*)
	];
	If[OptionValue[PermSq]==={},
		(*in this case, we will search letters by relations between poles and letter ansatz. This approach is not preferred now*)
		Do[
			polepath=poles[[j,2]]//SortBy[#,Length]&;
			(*polepath=RefinedPolePath[glist,poles[[j]],poles,krep]//SortBy[#,Length]&;*)(*some poles are special and they can appear in subsectors*)
			If[AllTrue[polepath,ContainsAny[#,OptionValue[Sector]]&],(*if there is any element in the path, it means this pole only appear in subsectors, so if we only want pole appearing top sector they will be discarded*)
				Continue[],
				If[AnyTrue[polepath,(*MemberQ[letterpath[[2]],#]&*)(MemberQ[Join[letterpath[[1]],letterpath[[2]]],#])&],
				(*If[!AnyTrue[path,IntersectingQ[#,polepath]&],Continue[]];*)
					If[!AllTrue[path,CompatiblePathQ[#,polepath,TopOnly->True]&],Continue[]];
					AppendTo[nletters,{letters[[i,1]],poles[[j,1]],letters[[i,4]]}];
				]
			]
		,{j,1,Length[poles]}],
		(*In this case, we use information of square roots to constrain the form*)
		Do[
			If[!AdmissiblePoleQ[glist,poles[[j,1]],OptionValue[PermSq],krep],
				Continue[],
				If[!OptionValue[NoPathInfo],
					(*if we put on the information of pole's path*)
					polepath=poles[[j,2]]//SortBy[#,Length]&;
					If[!AllTrue[path,CompatiblePathQ[#,polepath]&],Continue[]]
				];
				AppendTo[nletters,{letters[[i,1]],poles[[j,1]],letters[[i,4]]}]
			]
		,{j,1,Length[poles]}]
	];
	If[OptionValue[deBug],Print[SessionTime[]-start]]
,{i,1,Length[letters]}];
Return[{nletters,kletters,sletters}];
];


GetQuadMatrix[exp_]:=Module[{xl,tem,d,Q},
xl=Cases[exp,Subscript[x,_],Infinity]//DeleteDuplicates//Sort;
xl=Append[xl,Subscript[x, 0]];
d=Length[xl];
tem=Numerator[exp/.{Subscript[x,i_]:>Subscript[x,i]/Subscript[x, 0]}//Together];
Q=Table[If[i==j,Coefficient[tem,xl[[i]],2]//Factor,Coefficient[tem,xl[[i]]*xl[[j]]]/2//Factor],{i,1,d},{j,1,d}];
Return[{Q,xl}/.{Subscript[x, 0]->1}];(*here x0 is just an auxiliary variable*)
];


MixList[list1_,list2_]:=Module[{l,r},
l=Min[Length[list1],Length[list2]];
r=Flatten/@Partition[Riffle[Take[list1,l],Take[list2,l]],2];(*make sure every element is a replacement list*)
If[Length[list2]>l,r=Join[r,Drop[list2,l]],r=Join[r,Drop[list1,l]]];
If[FreeQ[r,x]||FreeQ[r,y],r=r//Flatten];
Return[r//DeleteCases[#,{}]&];
];


Options[ApplyPoleToQM]={deBug->False,AdRep->{}};
ApplyPoleToQM[poles_,Gram_,krep_,OptionsPattern[]]:=Module[{QM,tem,ntem,rep,rep1,isp,xll,pole,result={},rpolesx,rpolesy,pos},
If[Length[poles]<2,Message[ApplyPoleToQM::err,poles];Return[$Failed]];

rpolesx=Table[{},{i,1,Length[poles]}];
rpolesy=Table[{},{i,1,Length[poles]}];
QM=GetQuadMatrix[Det[Gram/.{G[a_,b_]:>GramMat[a,b,krep]}]];
If[OptionValue[deBug],Print["QM: ",QM]];
xll=Table[QM[[2]]/.poles[[i,1]],{i,1,Length[poles]}];
Do[
	pole=poles[[j]];
	Do[
		If[!FreeQ[pole[[k]],Infinity],
			rep=Split[pole[[k]],(Last[#1]===Last[#2])&];
			ntem=Table[If[FreeQ[rep[[i]],Infinity],rep[[i]],{rep[[i]],Subscript[x, 0]->0}/.{Rule[Subscript[x,a_],Infinity]:>Rule[Subscript[x,a],Subscript[x,a]/Subscript[x,0]]}],{i,1,Length[rep]}]//Flatten[#,1]&;
			isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;(*keep track of variables taken to infinity*)
			(*Do[
				pos=Position[Keys[ntem],isp[[i]]];
				If[Length[pos]>1,Continue[],If[ntem[[pos[[1,1]]+1,1]]=!=Subscript[x, 0],ntem=Insert[ntem,Subscript[x, 0]->0,pos[[1,1]]+1]]];
			,{i,1,Length[isp]}];(*it is different whether taking variables to infinity one by one or taking them to infinity at one time. *)*)
			If[OptionValue[deBug],Print["ntem: ",ntem];Print["poles: ",pole[[k]]];Print["isp: ",isp]];
			(*rep=Thread@Rule[isp,isp/Subscript[x, 0]];
			ntem=DeleteCases[pole[[k]],Rule[_,Infinity]];*)
			rep={};
			rep1=If[ntem=!={},If[MemberQ[isp,ntem[[-1,1]]],Drop[ntem,-1],Print["not projective curve!: ",ntem];ntem],{}]
			(*when projecting to infinity space, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*),
			
			pos=Position[pole[[k]],{__},1];
			If[pos==={},rep=pole[[k]];rep1={Subscript[x, 0]->1},rep={pole[[k]][[pos[[1,1]]]],Subscript[x, 0]->0};rep1=Delete[pole[[k]],pos[[1,1]]]];
		];
		rpolesx[[j]]=Join[rpolesx[[j]],rep,rep1];
		rpolesy[[j]]=rpolesx[[j]]/.{x->y};
		(*xll[[j]]=xll[[j]]//.rep//.rep1*)
	,{k,2,Length[pole]}]
,{j,1,Length[poles]}];
(*Now we get all the possible xll*)
(*tem=Partition[Riffle[xll,poles],2]//DeleteDuplicatesBy[#,First]&;
If[OptionValue[deBug],Print["tem: ",tem]];
xll=tem[[All,1]];
rpoles=tem[[All,2]];*)
Do[
	Do[
		tem=SpecialSimplify[{xll[[i]] . QM[[1]] . (xll[[j]]/.{x->y}),(xll[[i]] . QM[[1]] . xll[[i]])*((xll[[j]]/.{x->y}) . QM[[1]] . (xll[[j]]/.{x->y}))},AdRep->Join[MixList[rpolesx[[i]],rpolesy[[j]]],OptionValue[AdRep]]];
		If[((tem[[1]]-Sqrt[tem[[2]]]//Factor)===0)||((tem[[1]]+Sqrt[tem[[2]]]//Factor)===0),Continue[]];
		AppendTo[result,{Log[(tem[[1]]+Sqrt[tem[[2]]])/(tem[[1]]-Sqrt[tem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&],Gram,{poles[[i]],poles[[j]]}}]
	,{j,i+1,Length[xll]}]
,{i,1,Length[xll]}];
Return[result];
];


Options[ApplyPolesToAlgLetter2]={AdRep->{},deBug->False,PermSq->{}};
ApplyPolesToAlgLetter2[opoles_,ogram_,reform_,krep_,OptionsPattern[]]:=Module[{start,gram,polepath,poles,tem,ntem,var,glist,Glist,ll,nletters={},pos,path,adpole={},assoc},
start=SessionTime[];
If[OptionValue[deBug],Print["Calculation starts!"]];

If[ogram==={},gram=reform[[All,1]]//DeleteDuplicates[#,EquivalentGramQ[#1,#2,krep]&]&,gram=ogram//DeleteDuplicates[#,EquivalentGramQ[#1,#2,krep]&]&];
(*when ogram is {}, we directly calculate all possible letter from reform. We will delete those equivalent Grams*)
ntem=Reap[
	Do[
		tem=gram[[i]]/.{G[a_,b_]:>GramMat[a,b,krep]};
		If[FreeQ[tem,x],Continue[]];
		If[Exponent[Det[tem/.{Subscript[x,e_]:>Subscript[x,e]*R}],R]>2,Continue[]];(*polynomial not a quadratic cannot be used in the second ansatz*)
		Sow[gram[[i]],a];(*keep those grams satisfying constraints*)
		var=Cases[tem,Subscript[x,_],Infinity]//DeleteDuplicates//Sort;
		Sow[var,b](*the Baikov variables of corresponding gram*)
	,{i,1,Length[gram]}]
][[2]];
If[ntem=!={},gram=ntem[[1]];var=ntem[[2]]//DeleteDuplicates];
assoc=Association[Table[Rule[var[[i]],FindPoleMaps[opoles,var[[i]]]],{i,1,Length[var]}]];(*the association that associate a baikov variable list to a list of pole maps*)
If[OptionValue[deBug],Print["var: ",Short[Keys[assoc],10]," length of grams: ",Length[gram]]];
Do[
	adpole={};
	glist=Cases[{gram[[i]]},_G,Infinity]//DeleteDuplicates;
	var=Cases[gram[[i]]/.{G[a_,b_]:>GramMat[a,b,krep]},Subscript[x,_],Infinity]//DeleteDuplicates//Sort;
	(*If[OptionValue[deBug],Print["glist: ",glist]];*)
	(*If[OptionValue[deBug],Print["gram: ",glist[[1]]]];*)
	If[OptionValue[PermSq]==={},
		ll=(Length[#[[1]]]&/@glist);
		Glist=Select[reform,MemberQ[ll,Length[#[[1,1]]]]&];
		pos=FindGram[glist,Glist[[All,1]],krep];
		path=(Union@@(Glist[[#,3]]))&/@pos;(*paths for grams in the square root*)
		Do[
			polepath=Complement[xl,Keys[Flatten[opoles[[1,j]]]]]//ReverseSort;
			If[IntersectingQ[path[[1]],polepath],AppendTo[adpole,opoles[[1,j]]]];
			(*get all admissible poles that can be substituted into the gram matrix*)
		,{j,1,Length[opoles[[1]]]}],
		(*In this case, we use additional information to constrain the form*)
		pos=assoc[var];
		poles=Table[If[Head[pos[[i]]]===List,{MergePoles[opoles[[1]][[pos[[i]]]]],pos[[i]]},{MergePoles[{opoles[[1,pos[[i]]]]}],{pos[[i]]}}],{i,1,Length[pos]}];(*the pole maps needed*)
		Do[
			If[AdmissiblePoleQ[glist,poles[[j,1]],OptionValue[PermSq],krep],AppendTo[adpole,poles[[j,1]]]]
		,{j,1,Length[poles]}]
	];
	(*If[OptionValue[deBug],Print["adpole: ",adpole]];*)
	If[Length[adpole]<2,Continue[]];(*if there are less than 2 admissible poles, then we can't construct from this gram*)
	AppendTo[nletters,ApplyPoleToQM[adpole,gram[[i]],krep]//DeleteDuplicatesBy[#,First]&//DeleteCases[#,_?(FreeQ[#,Power[_,1/2]]&)]&];
,{i,1,Length[gram]}];
Return[nletters//Flatten[#,1]&//DeleteDuplicatesBy[#,First]&];
];


Options[AllAlgLetters]={deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True};
AllAlgLetters[poles_,algletter_,reform_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,eresult2,a,b,permsq},
start=SessionTime[];
l=Length[algletter];
permsq=OptionValue[PermSq];
permsq=Join[permsq,Table[Times@@(PerfectSquareSplit[permsq[[i]]*permsq[[j]]][[2]]),{i,1,Length[permsq]},{j,i+1,Length[permsq]}]//Flatten]//DeleteDuplicates;
Print["totally ",l," sectors need to be analyzed"];
Print["analyzing first type construction..."];
Monitor[
	Do[
		If[OptionValue[deBug],PrintTemporary["subset: ",algletter[[a,2]]," session time: ",SessionTime[]-start]];
		tem=ApplyPolesToAlgLetter1[poles,algletter[[a,1]],reform,krep,PermSq->permsq,PathDis->OptionValue[PathDis],KinePath->OptionValue[KinePath],LoopPath->OptionValue[LoopPath](*,deBug->If[algletter[[a,2]]==={3,4,5,6,7},True,False]*)];
		ntem=Flatten[tem[[{1,2}]],1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;
		(*delete duplicates according to the first (Letter form) and second (poles info) term of the unit.*)
		AppendTo[result,ntem];
		If[OptionValue[deBug],AppendTo[spresult,tem[[3]]]]
	,{a,1,Length[algletter]}]
,ProgressIndicator[a,{1,Length[algletter]}]];
Print["session time: ",SessionTime[]-start];
Print["Substituting poles into expressions..."];
tem=Flatten[result,1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;(*remove duplicate expression again*)
eresult=Monitor[Table[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["session time: ",SessionTime[]-start];
Print["analyzing second type construction..."];
eresult2=ApplyPolesToAlgLetter2[poles,{},reform,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{eresult,eresult2,spresult}]];
Return[{eresult,eresult2}];
];


Options[AllAlgLettersPL]={deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True};
AllAlgLettersPL[poles_,algletter_,reform_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,eresult2,part,pathdis,permsq,looppath,kinepath,b},
start=SessionTime[];
l=Length[algletter];
permsq=OptionValue[PermSq];
permsq=Join[permsq,Table[Times@@(PerfectSquareSplit[permsq[[i]]*permsq[[j]]][[2]]),{i,1,Length[permsq]},{j,i+1,Length[permsq]}]//Flatten]//DeleteDuplicates;
Print["totally ",l," sectors need to be analyzed"];
Print["analyzing first type construction..."];
(*Export[OptionValue[tmpDir]<>"input.mx",{poles,reform,krep}];
part=Partition[algletter,UpTo[Quotient[Length[algletter],OptionValue[NThreads]]+1]];
Table[Export[OptionValue[tmpDir]<>"list"<>ToString[i]<>".mx",part[[i]]],{i,1,Length[part]}];
(*export the data to parallelly run them in terminal*)*)
pathdis=OptionValue[PathDis];
looppath=OptionValue[LoopPath];
kinepath=OptionValue[KinePath];
DistributeDefinitions[poles,algletter,reform,krep,permsq,pathdis,looppath,kinepath,ApplyPolesToAlgLetter1,AllAlgLettersPL];
SetSharedVariable[result,spresult];
ParallelDo[
(*delete duplicates according to the first (Letter form) and second (poles info) term of the unit.*)
	AppendTo[result,ApplyPolesToAlgLetter1[poles,algletter[[i,1]],reform,krep,PermSq->permsq,PathDis->pathdis,KinePath->kinepath,LoopPath->looppath]//Flatten[#[[{1,2}]],1]&//DeleteDuplicatesBy[#,#[[{1,2}]]&]&];
	If[OptionValue[deBug],AppendTo[spresult,tem[[3]]]]
,{i,1,Length[algletter]}];
Print["session time: ",SessionTime[]-start];
Print["Substituting poles into expressions..."];
tem=Flatten[result,1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;(*remove duplicate expression again*)
(*DistributeDefinitions[tem,ApplyPoleToLetter];*)
(*Return[tem];*)
eresult=Monitor[Table[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
(*eresult=ParallelTable[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]},Method->"CoarsestGrained"];*)
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["session time: ",SessionTime[]-start];
Print["analyzing second type construction..."];
eresult2=ApplyPolesToAlgLetter2[poles,{},reform,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{eresult,eresult2,spresult}]];
Return[{eresult,eresult2}];
];


Options[SameAlgLetterQ]={deBug->False};
SameAlgLetterQ[a1_,a2_,OptionsPattern[]]:=Catch@Module[{start,sq1,sq2,r1,r2,var1,var2,numrep},
start=SessionTime[];
var1=Variables[a1[[1]]]//Sort;
var2=Variables[a2[[1]]]//Sort;
If[OptionValue[deBug],Print["var1,var2: ",{var1,var2}]];
If[var1=!=var2,Throw[False]];
sq1=Cases[a1,Power[_,1/2],Infinity]//DeleteDuplicates;
sq2=Cases[a2,Power[_,1/2],Infinity]//DeleteDuplicates;
If[OptionValue[deBug],Print[SessionTime[]-start]];
If[sq1=!=sq2,Throw[False]];(*if square root is not the same, then it is not the same*)
(*If[sq1==={}||sq2==={},Throw[False]];
If[var1==={},Throw[False]];*)
Do[
numrep=If[var1=!={},Thread@Rule[var1,RandomPrime[{100,1000},Length[var1]]],{}];
r1=Times@@NumeratorDenominator[a1[[1]]]/.numrep//Expand;
r2=Times@@NumeratorDenominator[a2[[1]]]/.numrep//Expand;
If[OptionValue[deBug],Print["r1,r2: ",{r1,r2}]];
If[(r1+r2//Factor)=!=0&&(r1-r2//Factor)=!=0,Throw[False]]
,{i,1,3}];
If[OptionValue[deBug],Print[SessionTime[]-start]];
Throw[True];
];

GetCoord[exp_]:=Module[{sq,r},
sq=Cases[exp,Power[_,1/2],Infinity]//DeleteDuplicates;
r=Times@@NumeratorDenominator[exp]//Expand;
Return[{sq,r}];
];


Options[DeleteSameAlgLetter]={deBug->False};
DeleteSameAlgLetter[alglist_,OptionsPattern[]]:=Module[{list,numrep,var,nlist,coord,den,dennum,pos},
list=alglist/.{Log[z_]:>z};
var=Variables[list]//Sort;
numrep=Table[Thread@Rule[var,Table[Prime[300*i+3*j],{j,1,Length[var]}]],{i,1,3}];
den=Denominator/@list;
dennum=Table[den/.numrep[[i]]//Factor,{i,1,3}]//Transpose;
pos=Position[dennum,{0,0,0}];
If[OptionValue[deBug],Print["pos: ",pos]];
list=Delete[list,pos];(*remove terms that are actually singular, this can happen in degenerate case which is common in massless one-loop n-gon*)
If[OptionValue[deBug],Print["numrep: ",numrep]];
coord=Table[GetCoord/@(list/.numrep[[i]]//Factor),{i,1,3}]//Transpose//Quiet;
nlist=Partition[Riffle[coord,list],2];
nlist=DeleteDuplicatesBy[nlist,First];
Return[Log/@(nlist[[All,2]])];
];


ToRational[alg_]:=Times@@NumeratorDenominator[alg[[1]]]//Expand//Factor;

GenerateNumReal[sq_,var_]:=Module[{rep,tem,k=1},
rep=Thread@Rule[var,RandomChoice[{1,-1}]*RandomPrime[{1000,10000},Length[var]]];
tem=sq/.{Power[a_,1/2]:>a};
While[Length[Values[rep]//Union]<Length[var]||((tem/.rep)<0&&k<1000),
rep=Thread@Rule[var,RandomChoice[{1,-1}]*RandomPrime[{1000,10000},Length[var]]];
k=k+1;
];
Return[rep];
];

SearchIndepLetterNum[letter_,repnum_]:=Module[{lnum,ind={},tem},
lnum=letter/.{Log[a_]:>Log[Abs[a/.repnum]]};
AppendTo[ind,1];
Do[
tem=Quiet[FindIntegerNullVector[Append[lnum[[ind]],lnum[[i]]],30]];
If[Head[tem]=!=List,AppendTo[ind,i],Continue[]];
,{i,2,Length[letter]}];
Return[ind];
];


Options[GetAllAlgIndepLetter]={deBug->False};
GetAllAlgIndepLetter[algletter_,OptionsPattern[]]:=Module[{gr,var,sq,rep,pos,tem,ind={}},
gr=GatherBy[algletter,(Cases[#,Power[_,1/2],Infinity]//DeleteDuplicates)&];
Print["totally ",Length[gr]," different square roots"];
var=algletter/.{Log[z_]:>z}//Variables;
Do[
sq=Cases[gr[[i]],Power[_,1/2],Infinity]//DeleteDuplicates;
If[sq==={},Continue[]];
tem=gr[[i]]//SortBy[#,LeafCount]&;
If[OptionValue[deBug],Print["square root: ",Short[sq,20]];Print["length: ",Length[tem]]];
rep=Table[GenerateNumReal[sq[[1]],var],{i,1,4}];
pos=Commonest[SearchIndepLetterNum[tem,#]&/@rep][[1]];
AppendTo[ind,tem[[pos]]]
,{i,1,Length[gr]}];
Return[ind];
];

Options[FindLetterLinearRelation]={deBug->False};
FindLetterLinearRelation[basis_,target_,OptionsPattern[]]:=Module[{gr,var,sq,rep,rel},
gr=Append[basis,target];
var=gr/.{Log[z_]:>z}//Variables;

sq=Cases[{gr[[-1]]},Power[_,1/2],Infinity]//DeleteDuplicates;
If[OptionValue[deBug],Print["square root: ",Short[sq,20]]];
rep=Table[GenerateNumReal[sq[[1]],var],{i,1,3}];
rel=Commonest[Quiet[FindIntegerNullVector[ReplaceAll[gr/.{Log[z_]:>Log[Abs[z]]},#],30]]&/@rep];
Return[rel];
];


LetterInfo[letter_,algresult_,polestructure_]:=Module[{nl,pos,tem},
If[Head[letter]=!=Log,nl=Log[letter],nl=letter];
If[FreeQ[nl,Power[_,1/2]],
(*this is a rational letter*)
pos=Position[polestructure[[All,3]],nl[[1]]][[All,1]];
Print["This letter can be generated from following sectors: ",polestructure[[pos,-1]]];
Print["The corresponding path is: ",polestructure[[pos,1,All,-1]]];
Print["To see more information, search the following position of output of PolesAnalyze[]: ",pos],
(*then this is an algebraic letter*)
pos=Position[algresult,nl];
tem=pos[[All,1]]//Union;
Switch[tem,
{1},Print["This letter comes from the first type construction."],
{2},Print["This letter comes from the second type construction."],
{1,2},Print["This letter can either come from first type or second type construction."],
{},Print["This letter doesn't come from  this construction!"]];
pos=Append[#,{2,3}]&/@(Take[#,2]&/@pos);
tem=Table[algresult[[Sequence@@pos[[i]]]],{i,1,Length[pos]}]//DeleteDuplicates;
Print["It can come from the following gram and poles: ",Short[tem,100]];
];
];


FindGramFromPoly[poly_,poles_,reform_,krep_]:=Module[{tem,pos={},var,numrep,nkrep},
Do[
	tem=reform[[i,1]]/.{G[a_,b_]:>GramMat[a,b,krep]}//Factor;
	If[FreeQ[tem,Subscript[x,_]],If[NumericQ[Det[tem]/poly//Factor],AppendTo[pos,{{},reform[[i,1]]}]];Continue[]];
	Do[
		tem=(ApplyPoleToGram[poles[[j,1]],{reform[[i,1]]},krep]/.{Subscript[x, _]->1});
		If[tem=!=0&&NumericQ[tem/poly//Factor],AppendTo[pos,{poles[[j,1]],reform[[i,1]]}]];
	,{j,1,Length[poles]}]
,{i,1,Length[reform]}];
var=Variables[Values[krep]];
numrep=Thread@Rule[var,Prime[1000+#]&/@Range[Length[var]]];
nkrep=Thread@Rule[Keys[krep],Values[krep]/.numrep];
pos=pos//GatherBy[#,Gram2Poly[Last[#],nkrep]&]&;
Return[Table[{pos[[i,All,1]]//DeleteDuplicates,pos[[i,1,-1]]},{i,1,Length[pos]}]];
];


RegularizeSquareRoots[list_]:=(If[FreeQ[#,Complex],#/.{Power[a_,1/2]:>Power[a//Expand//Factor,1/2]},#/.{Power[a_,1/2]:>I*Power[-a//Expand//Factor,1/2]}]&)/@list;


CheckPosition[posl_]:=Module[{g},
Print["True: all in one column and one row;\ False: not in the preferred configuration"];
Table[
g=GatherBy[posl[[i]],First]//ReverseSortBy[#,Length]&;
g=Drop[g,1]//Flatten[#,1]&//GatherBy[#,Last]&;
If[Length[g]>1,False,True]
,{i,1,Length[posl]}]
];


End[];


EndPackage[];

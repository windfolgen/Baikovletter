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

PolyDim::usage="PolyDim[poly] returns the mass dimension of a homogeneous polynomial.";
RemoveCoeff::usage="Remove the constant numeric factor from expression."
TakeSquareRoot::usage="TakeSquareRoot[exp] take the square roots of an expression.";
EnhancedFactor::usage="EnhancedFactor[list,z] factor the expression with respect to z.";

ResolveSingleVariable::usage="ResolveSingleVariable[intset,hintset,z] resolves singularities of z from the integer power set: intset and half-integer power set: hintset.";
ResolveSingleVariable::warning="Odd power under square: `1` !";
ResolveSingleVariable::err="In ResolveSingleVariable[intset,hintset,z], intset and hintset can't both be {}";
ResolveSingleVariable::elliptic="Elliptic case encountered! `1`";
ResolveSingleVariable::additional="There may be additional singularities: {poly, poly list under square root} `1`.";

ResolveSingularities::usage="ResolveSingularities[intset,hintset] resolves singularities for all the variables in intset and hintset. It is the multivariate version of ResolveSingleVariable[].";
ResolveSingularities::err="There is something wrong with the arguments. They are all irrelevant to Baikov variables.";
ResolveSingularities::warning="In this case it is equivalent to integrating this variable out, so we won't consider its pole here but will include relevant Landau varieties. `1`.";
ResolveSingularities::oddpower="There is an odd power higher than 1 of variable under square root or it's a linear expression but with no other poles in the denominator! The infinity is a branch cut. Its expression is `1` and the input is `2`.";
ResolveSingularities::elliptic="Elliptic or hyperelliptic case encountered! `1`";
ResolveSingularities::additional="There may be additional singularities: {poly, poly list under square root} `1`.";
RemoveSquareRootsSol::usage="RemoveSquareRootsSol[sol] remove the solutions which contains irreducible square roots.";

Trans2Inf::usage="Trans2Inf[intset,hintset,xl] transforms the polynomial to infinity plane which corresponds to second Landau singularities.";


CheckValidity::usage="CheckValidity[brep,subset] checks whether a reprentation is valid under the maximal cut, if there is no valid representation, it will return {}.";


ExistRelationQ::usage="ExistRelationQ[poly,sqpolylist] detects whether poly has some relation with a list of other polynomials under square roots sqpolylist.";


FindCover::usage="FindCover[list] find a cover of the set of all variables present from a set of sets of these variables.";


RFindInstance::usage="Refined version of FindInstance[]. Same usage as FindInstance[poly,var,domain,n] except that the first argument should be a polynomial.";
TFindInstance::usage="Refined version of FindInstance[] for quadratic polynomials. Internal use only.";


PolesAnalyze::usage="PolesAnalyze[result,topsector,krep,n] gives the possible poles, their paths and leading singularities related. topsector is specified by a number list like {1,2,3,4} where 1-st, 2-nd, 3-rd and 4-th propagators are present in this sector. n is the total number of Baikov variables in a family. Actually, this serves as a limited case of dlog construction. The output will be a list of forms like {{{pole,poly,path},...},polylist1,polylist2,sector}. poles are classified by their sectors.";
PolesAnalyze::misp = "We need consider at least next-to-next-to-maximal cut for this sector: `1` and only `2` of the propagators can be set to 0.";


RIntersection::usage="RIntersection[list1,list2] calculates the intersection of two lists where terms with different signs are taken as the same.";
AllRuleQ::usage="AllRuleQ[list] returns True if all elements in the list are rules.";


ApplyPoleToVar::usage="ApplyPoleToVar[poles,var] applies poles to a list of Baikov variables var.";
ExtractPoleInfo::usage="ExtractPoleInfo[exp,n] extracts the information of pole from the output of PolesAnalyze[] function. n is the total number of Baikov variables in an integral family. There is an option OutputLevel, it is set to 1 by default. When it is 2, the output will include some second type of poles.";

CompatiblePoleQ::usage="CompatiblePoleQ[p1,p2] checks whether two poles p1 and p2 are compatible.";
DecendentQ::usage="DecendentQ[basis,pole] checks whether a pole can be generated from the basis.";
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

PolyMemberQ::usage="PolyMemberQ[polylist,poly] checks whether a polynomial poly is in the list of polylist. Note that they can differ by a numeric constant.";
AdmissiblePoleQ::usage="AdmissiblePoleQ[glist,pole,permsq,krep] decides whether a pole can be substituted into glist. permsq is a list of polynomials the square root of which are leading singularities appearing in the construction of UT basis.";
PathsDistance::usage="PathsDistance[p1,p2] calculate the minimum distance of two path list.";


SelectAlgLetter::usage="SelectAlgLetter[alg,rl] selects algebraic letters in the list alg according to whether this algebraic letter is related to some rational letter in list rl.";


GramPCQ::err="There must be two grams in the arguments";
GramPCQ::usage="GramPCQ[glist] finds whether two grams in glist are parent Gram and child Gram. For example, GramPCQ[{G[{1,3,4}],G[{3,4}]}] will be True";


ApplyPolesToAlgLetter1::usage="ApplyPolesToAlgLetter1[poles,letters,reform,krep] apply all possible poles to the letters constructed. reform is the output of ReformRep[result]. This will give all first type letters. The output consists of three part: nletters (which is the normal letter), kletters (which consists of only pure kinematics grams) and sletters (which are not permissible because the distance between the selection rule we put on)";


ApplyPolesToAlgLetter1::permsq="The option \"PermSq\" can not be empty!";


ApplyPolesToAlgLetter2::permsq="The option \"PermSq\" can not be empty!";


GetQuadMatrix::usage="GetQuadMatrix[exp] get the quadratic matrix from a one-loop baikov polynomial.";
ApplyPoleToQM::usage="ApplyPoleToQM[poles,G,krep] apply two poles to one Q matrix. G is the Gram G[]. The ouput will be {{Log[...],G,{pole1,pole2}},...}.";
ApplyPolesToAlgLetter2::usage="ApplyPolesToAlgLetter2[poles,gram,reform,krep] apply poles to one Gram to construct the second type of letter. When gram={}, reform will give all possible gram. The output will be in the form {{{Log[],Gram1,{pole1,pole2}},{Log[],Gram1,{pole3,pole4}},...},{{Log[],Gram2,{pole1,pole2}},{Log[],Gram2,{pole3,pole4}},...},...}";


AllAlgLetters::usage="AllAlgLetters[poles,algletter,reform,krep] get all possible algbraic letter for a family. algletter is the output of ExtractAlgLetter[]. reform is the output of ReformRep[result]";
AllAlgLettersPL::usage="The parallel version of AllAlgLetters[poles,algletter,reform,krep]";


ApplyPolesToAlgLetter3::usage="ApplyPolesToAlgLetter3[poles,letters,krep] apply the poles with integrated variables set to 0 to letters.";
AllAlgLettersSupplement::usage="AllAlgLettersSupplement[poles,algletter,krep] gives the supplement set of algebraic letters calculated by ApplyPolesToAlgLetter3[];";
AllAlgLettersSupplementPL::usage="AllAlgLettersSupplementPL[poles,algletter,krep] gives the supplement set of algebraic letters calculated by ApplyPolesToAlgLetter3[]. The parallel version of AllAlgLettersSupplement.";


SameAlgLetterQ::usage="SameAlgLetterQ[a1,a2] determines whether a1 and a2 are two equivalent algebraic letters.";
SameGramQ::usage="SameGramQ[g1,g2] decides whether two monomials of Grams equal to each other (or only differ by a minus sign, we take this as the same too).";
SameAlgLetterGramQ::usage="Internal usage! SameAlgLetterGramQ[letter1,letter2] decides whether two algbraic letters with abstract expression of G equals to each other. letter1 and letter2 should be in the form {Log[a+Sqrt[b]/(a-Sqrt[b])],a^2-b,...}.";
GetCoord::usage="GetCoord[exp] gives {sq,r} from an algebraic letter (b+Sqrt[sq])/(b-Sqrt[sq]). r= b^2-sq.";
DeleteSameAlgLetter::usage="DeleteSameAlgLetter[alglist] delete the same algletter in a list.";
RemoveDegenerateAlgLetter::usage="RemoveDegenerateAlgLetter[alglist,krep] delete degenerate algletter under special kinematics which is represented by a replacement rule for kinematic variables, krep."


SearchIndepLetterNum::usage="SearchIndepLetterNum[letter,repnum] finds the independent letter with all its variables have been replaced by random num.";
GetAllAlgIndepLetter::usage="GetAllAlgIndepLetter[alglist] gets all independent algebraic letters in alglist. The result is classified by the square roots appearing in the letters.";


LetterInfo::usage="LetterInfo[letter,algresult,polestructure] show information about all this letter is constructed and where it comes from.";


FindLetterLinearRelation::usage="FindLetterLinearRelation[basis,target] expand the target letter on a set of letter basis.";


FindGramFromPoly::usage="FindGramFromPoly[poly,poles,reform,krep] finds the poles and Grams corresponding to the poly.";


RegularizeSquareRoots::usage="RegularizeSquareRoots[list] makes the square roots in the algebraic letter appearing in a canonical form which is defined by mathematica function Factor[]. This avoids the same square root being taken as different ones due to their different forms caused by the simplification process.";


GetIndepAlgLetters::usage="GetIndepAlgLetter[result_,rationalset] gives all the independent algebraic letters from the result which is the output of AllAlgLettersPL[] or AllAlgLetters[]. rationalset if the set of leading singularities (candidates for rational letters which are used to filter spurious algebraic letters.) ";


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
If[Not@FreeQ[expr,R],(*if the expression contains abbreviated R, then we will not factorize it*)
	If[Head[expr]===Times,list=List@@expr,Return[{{},{expr}}]];
	Do[
		If[MatchQ[list[[i]],Power[_,_?EvenQ]],AppendTo[sl,list[[i]]/.{Power[z_,n_]:>z}],AppendTo[nsl,list[[i]]]]
	,{i,1,Length[list]}];
	Return[{sl,nsl}]
];
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


Options[TakeSquareRoot]={"test"->False};
TakeSquareRoot[exp_,OptionsPattern[]]:=Module[{list,a,s,id,pos},
	list=FactorList[exp];
	If[MemberQ[list[[All,1]],0],Return[0]];(*the square root is actually 0*)
	(*Print["list: ",list];*)
	a=Product[Power[list[[i,1]],Quotient[list[[i,2]],2]],{i,1,Length[list]}];
	s=Product[Power[list[[i,1]],Mod[list[[i,2]],2]],{i,1,Length[list]}];
	If[s=!=1,
		If[!FreeQ[s,x],Return[-1]];
		If[OptionValue["test"],Return[1]];
		If[!MemberQ[Rlist,s],
		Unprotect[R,Rlist];
		AppendTo[Rlist,s];
		id=Length[Rlist];
		R/:R[id]*R[id]=Rlist[[-1]];
		R/:Power[R[id],n_/;(EvenQ[n])]=Power[Rlist[[-1]],n/2];
		R/:Power[R[id],n_/;(OddQ[n])]=Power[Rlist[[-1]],Quotient[n,2]]*R[id];
		Protect[R,Rlist];
		,
		id=Position[Rlist,s][[1,1]]
		];
		s=R[id];
	];
	Return[a*s];
];


EnhancedFactor[intset_,z_]:=Module[{tem,result,pow,sol,sq,rep,var,numrep},
	(*Return[intset];*)
	result=Reap[
	Do[
		If[FreeQ[intset[[i]],z],
			Sow[intset[[i]]];Continue[],
			pow=Exponent[intset[[i]],z];
			If[pow==1||pow>2,
				Sow[intset[[i]]];Continue[],
				var=Variables[intset[[i]]]//DeleteCases[#,_?(!FreeQ[#,x]&)]&;
				If[var=!={},(*we first use a set of numerical value to test which may avoid unnecessary calculation*)
					numrep=Thread@Rule[var,RandomPrime[10000,Length[var]]];
					sol=Solve[(intset[[i]]/.numrep)==0,z];
					sq=Cases[sol,Power[_,1/2],Infinity]//DeleteDuplicates;
					tem=TakeSquareRoot[#,"test"->True]&/@(sq/.{Power[a_,1/2]:>a});
					If[MemberQ[tem,-1],Sow[intset[[i]]];Continue[]];(*in this case, there are still Baikov variables under square roots, we don't consider it*)
				];
				sol=Solve[intset[[i]]==0,z];
				sq=Cases[sol,Power[_,1/2],Infinity]//DeleteDuplicates;
				If[sq==={},Sow[(Power@@@(FactorList[intset[[i]]]))//DeleteCases[#,_?NumericQ]&];Continue[]];
				tem=TakeSquareRoot/@(sq/.{Power[a_,1/2]:>a});
				If[MemberQ[tem,-1],Sow[intset[[i]]];Continue[]];(*in this case, there are still Baikov variables under square roots, we don't consider it*)
				rep=Thread@Rule[sq,tem];
				tem=Times@@(z-(z/.(sol/.rep))//Factor//Numerator//DeleteDuplicates);
				Sow[tem];
				(*Do[Sow[tem[[j]]],{j,1,Length[tem]}];*)
			]
		]
	,{i,1,Length[intset]}]
	][[2]];
	If[result=!={},Return[result[[1]]//Flatten],Return[result]];
];


Options[ResolveSingleVariable]={deBug->False,LastVar->False,HigherPower->False,DLog->False,AdInfo->True,SelectQ->True,"GeneralPinch"->True,"SpecialKinematics"->{}};
ResolveSingleVariable[iintset_,ihintset_,z_,OptionsPattern[]]:=Module[{intset,hintset,p,tem,tem1,tem2,tem3,coe,irintset={},rintset={},irhintset={},rhintset={},result={},sresult={},sol,temsol,hpintset={},rintsetp,dlog={},a1,b1,i,j,k,l},
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


(*resolve singularities relevant to square roots*)
If[rhintset=!={},
	tem=Times@@rhintset;(*expression under square root*)
	tem1=Exponent[tem,z];
	(*we will include the 'intersection' from two polynomials under square root in elliptic or quadratic case*)
	If[OptionValue[LastVar]&&OptionValue["GeneralPinch"],
		If[Length[rhintset]>1&&tem1>2,
			Table[irintset=Join[irintset,{Resultant[rhintset[[i]],rhintset[[j]],z]//Factor}];,{i,1,Length[rhintset]},{j,i+1,Length[rhintset]}]];
		If[tem1==2,
			irintset=Join[irintset,{Discriminant[tem,z]//Factor}]
		]
	];
	If[OptionValue[LastVar]&&OptionValue["SpecialKinematics"]=!={}&&Not@OptionValue["GeneralPinch"],(*in this case, we consider the degeneration of polynomials under special kinematics*)
		Do[
			If[Exponent[rhintset[[i]],z]>2,Continue[]];(*powers too high*)
			tem2=Discriminant[rhintset[[i]]/.OptionValue["SpecialKinematics"],z]//Factor;
			If[tem2===0,
				Print["Special kinematics work. We will take pinch for ",Short[rhintset[[i]],40]];
				irintset=Join[irintset,{Discriminant[rhintset[[i]],z]//Factor}]
			]
		,{i,1,Length[rhintset]}];
	];
	If[OddQ[tem1],(*if polynomial under square root is of odd power*)
		If[tem1>1,
			If[OptionValue[LastVar],Message[ResolveSingularities::elliptic,tem]];(*when there is only one integration variable, it is elliptic or hyperelliptic*)
			coe=0,(*we don't consider this case since it is hardly related to MPL*)
			coe=Prime[2024](*when this is a linear polynomial it can be related to dlog like Sqrt[b*a-c]/(z-a)/Sqrt[b*z-c], but we need to keep in mind that Sqrt[b]/Sqrt[b*z-c] is not a dlog*)
			(*here we use coe=Prime[2024] which is chosen randomly to distinguish it from the following case Sqrt[a]/Sqrt[az^2+bz+c] which is a dlog form*)
		],
		(*even power case*)
		If[tem1>2,
			If[OptionValue[LastVar],Message[ResolveSingularities::elliptic,tem]];(*when there is only one integration variable, it is elliptic or hyperelliptic*)
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
		If[OptionValue[SelectQ], (*now this is turned off by default*)
			tem1={{},tem1[[2]]};(*--------------------(^_^)--------------------*)
			(*here we remove those terms from perfect square roots under square root which we believe to be spurious letters*)
		];
		AppendTo[sresult,{{Sqrt[Times@@rhintset]},Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}];
		AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}];
	];
];

If[OptionValue[deBug],Print["coe: ",coe];Print["result: ",result]];

(*resolve singularities from rational part*)
rintsetp=EnhancedFactor[rintset,z];(*Enhanced factor can factor some quadratic polynomials and abbreviate some square roots of kinematic variables*)
p=Exponent[#,z]&/@rintset;(*power of rational polynomials w.r.t variable z*)
sol=Table[If[p[[i]]<3&&FreeQ[rintset[[i]],R],(Solve[rintsetp[[i]]==0,z]//DeleteDuplicates),Null],{i,1,Length[rintsetp]}];(*only solve linear polynomials and quadratic polynomials without square roots from the last-level calculation*)
tem2={};
Do[
	(*when there are polynomials with power of z higher than 2, we don't consider them since they may generate singularities with power other than 1/2*)
	If[p[[i]]>2,AppendTo[hpintset,rintset[[i]]];Continue[]];
	If[sol[[i]]===Null,Continue[]];(*if this polynomial has not been solved due to it involves square roots*)
	
			
	If[Length[sol[[i]]]==1,(*if this pole is a linear rational pole, we will keep this pole*)
		If[rhintset==={},
			(*if there are no square roots*)
			(*take pole at this solution*)
			If[!OptionValue[DLog],
				AppendTo[result,{sol[[i,1]],Join[irintset,(*{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)},*)(RemoveCoeff/@(NumeratorDenominator[(rintset/.sol[[i,1]]//Expand//Factor//DeleteCases[#,_?NumericQ]&)]//Flatten))]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,irhintset}],
				(*dlog type construction*)
				AppendTo[dlog,{Join[irintset,{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)}],irhintset}];(*a/(a*z-b) type dlog*)
				Table[If[k!=i,AppendTo[dlog,{Join[irintset,{RemoveCoeff[Times@@NumeratorDenominator[rintset[[k]]/.sol[[i,1]]//Expand//Factor]]}],irhintset}]],{k,1,Length[rintset]}](*(bc-ad)/(az-b)/(cz-d) type dlog*)
			]
			,
			(*if there are square roots*)
			If[OptionValue[LastVar]&&Exponent[Times@@rhintset,z]>2,(*Message[ResolveSingularities::elliptic,{z,rhintset}];*)Continue[]];(*when it is the last variable and the power is higher than 2 in square root, then this is elliptic or hyperelliptic function*)
			tem=Times@@@NumeratorDenominator[rhintset/.sol[[i,1]]//Factor];(*move denominator to numerator*)
			(*when the polynomials under square root vanish, then this is a branch cut, not a pole*)
			If[!FreeQ[tem,0],Continue[]];
			If[OptionValue[SelectQ],(*now this is turned off by default*)
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
		If[!FreeQ[sol[[i]],Power[_,1/2]],(*if this quadratic polynomial cannot be factorized by EnhancedFactor*)
			(*---------------------------when there is an irreducible quadratic polynomial in the denominator----------------------*)
			If[OptionValue[LastVar]&&Exponent[Times@@rhintset,z]>2,(*Message[ResolveSingularities::elliptic,{z,rhintset}];*)Continue[]];
			If[rhintset=!={}&&!ExistRelationQ[rintset[[i]],rhintset],
				If[OptionValue[AdInfo],Message[ResolveSingleVariable::additional,{rintset[[i]],rhintset}]];Continue[]
			];
			tem=Discriminant[rintset[[i]],z]*coe//Factor;
			If[tem===0,Continue[]];(*if tem is 0, this construction will stop here*)
			tem1=PerfectSquareSplit[tem];
			If[OptionValue[SelectQ],(*now this is turned off by default*)
				tem1={{},tem1[[2]]};(*---------------------------(^_^)-------------------------------*)
				(*here we remove those terms from perfect squares under square roots which we believe to be spurious letters*)
			];
			If[!OptionValue[DLog],
				AppendTo[sresult,{{rintset[[i]]}(*Solve[D[rintset[[i]],z]==0,z][[1]]*),Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}],
				(*a special kind of dlog, we slightly enlarge its possible realization.*)
				AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}]
			],
			(*if it can be factorized by EnhancedFactor*)
			Do[
				tem=Times@@((List@@(rintsetp[[i]]))/.sol[[i,l]]//Factor//DeleteCases[#,_?NumericQ]&);(*this quadratic polynomial will be handled separately*)
				tem1=Times@@(FactorList[tem][[All,1]]//DeleteCases[#,_?(FreeQ[#,R]&)]&);(*the genuine square root part*)
				tem=tem/tem1;
				If[!FreeQ[tem,R],Print["temR: ",tem]];
				If[rhintset==={},
					(*if there are no polynomials under square root*)
					(*take pole at this solution*)
					If[FreeQ[tem,x],(*if the square root has no chance to be 0, then we don't substitute it into other expressions for simplicity*)
						AppendTo[result,{sol[[i,l]],Join[irintset,(*{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)},*)(RemoveCoeff/@(NumeratorDenominator[Join[{},{tem}]]//Flatten//DeleteCases[#,_?NumericQ]&))]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,{tem1^2}]}],
						If[!OptionValue[DLog],
							AppendTo[result,{sol[[i,l]],Join[irintset,(*{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)},*)(RemoveCoeff/@(NumeratorDenominator[Join[(Delete[rintsetp,i]/.sol[[i,l]]//Factor),{tem}]]//Flatten//DeleteCases[#,_?NumericQ]&))]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,{tem1^2}]}],
							(*dlog type construction*)
							AppendTo[dlog,{Join[irintset,{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)}],irhintset}];(*a/(a*z-b) type dlog*)
							Table[If[k!=i,AppendTo[dlog,{Join[irintset,{RemoveCoeff[Times@@NumeratorDenominator[rintset[[k]]/.sol[[i,l]]//Factor]]}],irhintset}]],{k,1,Length[rintset]}](*(bc-ad)/(az-b)/(cz-d) type dlog*)
						]
					]
					,
					(*if there are square roots*)
					If[OptionValue[LastVar]&&Exponent[Times@@rhintset,z]>2,(*Message[ResolveSingularities::elliptic,{z,rhintset}];*)Continue[]];(*when it is the last variable and the power is higher than 2 in square root, then this is elliptic or hyperelliptic function*)
					If[FreeQ[tem,x],(*if the square root in solution has no chance to be 0, then we don't substitute it into other expressions for simplicity*)
						If[!ExistRelationQ[rintset[[i]],rhintset],If[OptionValue[AdInfo],Message[ResolveSingleVariable::additional,{rintset[[i]],rhintset}]];Continue[]];
						tem3=Discriminant[rintset[[i]],z]*coe//Factor;
						If[tem3===0,Continue[]];(*if tem is 0, this construction will stop here*)
						tem1=PerfectSquareSplit[tem3];
						If[OptionValue[SelectQ],(*now this is turned off by default*)
							tem1={{},tem1[[2]]};(*---------------------------(^_^)-------------------------------*)
							(*here we remove those terms from perfect squares under square roots which we believe to be spurious letters*)
						];
						If[!OptionValue[DLog],
							AppendTo[result,{sol[[i,l]],Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}],
							(*a special kind of dlog, we slightly enlarge its possible realization.*)
							AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]]//DeleteCases[#,_?NumericQ]&}]
						],
						tem3=Times@@@NumeratorDenominator[rhintset/.sol[[i,l]]//Factor];(*move denominator to numerator*)
						(*when the polynomials under square root vanish, then this is a branch cut, not a pole*)
						If[!FreeQ[tem3,0],Continue[]];
						If[OptionValue[SelectQ],(*now this is turned off by default*)
							tem3=tem3//DeleteCases[#,_?PerfectSquareQ]&;(*------------------(^_^)-------------------*)
							(*remove perfect squares under square roots, they are conjectured to be spurious*)
						];
						tem3=PerfectSquareSplit[Times@@tem3];
						AppendTo[tem2,{i,Times@@@tem3}];(*for later construction*)
						If[!OptionValue[DLog],
							AppendTo[result,{sol[[i,l]],Join[irintset,(RemoveCoeff/@(NumeratorDenominator[Join[(Delete[rintsetp,i]/.sol[[i,l]]//Factor),tem3[[1]],{tem}]]//Flatten//DeleteCases[#,_?NumericQ]&))]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,{tem1^2},tem3[[2]]]}],
							(*dlog type construction*)
							AppendTo[dlog,{Join[irintset,{Times@@Power@@@(FactorList[Coefficient[rintset[[i]],z,p[[i]]]]//DeleteCases[#,_?(NumericQ[First[#]]&)]&)}],irhintset}];(*a/(a*z-b) type dlog*)
							Table[If[k!=i,AppendTo[dlog,{Join[irintset,{RemoveCoeff[Times@@NumeratorDenominator[rintset[[k]]/.sol[[i,l]]//Factor]]}],irhintset}]],{k,1,Length[rintset]}](*(bc-ad)/(az-b)/(cz-d) type dlog*)
						]
					];
					(*tem=Times@@@NumeratorDenominator[rhintset/.sol[[i,l]]//Factor];(*move denominator to numerator*)
					(*when the polynomials under square root vanish, then this is a branch cut, not a pole*)
					If[!FreeQ[tem,0],Continue[]];
					If[OptionValue[SelectQ],(*now this is turned off by default*)
						tem=tem//DeleteCases[#,_?PerfectSquareQ]&;(*------------------(^_^)-------------------*)
						(*remove perfect squares under square roots, they are conjectured to be spurious*)
					];
					tem1=PerfectSquareSplit[Times@@tem];
					AppendTo[tem2,{i,Times@@@tem1}];(*for later construction*)
					If[!OptionValue[DLog],
						(*AppendTo[result,{sol[[i,1]],Join[irintset,tem1[[1]]]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}]*)
						AppendTo[result,{sol[[i,l]],Join[irintset,tem1[[1]],(Join[(((rintset(*//DeleteCases[#,_?(Exponent[#,z]>1&)]&*))/.sol[[i,l]]//Factor)//Numerator),{}])//Flatten]//DeleteCases[#,_?NumericQ]&//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}],
						(*dlog of type (Sqrt[ad^2+bd+c]/(z-d)/(Sqrt[az^2+bz+c]))*)
						(*Note that when two polynomials are rational, we can always perform the variable transformation, so we can set some poly in rintset at this pole in the same time*)
						AppendTo[dlog,{Join[irintset,tem1[[1]]//DeleteCases[#,_?NumericQ]&]//DeleteDuplicates,Join[irhintset,tem1[[2]]//DeleteCases[#,_?NumericQ]&]}]
					]*)
				]
			,{l,1,1}](*one of the two solutions is enough*)
		]
	]
,{i,1,Length[rintset]}];

(*a sepcial kind of dlog, they are actually combinations of above dlogs already constructed, but they may give new varieties*)
If[Length[tem2]>1,
	(*we search for the possible combinations of dlogs, they may give new varieties.*)
	(*Print["special construction: ",tem2];*)
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


Options[RemoveSquareRootsSol]={deBug->False,"infty"->False};
RemoveSquareRootsSol[sol_,OptionsPattern[]]:=Module[{keys,values,den,num,tem,temq,temsol,result,i,j,flag,t=1},
	result=Reap[
		Do[
			If[FreeQ[sol[[i]],R],Sow[sol[[i]]];Continue[]];(*we only want to remove those terms that involve R*)
			keys=Keys[sol[[i]]];
			values=Values[sol[[i]]];
			If[OptionValue["infty"],t=2];(*if the pole will map variables to infinity plane, then the last replacement should be considered separately, it can not be taken as a replacement*)
			temsol=Take[sol[[i]],{-t}];
			flag=1;
			Do[
				den=Denominator[values[[-j]]]/.temsol//Expand//Factor;
				If[!FreeQ[den,R],flag=0;Break[]];
				If[OptionValue[deBug],Print["den: ",den]];
				num=Numerator[values[[-j]]]/.temsol//Expand//Factor;
				If[!FreeQ[num,R],flag=0;Break[]];
				If[OptionValue[deBug],Print["num: ",num]];
				If[den===0,
					PrependTo[temsol,sol[[i,-j]]],
					PrependTo[temsol,Rule[keys[[-j]],num/den//Factor]]
				];
			,{j,t+1,Length[sol[[i]]]}];
			If[flag===0,
				Continue[],
				If[OptionValue["infty"],Sow[Append[temsol,sol[[i,-1]]]],Sow[temsol]];
			]
		,{i,1,Length[sol]}]
	][[2]];
	If[result=!={},Return[result[[1]]],Return[{}]];
];


(*ConjugatePairQ[p1_,p2_]:=Module[{tem,sq},
	sq=Cases[p1,R[_],Infinity]//DeleteDuplicates;
	If[sq==={},Return[True]];
	If[Length[sq]>1,Return[False]];
	tem={(p1/.{R[a_]:>-R[a]})-p2,(p1/.{R[a_]:>-R[a]})+p2}//Factor;
	If[MemberQ[tem,0],Return[True],Return[False]];
];
RemoveSquareRootsPoly[list_]:=Module[{tem,mark={},m,l},
	If[FreeQ[list,R],Return[list]];
	tem=Reap[
		Do[
			If[FreeQ[list[[l]],R],Sow[list[[l]]];Continue[]];
			If[MemberQ[mark,l],Continue[]];
			Do[
				If[ConjugatePairQ[list[[l]],list[[m]]],Sow[list[[l]]*list[[m]]//Expand];AppendTo[mark,m];Break[]]
			,{m,l+1,Length[list]}]
		,{l,1,Length[list]}]
	][[2]];
	If[tem==={},Return[{}],Return[tem[[1]]//DeleteDuplicates]];
];*)


Options[ResolveSingularities]={deBug->False,SortQ->True,AdInfo->True,SelectQ->False,RemoveSquareRoots->True,"GeneralPinch"->False,"SpecialKinematics"->{},"infty"->False};
ResolveSingularities[iintset_,ihintset_,OptionsPattern[]]:=Module[{intset,hintset,xl,xla,nl,lv={},sq={},sol={},tem,p,tem1,tem2,sol1={},ssol={},result,a1,b1,i,j,k,supv},
intset=iintset//SortBy[#,LeafCount]&;hintset=ihintset//SortBy[#,LeafCount]&;
xl=Cases[{intset,hintset},Subscript[x,_],Infinity]//DeleteDuplicates;(*number of remaining ISPs*)
If[xl==={},Message[ResolveSingularities::err];Return[{{},intset,hintset}]];

(*sort the variable by their appearance in polynomials*)
If[OptionValue[SortQ],
	nl=Table[{Count[FreeQ[#,xl[[i]]]&/@Join[(*intset,*)hintset],False],(Exponent[#,xl[[i]]]&/@intset)//DeleteCases[#,0]&//Min,Count[FreeQ[#,xl[[i]]]&/@intset,True],FirstPosition[intset,xl[[i]],Missing["NotFound"],Infinity][[1]],xl[[i]]},{i,1,Length[xl]}]//SortBy[#,{First,#[[2]]&,#[[3]]&,#[[4]]&}]&;
	xl=nl[[All,-1]];
	(*sort the variables by their appearing frequences*)
	(*here we set the rule for the order. In principle, we should consider every order. First, we choose those that appear less in the square root with most priority. Then we choose those with less power outside the square root. At last, we choose those appear more outside the square root *)
];
If[OptionValue[deBug],Print["xl: ",xl]];

(*we consider the single ISP case separately*)
If[Length[xl]==1,(*if there is only one ISP remaining*)
	Do[(*this will give all possible poles given by linear and quadratic equations*)
		If[Exponent[intset[[i]],xl[[1]]]<2,
			(*linear pole*)
			AppendTo[sol,Solve[intset[[i]]==0,xl][[1]]];
			AppendTo[lv,Coefficient[intset[[i]],xl[[1]],1]];
			If[OptionValue[deBug],Print["lv: ",lv," intset:",intset[[i]]]],
			(*quadratic pole or higher, for the quardratic polynomial, we just solve D[poly,x]==0, for higher polynomial, we keep its form*)
			AppendTo[ssol,{intset[[i]]}(*If[Exponent[intset[[i]],xl[[1]]]==2,Solve[Equal[D[intset[[i]],xl[[1]]],0],xl][[1]],{intset[[i]]}]*)];
			If[hintset=!={}&&!ExistRelationQ[intset[[i]],hintset],
				If[OptionValue[AdInfo],Message[ResolveSingularities::additional,{intset[[i]],hintset}];Print[{intset[[i]],hintset}]],
				If[OptionValue[SelectQ],
					(*this option has been turned off in current version*)
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
		lv=Join[lv,((Table[RemoveCoeff/@(NumeratorDenominator[intset/.sol[[i]]//Factor]//Flatten//DeleteDuplicates),{i,1,Length[sol]}]//Flatten))//DeleteCases[#,_?NumericQ]&//DeleteDuplicates];
		(*here, it need to stress that (a*b)/(z*(z+a)(z-b)) can not be a dlog, only two linear term can be combined to a dlog form a/(z(z-a)) or ((a-b)/((z-a)(z-b))). but familiar construction can work for multivariate case because we can perform a change of variables*)
		If[OptionValue[deBug],Print["lv: ",lv," intset:",intset]];
		AppendTo[sq,1];(*this corresponds to the case 1/x, that is one single pole of sol[[i]], it will multiply a factor from coefficients*)
		sol1=sol;
		(*AppendTo[sol1,{xl[[1]]->Infinity,xl[[1]]->0}](*in this case Infinity pole is present, not a branch cut*)*),
		
		(*when there are square roots in the denominator*)
		(*here we only concern square roots in dlog, so we won't consider the combinations of dlog any more since this won't give us new square root*)
		tem=Times@@hintset;
		p=Exponent[tem,xl[[1]]];
		If[p>2,(*elliptic or hyperelliptic case*)
			Message[ResolveSingularities::elliptic,tem];
			If[OptionValue["GeneralPinch"],
				If[Length[hintset]>1,Table[lv=Join[lv,{Resultant[hintset[[i]],hintset[[j]],xl[[1]]]//Factor}],{i,1,Length[hintset]},{j,i+1,Length[hintset]}]]
			],
			(*quadratic case*)
			If[p==2,
				If[OptionValue["GeneralPinch"],
					lv=Join[lv,{Discriminant[tem,xl[[1]]]//Factor}]
				]
			]
		];
		If[OptionValue["SpecialKinematics"]=!={}&&Not@OptionValue["GeneralPinch"],(*in this case, we consider the degeneration of polynomials under special kinematics*)
			Do[
				If[Exponent[hintset[[i]],xl[[1]]]>2,Continue[]];(*powers too high*)
				tem1=Discriminant[hintset[[i]]/.OptionValue["SpecialKinematics"],xl[[1]]]//Factor;
				If[tem1===0,
					Print["Special kinematics work. We will take more general pinch for ",Short[hintset[[i]],40]];
					lv=Join[lv,{Discriminant[hintset[[i]],xl[[1]]]//Factor}]
				]
			,{i,1,Length[hintset]}];
		];
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
					If[OptionValue[SelectQ],(*now this is turned off by default*)
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
		
		(*a sepcial kind of ls (dlog), they are actually combinations of above ls (dlogs) already got, but they may give cancel some square roots in the expression and give a rational expression*)
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
	xla=xl;(*we may adjust the order in the intermediate step*)
	supv=TimeConstrained[Do[
		(*if one element of tem doesn't contain the variable to be studied at all, then this is not what we want*)
		tem[[1]]=DeleteCases[tem[[1]],_?(FreeQ[#[[{2,3}]],xla[[i]]]&)];
		tem[[2]]=DeleteCases[tem[[2]],_?(FreeQ[#[[{2,3}]],xla[[i]]]&)];
		tem1={{},{}};(*structure: {{first type},{second type}}, element:{{},intset,hintset}*)
		Do[
			tem2=ResolveSingleVariable[tem[[1,j,2]],tem[[1,j,3]],xla[[i]],LastVar->If[i==Length[xl],True,False],SelectQ->OptionValue[SelectQ],"GeneralPinch"->OptionValue["GeneralPinch"],"SpecialKinematics"->OptionValue["SpecialKinematics"]];
			If[tem2[[1]]=!={},tem1[[1]]=Join[tem1[[1]],MapAt[Join[tem[[1,j,1]],#]&,tem2[[1]],{All,1}]]];
			If[tem2[[2]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[1,j,1]],#]&,tem2[[2]],{All,1}]]]
		,{j,1,Length[tem[[1]]]}];(*Analyze the first part of tem*)
		Do[
			tem2=ResolveSingleVariable[tem[[2,j,2]],tem[[2,j,3]],xla[[i]],LastVar->If[i==Length[xl],True,False],SelectQ->OptionValue[SelectQ],"GeneralPinch"->OptionValue["GeneralPinch"],"SpecialKinematics"->OptionValue["SpecialKinematics"]];
			If[tem2[[1]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[2,j,1]],#]&,tem2[[1]],{All,1}]]];
			If[tem2[[2]]=!={},tem1[[2]]=Join[tem1[[2]],MapAt[Join[tem[[2,j,1]],#]&,tem2[[2]],{All,1}]]]
		,{j,1,Length[tem[[2]]]}];(*Analyze the second part of tem, these will all be classified as second type*)
		tem=tem1;(*renew the original list*)
		If[tem1[[1]]=!={},(*we adjust the order according to the simplicity of variables appearing in the remaining expressions*)
			nl=Table[{(Exponent[#,xla[[j]]]&/@tem1[[1,1,2]])//DeleteCases[#,0]&//Min,FirstPosition[tem1[[1,1,2]]//SortBy[#,LeafCount]&,xla[[j]],Missing["NotFound"],Infinity][[1]],xla[[j]]},{j,i+1,Length[xl]}]//SortBy[#,{First,#[[2]]&}]&;
			If[nl=!={},xla=Join[Take[xla,i],nl[[All,-1]]]];
			If[OptionValue[deBug],Print["xl order renewed: ",xla];Print["tem: ",tem]];
		];
	,{i,1,Length[xl]}],500,"overtime"];(*we set time constraints in case some calculation stuck here*)
	If[supv==="overtime",Print["ResolveSingularities[] hasn't finished in 500 seconds. Baikov variables: ",xla];Return["overtime"]];
	If[tem[[1]]=!={},
		sol1=tem[[1,All,1]];(*form like {{a->1,b->2},{a->2,b->1}}*)
		If[OptionValue[deBug],Print["sol1: ",sol1]];
		lv=Join[lv,tem[[1,All,{2,3}]]//Flatten];
		tem2={};
		Do[(*move the abbreviate R[i] from rational letter list to square root list*)
			tem1=Times@@Cases[(FactorList[#][[All,1]]&/@tem[[1,s,2]])//Flatten,R[_]]/.{R[a_]:>Rlist[[a]]};
			AppendTo[tem2,tem1];
		,{s,1,Length[tem[[1]]]}];
		sq=Join[sq,Table[Times@@(Join[tem[[1,s,{3}]],{tem2[[s]]}]//Flatten),{s,1,Length[tem[[1]]]}]];
		If[OptionValue[RemoveSquareRoots],(*remove those solution which will result in square roots in the end*)
			sol1=RemoveSquareRootsSol[sol1,"infty"->OptionValue["infty"]];
			lv=lv//DeleteCases[#,_?(!FreeQ[#,R]&)]&;
			sq=sq//DeleteCases[#,_?(!FreeQ[#,R]&)]&;
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
			lv=lv//DeleteCases[#,_?(!FreeQ[#,R]&)]&;
			sq=sq//DeleteCases[#,_?(!FreeQ[#,R]&)]&;
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
Trans2Inf[intset_,hintset_,xl_,OptionsPattern[]]:=Module[{tem,tem1,nintset={},nhintset={},LVlocal={},sqlocal={},k},
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

CheckValidity[brep_,subset_,krep_]:=Catch@Module[{cut,tem,j},
cut=Thread@Rule[Subscript[x,#]&/@subset,0];
Do[
	tem=(brep[[j,2,1,All,1]]//Gram2Poly[#,krep]&)/.cut//Factor;
	If[!FreeQ[tem,0],Throw[{}]](*if one in the brep is 0 under cut, then we will return {}. This is designed to keep all the representations equal*)
,{j,1,Length[brep]}];
Throw[brep];
];


Options[FindCover]={deBug->False};
FindCover[list_,OptionsPattern[]]:=Module[{var,subset,flag=0,L=5,tem,temset,pos,result={},i,k=1},
	var=list//Flatten//Union;
	If[OptionValue[deBug],Print["var: ",var]];
	If[Length[var]<=6,Return[{{var,{}}}]];
	While[flag==0&&k<20,(*we will find a set of subsets which can cover the original set but with length smaller than the total set*)
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
		While[flag==0&&k<20,
			tem=Length[Complement[#,temset]]&/@subset;
			If[Max[tem]==0,flag=1;Break[]];
			pos=FirstPosition[tem,Max[tem]][[1]](*PositionLargest[tem][[1]]*);
			temset=Join[temset,subset[[pos]]]//Union;
			AppendTo[result,{subset[[pos]],Thread@Rule[#,Table[Prime[20+i*3],{i,1,Length[#]}]]&[Complement[var,subset[[pos]]]]}];
			k=k+1;
		];
		k=k+1
	];
	Return[result];
];


PolyDim[exp_]:=Module[{var,rep},
	(*return the dimension of one homogeneous multi-variate polynomial. If this is not a homogeous one, then it returns the maximal dimension of some monomials in this polynomial*)
	var=Variables[exp]//Flatten;
	rep=Thread@Rule[var,var*R];
	Return[Exponent[exp/.rep,R]];
]


TFindInstance[poly_,var_,domain_]:=Module[{pow,pos,sol,sq,sol1},
	pow=Exponent[poly,#]&/@var;
	pos=FirstPosition[pow,2];(*if this is a quadratic polynomial of some variable*)
	If[pos===Missing["NotFound"],Return[FindInstance[{(poly)==0(*,var[[1]]!=0*)},var,domain]//Quiet]];
	sol=Solve[poly==0,var[[pos]]][[1]];
	sq=Cases[sol,Power[a_,1/2]->a,Infinity]//DeleteDuplicates//DeleteCases[#,_?NumericQ]&;
	sol1=FindInstance[sq[[1]]==0,Complement[var,var[[pos]]],domain]//Quiet;
	If[sol1==={}||Head[sol1]===FindInstance,Return[FindInstance[{(poly)==0(*,var[[1]]!=0*)},var,domain]//Quiet]];
	Return[{Join[sol/.sol1[[1]],sol1[[1]]]}];
];


Options[RFindInstance]={deBug->False};
RFindInstance[poly_,ovar_,domain_,n_,OptionsPattern[]]:=Module[{var,pow,list,rep,sol,result,k,a,m=1,flag,count=1},
	var=RandomSample[ovar];
	pow=Exponent[poly,#]&/@var;
	list=(Partition[Riffle[pow,var],2]//SortBy[#,First]&)[[All,2]];(*sort the variables by its power*)
	(*Print["list: ",list];*)
	If[Length[var]<=2,Return[FindInstance[{poly==0,var[[1]]!=0},var,domain,n]]];
	pow=Exponent[poly/.Thread@Rule[var,var*a],a];
	If[pow<=2,(*if the total power of polynomial is less or equal to 2, it usually can be handled by FindInstance[]*)
		result=TimeConstrained[FindInstance[{poly==0,(Times@@var)!=0},var,domain,n],10,$Failed];
		If[result=!=$Failed&&result=!={}&&Head[result]=!=FindInstance,Return[result]];
	];
	result=Reap[
		flag=0;
		Do[
			If[count>n,Break[]];
			var=RandomSample[var];
			rep=Thread@Rule[Drop[var,m+1],RandomInteger[{1,5000},Length[var]-m-1]];
			If[OptionValue[deBug],Print["rep: ",rep]];
			sol=TFindInstance[{(poly/.rep)(*,var[[1]]!=0*)},Take[var,m+1],domain]//Quiet;
			If[OptionValue[deBug],Print["sol: ",sol]];
			If[sol==={}||Head[sol]===FindInstance,
				If[flag<10,flag=flag+1;Continue[],m=m+1;If[m==Length[var],Break[]]],
				Sow[Join[sol[[1]],rep]];count=count+1;
			]
		,{k,1,100}];(*maximal 100 times*)
	][[2]];
	If[Not@FreeQ[result[[1]],Equal],Print["poly, var: ",{poly,var}]];
	If[result==={},Return[{}],Return[result[[1]]]];
];


Options[PolesAnalyze]={deBug->False,AugAna->True,SelectQ->False,SelectAllQ->False,AdInfo->False,RemoveSquareRoots->True,"GeneralPinch"->False,"SpecialKinematics"->{},"samplingNum"->2};
PolesAnalyze[result_,topsector_,krep_,n_,OptionsPattern[]]:=Module[{start,subset,zerolist,brep,supersec,xl,cutr,cut,adcut,singular={},sol={},LVlocal,sqlocal,sqt,intset,hintset,tem,tem1,tem2,tem3,rl,tsingular={},cc,len,var,groebner,pos,pow,flag,brepr,benchlist,bench,u,tsol,vrep,cover,dim},
start=SessionTime[];
subset=Subsets[topsector]//ReverseSortBy[#,Length]&//DeleteCases[#,{}]&;
zerolist=GetMatZeroSector[result,n,Complement[Range[n],topsector]];(*all zero sectors*)
subset=DeleteCases[subset,_?(MemberQ[zerolist,Sector2Digits[#]]&)];(*remove all zero sectors*)
Print["totally ",Length[subset]," sectors need to be analyzed!"];

(*subset={{1,2,3,4,6,7,8}};*)
Monitor[Do[(*analyze sector by sector*)
	If[OptionValue[deBug],Print["subset: ",subset[[a]]]];(*//////////////////////////////////*)
	singular={};
	brep=GetBaikovMatRep[result,subset[[a]],n,"looporder"->ExtractLoopOrder[krep]];
	(*If[OptionValue[deBug],Print["path: ",brep]];*)
	If[OptionValue[AugAna],
		If[CheckValidity[brep,{},krep]==={},Continue[]];(*if some Grams in the representation equal to 0 before cut, then this sector is actually reducible already. We don't consider it.*)
		tem={CheckValidity[brep,subset[[a]],krep],Thread@Rule[Subscript[x,#]&/@subset[[a]],0]};(*check validity of the representations*)
		cc=Length[subset[[a]]];
		While[tem[[1]]==={},
			(*when maximal cut is 0, we consider the next-to-maximal cut and iterate this procedure*)
			(*we assume that there is at least one next-to-maximal cut is not 0, when it is not the case, a warning will be generated*)
			cc=cc-1;
			If[cc<=Length[subset[[a]]]-2,If[OptionValue[AdInfo],Message[PolesAnalyze::misp,subset[[a]],cc]]];
			cutr=Subsets[subset[[a]],{cc}];
			Do[
				tem={CheckValidity[brep,cutr[[j]],krep],Thread@Rule[Subscript[x,#]&/@cutr[[j]],0]};
				If[tem[[1]]=!={},Break[]]
			,{j,1,Length[cutr]}];
		];
		brep=tem[[1]];
		cutr=tem[[2]]
		,
		brep=GetBaikovMatRep[result,subset[[a]],n];
		cutr=Thread@Rule[Subscript[x,#]&/@subset[[a]],0];
	];
	If[OptionValue[deBug],Print["sec rep: ",brep[[All,1]]]];(*//////////////////////////////////*)
	If[OptionValue[SelectQ],If[OptionValue[SelectAllQ],len=2,len=Length[brep](*the number of minimal Baikov representations*)],len=1];
	pos={};(*keep track of representations which hasn't finished in given time*)
	Do[(*iteration for different representations of one sector*)
		intset={};
		hintset={};
		sol={};
		LVlocal={};
		sqlocal={};
		adcut=Complement[Subscript[x,#]&/@subset[[a]],Keys[cutr]];
		(*If[OptionValue[deBug],Print["adcut: ",adcut]];*)
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
		(*If[OptionValue[deBug],Print["intset,hintset: ",{intset,hintset}]];*)
		tem=ResolveSingularities[intset,hintset,SortQ->True,SelectQ->If[len>1,True,False],AdInfo->OptionValue[AdInfo],RemoveSquareRoots->OptionValue[RemoveSquareRoots],"GeneralPinch"->OptionValue["GeneralPinch"],"SpecialKinematics"->OptionValue["SpecialKinematics"]];(*when len>1, it means that we should always apply our selection rule*)
		(*If[OptionValue[deBug],Print["finished "]];*)
		(*If[OptionValue[deBug],Print["tem: ",tem];Print["sqlocal: ",sqlocal]];*)
		If[tem==="overtime",Print["section: ",subset[[a]]," repr: ",brep[[j]]," hasn't finished in given time! we give up the calculation of this representation."];AppendTo[pos,{j}];Continue[]];
		tem1=tem[[1]];
		sol=Join[sol,tem1//Flatten[#,1]&];
		LVlocal=Join[LVlocal,tem[[2]]];
		sqlocal=Join[{},sqt*tem[[3]]];
		(*the square roots from the coefficients should be multiplied to those from integrands*)

		(*If[OptionValue[deBug],Print["infinity plane calculation"]];*)
		(*now we consider sending all the isps to infinity plane which corresponds to the second Landau singularity*)
		xl=Cases[{intset,hintset},Subscript[x,_],Infinity]//DeleteDuplicates;
		If[Length[xl]>0,(*Length[xl]=1 case has been include before*)
			tem=Trans2Inf[intset,hintset,xl];
			intset=tem[[1]];
			hintset=tem[[2]];
			LVlocal=Join[LVlocal,tem[[3]]];
			sqt=sqt*Times@@(tem[[4]]);
			(*If[OptionValue[deBug],Print["infty intset,hintset: ",{intset,hintset}]];*)
			tem1=ResolveSingularities[intset,hintset,SortQ->True,SelectQ->If[len>1,True,False],AdInfo->OptionValue[AdInfo],RemoveSquareRoots->OptionValue[RemoveSquareRoots],"GeneralPinch"->OptionValue["GeneralPinch"],"SpecialKinematics"->OptionValue["SpecialKinematics"],"infty"->True];
			(*If[OptionValue[deBug],Print["tem1: ",tem1]];*)
			If[tem1==="overtime",Print["section: ",subset[[a]]," repr: ",brep[[j]]," hasn't finished in given time! we give up the calculation of this representation."];AppendTo[pos,{j}];Continue[]];
			tem=Join[Thread@Rule[xl,Infinity],#]&/@((tem1[[1]]//Flatten[#,1]&));
			sol=Join[sol,tem(*//Flatten[#,1]&*)];
			LVlocal=Join[LVlocal,tem1[[2]]];
			sqlocal=Join[sqlocal,sqt*tem1[[3]]]
		];
		(*If[OptionValue[deBug],Print["infinity plane calculation finished"]];*)
		AppendTo[singular,{sol//DeleteDuplicates,LVlocal//DeleteDuplicates,sqlocal//DeleteDuplicates,brep[[j,1]]}]
	,{j,1,Length[brep]}];
	If[singular=!={},tem=singular[[All,2]],Print["The calculation of subset ",subset[[a]]," hasn't finished in given time. It has been skipped!"];Continue[]];
	If[pos=!={},brep=Delete[brep,pos]];(*delete those representations whose calculations are not finished in given time*)
	If[!OptionValue[SelectQ]||OptionValue[SelectQ],
		(*now we try to remove all spurious letters in a more rigorous way*)
		tem3=tem;(*Print["tem: ",tem];*)
		(*among different representations, not all analysis of them are reliable, especially when the number of remaining Baikov variables are larger than 3*)
		(*now we will distinguish different representations according to the number of remaining Baikov variables*)
		rl=n-Length[subset[[a]]]-(Length/@brep[[All,1]]);
		(*If[OptionValue[deBug],Print["rl: ",rl]];*)
		If[Length[rl]==1,
			tem1=(FactorList[#][[All,1]]&/@(tem3[[1]]))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates;
			tem2=tem1;
			tem3={},
			(*if these representations are equal in the remaining number of Baikov variables, all of them will be taken*)
			tem=Table[(FactorList[#][[All,1]]&/@tem[[k]])//Flatten//DeleteCases[#,_?NumericQ]&,{k,1,Length[tem]}];
			tem2=If[tem=!={},Intersection[Sequence@@(tem)],{}];
			tem1=tem//Flatten//DeleteDuplicates;
			tem3=Complement[tem1,tem2];
		];
		(*Print["tem1: ",tem1];Print["tem3: ",tem3];*)
		(*If[subset[[a]]==={1,2,3,4,5,6,8},Print[" tem3: ",tem3];Print["cutr: ",cutr]];*)
		benchlist={};
		Do[(*calculating the bench result*)
			u=(Times@@Power@@@(((brep[[l,2,1]])//Gram2Poly[#,krep]&)/.cutr//Factor))//PowerExpand;
			If[Head[u]===Times,u=(List@@@((List@@u)//DeleteCases[#,_?NumericQ]&))/.adcut//Factor//DeleteCases[#,_?(FreeQ[First[#],x]&)]&,u={List@@u}/.adcut//Factor//DeleteCases[#,_?(FreeQ[First[#],x]&)]&];
			(*Print["sol: ",sol];(*/////////////////////////////////////////////*)*)
			xl=Cases[u,Subscript[x,_],Infinity]//DeleteDuplicates;
			If[xl==={},Continue[]];
			u=Times@@Power@@@u;
			AppendTo[benchlist,{u,xl,GetDimension[u,xl,Time->300]}]
		,{l,1,Length[brep]}];
		(*If[OptionValue[deBug],Print["bench dimension: ",benchlist]];*)
		Do[
			flag=0;
			pos=Position[singular[[All,2]],_?(MemberQ[#,tem3[[k]]]&),1];
			brepr=(Delete[brep,pos]);
			var=Variables[tem3[[k]]];
			pow=Exponent[tem3[[k]],#]&/@var;
			pos=FirstPosition[pow,1];(*if there is one variable linear in the equation, then the rational number solution of this system is easy to get*)
			If[pos===Missing["NotFound"],
				(*If[subset[[a]]==={2,3,4,5,6,8},Print["tem3,var: ",{tem3,var}]];(*//////////////////////////////////*)*)
				sol=TimeConstrained[RFindInstance[tem3[[k]](*Join[{tem3[[k]]==0,var[[1]]>0},Thread@Unequal[var,0]]*),var,Rationals,Max[OptionValue["samplingNum"],2]],120,0];(*find two numeric solutions to check*)
				If[sol===0,Print["haven't find solutions in given time. add this singularity: ",tem3[[k]]];(*tem1=Complement[tem1,{tem3[[k]]}];*)Continue[]];
				If[sol==={},Print["this polynomial has no non-zero solution: ",tem3[[k]]," sec: ",a];sol=FindInstance[{tem3[[k]]==0},var]],
				tsol=Solve[tem3[[k]]==0,var[[pos]]][[1]];
				sol=Table[vrep=Thread@Rule[Delete[var,pos],Table[Prime[2024+l*m],{m,1,Length[var]-1}]];Join[Thread@Rule[var[[pos]],var[[pos]]/.tsol/.vrep],vrep],{l,1,2}]//DeleteDuplicates
			];
			(*Print["dimension calculation"];(*//////////////////////////////////*)*)
			(*Print["tem3: ",tem3[[k]]];*)
			Do[
				pos=Position[brep,brepr[[l,1]]][[1,1]];
				bench=benchlist[[pos]];
				If[bench[[3]]===$Failed,Continue[]];(*overtime or degenerate, we will pass this case*)
				u=bench[[1]];
				xl=bench[[2]];
				If[Max[Table[GetDimension[u/.sol[[m]],xl,Time->300],{m,1,Length[sol]}]]>=bench[[3]],flag=1;Break[]]
			,{l,1,Length[brepr]}];
			(*Print["dimension calculation finished!"];(*//////////////////////////////////*)*)
			If[flag==1,tem1=Complement[tem1,{tem3[[k]]}]]
		,{k,1,Length[tem3]}];
		(*Now we also check whether element in tem2 is spurious*)
		(*Do[
			flag=0;
			var=Variables[tem2[[k]]];
			pow=Exponent[tem2[[k]],#]&/@var;
			pos=FirstPosition[pow,1];(*if there is one variable linear in the equation, then the rational number solution of this system is easy to get*)
			If[pos===Missing["NotFound"],
				(*If[subset[[a]]==={2,3,4,5,6,8},Print["tem3,var: ",{tem3,var}]];(*//////////////////////////////////*)*)
				sol=TimeConstrained[RFindInstance[tem2[[k]](*Join[{tem3[[k]]==0,var[[1]]>0},Thread@Unequal[var,0]]*),var,Rationals,2],120,0];(*find two numeric solutions to check*)
				If[sol===0,Print["haven't find solutions in given time. add this singularity: ",tem2[[k]]];(*tem1=Complement[tem1,{tem3[[k]]}];*)Continue[]];
				If[sol==={},Print["this polynomial has no non-zero solution: ",tem2[[k]]," sec: ",a];sol=FindInstance[{tem2[[k]]==0},var]],
				tsol=Solve[tem2[[k]]==0,var[[pos]]][[1]];
				sol=Table[vrep=Thread@Rule[Delete[var,pos],Table[Prime[2024+l*m],{m,1,Length[var]-1}]];Join[Thread@Rule[var[[pos]],var[[pos]]/.tsol/.vrep],vrep],{l,1,2}]//DeleteDuplicates
			];
			(*Print["dimension calculation"];(*//////////////////////////////////*)*)
			(*Print["tem3: ",tem3[[k]]];*)
			Do[
				bench=benchlist[[l]];
				If[bench[[3]]===$Failed,Continue[]];(*overtime or degenerate, we will pass this case*)
				u=bench[[1]];
				xl=bench[[2]];
				If[Max[Table[GetDimension[u/.sol[[m]],xl],{m,1,Length[sol]}]]>=bench[[3]],flag=1;Break[]]
			,{l,1,Length[brep]}];
			(*Print["dimension calculation finished!"];(*//////////////////////////////////*)*)
			If[flag==1,tem1=Complement[tem1,{tem2[[k]]}]]
		,{k,1,Length[tem2]}]*)
		,
		tem=Table[(FactorList[#][[All,1]]&/@tem[[k]])//Flatten//DeleteCases[#,_?NumericQ]&,{k,1,Length[tem]}];
		tem2=If[tem=!={},Intersection[Sequence@@(tem)],{}];
		tem1=If[tem=!={},Union[Sequence@@(tem)],{}](*our old way to remove spurious letters which is experimental*)
	];
	tsingular=AppendTo[tsingular,{singular[[All,{1,2,3,4}]],tem2,tem1,singular[[All,3]]//Flatten//DeleteDuplicates,subset[[a]]}]
,{a,1,Length[subset]}],ProgressIndicator[a,{1,Length[subset]}]];
Print["time consuming: ",SessionTime[]-start];
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


Options[CompatiblePoleQ]={deBug->False,CutConstraint->True};
CompatiblePoleQ[p1_,p2_,OptionsPattern[]]:=Module[{pt1,pt2,tem,xint},
	If[p1==={}||p2==={},Return[True]];
	(*If[Xor[FreeQ[p1,Rule[_,Infinity]],FreeQ[p2,Rule[_,Infinity]]],
		If[Not@ContainsOnly[Values[DeleteCases[{p1[[2]],p2[[2]]},_?(Not@FreeQ[#,Infinity]&)][[1]]],{0}],Return[False]]
	];(*poles with infinity can not be mixed with poles without infinity except for poles taken to be 0 since otherwise there may be ambiguities when taking the limits*)*)
	If[OptionValue[CutConstraint],(*poles are compatible only when they are in the same sector*)
		xint={Cases[p1[[1]],Subscript[x,_],Infinity]//DeleteDuplicates,Cases[p2[[1]],Subscript[x,_],Infinity]//DeleteDuplicates}//SortBy[#,Length]&;
		(* x variables which are cut variables*)
		If[!ContainsExactly[xint[[2]],xint[[1]]],Return[False]];
	];
	If[p1[[2]]==={}||p2[[2]]==={},Return[True]];(*the non-cut part is empty, so they must be compatible*)
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


Options[CompatiblePoleGraph]={deBug->False,OutputLevel->2,Info->True};
CompatiblePoleGraph[poles_,OptionsPattern[]]:=Module[{remain,graph,tem},
	If[poles==={},Return[{{},{}}]];
	If[Length[poles]==1,Return[{{1},{{1}}}]];(*in this case, one pole is compatible with itself*)
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
	][[2]];
	If[graph=!={},graph=graph[[1]]];
	If[OptionValue[OutputLevel]==1,Return[{graph}]];
	(*edge=EdgeList[graph]/.{UndirectedEdge->List};*)
	(*If[OptionValue[deBug],Print[FindClique[graph,{3},All]]];*)
	tem=FindClique[graph,Infinity,All]//Reverse;(*we find complete subgraphs formed by compatible poles*)
	If[tem=!={}&&OptionValue[Info],Print["maximal clique: ",Length[tem[[-1]]],"; # of cliques: ",Length[tem]]];
	If[OptionValue[OutputLevel]==2,Return[{graph,tem}]];
];


(*Options[MaxPole]={deBug->False};
MaxPole[p1_,p2_]:=Module[{a,b,c,d},
	If[!CompatiblePoleQ[p1,p2],Return[0]];
	c=Cases[p1[[2]],Subscript[x,_],Infinity]//DeleteDuplicates//Sort;
	d=Cases[p2[[2]],Subscript[x,_],Infinity]//DeleteDuplicates//Sort;
	If[c=!=d,Return[0]];(*when the ISPs are not the same and they state the pole for one common variable, we take this two poles as intrinsically different.*)
	(*a=Cases[p1,Subscript[x,_],Infinity]//DeleteDuplicates;
	b=Cases[p2,Subscript[x,_],Infinity]//DeleteDuplicates;
	If[ContainsAll[a,b],Return[1]];
	If[ContainsAll[b,a],Return[2]];*)
	Return[0];
];*)


(*Options[RemoveSPole]={deBug->False};
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
];*)
(*Options[DecendentQ]={deBug->False};
DecendentQ[basis_,pole_]:=Module[{tem,graph,xl,rep,flag=0},
	tem=Select[basis,CompatiblePoleQ[#,pole]&];(*First select all poles in basis that are compatible with the pole we want to check*)
	graph=CompatiblePoleGraph[tem,Info->False];(*generate the compatible pole graph which we will use in the next step*)
	If[OptionValue[deBug],Print["graph",graph]];
	xl=Cases[pole, Subscript[x,_],Infinity]//DeleteDuplicates;
	Do[
		rep=tem[[graph[[2,i]]]]//Flatten//DeleteDuplicates//Keys;
		If[ContainsAll[rep,xl],flag=1;Break[]];
	,{i,1,Length[graph[[2]]]}];
	If[flag==1,Return[True],Return[False]];
];

Options[RemoveSPole]={deBug->False};
RemoveSPole[poles_,OptionsPattern[]]:=Module[{tem,tem1,basis,i},
	tem=poles//SortBy[#,Length[(Cases[#,Subscript[x,_],Infinity]//DeleteDuplicates)]&]&;(*sort the poles by their length*)
	(*the problem is to use a small set of poles in the list to generate all possible poles we get.*)
	If[tem==={},Return[tem]];
	basis=Take[tem,1];
	Do[
		(*for every pole, we determine whether it can be covered by a set of smaller sets*)
		If[OptionValue[deBug]&&i<10,Print["basis",basis];Print["poles: ",tem[[i]]]];
		If[DecendentQ[basis,tem[[i]]],Continue[],AppendTo[basis,tem[[i]]]]
	,{i,2,Length[tem]}];
	Return[basis];
];*)


(*extension of RemoveSPole*)
Options[RemoveSPoleBeta]={deBug->False};
RemoveSPoleBeta[poles_,OptionsPattern[]]:=Module[{flag,tem,tem1,den,num,pos,pos1,var,sol,isprep,inf,infr,i,k},
	tem=poles;
	(*This function simplifies the poles which may be very simple if the former variable depends on the latter variable*)
	tem1=Reap[
		Do[
			If[!FreeQ[tem[[i]],Infinity],
				inf=DeleteCases[tem[[i,2]],_?(FreeQ[#,Infinity]&)];
				infr=Drop[DeleteCases[tem[[i,2]],_?(!FreeQ[#,Infinity]&)],-1];
				den=Denominator[Values[infr]]//.infr//Factor//Quiet;
				If[!FreeQ[den,0],
					Sow[tem[[i]]];
					pos=Position[den,0,1]//Flatten;
					(*check if the numerator is also 0 or not, if it is an indeterminant, then it usually means that a homogeneous polynomial exist, we need to be careful that the solution can be determined in a more careful way*)
					num=Numerator[Values[infr]]//.infr//Factor//Quiet;
					pos1=Position[num,0,1]//Flatten;
					If[Not@ContainsAll[pos1,pos],Continue[]];
					(*only when both the numerator and denominator are 0, then it is an indeterminant*)
					flag=0;
					Do[
						var=Cases[{Denominator[infr[[pos[[-k]],2]]]},Subscript[x,_],Infinity]//DeleteDuplicates;
						sol=Solve[Thread@Equal[(D[(infr[[pos[[-k]],1]]-infr[[pos[[-k]],2]])//Together//Numerator,#]&/@var)//.(Drop[infr,pos[[-k]]]),0],infr[[pos[[-k]],1]]];
						If[sol==={}||sol==={{}},flag=1;Break[],infr[[pos[[-k]]]]=(sol[[1,1]])]
					,{k,1,Length[pos]}];
					If[flag===1,Continue[]];
					Sow[{tem[[i,1]],Join[inf,infr,Take[tem[[i,2]],-1]]}];
					If[OptionValue[deBug],Print["indeterminant calculated: ",infr," from pole: ",tem[[i]]]];
					Continue[]
				];(*avoid the situation where some variable becomes infinity after substituting values*)
				isprep=Thread@Rule[Keys[infr],Values[infr]//.infr//Factor](*//Sort*);
				Sow[{tem[[i,1]],Join[inf,isprep,Take[tem[[i,2]],-1]]}]
				,
				den=Denominator[Values[tem[[i,2]]]]//.tem[[i,2]]//Factor//Quiet;
				If[!FreeQ[den,0],
					Sow[tem[[i]]];
					pos=Position[den,0,1]//Flatten;
					infr=tem[[i,2]];
					(*check if the numerator is also 0 or not, if it is an indeterminant, then it usually means that a homogeneous polynomial exist, we need to be careful that the solution can be determined in a more careful way*)
					num=Numerator[Values[infr]]//.infr//Factor//Quiet;
					pos1=Position[num,0,1]//Flatten;
					If[Not@ContainsAll[pos1,pos],Continue[]];
					(*only when both the numerator and denominator are 0, then it is an indeterminant*)
					flag=0;
					Do[
						var=Cases[{Denominator[infr[[pos[[-k]],2]]]},Subscript[x,_],Infinity]//DeleteDuplicates;
						sol=Solve[Thread@Equal[(D[(infr[[pos[[-k]],1]]-infr[[pos[[-k]],2]])//Together//Numerator,#]&/@var)//.(Drop[infr,pos[[-k]]]),0],infr[[pos[[-k]],1]]];
						If[sol==={}||sol==={{}},flag=1;Break[],infr[[pos[[-k]]]]=(sol[[1,1]])]
					,{k,1,Length[pos]}];
					If[flag===1,Continue[]];
					Sow[{tem[[i,1]],infr}];
					If[OptionValue[deBug],Print["indeterminant calculated: ",infr," from pole: ",tem[[i]]]];
					Continue[]
				];(*avoid the situation where some variable becomes infinity after substituting values*)
				isprep=Thread@Rule[Keys[tem[[i,2]]],Values[tem[[i,2]]]//.tem[[i,2]]//Factor](*//Sort*);(*the second part of tem[[i]] is for isp*)
				Sow[{tem[[i,1]],isprep}]
			];
		,{i,1,Length[tem]}];
	][[2]];
	If[tem1=!={},tem1=DeleteDuplicates[tem1[[1]]//Expand//Factor]];
	Return[tem1];
];


Options[ExtractPoleInfo]={deBug->False,OutputLevel->1,Simp->True,Enhanced->True,RemoveSquareRoots->True};
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
If[tem1==={},tem1={{Thread@Rule[xl,0],{}}}];(*if it is empty set, this is a one-loop problem*)
Print["original total ",Length[tem1]," poles"];
If[!OptionValue[Simp],Return[tem1]];
poles=tem1;
(*poles=RemoveSPole[tem1];(*remove poles that are already contained in other poles*)*)
If[OptionValue[Enhanced],poles=RemoveSPoleBeta[poles]];(*in enhanced version, we further remove some poles which can be identified to each other by substitute the values of the variables. For example, {x1->s-x2,x2->s} is equivalent to {x1->0,x2->s}*)
If[OptionValue[RemoveSquareRoots],poles=poles//DeleteCases[#,_?(Not@FreeQ[#,R]&)]&];(*remove poles containing R[i], that is, containing square roots*)
(*tem=CompatiblePoleMerge[poles];(*merge compatible poles*)
poles=RemoveSPole[tem];*)
Print["final total ",Length[poles]," poles"];
Return[{poles,Sequence@@CompatiblePoleGraph[poles]}];
];


Options[AllRationalLetters]={OutputLevel->2};
AllRationalLetters[polestructure_,OptionsPattern[]]:=
If[OptionValue[OutputLevel]==1,
Return[(FactorList[#][[All,1]]&/@(polestructure[[All,2]]//Flatten//DeleteDuplicates))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates],
Return[(FactorList[#][[All,1]]&/@(polestructure[[All,3]]//Flatten//DeleteDuplicates))//Flatten//DeleteCases[#,_?NumericQ]&//DeleteDuplicates]
];


ExtractSquareRoots[exp_]:=Module[{tem,rllist,result,i,tem1,tem2},
	tem=((Times@@Power@@@(MapAt[Mod[#,2]&,DeleteCases[FactorList[#],{_?NumericQ,_}],{All,2}]))&/@(exp[[All,4]]//Flatten//DeleteCases[#,_?NumericQ]&))//DeleteCases[#,1]&//DeleteDuplicates[#,((#1-#2//Expand)===0||(#1+#2//Expand)===0)&]&;
	(*then we remove terms that are spurious, that is, some factor is not in rllist*)
	rllist=AllRationalLetters[exp];
	result=Reap[
		Do[
			tem1=FactorList[tem[[i]]][[All,1]]//DeleteCases[#,_?NumericQ]&;
			tem2=RIntersection[rllist,tem1];
			tem1=Complement[tem1,tem2[[2]]];
			If[tem1==={},Sow[tem[[i]]]]
		,{i,1,Length[tem]}];
	][[2]];
	If[result=!={},Return[result[[1]]],Return[{}]];
];


Options[ConstructLetter]={deBug->False};
ConstructLetter[Gl1_,Gl2_,type_,numrep_,OptionsPattern[]]:=Module[{G1,G2,l1,l2,od,mm,MM,var,var1,var2,s1,s2,tem,tem1,sys,sol,mat,mat1,t,t1,pos,r},
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
		mat=(mat/.sol)//RowReduce;(*note that there are several set solutions, so this is a matrix*)
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
		Return[{dlogForm[G[Sequence@@od],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]],-G[mm,mm]*G[MM,MM],{Abs@Sum[Power[10,i]*(G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}],Abs@Sum[Power[10,i]*(-G[mm,mm]*G[MM,MM]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}]},1}],
		Return[{dlogForm[G[Sequence@@od],-G[mm,mm]*G[MM,MM]],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]],{Abs@Sum[Power[10,i]*(-G[mm,mm]*G[MM,MM]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}],Abs@Sum[Power[10,i]*(G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}]},1}]
	],
	If[Abs[l1-l2]==0,
		Return[{dlogForm[G[Sequence@@od],-G[mm,mm]*G[MM,MM]],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]],{Abs@Sum[Power[10,i]*(-G[mm,mm]*G[MM,MM]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}],Abs@Sum[Power[10,i]*(G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}]},2}],
		Return[{dlogForm[G[Sequence@@od],G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]],-G[mm,mm]*G[MM,MM],{Abs@Sum[Power[10,i]*(G[od[[1]],od[[1]]]*G[od[[2]],od[[2]]]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}],Abs@Sum[Power[10,i]*(-G[mm,mm]*G[MM,MM]/.{G[{},{}]->1}//Gram2Poly[#,numrep[[i]]]&),{i,1,3}]},2}]
	]
];
];

ConstructFromGram[gramlist_,numrep_]:=Table[ConstructLetter[gramlist[[i,1]],gramlist[[i,j]],If[IntegerQ[gramlist[[i,1,2]]/.{\[Epsilon]->0}],2,1],numrep],{i,1,Length[gramlist]},{j,2,Length[gramlist[[i]]]}]//Flatten[#,1]&//DeleteDuplicatesBy[#,({First[#],Last[#]}&)]&;

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


SameGramQ[gp1_,gp2_]:=Module[{var,numrep,rep,i},
	var=(Cases[{gp1,gp2},_G,Infinity]//DeleteDuplicates)/.{G[a_,b_]:>GramMat[a,b,{}]}//Variables;
	numrep=Thread@Rule[var,Table[Prime[100+23*i],{i,1,Length[var]}]];
	If[Not@FreeQ[{gp1-gp2,gp1+gp2}//Gram2Poly[#,numrep]&,0],Return[True],Return[False]];
];


Options[SameAlgLetterGramQ]={deBug->False};
SameAlgLetterGramQ[letter1_,letter2_,OptionsPattern[]]:=Module[{sq1,sq2},
	sq1=Cases[letter1,Power[a_,1/2]->a,Infinity]//DeleteDuplicates;
	sq2=Cases[letter2,Power[a_,1/2]->a,Infinity]//DeleteDuplicates;
	If[SameGramQ[sq1,sq2]&&SameGramQ[letter1[[2]],letter2[[2]]],Return[True],Return[False]];
];


Options[ExtractAlgLetter]={deBug->False,KineVar->{},KineRestrict->False,LoopRestrict->False,OneLoop->False};
ExtractAlgLetter[result_,topsector_,n_,OptionsPattern[]]:=Module[{subset,zerolist,reform,var,numrep,rep,glist,tem,cmGram,letter,letterref,summary={},i,a},
reform=ReformRep[result];(*Note that the structure of reform is {{g,c,{paths}},{g,c,{paths}},...}*)
subset=Subsets[topsector]//ReverseSortBy[#,Length]&;
zerolist=GetMatZeroSector[result,n,Complement[Range[n],topsector]];(*all zero sectors*)
subset=DeleteCases[subset,_?(MemberQ[zerolist,Sector2Digits[#]]&)];(*remove all zero sectors*)
Print["totally ",Length[subset]," sectors need to be analyzed!"];

(*get all variables and take numeric values of them to give every letter we construct an identification code*)
var=Variables[Cases[result[[1]],_G,Infinity]/.{G[o_,z_]:>GramMat[o,z,{}]}];
numrep=Table[Thread@Rule[var,RandomPrime[{100,10000},Length[var]]],{i,1,3}];

letterref={};(*keep track of algbraic letters so that repeated ones will not included in the result*)
Monitor[Do[
	If[OptionValue[deBug],Print["subset: ",subset[[a]]]];
	rep=GetBaikovMatRep[result,subset[[a]],n]//Simplify;
	cmGram=ExtractRelevantGramSim[rep,reform,n,KineVar->OptionValue[KineVar],KineRestrict->OptionValue[KineRestrict],LoopRestrict->OptionValue[LoopRestrict],OneLoop->OptionValue[OneLoop](*,deBug->If[subset[[a]]==={2,3,4,6},True,False]*)];
	(*If[OptionValue[deBug],Print["cmGram: ",cmGram]];*)
	letter=ConstructFromGram[cmGram,numrep];(*//DeleteCases[#,_?(!AdmissibleLetterQ[#]&)]&*);(*element: {Log[_],_,1|2}*)
	letter=Complement[letter,letterref];
	letterref=Join[letterref,letter]//DeleteDuplicatesBy[#,#[[3]]&]&;(*Delete duplicates by the first term.*)
	letter=Intersection[letter,letterref];(*take only unique terms that haven't appeared so far*)
	(*letterref=LetterPathRefined[reform,letter,krep];*)
	AppendTo[summary,{letter,subset[[a]]}];
	(*If[a==3,Break[]]*)
	,{a,1,Length[subset]}],
ProgressIndicator[a,{1,Length[subset]}]];
Return[summary];
];


Options[SpecialSimplify]={AdRep->{},deBug->False};
SpecialSimplify[{gr_,ga_},OptionsPattern[]]:=Catch@Module[{rep,rp,tem,sq,numsq,densq,R},
sq=ga//Together;
If[sq===0,Throw[{0,0}]];
rp=gr//Together;
If[rp===0,Throw[{0,0}]];(*If there is one zero, then this letter will be degenerate*)
rep=OptionValue[AdRep];
(*we apply the replacement to the expression one-by-one*)
If[rep==={},
	If[Head[sq]===Times,sq=List@@sq,sq={sq}];
	tem=Table[If[MatchQ[sq[[i]],Power[_,_?EvenQ]],Replace[sq[[i]],Power[z_,a_]:>Power[z,a/2]],If[MatchQ[sq[[i]],Power[_,_?OddQ]],Replace[sq[[i]],Power[z_,a_]:>Power[R*z,1/2]*Power[z,Quotient[a,2]]],Power[sq[[i]],1/2]]],{i,1,Length[sq]}];(*take the square root, RR is introduced to avoid things like Power[a,3/2]*)
	tem=((Times@@tem)/(rp))//Factor;(*cancel common parts between gr and Sqrt[ga]*)
	If[tem===0,Throw[{0,0}]];
	numsq=Cases[{Numerator[tem]},Power[z_,1/2]->z,Infinity];(*treat the square root in numerator and denomiantor separately*)
	densq=Cases[{Denominator[tem]},Power[z_,1/2]->z,Infinity];
	sq=Join[numsq,densq];
	densq=Times@@densq;
	If[OptionValue[deBug],Print["tem ",tem];Print["sq ",sq];Print["densq: ",densq]];
	tem=tem/.{Power[_,1/2|-1/2]->1}//Factor;(*remove square root part*)
	Throw[{densq/(tem)/.{R->1}//Factor,(Times@@(sq(*/.{Power[z_,1/2]:>z}*)/.{R->1}))//Factor}],
	
	(*we add time constraints in case the Factor[] take too long time*)
	Do[
		sq=TimeConstrained[sq/.rep[[i]]//Factor,60,sq/.rep[[i]]//Together];
		If[sq===0,Throw[{0,0}]];
		If[(Denominator[rp]/.rep[[i]]//Factor)===0,rp=0;sq=0;Break[]];(*when the rational part becomes infinity after replacement*)
		rp=TimeConstrained[rp/.rep[[i]]//Factor,60,rp/.rep[[i]]//Together];
		If[rp===0,Throw[{0,0}]];(*when the rational part or square root part become 0 after replacement*)
		If[Head[sq]===Times,sq=List@@sq,sq={sq}];
		tem=Table[If[MatchQ[sq[[i]],Power[_,_?EvenQ]],Replace[sq[[i]],Power[z_,a_]:>Power[z,a/2]],If[MatchQ[sq[[i]],Power[_,_?OddQ]],Replace[sq[[i]],Power[z_,a_]:>Power[R*z,1/2]*Power[z,Quotient[a,2]]],Power[R*sq[[i]],1/2]]],{i,1,Length[sq]}];
		If[OptionValue[deBug],Print["sq: ",sq];Print["rp: ",rp];Print["tem: ",tem]];
		tem=TimeConstrained[((Times@@tem)/(rp))//Factor,60,((Times@@tem)/(rp))//Together//Cancel];
		If[tem===0,rp=0;sq=0;Break[]];
		numsq=Cases[{Numerator[tem]},Power[z_,1/2]->z,Infinity];(*treat the square root in numerator and denomiantor separately*)
		densq=Cases[{Denominator[tem]},Power[z_,1/2]->z,Infinity];
		sq=Join[numsq,densq];
		densq=Times@@densq;
		If[OptionValue[deBug],Print["tem ",tem];Print["sq ",sq];Print["densq: ",densq]];
		tem=TimeConstrained[tem/.{Power[_,1/2|-1/2]->1}//Factor,60,tem/.{Power[_,1/2|-1/2]->1}];
		rp=TimeConstrained[densq/(tem)//Factor,60,densq/(tem)//Together//Cancel];
		sq=TimeConstrained[(Times@@(sq(*/.{Power[z_,1/2]:>z}*)))//Factor,60,(Times@@(sq(*/.{Power[z_,1/2]:>z}*)))];
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
			If[!FreeQ[pole[[k]],Subscript[x, 0]]||!FreeQ[Values[pole[[k]]],x],(*decides whether we can substitute the variables once and for all, or we need to perform step-by-step, when there is infinity, this always needs to be done step by step.*)
				rep=pole[[k]];rep1={},
				If[({Numerator[tem[[1]]],Denominator[tem[[1]]]^2*tem[[2]]}/.pole[[k]]//Expand)=!={0,0},rep={pole[[k]]};rep1={},rep=pole[[k]];rep1={}]
			];
			(*rep=pole[[k]];rep1={}*)
		];
		If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1];Print["tem: ",tem]];
		tem=SpecialSimplify[tem,AdRep->Join[rep,rep1,{OptionValue[AdRep]}]];
	,{k,2,Length[pole]}],
	(*if pole==={}, then this only involves pure kinematics grams*)
	tem=SpecialSimplify[({gr,ga}//Gram2Poly[#,krep]&)]
];
If[((tem[[1]]^2-tem[[2]]//Together)===0),Return[0]];
Return[Log[(tem[[1]]+Sqrt[tem[[2]]])/(tem[[1]]-Sqrt[tem[[2]]])//Factor//Map[Collect[#,Power[_,1/2],Factor]&,#,2]&]]
];

Options[ApplyPoleToGram]={deBug->False,AdRep->{},"RemovePerfectSquare"->True};
ApplyPoleToGram[pole_,gram_,krep_,OptionsPattern[]]:=Module[{tem,tem1,ntem,fac,rep,rep1,isp,pos,k,intexp,i},
	tem=(Det/@((gram/.{G[a_,b_]:>GramMat[a,b,krep]})/.pole[[1]]))/.OptionValue[AdRep];
	intexp=(FactorList/@tem)//Flatten[#,1]&//DeleteCases[#,{_,_?(EvenQ)}]&;
	ntem=(intexp[[All,1]]);
	If[OptionValue[deBug],Print["ntem: ",ntem];];
	If[Not@OptionValue["RemovePerfectSquare"],(*if we don't remove the perfect square in the Gram*)
		Do[
		If[!FreeQ[pole[[k]],Infinity],
			isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
			If[OptionValue[deBug],Print["isp: ",isp]];
			rep=Thread@Rule[isp,isp/Subscript[x, 0]];
			tem1=DeleteCases[pole[[k]],Rule[_,Infinity]];
			rep1=If[tem1=!={},If[MemberQ[isp,tem1[[-1,1]]],Drop[tem1,-1],tem1],{}]
			(*when projecting to infinity plane, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*)
			,
			pos=Position[pole[[k]],{__},1];
			If[pos==={},rep=pole[[k]];rep1={},rep=Take[pole[[k]],pos[[1,1]]]//Flatten;rep1=Drop[pole[[k]],pos[[1,1]]]];
		];
		(*tem=SpecialSimplify[{1,Times@@tem},AdRep->Join[rep,rep1,{OptionValue[AdRep]}]];*)
		If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
		
		If[!FreeQ[rep,Subscript[x, 0]],
			ntem=Numerator[ntem/.rep//Together];
			ntem=ntem/.{Subscript[x, 0]->0}//Factor;(*here we take the variables to infinity*)
			If[OptionValue[deBug],Print["ntem: ",ntem]];
			Do[
				If[Keys[rep1[[i]]]===Subscript[x, 0],Continue[]];(*Subscript[x, 0] has been removed above*)
				ntem=NumeratorDenominator[ntem/.rep1[[i]]//Together]//Flatten;
			,{i,1,Length[rep1]}],
			If[FreeQ[Values[rep],x],(*the replacement is free of order*)
				ntem=ntem/.rep//Together,
				Do[
				    ntem=NumeratorDenominator[ntem/.rep[[i]]//Together]//Flatten;
				,{i,1,Length[rep]}];
			];
		];
		,{k,2,Length[pole]}];
		Return[Times@@ntem//Factor];
	];
	
	Do[
		If[!FreeQ[pole[[k]],Infinity],
			isp=Cases[pole[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
			If[OptionValue[deBug],Print["isp: ",isp]];
			rep=Thread@Rule[isp,isp/Subscript[x, 0]];
			tem1=DeleteCases[pole[[k]],Rule[_,Infinity]];
			rep1=If[tem1=!={},If[MemberQ[isp,tem1[[-1,1]]],Drop[tem1,-1],tem1],{}]
			(*when projecting to infinity plane, the polynomial is homogeneous so the last variable will certainly be cancelled from expression, we don't need its replacement*)
			,
			pos=Position[pole[[k]],{__},1];
			If[pos==={},rep=pole[[k]],rep=Take[pole[[k]],pos[[1,1]]]//Flatten;rep1=Drop[pole[[k]],pos[[1,1]]]];
		];
		(*tem=SpecialSimplify[{1,Times@@tem},AdRep->Join[rep,rep1,{OptionValue[AdRep]}]];*)
		If[OptionValue[deBug],Print["rep: ",rep];Print["rep1: ",rep1]];
		
		
		If[!FreeQ[rep,Subscript[x, 0]],
			ntem=NumeratorDenominator[ntem/.rep//Together]//Flatten;
			intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;(*remove perfect square*)
			ntem=(Tally[intexp]//DeleteCases[#,{_,_?EvenQ}]&)[[All,1]];(*remove perfect square*)
			ntem=ntem/.{Subscript[x, 0]->0};(*here we take the variables to infinity*)
			intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;
			ntem=(intexp);
			If[OptionValue[deBug],Print["ntem: ",ntem]];
			Do[
				If[Keys[rep1[[i]]]===Subscript[x, 0],Continue[]];(*Subscript[x, 0] has been removed above*)
				ntem=NumeratorDenominator[ntem/.rep1[[i]]//Together]//Flatten;
				intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;
				ntem=(intexp)(*remove perfect square*)
			,{i,1,Length[rep1]}],
			If[FreeQ[Values[rep],x],(*the replacement is free of order*)
				ntem=NumeratorDenominator[ntem/.rep//Together]//Flatten;
				(*If[OptionValue[deBug],Print["ntem: ",ntem]];*)
				intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;(*remove perfect square*)
				intexp=(Tally[intexp]//DeleteCases[#,{_,_?EvenQ}]&)[[All,1]];(*remove perfect square*)
				If[Not@FreeQ[intexp,0],
				   (*when the result is 0, we need to be careful since intermediate step may generate an expression which is a perfect square and it is 0 at last, then this 0 is not what we want and we must apply the replacement step by step*)
				   Do[
				       ntem=NumeratorDenominator[ntem/.rep[[i]]//Together]//Flatten;
					   intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;
					   ntem=(intexp)(*remove perfect square*)
				   ,{i,1,Length[rep]}],
				  ntem=intexp
				],
				Do[
				    ntem=NumeratorDenominator[ntem/.rep[[i]]//Together]//Flatten;
				    (*If[OptionValue[deBug],Print["ntem: ",ntem]];*)
					intexp=Table[TimeConstrained[PerfectSquareSplit[ntem[[i]]][[2]],30,ntem[[i]]],{i,1,Length[ntem]}]//Flatten//DeleteCases[#,1|-1]&;
					ntem=(intexp)(*remove perfect square*)
				,{i,1,Length[rep]}];
			];
		];
	,{k,2,Length[pole]}];
	ntem=(Tally[ntem]//DeleteCases[#,{_,_?EvenQ}]&)[[All,1]];
	tem=Times@@ntem;
Return[tem];
];


Options[PolyMemberQ]={deBug->False};
PolyMemberQ[polylist_,poly_,OptionsPattern[]]:=Catch@Module[{var,numrep,i,tem,tem1,check={}},
	(*check whether poly is an item in polylist, they may differ by a const number*)
	var=Variables[{poly,polylist}];
	numrep=Reap[Do[tem=Thread@Rule[var,RandomPrime[{1000,10000},Length[var]]];If[(poly/.tem)=!=0,Sow[tem]],{i,1,3}]][[2]];
	If[numrep=!={},numrep=numrep[[1]],Print["this polynomial is singular: ",poly]];
	tem=Table[polylist/.numrep[[i]],{i,1,Length[numrep]}];
	tem1=DiagonalMatrix[poly/.numrep]//Inverse;
	check=Length[Union[#]]&/@(tem1 . tem//Transpose);
	(*If[OptionValue[deBug],Print["tem: ",tem];Print["check: ",check]];*)
	If[MemberQ[check,1],Throw[True],Throw[False]];
];


Options[AdmissiblePoleQ]={deBug->False};
AdmissiblePoleQ[glist_,pole_,permsq_,krep_,OptionsPattern[]]:=Module[{tem,ntem(*,c,d*)},
If[pole==={},
	(*cases where glist is list of Grams of external kinematics*)
	tem=(Times@@glist)//Gram2Poly[#,krep]&//Factor,
	tem=ApplyPoleToGram[pole,glist,krep]/.{Subscript[x,_]->1};
];

If[OptionValue[deBug],Print["tem: ",tem];];
If[NumericQ[tem],Return[{False}]];(*when the polynomials under square root are all perfect square, it is actually a rational letter, we don't count them here*)
(*c=PolyMemberQ[permsq,tem];d=AnyTrue[permsq,NumericQ[#/tem//Cancel]&];
If[c=!=d,Print["exp: ",permsq,tem]];*)
If[PolyMemberQ[permsq,tem](*AnyTrue[permsq,NumericQ[#/tem//Cancel]&]*),Return[{True,tem}],Return[{False}]];
];


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

TDivisible[list_,targetlist_]:=Module[{intexp,i},
	intexp=IntegerExponent[Abs[list],Abs[#]]&/@targetlist;
	Do[
	 If[Length[intexp[[i]]//Union]!=1,Return[False]];
	,{i,1,Length[intexp]}];
	intexp=intexp[[All,1]];
	If[Length@Union[list/Product[Power[targetlist[[i]],intexp[[i]]],{i,1,Length[targetlist]}]]!=1,Return[False]];
	Return[True];
]

Options[SelectAlgLetter]={deBug->False};
SelectAlgLetter[alg_,rl_,OptionsPattern[]]:=Module[{tem,tem1,int,select={},var,numrep,rlnum,temnum,pos},
var=Variables[rl];
numrep=Table[Thread@Rule[var,RandomPrime[{10^5,10^6},Length[var]]],{i,1,4}];(*we will pre-select by numerical test*)
rlnum=Table[rl/.numrep[[i]],{i,1,4}]//Transpose;
Do[
	tem=(Numerator[#[[1]]]*Denominator[#[[1]]]&[alg[[i]]])//Expand;
	temnum=Table[tem/.numrep[[i]],{i,1,4}];
	tem1=Mod[temnum,#]&/@rlnum;
	If[OptionValue[deBug],Print["tem1: ",tem1]];
	If[Not@MemberQ[tem1,{0,0,0,0}],
		Continue[],(*if numeric test fails, continue*)
		pos=Position[tem1,{0,0,0,0},1]//Flatten;
		If[OptionValue[deBug],Print["pos: ",pos]];
		If[Not@TDivisible[temnum,rlnum[[pos]]],Continue[],AppendTo[select,alg[[i]]]]
	];
	(*tem=FactorList[tem][[All,1]]//DeleteCases[#,_?NumericQ]&;
	int=RIntersection[tem,rl];
	If[OptionValue[deBug]==True,Print["int: ",int];Print["tem: ",tem]];
	If[Complement[tem,int[[1]]]=!={},Continue[],AppendTo[select,alg[[i]]]];*)
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
MergePoles[poleslist_,OptionsPattern[]]:=Module[{rep1,rep2,tem1,tem2,ntem2,isp,int,intv,tuple,k},
	(*poleslist should be in the form {{{Subscript[x, 1]->0,Subscript[x, 2]->0,...},{Subscript[x, 6]->Infinity,...}},{...}}*)
	rep1=poleslist[[All,1]]//Flatten//DeleteDuplicates//Sort;(*poles can be merged must have the same first part since they come from the same sector*)
	tem1=Cases[poleslist[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//DeleteDuplicates;(*Note that the second part has order so we can't sort them*)
	tem2=DeleteCases[poleslist[[All,2]],_?(FreeQ[#,Rule[_,Infinity]]&)]//DeleteDuplicates;(*poles with infinity will be pre-processed*)
	If[Length[tem2]!=0,
		(*if there are infinity poles*)
		tem2=Table[If[Length[tem2[[k]]]==1,Append[tem2[[k]],Rule[tem2[[k,1,1]],0]],tem2[[k]]],{k,1,Length[tem2]}];(*complete some single poles like {x->Infinity} to {x->Infinity,x->0}*)
		ntem2={};
		Do[
			isp=Cases[tem2[[k]],Rule[Subscript[x,a_],Infinity]->Subscript[x,a],Infinity]//DeleteDuplicates;
			AppendTo[ntem2,Join[{Thread@Rule[isp,isp/Subscript[x, 0]],Subscript[x, 0]->0},DeleteCases[tem2[[k]],Rule[_,Infinity]]//Delete[#,-1]&]]
		,{k,1,Length[tem2]}];
		tem1=Join[tem1,ntem2]
	];
	(*due to that the order of replacement is important, different orders may have different results, we give all possible order of two poles that can be merged*)
	tuple=Permutations[Range[Length[tem1]]];
	rep2=tem1[[#]]&/@tuple;(*all possible order*)
	Return[Table[{rep1,rep2[[k]]//Flatten[#,1]&},{k,1,Length[rep2]}]];
	(*the final data form will be like {{{Subscript[x, 1]->0,Subscript[x, 2]->0},{Subscript[x, 3]->s,Subscript[x, 4]->t}},{{Subscript[x, 1]->0,Subscript[x, 2]->0},{{Subscript[x, 3]->Subscript[x, 3]/Subscript[x, 0],Subscript[x, 4]->Subscript[x, 4]/Subscript[x, 0]},Subscript[x, 3]->s,Subscript[x, 4]->t}},...}*)
];


Options[ApplyPoleToVar]={deBug->False};
ApplyPoleToVar[poles_,var_]:=Catch@Module[{isp,tem,i},
	isp=Keys/@poles[[All,2]];
	Do[If[(Not@ContainsAll[var,isp[[i]]])&&(Not@ContainsNone[var,isp[[i]]]),Throw["forbidden"]],{i,1,Length[isp]}];
	tem=poles//Flatten//DeleteDuplicates//GatherBy[#,(Last[#]===Infinity)&]&;
	If[Length[tem]==1,Throw[var/.tem[[1]]//Factor]];
	Throw[Table[If[MemberQ[Keys[tem[[2]]],var[[i]]],g[var[[i]]/.tem[[1]]],var[[i]]/.tem[[1]]]//Factor,{i,1,Length[var]}]];
];


Options[FindPoleMaps]={deBug->False};
FindPoleMaps[poles_,var_,sector_:{},OptionsPattern[]]:=Module[{alist={},xl,tem,atem={},alisttem,l=2,max,subset,n,pos},
	(*poles consists of three parts: poles, graph, cliques of the graph.*)
	(*we find the minimal pole maps that can cover var*)
	(*minimal pole map means this map only involve the variables in var with extra variables as less as possible*)
	(*two poles can be merged only when they share the same cut*)
	If[var==={},Return[{}]];
	xl={Thread@Rule[Subscript[x,#]&/@sector,0],{}};
	If[xl[[1]]==={},(*if the sector has not been specified*)
		tem=poles[[1]];
		alist=Reap[
				Do[
					If[ContainsAll[Keys[tem[[i]]//Flatten],var],Sow[i]]
				,{i,1,Length[tem]}]
		][[2]];(*all the single pole maps that can cover var*)
		,
		tem=poles[[1]];
		alist=Reap[
				Do[
					If[ContainsAll[Keys[tem[[i]]//Flatten],var]&&CompatiblePoleQ[tem[[i]],xl],Sow[i]]
				,{i,1,Length[tem]}]
		][[2]];(*all the single pole maps that can cover var*)
	];
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
	If[OptionValue[deBug],Print["atem: ",atem]];
	(*then collect all different maps and remove duplicate cases*)
	atem=GatherBy[atem//DeleteCases[#,{"forbidden",_}]&,First];
	Return[atem[[All,1,2]]];(*for each group, we only need the first element, the remaining map will give the same results*)

	Return[alist];
];


Options[ApplyPolesToAlgLetter1]={Sector->{},AdRep->{},deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True,NoPathInfo->True};
ApplyPolesToAlgLetter1[opoles_,letters_,reform_,krep_,OptionsPattern[]]:=Catch@Module[{start,poles,polepath,letterpath,kletters={},tem,ntem,glist,glistd,Glist,ll,nletters={},sletters={},var,assoc,rep,rep1,gr,ga,isp,pos,path,dis},
If[letters==={},Return[{nletters,kletters,sletters}]];(*when the input is null, return empty set*)
start=SessionTime[];
(*first find all possible pole maps for the algebraic letters*)
tem=Reap[
	Do[
		Sow[(letters[[i,1]]//Cases[#,_G,Infinity]&//DeleteDuplicates)/.{G[a_,b_]:>GramMat[a,b,krep]}//Cases[#,Subscript[x,_],Infinity]&//DeleteDuplicates//Sort];(*all Baikov variables involved in an algebraic letter*)		
	,{i,1,Length[letters]}]
][[2]];
If[tem=!={},var=tem[[1]]//DeleteDuplicates,var={}];
assoc=Association[Table[Rule[var[[i]],FindPoleMaps[opoles,var[[i]]]],{i,1,Length[var]}]];(*the association that associate a baikov variable list to a list of pole maps*)
If[OptionValue[deBug],PrintTemporary["all poles map found ",SessionTime[]-start];PrintTemporary[Keys[assoc]]];

Do[
	If[OptionValue[deBug],PrintTemporary[SessionTime[]-start]];
	(*letterpath=letters[[i,3]];*)
	ga=((Cases[letters[[i,1]],Power[z_,1/2]->z,Infinity]//DeleteDuplicates)[[1]])/.{G[{},{}]->1};(*algebraic part of letter*)
	gr=(letters[[i,1]][[1]]//Numerator)/.{Power[_,1/2]->0};(*rational part of letter*)
	tem=({gr,ga}//Cases[#,_G,Infinity]&)/.{G[a_,b_]:>GramMat[a,b,krep]};
	var=Cases[tem,Subscript[x,_],Infinity]//DeleteDuplicates//Sort;(*get all Baikov variables involved in this letter*)
	glist=Cases[{ga},_G,Infinity]//DeleteDuplicates;(*grams under square root*)
	glistd=Cases[{letters[[i,2]]},_G,Infinity]//DeleteDuplicates;(*grams corresponding to the multiplication of numerator and denominator of algebraic letter*)
	If[OptionValue[deBug],PrintTemporary["ga: ",ga];PrintTemporary["gr: ",gr];PrintTemporary["var: ",var]];
		
	If[OptionValue[deBug],PrintTemporary[SessionTime[]-start]];
	If[var==={},
		(*when this letter consists only of pure kinematics gram*)		
		If[!(AdmissiblePoleQ[glist,{},OptionValue[PermSq],krep][[1]]),Continue[],AppendTo[kletters,{letters[[i,1]],{},letters[[i,4]]}]];
		Continue[],
		
		pos=assoc[var];
		poles=Table[If[Head[pos[[i]]]===List,{#,pos[[i]]}&/@MergePoles[opoles[[1]][[pos[[i]]]]],{#,{pos[[i]]}}&/@MergePoles[{opoles[[1,pos[[i]]]]}]],{i,1,Length[pos]}]//Flatten[#,1]&;(*the pole maps needed*)
	];
	(*Print["var: ",var];
	If[var===(Subscript[x,#]&/@Range[11]),Print[glist];Throw[poles]];(*deBug use only//////////////////////////////////////////////*)*)
	If[OptionValue[PermSq]==={},
		(*in this case, we will give a warning messge. This approach is not used now*)
		Message[ApplyPolesToAlgLetter1::permsq];
		Throw[{{},{},{}}]
		,
		(*In this case, we use information of square roots to constrain the form*)
		ntem={};
		Do[
			tem=AdmissiblePoleQ[glist,poles[[j,1]],OptionValue[PermSq],krep];
			If[!(tem[[1]]),
				Continue[],
				If[MemberQ[ntem,{poles[[j,1]]//Flatten//Sort,tem[[2]]}],
					Continue[](*cases that has been considered*),
					AppendTo[ntem,{poles[[j,1]]//Flatten//Sort,tem[[2]]}];
				];
				AppendTo[nletters,{letters[[i,1]],poles[[j,1]],letters[[i,4]]}]
			]
		,{j,1,Length[poles]}]
	];
	If[OptionValue[deBug],PrintTemporary[SessionTime[]-start]]
,{i,1,Length[letters]}];
Throw[{nletters,kletters,sletters}];
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
		If[Not@ContainsAll[poles[[i,1]],poles[[j,1]]]&&Not@ContainsAll[poles[[j,1]],poles[[i,1]]],Continue[]];
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
		Message[ApplyPolesToAlgLetter2::permsq];
		Throw[{}]
		,
		(*In this case, we use additional information to constrain the form*)
		pos=assoc[var];
		poles=Table[If[Head[pos[[i]]]===List,{#,pos[[i]]}&/@MergePoles[opoles[[1]][[pos[[i]]]]],{#,{pos[[i]]}}&/@MergePoles[{opoles[[1,pos[[i]]]]}]],{i,1,Length[pos]}]//Flatten[#,1]&;(*the pole maps needed*)
		ntem={};
		Do[
			tem=AdmissiblePoleQ[glist,poles[[j,1]],OptionValue[PermSq],krep];
			If[!(tem[[1]]),
				Continue[],
				If[MemberQ[ntem,{poles[[j,1]]//Flatten//Sort,tem[[2]]}],
					Continue[](*cases that has been considered*),
					AppendTo[ntem,{poles[[j,1]]//Flatten//Sort,tem[[2]]}];
				];
				AppendTo[adpole,poles[[j,1]]]
			]
		,{j,1,Length[poles]}]
	];
	(*If[OptionValue[deBug],Print["adpole: ",adpole]];*)
	If[Length[adpole]<2,Continue[]];(*if there are less than 2 admissible poles, then we can't construct from this gram*)
	AppendTo[nletters,ApplyPoleToQM[adpole,gram[[i]],krep]//DeleteDuplicatesBy[#,First]&//DeleteCases[#,_?(FreeQ[#,Power[_,1/2]]&)]&];
,{i,1,Length[gram]}];
Return[nletters//Flatten[#,1]&//DeleteDuplicatesBy[#,First]&];
];


Options[ApplyPolesToAlgLetter3]={Sector->{},AdRep->{},deBug->False,PermSq->{}};
ApplyPolesToAlgLetter3[opoles_,letters_,krep_,OptionsPattern[]]:=Catch@Module[{start,poles,tempole,temvar,tem,ntem,glist,glistd,Glist,ll,nletters={},var,rep,rep1,gr,ga,isp,pos,path,dis},
If[letters==={},Return[{nletters,{},{}}]];(*when the input is null, return empty set*)
start=SessionTime[];
(*first find all possible pole maps for the algebraic letters*)
poles=opoles[[1]];
Do[
	var=(letters[[i,1]]//Cases[#,_G,Infinity]&//DeleteDuplicates)/.{G[a_,b_]:>GramMat[a,b,krep]}//Cases[#,Subscript[x,_],Infinity]&//DeleteDuplicates//Sort;(*all Baikov variables involved in an algebraic letter*)
	If[var==={},Continue[]];(*this case has been considered already*)
	Do[
		temvar=Cases[poles[[k]],Subscript[x,_],Infinity]//DeleteDuplicates;
		If[ContainsAll[temvar,var],Continue[]];(*this case has been considered already*)
		temvar=Complement[var,temvar];(*the additional variables present*)
		If[Length[temvar]>2,Continue[]];(*no more than 2 additional variables*)
			
		If[OptionValue[deBug],PrintTemporary[SessionTime[]-start];Print["temvar: ",temvar]];
		(*letterpath=letters[[i,3]];*)
		ga=((Cases[letters[[i,1]],Power[z_,1/2]->z,Infinity]//DeleteDuplicates)[[1]])/.{G[{},{}]->1};(*algebraic part of letter*)
		gr=(letters[[i,1]][[1]]//Numerator)/.{Power[_,1/2]->0};(*rational part of letter*)
		glist=Cases[{ga},_G,Infinity]//DeleteDuplicates;(*grams under square root*)
		glistd=Cases[{letters[[i,2]]},_G,Infinity]//DeleteDuplicates;(*grams corresponding to the multiplication of numerator and denominator of algebraic letter*)
		If[OptionValue[deBug],PrintTemporary["ga: ",ga];PrintTemporary["gr: ",gr];PrintTemporary["var: ",var]];
			
		If[OptionValue[deBug],PrintTemporary[SessionTime[]-start]];
		tempole=MergePoles[{Join[{Join[poles[[k,1]],Thread@Rule[temvar,0]]},Drop[poles[[k]],1]]}][[1]];
		If[OptionValue[deBug],Print["tempole: ",tempole]];
		tem=AdmissiblePoleQ[glist,tempole,OptionValue[PermSq],krep];
		If[!(tem[[1]]),
			Continue[],
			AppendTo[nletters,{letters[[i,1]],tempole,letters[[i,4]]}]
		]
	,{k,1,Length[poles]}]		
,{i,1,Length[letters]}];

Throw[{nletters,{},{}}];
];


Options[AllAlgLettersSupplement]={deBug->False,PermSq->{}};
AllAlgLettersSupplement[poles_,algletter_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,part,pathdis,permsq,looppath,kinepath,b},
start=SessionTime[];
l=Length[algletter];
permsq=OptionValue[PermSq];
permsq=Join[permsq,Table[Times@@(PerfectSquareSplit[permsq[[i]]*permsq[[j]]][[2]]),{i,1,Length[permsq]},{j,i+1,Length[permsq]}]//Flatten]//DeleteDuplicates;
(*Export[OptionValue[tmpDir]<>"input.mx",{poles,reform,krep}];
part=Partition[algletter,UpTo[Quotient[Length[algletter],OptionValue[NThreads]]+1]];
Table[Export[OptionValue[tmpDir]<>"list"<>ToString[i]<>".mx",part[[i]]],{i,1,Length[part]}];
(*export the data to parallelly run them in terminal*)*)
Monitor[Do[
(*delete duplicates according to the first (Letter form) and second (poles info) term of the unit.*)
	(*If[algletter[[i,-1]]=!={1,2,3,4,6,8},Continue[]];*)
	AppendTo[result,ApplyPolesToAlgLetter3[poles,algletter[[b,1]],krep,PermSq->permsq]//Flatten[#[[{1,2}]],1]&//DeleteDuplicatesBy[#,#[[{1,2}]]&]&];
,{b,1,Length[algletter]}],ProgressIndicator[b,{1,Length[algletter]}]];
Print["        session time: ",SessionTime[]-start];
Print["        Substituting poles into expressions..."];
tem=Flatten[result,1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;(*remove duplicate expression again*)
(*DistributeDefinitions[tem,ApplyPoleToLetter];*)
(*Return[tem];*)
eresult=Monitor[Table[{TimeConstrained[ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],1000,$Failed],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
(*eresult=ParallelTable[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]},Method->"CoarsestGrained"];*)
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["        session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{eresult,spresult}]];
Return[{eresult}];
];


Options[AllAlgLettersSupplementPL]={deBug->False,PermSq->{}};
AllAlgLettersSupplementPL[poles_,algletter_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,part,permsq,b},
start=SessionTime[];
l=Length[algletter];
permsq=OptionValue[PermSq];
permsq=Join[permsq,Table[Times@@(PerfectSquareSplit[permsq[[i]]*permsq[[j]]][[2]]),{i,1,Length[permsq]},{j,i+1,Length[permsq]}]//Flatten]//DeleteDuplicates;
(*Export[OptionValue[tmpDir]<>"input.mx",{poles,reform,krep}];
part=Partition[algletter,UpTo[Quotient[Length[algletter],OptionValue[NThreads]]+1]];
Table[Export[OptionValue[tmpDir]<>"list"<>ToString[i]<>".mx",part[[i]]],{i,1,Length[part]}];
(*export the data to parallelly run them in terminal*)*)
DistributeDefinitions[poles,algletter,krep,permsq,ApplyPolesToAlgLetter3,AllAlgLettersSupplementPL];
SetSharedVariable[result,spresult];
ParallelDo[
(*delete duplicates according to the first (Letter form) and second (poles info) term of the unit.*)
	(*If[algletter[[i,-1]]=!={1,2,3,4,6,8},Continue[]];*)
	AppendTo[result,ApplyPolesToAlgLetter3[poles,algletter[[i,1]],krep,PermSq->permsq]//Flatten[#[[{1,2}]],1]&//DeleteDuplicatesBy[#,#[[{1,2}]]&]&];
,{i,1,Length[algletter]}];
Print["        session time: ",SessionTime[]-start];
Print["        Substituting poles into expressions..."];
tem=Flatten[result,1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;(*remove duplicate expression again*)
(*DistributeDefinitions[tem,ApplyPoleToLetter];*)
(*Return[tem];*)
eresult=Monitor[Table[{TimeConstrained[ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],1000,$Failed],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
(*eresult=ParallelTable[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]},Method->"CoarsestGrained"];*)
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["        session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{eresult,spresult}]];
Return[{eresult}];
];


Options[AllAlgLetters]={deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True};
AllAlgLetters[poles_,algletter_,reform_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,eresult2,eresult3,a,b,permsq},
start=SessionTime[];
l=Length[algletter];
permsq=OptionValue[PermSq];
permsq=Join[permsq,Table[Times@@(PerfectSquareSplit[permsq[[i]]*permsq[[j]]][[2]]),{i,1,Length[permsq]},{j,i+1,Length[permsq]}]//Flatten]//DeleteDuplicates;
Print["totally ",l," sectors need to be analyzed"];
Print["analyzing first type construction..."];
Monitor[
	Do[
		If[OptionValue[deBug],PrintTemporary["subset: ",algletter[[a,2]]," session time: ",SessionTime[]-start]];
		tem=ApplyPolesToAlgLetter1[poles,algletter[[a,1]],reform,krep,PermSq->permsq,PathDis->OptionValue[PathDis],KinePath->OptionValue[KinePath],LoopPath->OptionValue[LoopPath]];
		ntem=Flatten[tem[[{1,2}]],1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;
		(*delete duplicates according to the first (Letter form) and second (poles info) term of the unit.*)
		AppendTo[result,ntem];
		If[OptionValue[deBug],AppendTo[spresult,tem[[3]]]]
	,{a,1,Length[algletter]}]
,ProgressIndicator[a,{1,Length[algletter]}]];
Print["session time: ",SessionTime[]-start];
Print["Substituting poles into expressions..."];
tem=Flatten[result,1]//DeleteDuplicatesBy[#,#[[{1,2}]]&]&;(*remove duplicate expression again*)
eresult=Monitor[Table[{TimeConstrained[ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],1000,$Failed],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["session time: ",SessionTime[]-start];
Print["analyzing second type construction..."];
eresult2=ApplyPolesToAlgLetter2[poles,{},reform,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
Print["analyzing supplemental poles..."];
eresult3=AllAlgLettersSupplement[poles,algletter,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{Join[eresult,eresult3[[1]]],eresult2,spresult}]];
Return[{Join[eresult,eresult3[[1]]],eresult2}];
];


Options[AllAlgLettersPL]={deBug->False,PermSq->{},PathDis->False,KinePath->False,LoopPath->True};
AllAlgLettersPL[poles_,algletter_,reform_,krep_,OptionsPattern[]]:=Module[{start,l,tem,ntem,result={},spresult={},eresult,eresult2,eresult3,part,pathdis,permsq,looppath,kinepath,b},
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
eresult=Monitor[Table[{TimeConstrained[ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],1000,$Failed],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]}],ProgressIndicator[b,{1,Length[tem]}]];
(*eresult=ParallelTable[{ApplyPoleToLetter[tem[[b,2]],tem[[b,1]],krep],tem[[b,1]],tem[[b,2]]},{b,1,Length[tem]},Method->"CoarsestGrained"];*)
eresult=DeleteCases[eresult,_?(FreeQ[First[#],Power[_,1/2]]&)];(*remove expression without square roots*)
Print["session time: ",SessionTime[]-start];
Print["analyzing second type construction..."];
eresult2=ApplyPolesToAlgLetter2[poles,{},reform,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
Print["analyzing supplemental poles..."];
eresult3=AllAlgLettersSupplementPL[poles,algletter,krep,PermSq->permsq];
Print["session time: ",SessionTime[]-start];
If[OptionValue[deBug],Return[{Join[eresult,eresult3[[1]]],eresult2,spresult}]];
Return[{Join[eresult,eresult3[[1]]],eresult2}];
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


Options[RemoveDegenerateAlgLetter]={deBug->False};
RemoveDegenerateAlgLetter[alglist_,krep_,OptionsPattern[]]:=Module[{list,numrep,var,nvar,nlist,coord,num,numnum,den,dennum,pos},
list=alglist/.{Log[z_]:>z};
var=Variables[list]//Sort;
nvar=Variables[var/.krep]//Sort;
numrep=Table[Thread@Rule[var,var/.krep/.(Thread@Rule[nvar,Table[Prime[300*i+3*j],{j,1,Length[nvar]}]])],{i,1,3}];
den=(Denominator/@list);
dennum=Table[den/.numrep[[i]]//Factor,{i,1,3}]//Transpose;
pos=Position[dennum,{0,0,0}];
num=(Numerator/@list);
numnum=Table[num/.numrep[[i]]//Factor,{i,1,3}]//Transpose;
pos=Join[pos,Position[numnum,{0,0,0}]];
If[OptionValue[deBug],Print["pos: ",pos]];
list=Delete[list,pos];(*remove terms that are actually singular, this can happen in degenerate case which is common in massless one-loop n-gon or 4d kinematics*)
If[OptionValue[deBug],Print["numrep: ",numrep]];
Return[Log/@(list)];
];


ToRational[alg_]:=Times@@NumeratorDenominator[alg[[1]]]//Expand//Factor;

Options[GenerateNumReal]={"specialkinematics"->{}};
GenerateNumReal[sq_,ovar_,OptionsPattern[]]:=Module[{var,rep,tem,k=1},
var=ovar/.OptionValue["specialkinematics"]//Variables;
rep=Thread@Rule[var,RandomChoice[{1,-1}]*RandomPrime[{10^5,10^6},Length[var]]];
tem=sq/.{Power[a_,1/2]:>a}/.OptionValue["specialkinematics"];
While[Length[Values[rep]//Union]<Length[var]||((tem/.rep)<0&&k<1000),
rep=Thread@Rule[var,RandomChoice[{1,-1}]*RandomPrime[{10^5,10^6},Length[var]]];
k=k+1;
];
Return[Thread@Rule[ovar,ovar/.OptionValue["specialkinematics"]/.rep]];
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


Options[GetAllAlgIndepLetter]={deBug->False,"specialkinematics"->{}};
GetAllAlgIndepLetter[algletter_,OptionsPattern[]]:=Module[{gr,var,sq,rep,pos,tem,ind={}},
gr=GatherBy[algletter,(Cases[#,Power[_,1/2],Infinity]//DeleteDuplicates)&];
Print["totally ",Length[gr]," different square roots"];
var=algletter/.{Log[z_]:>z}//Variables;
Do[
sq=Cases[gr[[i]],Power[_,1/2],Infinity]//DeleteDuplicates;
If[sq==={},Continue[]];
tem=gr[[i]]//SortBy[#,LeafCount]&;
If[OptionValue[deBug],Print["square root: ",Short[sq,20]];Print["length: ",Length[tem]]];
rep=Table[GenerateNumReal[sq[[1]],var,"specialkinematics"->OptionValue["specialkinematics"]],{i,1,5}];
pos=Commonest[SearchIndepLetterNum[tem,#]&/@rep][[1]];
AppendTo[ind,tem[[pos]]]
,{i,1,Length[gr]}];
Return[ind];
];

Options[FindLetterLinearRelation]={deBug->False,"specialkinematics"->{}};
FindLetterLinearRelation[basis_,target_,OptionsPattern[]]:=Module[{gr,var,sq,rep,rel},
gr=Append[basis,target];
var=gr/.{Log[z_]:>z}//Variables;

sq=Cases[{gr[[-1]]},Power[_,1/2],Infinity]//DeleteDuplicates;
If[OptionValue[deBug],Print["square root: ",Short[sq,20]]];
rep=Table[GenerateNumReal[sq[[1]],var,"specialkinematics"->OptionValue["specialkinematics"]],{i,1,3}];
rel=Commonest[Quiet[FindIntegerNullVector[ReplaceAll[gr/.{Log[z_]:>Log[Abs[z]]},#],30]]&/@rep];
Return[rel];
];


LetterInfo[letter_,algresult_,polestructure_]:=Module[{nl,pos,tem,tempos,apos,i},
If[Head[letter]=!=Log,nl=Log[letter],nl=letter];
If[FreeQ[nl,Power[_,1/2]],
(*this is a rational letter*)
pos=Position[polestructure[[All,3]],nl[[1]]][[All,1]];
Print["This letter can be generated from following sectors: ",polestructure[[pos,-1]]];
tem=polestructure[[pos,1]]//Flatten[#,1]&;
apos=Reap[
	Do[
		tempos=Join[Position[tem[[i]],nl[[1]]],Position[tem[[i]],-nl[[1]]//Factor]];
		If[tempos=!={},Sow[tem[[i,-1]]]]
	,{i,1,Length[tem]}]
][[2]];
If[apos==={},Print["The explicit path is not determined, loop up the following position of output of PolesAnalyze[]: ",pos],apos=apos[[1]]];
Print["The corresponding path is: ",apos];
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


GetIndepAlgLetters[result_,rllist_]:=Module[{algcand,algselect,algind},
	algcand=Cases[result,Log[_?(FreeQ[#,G]&)],Infinity]//DeleteCases[#,_?(FreeQ[#,Power[_,1/2]]&)]&//DeleteSameAlgLetter;
	algselect=SelectAlgLetter[algcand//RegularizeSquareRoots,rllist];
	algind=GetAllAlgIndepLetter[algselect];
	Return[algind];
];


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

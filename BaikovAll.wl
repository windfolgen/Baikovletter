(* ::Package:: *)

(*
Package Name: BaikovAll
Version: 1.0.1
Description: Generate all Baikov representation in an integral family by analyzing the Gram matrix and integrate variables one-by-one. It also provides some tools to analyze the Baikov representation.
The output will be in the form:
{{{variables integrated out},{{Gram1,power1},{Gram2, power2},const}},...} The results will be graded by the number of variables which have been integrated out.
The corresponding Baikov representation will be like const*Power[Gram1,power1]*Power[Gram2,power2]*..., 
*)
(*
If you find any bug or suggest for this program, please contact: xhjiang@itp.ac.cn
*)


BeginPackage["Baikov`"];


x::usage="the variable for Baikov representation";
y::usage="the variable for Baikov representation";
G::usage="The function for Feynman integrals";
j::usage="The function for Feynman integrals in LiteRed";
\[Epsilon]::usage="the parameter for dimentional regularization";
R::usage="Possible alias for square roots.";
Rlist::usage="Possible list for square roots";
Unprotect[Rlist];
Rlist={};
Protect[x,y,G,j,\[Epsilon],R,Rlist];

(*some global path*)
$SingularFilePath::usage="Path for temporary files of Singular.";
$SingularPath::usage="Path for Singular excutable file.";

(*dimension calculation*)
GetDimension::uerr="u should in the form p1^a1...pn^an.`1`";
GetDimension::oerr="Omega should be a matrix.";
GetDimension::time="Computation not completed within given time. You can adjust the Option of Time to a larger value. Default value is `1`s.";
GetDimension::usage="GetDimension[u_or_Om,z] calculates the dimension of a system from the u function or Omega matrix with respect to z. There are two method. One is using Singular to calculate the proper critical points from u, in which case z is a list of variables. Another is using Omega from intersection calculation in which case z is the last variable when performing intersection number calculation in some order.";

SearchDimension::usage="SearchDimension[list,zl,replist,Eps,sw:0] searches the dimension for a given integral family u, list is dlog u for all variables zl, replist is the replacement rule for irrelevant kinematics variables. Eps is the precision constraint of the numerical value. sw controls whether solving the system by Solve or by NSolve. ";

MPermute::usage="MPermute[exp,sym] is basically the same with Permute[] except that it can permute list whose length is smaller than the length sym requires.";
(*polynomial operation*)
SymGroupStabilizer::usage=="SymGroupStabilizer[l] generates the permutation group which keeps the list l invariant.";

PolySym::usage="PolySym[p,v,s] searches for the permutation symmetries of v which keeps polynomial p invariant. The first s elements of v have different weight with the rest which means that we don't permute these two sets of elements with each other. 
PolySym[{p1,v1,s1},{p2,v2,s2}] finds the map of v1 to v2 that takes p1 to p2.";

PolySymF::usage="PolySymF[{p1,v1,s1},{p2,v2,s2}] finds the map of v1 to v2 that takes p1 to p2.";

PolySymCheck::usage="PolySymCheck[p,q,rule] checks whether two polynomials p and q are equivalent under map rule";

PolyOrder::usage="PolyOrder[mon1,mon2] gives the 'larger' polynomial according to the order we define. Here we use the convention that the term with larger corresponding binary number will be kept. For example, x1*x3 corresponds to 101, x2*x4 corresponds to 0101, so the latter is larger.";

PolyEqual::usage="PolyEqual[monomial,sym] gives the equivalent monomial according to the symmetry given. The convention about which term to keep is the same with PolyFold[]. sym is a list of replace rules.";

PolyFold::usage="PolyFold[poly,sym] folds a polynomial according to the symmetry given (sym should be a list of replacement). However, we should specify the way which term should keep. For example, if x1*x3 is symmetric to x2*x4 we should specify which term we will keep. Here we use the convention that the term with larger corresponding binary number will be kept. For example, x1*x3 corresponds to 101, x2*x4 corresponds to 0101, so the latter is larger.";
(*baikov representation *)
SProd::usage="SProd[p1,p2] defines a scalar product p1\[CenterDot]p2";

ExpandSP::usage="ExpandSP[exp] expand the arguments of SProd[] function.";

CayleyMengerForm::usage="CayleyMengerForm[list,rep] gives the Cayley-Menger form of a gram determinant of momenta list and rep is the variable replacement.";

GramDet::usage="GramDet[{p1,p2,p3}] calculate the determination of matrix {{p1.p1,p1.p2,p1.p3},{p2.p1,p2.p2,p2.p3},{p3.p1,p3.p2,p3.p3}}. And GramDet[{p1,p2},{q1,q2}] calculates {{p1.p1,p1.p2},{p2.p1,p2.p2}}.";

PD::usage="PD[loopm,extm,f] is used to represent a Feynman denominator. For example, PD[l1+l2,k1,-m2]=(l1+l2+k1)^2-m2";

LPD::usage="LPD[loopm,extm,f] is used to represent a linear Feynman denominator. For example, PD[l1,k1,-m2]=l1.k1-m2";

BaikovTrans::usage="BaikovTrans[dlist] gives the transformation between scalar product SProd[] and Baikov variables xi.";
BaikovTrans::marg="The freedom of denominator list `1` doesn't equal to `2`: the freedom of the integral measure.";

BaikovRep::usage="BaikovRep[dlist,llist,elist] calculates the Standard Baikov representation of a given denominator list.";
(*integrand in Baikov*)
ToCanonicalForm::usage="ToCanonicalForm[u,FI] transforms the FI which takes the form G[1,{...}] or j[tag,1,0,...] to the canonical form under Baikov representation. This canonical form will only have xi in the denominator has at most power 1 at the price of introducing some polynomials in the denominator.";

Trans2FI::"usage"="Trans2FI[num,list,plist,f] transfers an expression which consists of xi and kinematic variables to the form of FI such as G[{1,1,1,1,0,0}] with kinematic coefficient. list is a list of xis. plist is the power of xis in denominator. num is the numerator of the expression. f is the head of a function.";
(*generate all Baikov recursively*)
LocatePos::usage="LocatePos[matlist,v] locates v in a list of matrix and give its position and coefficients in the upper triangle part of the matrix.";

ReArrangeGram::usage="ReArrangeGram[g,pos,ratio,flag] rearrange the Gram matrix g according to pos, ratio and flag so the variable only appears in one position. flag=-1 means the same row and 1 means the same column.";
ReArrangeGram::err="flag can only be -1 or 1. flag: `1`";

ReduceMat::usage="ReduceMat[{g,power},pos,coef] reduces the Gram Matrix g by using recursive formula. coef is the coefficient of the variable to be reduced in the matrix.";

SimplifyGram::usage="SimplifyGram[gl,rep] simplifies the Gram matrix list provided. rep is the kinematics replacement rule.";

ZeroSectorMatQ::usage="ZeroSectorMatQ[gl,rep] decides whether this sector is zero no matter whether it has those propagators or not by its u.";
ZeroSectorMatQ::warning="This sector doesn't contain any baikov variables!";

ReduceRep::usage="ReduceRep[pl,v,rep] finds whether v in pl is reducible. If it is true, then return the results after reduction. If not, return False.";
ReduceRep::devoid="There is no variable `1` in the expression. There must be something wrong.";

AllSectorBaikovMat::err="there must be something wrong with path `1` and the variable `2`.";
AllSectorBaikovMat::usage="AllSectorBaikovMat[list,rep] generate all Baikov representation from the list provided, list is in the form {{{Gram1,power1},{Gram2,power2},...},const}.";

AllSectorBaikovMatC::usage="AllSectorBaikovMatC[resultmat,var,krep] continues the calculation of AllSectorBaikovMat[]. resultmat is the output of AllSectorBaikovMat[]. var is the list of all Baikov variables. krep is the replacement rule of scalar products including loop momenta. This is designed for cases where we need to continue the recursive derivation of resultmat.";

ExtractLoopOrder::usage="ExtractLoopOrder[krep] extracts the number of loops involved from the replacement rule krep.";
GetBaikovMatRep::usage="GetBaikovMatRep[result,var,n] gets the Baikov representation for given variables var. n is the number of Baikov variables. For example, var={1,2,3}. Option \"ExcVar\" gives the variables that must not be integrated out, that is, they will keep in the representation. It can be {9,10} for example. Option \"looporder\" is the loop number of the integral family, e.g. for two loop family it is 2.";
GetBaikovMatRep::zero="This is a zero sector.";
GetBaikovMatRep::err="All variables have been integrated out.";


AllSubSector::usage="AllZeroSector[zeroset] generates all sector from a set of sectors. Note that a subsector of a sector is again a sector. This is why we need to generate all.";
GetMatZeroSector::usage="GetMatZeroSector[allbaikovlist,n,isp] extracts the information about zero sectors from the output of AllSectorBaikov[]. n is the number of Baikov variables and isp is the list of isp's like {8,9}.";

GramMat::usage="GramMat[l1,l2,rep] gives the Gram matrix for G[l1,l2]. rep is kinematics replacement rule.";

Gram2Poly::usage="Gram2Poly[exp,krep] transforms all the G[l1,l2] expressions to polynomials. krep is the kinematics replacement rule.";


Begin["`Private`"]


Options[GetDimension] = {Method -> "Singular", deBug -> False, fileDir
     -> $SingularFilePath, fileName 
    -> "gen_dimension.sing", SingCommand -> $SingularPath,
     Time -> 600}; 

GetDimension[uorom_, z_, OptionsPattern[]] :=
    Module[{file, stream, pl, var, rep, zl, den, flag},
        If[OptionValue[Method] === "Det",
            Goto[detofOm]
        ];
        If[OptionValue[deBug],
            Print["Using Singular to get the dimension"]
        ];
        If[Head[uorom] =!= Times && Head[uorom] =!= Power,
            Message[GetDimension::uerr,uorom];
            Return[$Failed]
        ];
        If[(uorom/.{Power[0,_]->0})===0,Return[$Failed]];(*cases uorom contains Power[0,_], this means there is a singularity in the representation*)
        pl = (D[Log[uorom], #]& /@ z) // Together // Numerator;
        If[OptionValue[deBug],
            Echo[pl, "pl:"]
        ];
        var = Complement[Variables[{pl, uorom /. Power[y_, _] :> y}],
             z];
        If[OptionValue[deBug],
            Echo[var, "var:"]
        ];
        rep = Thread @ Rule[var, Table[Prime[3 i + 1], {i, 1, Length[
            var]}]];
        pl = (pl /. rep /. {Subscript[a_, b_] :> ToExpression[ToString[
            a] <> ToString[b]]}) // Factor;
        If[Head[uorom] === Times,
            den = uorom // Cases[#, Power[y_, _] :> y, {1}]&
            ,
            den = uorom // Cases[#, Power[y_, _] :> y, {0}]&
        ];
        den = Times @@ ((den /. rep /. {Subscript[a_, b_] :> ToExpression[
            ToString[a] <> ToString[b]]}) // Factor);
        zl = Prepend[z /. {Subscript[a_, b_] :> ToString[a] <> ToString[
            b]}, "x0"];
        file = OptionValue[fileDir] <> OptionValue[fileName];
        If[!FileExistsQ[file],
            CreateFile[file]
        ];
        stream = OpenWrite[file, CharacterEncoding -> "UTF-8"];
        Put[OutputForm["option(redSB);"], OutputForm["option(redTail);"
            ], stream];
        Put[OutputForm["ring r=0,(" <> StringJoin[Riffle[zl, ","]] <>
             "),dp;"], stream];
        Put[OutputForm["poly p=" <> ToString[den, InputForm] <> ";"],
             stream];
        Do[Put[OutputForm["poly g" <> ToString[i] <> "=" <> ToString[
            pl[[i]], InputForm] <> ";"], stream], {i, 1, Length[pl]}];
        Put[OutputForm["poly h=x0*p-1;"], stream];
        Put[OutputForm["ideal I=" <> StringJoin[Riffle[Table["g" <> ToString[
            i], {i, 1, Length[pl]}], ","]] <> ",h;"], stream];
        Put[OutputForm["int k=vdim(groebner(I));"], OutputForm["write(\":w "
             <> OptionValue[fileDir] <> "gen_result.txt\",k);"], OutputForm["exit;"
            ], stream];
        Close[stream];
        flag = TimeConstrained[RunProcess[{OptionValue[SingCommand], 
            file}], OptionValue[Time]];
        If[flag === $Aborted,
            Message[GetDimension::time,OptionValue[Time]];
            Return[$Failed]
        ];
        If[flag["StandardError"] =!= "",
            Return[flag["StandardError"]]
        ];
        Return[Get[OptionValue[fileDir] <> "gen_result.txt"]];
        Label[detofOm];
        If[OptionValue[deBug],
            Print["Using det of Omega to get the dimension"]
        ];
        If[!MatrixQ[uorom],
            Message[GetDimension::oerr];
            Return[$Failed]
        ];
        den = Det[uorom] // Together // Numerator;
        zl = {z} // Flatten;
        Return[Exponent[den, zl[[1]]]];
        
    ]; 

SearchDimension[list_, zl_, replist_, Eps_, sw_:0] :=
    Module[{nlist, dlist, v, rep, eqs, tem, s, ll = {}},
        nlist = Numerator[list // Together];
        dlist = Denominator[list // Together];
        v = Complement[Variables[nlist], zl];
        If[sw == 0 || sw == 1,
            rep = replist
            ,
            rep = Thread @ Rule[v, Table[Prime[i + 5], {i, 1, Length[
                v]}]]
        ];
        Print[rep];
        eqs = Thread @ Equal[nlist /. rep, Table[0, {i, 1, Length @ nlist
            }]];
        If[sw == 0,
            s = Solve[eqs, zl] // FullSimplify;
            tem = Length[s];
            dlist = Table[dlist /. rep /. s[[i]], {i, 1, tem}] // FullSimplify
                ;
            Print[dlist];
            Do[
                If[FreeQ[dlist[[i]], x_ /; (x == 0), 1],
                    AppendTo[ll, i]
                ]
                ,
                {i, tem}
            ]
            ,
            s = NSolve[eqs, zl];
            tem = Length[s];
            dlist = Table[dlist /. rep /. s[[i]], {i, 1, tem}];
            Do[
                If[FreeQ[dlist[[i]], x_ /; (Abs[x] < Eps)],
                    AppendTo[ll, i]
                ]
                ,
                {i, tem}
            ]
        ];
        (*Print[Compress[{dlist[[ll]],s[[ll]]}]];*)
        Print["length: ", Length[s], " ll: ", ll];
        Print[s[[ll]]];
        DeleteDuplicates[Floor[Values[s[[ll]]], 10^(-7)]] // N
    ]; 


Options[MPermute] = {ReVerse -> False}; 

MPermute[exp_, sym_, OptionsPattern[]] :=
    Module[{n, l, rep},
        n = Max[Cases[sym, _Integer, Infinity], {1}];
        l = Table[i, {i, 1, n}];
        rep =
            If[OptionValue[ReVerse],
                Thread @ Rule[Permute[l, sym], l]
                ,
                Thread @ Rule[l, Permute[l, sym]]
            ];
        exp /. rep
    ]; 


SymGroupStabilizer[l_] :=
    Module[{pi, permute},
        pi = Values[PositionIndex[l]] // DeleteCases[#, {_}]&;(*find 
            elements repeating in a list*)
        If[pi == {},
            Return[PermutationGroup[{}]]
        ];(*if not found, then this list cannot be permuted*)
        pi = Flatten[Partition[#, 2, 1]& /@ pi, 1];
        permute = Table[Cycles[{pi[[i]]}], {i, 1, Length[pi]}];
        PermutationGroup[permute]
    ]; 

Options[PolySym] = {deBug -> False}; 

PolySym[p_, v_, s_,OptionsPattern[]] :=
    Module[{pp, l, c, ml, tab, cha, stab, random, benchmark},
        pp = p // ExpandAll;
        l =
            If[Head[pp] =!= Plus,
                {pp}
                ,
                List @@ pp
            ];
        c = l /. (Thread @ Rule[v, 1]);
        (*Print["c: ",Short[c]];(*DeBug*)*)
        If[!NumberQ[c[[1]]],
            Print[pp];
            Abort[]
        ];
        ml = l / c;
        tab = Table[Exponent[ml[[i]], #]& /@ v, {i, 1, Length[ml]}];
        (*using a characteristic matrix 'tab' to describe the polynomial
            *)
        cha =
            Table[
                Abs[c] . tab[[All, i]] *
                    Power[
                        10
                        ,
                        If[i <= s,
                            0
                            ,
                            3
                        ]
                    ]
                ,
                {i, 1, Length[v]}
            ];(*get the characteristic value for each variable, using
                 'c' as the weight. first s elements in v have different weight with 
                the rest elements*)
        (*Return[{Length[ml],c,tab,cha}];*)
        (*stab=Monitor[SymGroupStabilizer[cha],"Finding the stabilizer..."
            ];*)
        stab = SymGroupStabilizer[cha];(*calculate the group elements
             which don't change the characteristic value of each variables in v*)
            
        If[OptionValue[deBug],
            Print[cha]
        ];
        If[stab == PermutationGroup[{}],
            Return[{}]
        ];
        random = Table[Prime[i], {i, 100, 100 + Length[v] - 1}];(*generate
             different prime value for each variable*)
        benchmark = p /. (Thread @ Rule[v, random]);(*calculate the benchmark
             value of the polynomial*)
        l = GroupOrbits[stab, {random}, Permute][[1]];
        If[l == {},
            Return[{}]
        ];
        c =
            Monitor[
                Table[
                    If[(p /. (Thread @ Rule[v, l[[i]]])) == benchmark,
                        
                        i
                        ,
                        0
                    ]
                    ,
                    {i, 1, Length[l]}
                ]
                ,
                "Total length: " <> ToString[Length[l]](*<>" layer: "
                    <>ToString[i]*) ] // DeleteCases[#, 0]&;(*search the permutations which
                     don't change the polynomial*)
        stab = GroupElements[stab][[c]] // DeleteCases[#, Cycles[{}]]&;
        {stab, DeleteCases[Thread @ Rule[v, #]& /@ (Permute[v, #]& /@
             stab), a_ -> a_, {2}]}
    ]; 

Options[PolySymF]={deBug->False};

PolySymF[{p1_, v1_, s1_}, {p2_, v2_, s2_}, OptionsPattern[]] :=
    Catch@Module[{pp, l1, l2, c, ml, tab, cha1, cha2, tran, stab, tranl, random, benchmark, result={}
        },
        If[Length[v1] != Length[v2] || s1 != s2,
            Throw[{}]
        ];
        pp = p1 // ExpandAll;
        l1 =
            If[Head[pp] =!= Plus,
                {pp}
                ,
                List @@ pp
            ];
        pp = p2 // ExpandAll;
        l2 =
            If[Head[pp] =!= Plus,
                {pp}
                ,
                List @@ pp
            ];
        If[Length[l1] != Length[l2],
            Throw[{}]
        ];
        c = l1 /. (Thread @ Rule[v1, 1]);
        ml = l1 / c;
        tab = Table[Exponent[ml[[i]], #]& /@ v1, {i, 1, Length[ml]}];
            
        (*using a characteristic matrix 'tab' to describe the polynomial
             p1*)
        cha1 =
            Table[
                c . tab[[All, i]] *
                    Power[
                        10
                        ,
                        If[i <= s1,
                            0
                            ,
                            3
                        ]
                    ]
                ,
                {i, 1, Length[v1]}
            ];(*get the characteristic value for each variable, using
                 'c' as the weight. first s elements in v have different weight with 
                the rest elements*)
        c = l2 /. (Thread @ Rule[v2, 1]);
        ml = l2 / c;
        tab = Table[Exponent[ml[[i]], #]& /@ v2, {i, 1, Length[ml]}];
            
        cha2 =
            Table[
                c . tab[[All, i]] *
                    Power[
                        10
                        ,
                        If[i <= s2,
                            0
                            ,
                            3
                        ]
                    ]
                ,
                {i, 1, Length[v2]}
            ];
        If[OptionValue[deBug],Print["cha1: ",cha1];Print["cha2: ",cha2]];
        If[Sort[cha1] != Sort[cha2],
            Throw[{}]
            ,
            tran = FindPermutation[cha1, cha2]
        ];
        stab = GroupElements@SymGroupStabilizer[cha2];
        If[OptionValue[deBug],Print["tran: ",tran];Print["stablizer: ",stab]];
        tranl = PermutationProduct[tran,#]&/@ stab;(*generate all possible permutation*)
        random = Table[Prime[i], {i, 100, 100 + Length[v1] - 1}];(*generate
             different prime value for each variable*)
        benchmark = p1 /. (Thread @ Rule[v1, random]); (*calculate the
             benchmark value of the polynomial*)
        If[OptionValue[deBug],Print["benchmark: ", benchmark]];
        Do[
           tran = tranl[[i]];
           l2 = Permute[random, tran];
           If[OptionValue[deBug],Print["num: ", (p2 /. (Thread @ Rule[v2, l2]))]];
           If[(p2 /. (Thread @ Rule[v2, l2])) == benchmark,
            AppendTo[result,{{tran}, {DeleteCases[Thread @ Rule[Permute[v1, tran
                ], v2], a_ -> a_, {2}]}}]
            ,
            Continue[]
            ];(*find all permutaion list until one satisfy the symmetry condition*)
        ,{i, 1, Length[tranl]}];
        Throw[result];
    ]; 

PolySymCheck[p_, q_, rule_] :=
    ParallelTable[
        If[(((p /. rule[[i]]) - q) // Factor) === 0,
            0
            ,
            1
        ]
        ,
        {i, 1, Length[rule]}
    ]; 

PolySymCheck[p_, rule_] :=
    ParallelTable[
        If[(((p /. rule[[i]]) - p) // Factor) === 0,
            0
            ,
            1
        ]
        ,
        {i, 1, Length[rule]}
    ]; 


Options[PolyOrder]={"order"->"normal"};

PolyOrder[mon1_, mon2_,OptionsPattern[]] :=
    Module[{term, xl, m, n},
        term = mon1 / mon2 // Factor;
        If[term === 1,
            Return[mon1]
        ]; (*the case two are equal*)
        xl = Cases[{Denominator[term]}, Subscript[x, i_] :> i, Infinity
            ] // DeleteDuplicates // Sort;
        m = Sum[Power[2, xl[[i]] - 1], {i, 1, Length[xl]}];
        xl = Cases[{Numerator[term]}, Subscript[x, i_] :> i, Infinity
            ] // DeleteDuplicates // Sort;
        n = Sum[Power[2, xl[[i]] - 1], {i, 1, Length[xl]}];
        If[m > n,
            If[OptionValue["order"]==="normal",Return[mon2],Return[mon1]]
            ,
            If[OptionValue["order"]==="normal",Return[mon1],Return[mon2]]
        ];
        
    ]; 

Options[PolyEqual]={"order"->"normal"};

PolyEqual[mono_, sym_,OptionsPattern[]] :=
    Module[{term = {}, max, tem},
        If[Head[mono] === Plus,
            If[FreeQ[mono, x],
                Return[mono]
                ,
                Message[PolyEqual::err];
                Return[$Failed]
            ]
        ];
        term = Prepend[Table[mono /. sym[[i]], {i, 1, Length[sym]}], 
            mono];
        max = term[[1]];
        Do[max = PolyOrder[max, term[[i]],"order"->OptionValue["order"]], {i, 2, Length[term]}];
        Return[max];
        
    ]; 

Options[PolyFold] = {deBug -> False, "order"->"normal"}; 

PolyFold[poly_, sym_, OptionsPattern[]] :=
    Module[{var, pl},
        var = Variables[poly] // DeleteCases[#, _?(FreeQ[#, x]&)]&;
        pl = MonomialList[poly, var];
        pl = PolyEqual[#, sym,"order"->OptionValue["order"]]& /@ pl;
        If[OptionValue[deBug],
            Echo[pl, "poly list:"]
        ];
        Return[Total[MonomialList[Total[pl], var]]];
        
    ]; 


Unprotect[SProd]; 

Options[SProd] = {ExpandQ -> True}; 

CenterDot[p_, 0, OptionsPattern[]] :=
    0; 

CenterDot[0, q_, OptionsPattern[]] :=
    0; 

SProd[p_, q_, OptionsPattern[]] :=
    If[OptionValue[ExpandQ],
        (CenterDot[p, q] //. {CenterDot[Plus[x_, y__], z_] :> SProd[x,
             z] + SProd[Plus[y], z], CenterDot[z_, Plus[x_, y__]] :> SProd[x, z] 
            + SProd[Plus[y], z], CenterDot[Times[a_Integer, b_], c_] :> a * SProd[
            b, c], CenterDot[c_, Times[a_Integer, b_]] :> a * SProd[b, c], CenterDot[Times[a_/;(NumericQ[a]), b_], c_] :> a * SProd[
            b, c], CenterDot[c_, Times[a_/;(NumericQ[a]), b_]] :> a * SProd[b, c]}) // Factor
            
        ,
        CenterDot[p, q]
    ]; 

SetAttributes[SProd, {Orderless, Listable, Protected}]; 

ExpandSP[exp_]:=Module[{rep},
rep={CenterDot[Plus[x_,y__],z_]:>SProd[x,z]+SProd[Plus[y],z],CenterDot[z_,Plus[x_,y__]]:>SProd[x,z]+SProd[Plus[y],z],CenterDot[Times[a_Integer,b_],c_]:>a*SProd[b,c],CenterDot[c_,Times[a_Integer,b_]]:>a*SProd[b,c]};
(exp//.rep)//Simplify
];

CayleyMengerForm[list_,rep_]:=Table[If[i==j,0,If[i==1||j==1,1,SProd[list[[i-1]]-list[[j-1]],list[[i-1]]-list[[j-1]],ExpandQ->True]]],{i,1,Length[list]+1},{j,1,Length[list]+1}]/.rep;


GramDet[list_] :=
    Module[{m, l},
        l = Length @ list;
        If[l == 0,
            Return[1]
        ];
        m = Table[SProd[list[[i]], list[[j]], ExpandQ -> False], {i, 
            1, l}, {j, 1, l}];
        (Det[m] /. {CenterDot[a_, b_] :> SProd[a, b, ExpandQ -> True]
            }) // Factor
    ]; 

GramDet[list1_, list2_] :=
    Module[{l, m},
        l = Length[list1];
        If[l == 0,
            Return[1]
        ];
        If[l != Length[list2],
            Message[GramDet::err];
            Return[$Failed]
        ];
        m = Table[SProd[list1[[i]], list2[[j]], ExpandQ -> False], {i,
             1, l}, {j, 1, l}];
        (Det[m] /. {CenterDot[a_, b_] :> SProd[a, b, ExpandQ -> True]
            }) // Factor
    ]; 

SyntaxInformation[PD] = {"ArgumentsPattern" -> {_, _, _}}; 

Options[PD] = {RepList -> {}}; 

PD[loopm_, extm_, f_, OptionsPattern[]] :=
    Module[{s, l},
        s = SProd[loopm, loopm] + 2 SProd[loopm, extm];
        l = SProd[extm, extm] /. OptionValue[RepList];
        {s, l + f}
    ]; 

SyntaxInformation[LPD] = {"ArgumentsPattern" -> {_, _, _}}; 

Options[LPD] = {RepList -> {}}; 

LPD[loopm_, extm_, f_, OptionsPattern[]] :=
    Module[{s, l},
        s = SProd[loopm, extm];
        l = f /. OptionValue[RepList];
        {s, l}
    ];

BaikovTrans[list_] :=
    Module[{l, v, m, xl, rep},
        l = Length @ list;
        xl = Table[Subscript[x, i], {i, 1, l}];
        v = Flatten[Variables /@ list[[All, 1]]] // DeleteDuplicates;
            
        m = Table[Coefficient[list[[i, 1]], #]& /@ v, {i, 1, l}];
        If[!SquareMatrixQ[m],
            Message[BaikovTrans::marg, Dimensions[m][[1]], Dimensions[
                m][[2]]];
            Return[0]
        ];
        rep = Thread @ Rule[v, LinearSolve[m, xl - list[[All, 2]]] //
             Factor];
        Return[{m, rep}];
        
    ]; 

Options[BaikovRep] = {Method -> "Standard", RepList -> {}, Abst -> True,
     D0 -> 4}; 

BaikovRep[dlist_, llist_, elist_, OptionsPattern[]] :=
    Module[{L, E, N, r, a, rep, K, Gt, c, l, P, int},
        L = Length @ llist;
        E = Length @ elist;
        N = L (L + 1) / 2 + L * E;
        r = BaikovTrans[dlist /. OptionValue[RepList]];
        If[r == 0,
            Return[0]
        ];
        a = Det[r[[1]]];
        rep = r[[2]];
        If[OptionValue[Method] === "Standard",
            K = L + E;
            Gt =
                If[OptionValue[Abst],
                    G[elist, elist]
                    ,
                    GramDet[elist]
                ];
            c =
                If[OptionValue[Abst],
                    Pi^((L - N) / 2) / (a * Product[Gamma[(D - K + i)
                         / 2], {i, 1, L}])
                    ,
                    Pi^((L - N) / 2) / (a * Product[Gamma[(D - K + i)
                         / 2], {i, 1, L}] * (Gt)^((D - E - 1) / 2))
                ];
            l = Join[llist, elist];
            P =
                If[OptionValue[Abst],
                    G[l, l]
                    ,
                    ((GramDet[l]) /. rep) // Factor
                ];
            
        ];
        int = P /. OptionValue[RepList] // Together;
        If[OptionValue[Abst],
            Return[{{{Gt, -(D - E - 1) / 2}, {P, (D - K - 1) / 2}}, c
                } /. {D -> OptionValue[D0] - 2 \[Epsilon]}]
        ];
        c = c * Power[Denominator[int], -(D - K - 1) / 2] /. OptionValue[
            RepList] /. {D -> OptionValue[D0] - 2 \[Epsilon]} // Simplify;
        int = Power[Numerator[int], (D - K - 1) / 2] /. {D -> OptionValue[
            D0] - 2 \[Epsilon]} // Factor;
        Return[{c, int}]
    ]; 


GramMat[l1_,l2_,rep_]:=Table[SProd[l1[[i]],l2[[j]]],{i,1,Length[l1]},{j,1,Length[l2]}]//.rep//Factor;
Gram2Poly[exp_,krep_]:=exp/.{G[x_,y_]:>(GramMat[x,y,krep]//Det//Factor)};


Options[ToCanonicalForm] = {Form -> 0, Cut -> False, deBug -> False};

ToCanonicalForm[u_, FI_, OptionsPattern[]] :=
    Module[{list, l, r},
        Switch[OptionValue[Form],
            0,(*the form j[family,1,0,1,...]*)list = Drop[FI, 1]
            ,
            1,(*the form G[1,{1,0,1,...}]*)list = FI[[2]]
            ,
            2, (*the form G[1,0,1,...]*)list = FI
        ];
        l = Length[list];
        r = u;
        Do[
            If[list[[i]] > 0,
                r = D[r, {Subscript[x, i], list[[i]] - 1}] / Factorial[
                    list[[i]] - 1]
            ]
            ,
            {i, l}
        ];
        r = r / u // Factor;
        Do[
            If[list[[i]] > 0,
                If[OptionValue[Cut],
                    r = r /. {Subscript[x, i] -> 0}
                    ,
                    r = r / Subscript[x, i]
                ]
            ]
            ,
            {i, l}
        ];
        Do[
            If[list[[i]] < 0,
                r = r * Power[Subscript[x, i], list[[i]]]
            ]
            ,
            {i, l}
        ];
        Return[r];
    ];


Options[Trans2FI] = {ReVerse -> Off};

Trans2FI[num_, list_, plist_, f_, OptionsPattern[]] :=
    Module[{l, cc, p, g, s, tt},
        If[OptionValue[ReVerse] === On,
            Goto[rev]
        ];
        l = List @@ (Expand[num]);
        cc = l /. Thread @ Rule[list, 1];
        p = plist - (Exponent[l, #]& /@ list) // Transpose;
        g = (Map[f[#]&, p, {1}]) /. {f[x_] :> f[family, Sequence @@ x
            ]};
        s = tt * (cc . g) // Collect[#, _f, Together]&;
        g =
            s //
            Cases[#, _f, Infinity]& //
            DeleteDuplicates;
        cc = Coefficient[s, #]& /@ g;
        Return[{cc /. {tt -> 1}, g}];
        Label[rev];
        l = Cases[num, _f, {0, Infinity}] // DeleteDuplicates;
        Print[l];
        cc =
            Table[
                p = l[[i]] /. {f[family, z__] :> {z}};
                Times @@ Table[Power[list[[j]], -p[[j]]], {j, 1, Length[
                    p]}]
                ,
                {i, 1, Length[l]}
            ];
        tt = Thread @ Rule[l, cc];
        Return[num /. tt];
    ];


LocatePos[matlist_, v_] :=
    Module[{pos},
        Return[
            Table[
                pos = Take[#, 2]& /@ Position[UpperTriangularize[matlist[[
                    i]]], v];
                If[pos === {},
                    Null
                    ,
                    Prepend[pos, i]
                ]
                ,
                {i, 1, Length[matlist]}
            ] // DeleteCases[#, Null]&
        ];
        
    ]; 

ReArrangeGram[g_, pos_, ratio_, flag_] :=
    Module[{row, col, pl},
        If[Length[ratio] == 1,
            Return[g]
        ];
        row = g[[1]];
        col = g[[2]];
        If[flag == -1,
            pl = pos[[All, 2]]
            ,
            If[flag == 1,
                pl = pos[[All, 1]]
            ]
            ,
            Message[ReArrangeGram::err, flag];
            Return[$Failed]
        ];
        If[flag == -1,
            Table[row[[pl[[i]]]] = row[[pl[[i]]]] - ratio[[i]] * row[[
                pl[[-1]]]], {i, 1, Length[pl] - 1}];
            Table[col[[pl[[i]]]] = col[[pl[[i]]]] - ratio[[i]] * col[[
                pl[[-1]]]], {i, 1, Length[pl] - 1}]
            ,
            Table[row[[pl[[i]]]] = row[[pl[[i]]]] - ratio[[i]] * row[[
                pl[[1]]]], {i, 2, Length[pl]}];
            Table[col[[pl[[i]]]] = col[[pl[[i]]]] - ratio[[i]] * col[[
                pl[[1]]]], {i, 2, Length[pl]}]
        ];
        (*Table[row[[pl[[i]]]]=row[[pl[[i]]]]-ratio[[i]]*row[[pl[[-1]
            ]]],{i,1,Length[pl]-1}];Table[col[[pl[[i]]]]=col[[pl[[i]]]]-ratio[[i]
            ]*col[[pl[[-1]]]],{i,1,Length[pl]-1}];*)
        Return[Head[g][row, col]];
        
    ]; 

ReduceMat[{g_, power_}, pos_, coef_] :=
    Module[{row, col, c, r1, c1, result = {}},
        row = g[[1]];
        col = g[[2]];
        c = Power[Abs[coef], -1] * Power[4, power + 1 / 2];(*last term
             comes from b^2-4a*c,first term is a rescale of integration variable*)
            
        c = c * (-((Gamma[power] Gamma[power + 1]) / (2 (2 power + 1)
             Gamma[2 power]))) // Simplify;
        r1 = Delete[row, {{pos[[1]]}, {pos[[2]]}}];
        c1 = Delete[col, {{pos[[1]]}, {pos[[2]]}}];
        If[r1 =!= {} && c1 =!= {},
            AppendTo[result, {-Head[g][r1, c1], -power - 1}]
            ,
            c = c * Power[-1, -power - 1]
        ];
        r1 = Delete[row, pos[[1]]];
        c1 = Delete[col, pos[[1]]];
        AppendTo[result, {Head[g][r1, c1], power + 1 / 2}];
        r1 = Delete[row, pos[[2]]];
        c1 = Delete[col, pos[[2]]];
        AppendTo[result, {Head[g][r1, c1], power + 1 / 2}];
        Return[AppendTo[result, c]];
        
    ]; 

Options[SimplifyGram] = {deBug -> False}; 

SimplifyGram[gl_, rep_, OptionsPattern[]] :=
    Module[{g, p, vl, nrep, ng, tl, temo = {}, tem = {}, temp = {}, temng
         = {}, result = {}, k = 1, flag = 1, sg = 1, t1, t2, c = 1},
        g = gl[[All, 1]] /. {G[pl1_, pl2_] :> GramMat[pl1, pl2, rep]}
            ;(*all gram determinant*)
        p = gl[[All, 2]];(*the power of them*)
        vl = Variables[g];
        nrep = Table[Thread @ Rule[vl, Table[RandomPrime[{10^3, 10^4}
            ], {i, 1, Length[vl]}]], {i, 1, 3}];(*generate three sets of random primes
            *)
        ng = Det /@ (g /. nrep[[1]]);
        temo = gl[[All, 1]];
        tem = g;
        temp = p;
        temng = ng;
        While[
            tem =!= {}
            ,
            flag = 0;
            Do[
                If[temng[[i]] == temng[[1]],
                    flag = i;
                    sg = 1;
                    Break[]
                ];
                If[temng[[i]] + temng[[1]] == 0,
                    flag = i;
                    sg = -1;
                    Break[]
                ]
                ,
                {i, 2, Length[temng]}
            ];
            If[OptionValue[deBug],
                Print["flag: ", flag]
            ];
            If[flag != 0,
                t1 = (Det[tem[[1]] /. nrep[[2]]] - sg * Det[tem[[flag
                    ]] /. nrep[[2]]]);
                t2 = (Det[tem[[1]] /. nrep[[3]]] - sg * Det[tem[[flag
                    ]] /. nrep[[3]]]);
                If[t1 != 0 || t2 != 0,
                    flag = 0
                ]
            ];(*if it doesn't pass another two numerical check, then 
                reset flag to 0*)
            If[OptionValue[deBug],
                Print["flag: ", flag]
            ];
            If[flag == 0,
                AppendTo[result, {temo[[1]], temp[[1]]}];
                temo = Drop[temo, 1];
                tem = Drop[tem, 1];
                temp = Drop[temp, 1];
                temng = Drop[temng, 1]
                ,
                If[(temp[[1]] + temp[[flag]] // Factor) =!= 0,
                    AppendTo[result, {temo[[1]], temp[[1]] + temp[[flag
                        ]] // Factor}]
                    ,
                    If[OptionValue[deBug],
                        Print["temo:", temo[[1]]];
                        Print["temp: ", temp[[1]]];
                        Print["temp[[flag]]", temp[[flag]]]
                    ]
                ];
                c = c * Power[sg, temp[[flag]]];
                temo = Delete[temo, {{1}, {flag}}];
                tem = Delete[tem, {{1}, {flag}}];
                temp = Delete[temp, {{1}, {flag}}];
                temng = Delete[temng, {{1}, {flag}}]
            ];
            
        ];
        Return[AppendTo[result, c]];
        
    ]; 

Options[ZeroSectorMatQ] = {deBug -> False}; 

ZeroSectorMatQ[gl_, rep_, OptionsPattern[]] :=
    Catch @
        Module[{ul, vl, numrep,(*pow,*)xl, set, pos, nrep, gr = {}, xg,
             r},
            ul = gl[[All, 1]] /. {G[pl1_, pl2_] :> GramMat[pl1, pl2, 
                rep]};(*all Baikov polynomial*)
            vl = Variables[ul] // DeleteCases[#, _?(!FreeQ[#, x]&)]&;
                
            numrep = Thread @ Rule[vl, Table[RandomPrime[{10^3, 10^4}
                ], {i, 1, Length[vl]}]];
            ul = Det /@ (ul /. numrep) // Factor;
            If[OptionValue[deBug],
                Print["ul: ", ul]
            ];
            ul = DeleteCases[ul, _?(FreeQ[#, x]&)];
            If[ul === {},
                Message[ZeroSectorMatQ::warning];
                Throw[True]
            ];
            xl = Table[Cases[ul[[i]], Subscript[x, _], {0, Infinity}]
                 // DeleteDuplicates, {i, 1, Length[ul]}];
            If[Length[xl] > 1,
                set = Subsets[Table[i, {i, 1, Length[xl]}], {2}];
                gr =
                    Table[
                        If[IntersectingQ[xl[[set[[i, 1]]]], xl[[set[[
                            i, 2]]]]],
                            UndirectedEdge[set[[i, 1]], set[[i, 2]]]
                            ,
                            Null
                        ]
                        ,
                        {i, 1, Length[set]}
                    ] // DeleteCases[#, Null]&;
                gr = ConnectedComponents[Graph[Table[i, {i, 1, Length[
                    xl]}], gr]]
                ,
                gr = {{1}}
            ];(*we get the connected component to find whether this group
                 of variables is homogeneous under rescaling*)
            If[OptionValue[deBug],
                Print["gr: ", gr]
            ];
            Do[
                xg = Union[Table[xl[[gr[[i, j]]]], {j, 1, Length[gr[[
                    i]]]}] // Flatten];
                If[xg === {},
                    Continue[]
                ];(*this part doesn't contain x variables, this probably
                     impossible if we have removed these polynomials*)
                nrep = Thread @ Rule[xg, R * xg];(*rescale this set of
                     variables*)
                r = (ul /. nrep) / ul // Factor;
                (*If[FreeQ[r,x],If[!FreeQ[Product[Power[r[[i]],pow[[i
                    ]]],{i,1,Length[ul]}]//Simplify,R],Throw[True]]];*)
                If[OptionValue[deBug],
                    Print["r: ", r]
                ];
                If[FreeQ[r, x],
                    Throw[True]
                ];
                (*If it is a homogeneous polynomial then we know it is
                     a zero sector*)
                ,
                {i, 1, Length[gr]}
            ];
            Throw[False];
            
        ]; 

Options[ReduceRep] = {deBug -> False}; 

ReduceRep[pl_, v_, rep_, OptionsPattern[]] :=
    Catch @
        Module[{gl, pos, i, int, int1, mat, elm, coef, ratio, flag = 
            0, tem, con, result},
            gl = pl[[All, 1]];
            pos = LocatePos[gl /. {G[pl1_, pl2_] :> GramMat[pl1, pl2,
                 rep]}, v];(*find position of v in the Gram det list*)
            If[OptionValue[deBug],
                Print["pos: ", pos]
            ];
            If[pos === {},
                Message[ReduceRep::devoid, v];
                Export["./debug.m",
                     {pl, v, rep}];
                Throw[$Failed]
            ];(*if v has not been found, then there must be something
                 wrong*)
            If[Length[pos] > 1,
                Throw[False]
            ];(*the reducible variable cannot appear in multiple Gram
                 determinants*)
            i = pos[[1, 1]];(*pick up the number of matrix which v is
                 in*)
            mat = gl[[i]] /. {G[pl1_, pl2_] :> GramMat[pl1, pl2, rep]
                };(*Extract the matrix*)
            pos = Drop[pos[[1]], 1] // Sort;
            If[OptionValue[deBug],
                Print["pos: ", pos]
            ];
            If[Length[pos] == 1 && pos[[1, 1]] == pos[[1, 2]],
                Throw[False]
            ];(*if v appears only in the diagonal position, it is linear
                 in polynomial*)
            int = Intersection[Sequence @@ Partition[pos[[All, 1]], 1
                ]];
            int1 = Intersection[Sequence @@ pos];(*in case v in the i
                -th column and i-th row*)
            If[int === {},
                If[Intersection[Sequence @@ Partition[pos[[All, 2]], 
                    1]] === {},
                    If[int1 === {},
                        Throw[False]
                        ,
                        flag = 0;
                        (*Print[v," dist fractured"]*)]
                    ,(*v is in one column and row but is fractured*)
                    flag = 1;
                    (*Print[v," dist in one column"]*)]
                ,(*v in one column*)
                flag = -1;
                (*Print[v," dist in one row"]*) ];(*v in one row*)
            (*if v appears in multiple rows and columns, it is not reducible
                . flag=-1 means in the same row and 1 means the same column*)
            elm = ((Extract[mat, pos] * R) /. {v -> v / R} // Factor)
                 /. {R -> 0};
            If[OptionValue[deBug],
                Print["elm: ", elm]
            ];
            If[flag != 0,(*two simpler case*)
                coef = elm[[flag]] / v;
                ratio = elm / elm[[flag]] // Factor
                ,
                (*the more complicated case*)
                coef = elm[[-1]] / v;
                ratio =
                    Table[
                        If[pos[[i, 1]] == pos[[i, 2]],
                            1 / 2 elm[[i]] / elm[[-1]]
                            ,
                            elm[[i]] / elm[[-1]]
                        ]
                        ,
                        {i, 1, Length[elm]}
                    ]
            ];
            (*if v in the same row, we eliminate others with the last
                 item. if in the same column, we eliminate with the first item*)
            If[flag == -1,
                If[pos[[1, 1]] == pos[[1, 2]],
                    ratio[[1]] = 1 / 2 ratio[[1]]
                ]
            ];
            If[flag == 1,
                If[pos[[-1, 1]] == pos[[-1, 2]],
                    ratio[[-1]] = 1 / 2 ratio[[-1]]
                ]
            ];
            If[flag == 0,
                pos = pos /. {{a_, int1[[1]]} :> {int1[[1]], a}};
                flag = -1
            ];
            (*if there is diag term, it should be treated separately,
                 note that we have folded the matrix into its uppertriangle form*)
            tem = ReArrangeGram[gl[[i]], pos, ratio, flag];(*rearrange
                 the matrix so that v only appear in one place*)
            If[OptionValue[deBug],
                Print["tem (matrix rearranged): ", tem]
            ];
            tem = ReduceMat[{tem, pl[[i, 2]]}, pos[[flag]], coef];
            con = tem[[-1]];(*extract the constant term after the recursion
                *)
            tem = Drop[tem, -1];
            tem = Join[Delete[pl, i], tem];(*Delete the origin matrix
                *)
            If[OptionValue[deBug],
                Print["tem: ", tem]
            ];
            result = SimplifyGram[tem, rep(*,deBug->True*)];
            con = con * result[[-1]];
            result = Drop[result, -1];
            If[ZeroSectorMatQ[result, rep],
                con = 0
            ];(*if this is a zero sector we drop its constant term since
                 the result will be 0*)
            Throw[{result, con}];
            
        ]; 

Options[AllSectorBaikovMat] = {Exc -> {}, deBug -> False}; 

AllSectorBaikovMat[list_, rep_, OptionsPattern[]] :=
    Catch @
        Module[{pl, glist, c, var, path = {}, flag = 1, k = 1, result
             = {}, pv, tem, intv, temr, tag = 0},
            pl = list[[1]];(*the Gram determinant and their power*)
            c = list[[2]];(*the constant*)
            glist = pl[[All, 1]];
            var = Complement[Variables[glist /. {G[l1_, l2_] :> GramMat[
                l1, l2, rep]}] // DeleteCases[#, _?(FreeQ[#, x]&)]&, OptionValue[Exc]
                ];(*variables that will be integrated out*)
            If[OptionValue[deBug],
                Print["var: ", var];
                (*Print["x: ", Context[x]];*)
                (*Print["var: ", Variables[glist /. {G[l1_, l2_] :> GramMat[
                l1, l2, rep]}]];
                Print["var: ", Variables[glist /. {G[l1_, l2_] :> GramMat[
                l1, l2, rep]}] // DeleteCases[#, _?(FreeQ[#, x]&)]&];*)
            ];
            result = {{{{}, list}}};(*store according to the layer involved
                *)
            While[
                flag == 1 && k < 100
                ,
                path = (ReverseSort /@ Subsets[var, {k}]) // ReverseSort
                    ;
                temr = {};
                Do[
                    pv = path[[i]];(*pick up a path*)
                    tag = 0;
                    Do[
                        intv = Complement[pv, result[[k]][[j, 1]]];
                        If[Length[intv] != 1,
                            Continue[]
                            ,
                            intv = intv[[1]]
                        ];(*find which variable to be integrated out 
                            based on the last layer result*)
                        (*If[OptionValue[deBug],Print["intv: ",intv]]
                            ;*)
                        If[result[[k]][[j, 2, 2]] === 0,
                            Break[]
                        ]; (*if it is a zero sector we won't reduce it
                             further*)
                        tem = ReduceRep[result[[k]][[j, 2, 1]], intv,
                             rep];
                                 (*If[tem===$Failed,Message[AllSectorBaikov
                            ::err,n,pv];Throw[$Failed]];
                        If[tem===Null,Break[]];(*If this is a zero sector
    , we break this loop*)*)
                        If[tem === $Failed,
                            Message[AllSectorBaikovMat::err, pv, intv
                                ];
                            Export["./debug_int.m",
                                 result];
                            Throw[{result[[k]][[j, 2, 1]], intv}]
                        ];
                        If[tem === False,
                            Continue[]
                        ];
                        tag = 1;
                        c = result[[k]][[j, 2, 2]] * tem[[-1]] // FullSimplify
                            ;
                        AppendTo[temr, {pv, {tem[[1]], c}}];
                        Break[]
                        ,
                        {j, 1, Length[result[[k]]]}
                    ];
                    If[tag != 0 && OptionValue[deBug],
                        Print["pv (effective path): ", pv]
                    ];
                    
                    ,
                    {i, 1, Length[path]}
                ];
                AppendTo[result, temr];
                k = k + 1;
                If[k >= Length[var],
                    flag = 0
                ]
            ];
            Throw[result];
            
        ]; 

Options[AllSectorBaikovMatC] = {deBug -> False}; 

AllSectorBaikovMatC[resultmat_,var_,rep_, OptionsPattern[]] :=
    Catch @
        Module[{ result, c, path = {}, flag = 1, k, pv, tem, intv, temr, tag = 0},
            result=DeleteCases[resultmat,{}];
            k=(DeleteCases[result,{}][[-1]][[1,1]]//Length)+1;
            While[
                flag == 1 && k < 100
                ,
                path = (ReverseSort /@ Subsets[var, {k}]) // ReverseSort
                    ;
                temr = {};
                Do[
                    pv = path[[i]];(*pick up a path*)
                    tag = 0;
                    Do[
                        intv = Complement[pv, result[[k]][[j, 1]]];
                        If[Length[intv] != 1,
                            Continue[]
                            ,
                            intv = intv[[1]]
                        ];(*find which variable to be integrated out 
                            based on the last layer result*)
                        (*If[OptionValue[deBug],Print["intv: ",intv]]
                            ;*)
                        If[result[[k]][[j, 2, 2]] === 0,
                            Break[]
                        ]; (*if it is a zero sector we won't reduce it
                             further*)
                        tem = ReduceRep[result[[k]][[j, 2, 1]], intv,
                             rep];
                                 (*If[tem===$Failed,Message[AllSectorBaikov
                            ::err,n,pv];Throw[$Failed]];
                        If[tem===Null,Break[]];(*If this is a zero sector
    , we break this loop*)*)
                        If[tem === $Failed,
                            Message[AllSectorBaikovMat::err, pv, intv
                                ];
                            Export["./debug_int.m",
                                 result];
                            Throw[{result[[k]][[j, 2, 1]], intv}]
                        ];
                        If[tem === False,
                            Continue[]
                        ];
                        tag = 1;
                        c = result[[k]][[j, 2, 2]] * tem[[-1]] // FullSimplify
                            ;
                        AppendTo[temr, {pv, {tem[[1]], c}}];
                        Break[]
                        ,
                        {j, 1, Length[result[[k]]]}
                    ];
                    If[tag != 0 && OptionValue[deBug],
                        Print["pv (effective path): ", pv]
                    ];
                    
                    ,
                    {i, 1, Length[path]}
                ];
                AppendTo[result, temr];
                k = k + 1;
                If[k >= Length[var],
                    flag = 0
                ]
            ];
            Throw[result];
            
         ]; 

ExtractLoopOrder[rep_]:=Module[{list,listk,pos,ext,loop},
	list=DeleteCases[rep,_?(FreeQ[#,x]&)];
	listk=Complement[rep,list];
	list=Keys[list]/.{SProd[a_,b_]:>{a,b},CenterDot[a_,b_]:>{a,b}}//Flatten//DeleteDuplicates;
	listk=Keys[listk]/.{SProd[a_,b_]:>{a,b},CenterDot[a_,b_]:>{a,b}}//Flatten//DeleteDuplicates;
	Return[Length[Complement[list,listk]]];
];


Options[GetBaikovMatRep] = {"ExcVar" -> {}, "looporder"-> 2, "ForceAdd" -> 0, deBug-> False}; 

GetBaikovMatRep[result_, var_, n_, OptionsPattern[]] :=
    Module[{intv, l, k = 1, tem, pos, temp, res = {},flag},
        intv = Complement[Table[i, {i, 1, n}], Join[var, OptionValue[
            "ExcVar"]]] // ReverseSort;
        intv = Subscript[x, #]& /@ intv;
        l = Length[intv];
        If[l >= Length[result],
            Message[GetBaikovMatRep::err]
        ];
        k = 1;
        flag = OptionValue["ForceAdd"];
        (*if the representation is not suitable, we need to force the representation to add one more isp*)
        While[
            k <= l
            ,
            tem = result[[l + 2 - k]];
            If[tem === {},
                If[OptionValue[deBug],Print["additional isp needed! ", k]];
                k = k + 1;
                Continue[]
            ];
            temp = tem[[All, 1]];
            Do[
                If[Length[Intersection[temp[[j]], intv]] > l - k,
                    If[res==={},(*this sector aims to find all representations that are independent*)
                        AppendTo[res, tem[[j]]],
                        If[flag!=0,
                            AppendTo[res, tem[[j]]],
	                        Do[
	                            If[OptionValue[deBug],Print["res: ",res[[All,1]]]];
	                            If[ContainsAll[res[[c,1]],temp[[j]]],
	                                Break[],
	                                AppendTo[res, tem[[j]]]
	                            ]
	                        ,{c,1,Length[res]}]
	                    ]
                    ];
                    If[tem[[j, 2, 2]] === 0,
                        Message[GetBaikovMatRep::zero];
                        Break[]
                    ]
                ]
                ,
                {j, 1, Length[temp]}
            ];
            If[res =!= {},
                If[flag > 0, 
                    flag = flag - 1, 
                    If[Length[res]<OptionValue["looporder"],k=k+1;Continue[],Break[]]
                ]
            ];
            If[OptionValue[deBug],Print["additional isp needed! ", k]];
            k = k + 1; 
        ];
        If[res == {},
            If[OptionValue[deBug],Print["Only the original representation can be its representation!"]];
            Return[result[[1]]]
        ];
        Return[res];
        
    ]; 


AllSubSector[zeroset_]:=Module[{l,nl,p,tp,result={}},
l=Length[zeroset];
Do[
nl=IntegerDigits[zeroset[[i]],2];
p=Position[nl,1]//Flatten;
tp=Tuples[{1,0},Length[p]];
result=Join[result,Table[ReplacePart[nl,Thread@Rule[p,tp[[i]]]],{i,1,Length[tp]}]];
,{i,1,l}];
result=result//DeleteDuplicates;
result=FromDigits[#,2]&/@result;
Return[result//DeleteDuplicates//Sort];
];

Options[GetMatZeroSector]={deBug->False};
GetMatZeroSector[list_,n_,isp_,OptionsPattern[]]:=Module[{pos,pl,nl,sl},
pos=Position[list,{_,0}];
pos=Drop[#,-1]&/@pos;
pl=Table[list[[Sequence@@(pos[[i]])]][[1]],{i,1,Length[pos]}]/.{Subscript[x,i_]:>i};
If[OptionValue[deBug],Print["pl: ",Short[pl,5]]];
nl=Table[i,{i,1,n}];
pl=Complement[nl,#]&/@pl;(*get those variables which haven't been integrated out*)
pl=Complement[#,isp]&/@pl;(*remove isps*)
If[OptionValue[deBug],Print["pl: ",Short[pl,5]]];
sl=Table[Sum[Power[2,pl[[i,j]]-1],{j,1,Length[pl[[i]]]}],{i,1,Length[pl]}]//DeleteDuplicates;
If[OptionValue[deBug],Print["sl: ",Short[sl,5]]];
Return[AllSubSector[sl]];
];


End[];


EndPackage[];

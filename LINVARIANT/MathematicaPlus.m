BeginPackage["LINVARIANT`MathematicaPlus`", {"LINVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
NOrderResponse              ::usage "NOrderResponse[eqn, var, n]"
GetVariationVar             ::usage "GetVariationVar[var]"
Complex2Exp                 ::usage "Complex2Exp[exp]"
Exp2Complex                 ::usage "Exp2Complex[exp]"
GetSubscriptInfo            ::usage "GetSubscriptInfo[atom]"
ReadNonEmptyLine            ::usage "ReadNonEmptyLine[stream]"
PermuteThrough              ::usage "PermuteThrough[tab]"
Var2Var                     ::usage "Var2Var[var1, var2, dir]"
ParseFortranNumber          ::usage "ParseFortranNumber[stream]"
MakeMatrixBlock             ::usage "MakeMatrixBlock[mat, dim]"
SimpsonIntegrate            ::usage "SimpsonIntegrate[f, x]"
Irrationalize               ::usage "Irrationalize[x, prec]"
RectangleIntPath            ::usage "ComplexIntegratePath[spt, ept, npt, ratio]"
partitionBy                 ::usage "partitionBy[l, p]"
CloneReshape2D              ::usage "CloneReshape2D[array0, array1]"
Model2P1                    ::usage "Model2P1[model]"
FixDoubleCounting           ::usage "FixDoubleCounting[expr]"
OnSiteExprQ                 ::usage "OnSiteExprQ[expr]"
VarBare                     ::usage "VarBare[var]"
VarHead                     ::usage "VarHead[var]"
VarInd                      ::usage "VarInd[var]"
VarSite                     ::usage "VarSite[var]"
VarGrid                     ::usage "VarGrid[var]"
VarGridShift                ::usage "VarGridShift[var, s]"
VarOnGrid                   ::usage "VarOnGrid[var, site]"
ExprOnGrid                  ::usage "ExprOnGrid[expr, site]"
Expr2Gamma                  ::usage "Expr2Gamma[model, grid]"
Model2Gamma                 ::usage "Model2Gamma[model]"
SortVarSub                  ::usage "SortVarSub[varsub, vars, uvars]"
ShiftVarSub                 ::usage "ShiftVarSub[varsub, shift]"
SwapVarSub                  ::usage "SwapVarSub[varsub]"
VarSubStar                  ::usage "VarSubStar[t5sub, t6sub, SymMat]"
MyTime                      ::usage "MyTime[inp, NewProjQ]"
ColorPots                   ::usage "ColorPots[v, n, tol]"
SimplifyTensorCommonFactor  ::usage "SimplifyTensorCommonFactor[mat]"
SimplifyTensor              ::usage "SimpifyTensor[T0, TT]"
DecomposeExpr               ::usage "DecomposeExpr[x]"
NumberCommonDivisor         ::usage "NumberCommonDivisor[NumList, prec_:10^-12]"
GetConstantFactor           ::usage "GetConstantFactor[expr]"
SimplifyCommonFactor        ::usage "SimplifyCommonFactor[expr, prec]"
SimplifyIrrationalNumber    ::usage "SimplifyIrrationalNumber[x, prec]"
Sub2List                    ::usage "Sub2List[sub]"
StrainInvQ                  ::usage = "StrainInvQ[x, e]"
XInvQ                       ::usage = "XInvQ[x, e]"
TSVars                      ::usage = "TSVars[ts]"
TSGrid                      ::usage = "TSGrid[ts]"
VarSubGrid                  ::usage = "VarSubGrid[vs]"
T5Vars                      ::usage = "T5Vars[t]"
T6Vars                      ::usage = "T6Vars[t]"
T5Data                      ::usage = "T5Data[t]"
T6Data                      ::usage = "T6Data[t]"
T56Data                     ::usage = "T56Data[t]"
OpenExtract                 ::usage = "OpenExtract[f, str, i, nline]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
CloneReshape2D[array0_, array1_] := Module[{clone, map, out},
  clone = Length[#] &/@ array0;
  map = {(Accumulate[clone]) - clone + 1, Accumulate[clone]}\[Transpose];
  out = Take[array1, #] &/@ map;
  Return[out]
]

Complex2Exp[exp_] := Module[{},
  Expand[exp /. z_?NumericQ :> Abs[z] Exp[I Arg[z]]]
]

Exp2Complex[exp_] := Module[{},
  Expand[ReIm[exp].{1, I}]
]

GetVariationVar[var_] := Module[{\[Delta]var, varnew},
  Which[ListQ[var], GetVariationVar[#] &/@ var,
  True,
  \[Delta]var = ToString[#] & /@ Level[var, 1];
  varnew = "\!\(\*SubscriptBox[\(" <> "\[Delta]" <> \[Delta]var[[1]] <> "\),\(" <> StringJoin[Riffle[\[Delta]var[[2 ;;]], ","]] <> "\)]\)";
  Return[varnew]]
]

NOrderResponse[eqn_, var_, n_] := Module[{\[Epsilon], exp, varnew},
  varnew = ToExpression[GetVariationVar[#]] &/@ var;
  exp = Normal[Series[eqn /. Thread[var -> (var + \[Epsilon] varnew)], {\[Epsilon], 0, n}]];
  Expand[(exp /. {\[Epsilon] -> 1}) - (exp /. {\[Epsilon] -> 0})] /.{0.->0}
]

GetSubscriptInfo[atom_] := Module[{},
  If[MatchQ[atom, _List], GetSubscriptInfo[#] & /@ atom, Which[AtomQ[atom], atom, MatchQ[atom, _Subscript], Level[atom, 1], MatchQ[atom, _Power], If[MatchQ[Level[atom, 1][[2]], _Integer], GetSubscriptInfo[#] & /@ ConstantArray[Level[atom, 1][[1]], Level[atom, 1][[2]]], ##&[]], True, GetSubscriptInfo[#] & /@ Level[atom, 1]]]
]

ReadNonEmptyLine[stream_] := Module[{templine},
  templine = ReadLine[stream];
  If[StringSplit[templine] === {}, ReadNonEmptyLine[stream], templine]
]

PermuteThrough[tab_] := Module[{depth, len, perm},
  {depth, len} = Dimensions[tab];
  perm = {};
  Do[perm = Permute[Prepend[perm, Range[len]]\[Transpose], FindPermutation[Range[len], tab[[i]]]]\[Transpose], {i, depth, 1, -1}];
  Return[perm]
]

Var2Var[var1_, var2_, dir_?IntegerQ, OptionsPattern[{"site"->False}]] := Module[{},
  If[OptionValue["site"],
     Which[dir == 1, 
           Thread[(Subscript@@(Join[GetSubscriptInfo[#],ToExpression["{ix_,iy_,iz_}"]])&/@var2)
                ->(Subscript@@(Join[GetSubscriptInfo[#],ToExpression["{ix,iy,iz}"]])&/@var1)], 
           dir == 2, 
           Thread[(Subscript@@(Join[GetSubscriptInfo[#],ToExpression["{ix_,iy_,iz_}"]])&/@var1)
                ->(Subscript@@(Join[GetSubscriptInfo[#],ToExpression["{ix,iy,iz}"]])&/@var2)]],
     Which[dir == 1, Thread[var2 -> var1], dir == 2, Thread[var1 -> var2]]
     ]
]

ParseFortranNumber[stream_] := Module[{out},
  If[ListQ[stream],
     ParseFortranNumber[#] &/@ stream,
     out = If[MatchQ[Head[stream], InputStream],
              ToExpression[StringSplit[StringReplace[ReadLine[stream], "e" | "E" -> "*^"]]],
              ToExpression[StringSplit[StringReplace[stream, "e" | "E" -> "*^"]]]];
     Return[If[Length[out] == 1, First@out, out]]]
]

MakeMatrixBlock[mat_, dim_] := Module[{MBlocked, col, row},
  {col, row} = dim;
  MBlocked = (Partition[#, col] & /@ ((Partition[#, row] & /@ mat)\[Transpose]))\[Transpose];
  Return[MBlocked]
]

SimpsonIntegrate[f_, x_] := Module[{npt, a, b, m},
  npt = Length[x];
  ParallelSum[a = x[[i]]; b = x[[i + 1]]; m = (a + b)/2;
   1/6 (b - a) (f[a] + 4 f[m] + f[b]), {i, npt - 1}]
]

RectangleIntPath[dir_, spt_, ept_, npt_, ratio_] := Module[{de, je, path},
  de = (ept - spt)/npt;
  je = ratio (ept - spt);
  path = {{spt + I #} & /@ Range[de/2, je, je/20], spt + # + I je & /@ Range[de/2, ept - spt, de], {ept + I #} & /@ Range[je - je/20/2, 0, -je/20]};
  Export[dir<>"/epath.dat", Join[{{Length[path]}}, {NumberForm[Re[#], {25, 15}], NumberForm[Im[#], {25, 15}]} & /@ path] // MatrixForm];
  Return[Flatten[path]]
]

partitionBy[l_, p_] := 
  MapThread[l[[# ;; #2]] &, {{0}~Join~Most@# + 1, #} &@Accumulate@p]

Model2P1[model_, OptionsPattern[{"round" -> 1.0*10^-16}]] := Module[{expr, terms, t, out, p},
  expr = Expand@Total[Times @@@ model];
  terms = If[Head[expr] === Plus, Level[Expand@expr, {1}], {expr}]; 
  out = Table[p = t /. {Subscript[__] -> 1}; 
              {Round[p, OptionValue["round"]], Rationalize[t/p]}, {t, terms}];
  Return[out]
]

FixDoubleCounting[expr_] := Module[{inv, n, v, i, ix, iy, iz, nsite, out, invlist},
  invlist = If[ListQ[expr], expr, {expr}];
  out = Table[inv = If[Head[invlist[[n]]] === Plus, First@Level[invlist[[n]], {1}], invlist[[n]]];
              nsite = Length@DeleteDuplicates[Cases[inv, Subscript[v_, i_, ix_, iy_, iz_] -> {ix, iy, iz}, Infinity]];
              Expand[invlist[[n]]/If[nsite==0, 1, nsite]], {n, Length@invlist}];
  Return[out]
]

Sub2List[sub_] := Transpose[{#1, #2} & @@@ sub]

VarGrid[var_] := Module[{v, i, ix, iy, iz, site, varinfo},
  If[ListQ[var], 
     VarGrid[#] &/@ var,
     varinfo = Level[var, Infinity];
     site = If[Length[varinfo] > 3, varinfo[[-3;;-1]], {0,0,0}]]
]

VarGridShift[var_, s_, grid_, CenterQ_:False] := Module[{v, i, ix, iy, iz, jx, jy, jz},
  If[ListQ[var],
     VarGridShift[#, s, grid, CenterQ] &/@ var,
     {v, i, ix, iy, iz} = Level[var, Infinity];
     {jx, jy, jz} = If[CenterQ, PbcDiff[{ix, iy, iz} + s, grid], Mod[{ix, iy, iz} + s, grid, 1]];
     Subscript@@{v, i, jx, jy, jz}]
]

VarOnGrid[var_, site_] := Module[{varlength},
  varlength = Length[Level[var, Infinity]];
  If[varlength == 2,
     Subscript @@ Join[Level[var, Infinity], site],
     var]
]

ExprOnGrid[expr_, site_] := Module[{vars, varsub},
  If[ListQ[expr],
     ExprOnGrid[#, site] &/@ expr,
     vars = Variables[expr];
     varsub = Thread[vars -> (VarOnGrid[#, site] &/@ vars)];
     Return[expr/.varsub]]
]

Expr2Gamma[model_, grid_] := Module[{ngx, ngy, ngz, ix, iy, iz, sites, s, out, m1, m2},
  If[ListQ[model], 
  Expr2Gamma[#, grid] &/@ model,
  {ngx, ngy, ngz} = grid;
  sites = Flatten[Table[PbcDiff[{ix, iy, iz} - {1, 1, 1}, grid], {iz, ngz}, {iy, ngy}, {ix, ngx}], 2];
  If[AtomQ[model]||Length@Level[Quiet@First@Variables[model], Infinity] == 5, model, Expand@Sum[ExprOnGrid[SimplifyCommonFactor[model, 10^-9], s], {s, sites}]]]
]

SortVarSub[varsub_] := Module[{grid, ngx, ngy, ngz, ix, iy, iz, sites, fullvars, supervars, v, out},
  grid = Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ (Sub2List[varsub][[1]])];
  {ngx, ngy, ngz} = grid;
  fullvars = DeleteDuplicates@VarBare[First@Sub2List[varsub]];
  sites = Flatten[Table[PbcDiff[{ix, iy, iz} - {1, 1, 1}, grid], {iz, ngz}, {iy, ngy}, {ix, ngx}], 2];
  supervars = Flatten[Table[VarOnGrid[v, #], {v, fullvars}] & /@ sites];
  out = Thread[supervars -> (supervars /. varsub)];
  Return[out]
]

ShiftVarSub[varsub_, shift_] := Module[{grid, ngx, ngy, ngz, ix, iy, iz, jx, jy, jz, sites, fullvars, supervars, v, s, out},
  grid = Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ (Sub2List[varsub][[1]])];
  {ngx, ngy, ngz} = grid;
  out = Table[{v, i, ix, iy, iz} = Level[First@s, Infinity];
              {jx, jy, jz} = PbcDiff[{ix, iy, iz} + shift, grid];
              (Subscript @@ {v, i, jx, jy, jz}) -> s[[2]], {s, varsub}];
  Return[SortVarSub[out]]
]

SwapVarSub[varsub_] := Module[{grid, ngx, ngy, ngz, ix, iy, iz, jx, jy, jz, sx, sy, sz, sites, fullvars, supervars, v, s, out, shift, tmp},
  grid = Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ (Sub2List[varsub][[1]])];
  {ngx, ngy, ngz} = grid;
  out = Table[shift = {sx, sy, sz} - {1, 1, 1};
              tmp = Table[{v, i, ix, iy, iz} = Level[First@s, Infinity];
                          {jx, jy, jz} = PbcDiff[{ix, iy, iz} + shift, grid];
                          (Subscript @@ {v, i, jx, jy, jz}) -> s[[2]], {s, varsub}];
                          SortVarSub[tmp], {sz, ngz}, {sy, ngy}, {sx, ngx}];
  Return[Flatten[out, 2]]
]

VarSubStar[t5sub_, t6sub_, SymMat_, OptionsPattern["absvalue"->True]] := Module[{grid, nvars0, vars, svars, t5, t6, out},
  vars = Sub2List[t5sub][[1]];
  svars = Sub2List[t6sub][[1]];
  grid = Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ vars];
  nvars0 = Length[vars]/Times @@ grid;
  out = Table[t5 = Thread[vars -> Flatten[(m[[1]] . #) & /@ Partition[Sub2List[t5sub][[2]], nvars0]]];
              t6 = Thread[svars -> m[[2]] . (Sub2List[t6sub][[2]])];
              If[OptionValue["absvalue"], Sign[Sub2List[Join[#, t6]][[2]]], Join[#, t6]] & /@ SwapVarSub[t5], {m, SymMat}];
  Return[DeleteDuplicates@Flatten[out\[Transpose], 1]]
]

OnSiteExprQ[expr_] := Module[{out, inv, n, ix, iy, iz, nsite},
  If[ListQ[expr], 
     OnSiteExprQ[#] &/@ expr, 
     inv = If[Head[expr] === Plus, First@Level[expr, {1}], expr];
     nsite = Length@DeleteDuplicates[Cases[inv, Subscript[v_, i_, ix_, iy_, iz_] -> {ix, iy, iz}, Infinity]];
     Return[nsite == 1]];
]

MyTime[inp_, NewProjQ_] := Module[{proj, txt, wedge, dir0, projects0, projects, items, pdata, name, AccumulateTime, TabData, ChartData, pltdata, chart, l, h, ip, timeused, t1, t0, timetag, WorkingQ, entry, entry0, AddNewQ, newnameQ, newproj, priority, PriorityPlt, spoly},
  dir0 = "/home/paulchern/Documents/Workshop/Working";
  {proj, txt} = inp;
  projects0 = Partition[ReadList[dir0 <> "/MyTime.dat", Word, RecordLists -> True, WordSeparators -> {"\n"}, RecordSeparators -> {"-----\n"}], 2];

  AddNewQ = Nor @@ Table[StringContainsQ[First@StringSplit[Last[projects0[[ip, 2]]]], "T"], {ip, Length@projects0}];
  newnameQ = Nor @@ Table[projects0[[ip, 1, 1]] == proj, {ip, Length@projects0}];

  projects = If[NewProjQ, If[newnameQ, Append[projects0, {{proj, "0.0"}, {DateString["ISODate"] <> " " <> RomanNumeral[0] <> " " <> txt}}], projects0],
    Table[h = Length[projects0[[ip, 2]]];
          name = StringTrim[projects0[[ip, 1, 1]]];
       timetag = First@StringSplit[Last[projects0[[ip, 2]]]];
      WorkingQ = proj == name && StringContainsQ[timetag, "T"];
            t1 = TimeObject@DateString["ISODateTime"];
      timeused = If[WorkingQ, QuantityMagnitude@UnitConvert[t1 - TimeObject@timetag, "hour"], 0];
        entry0 = DateString["ISODateTime"] <> " " <> RomanNumeral[h] <> " 0.0 " <> txt;
         entry = DateString["ISODate"] <> " " <> RomanNumeral[h - 1] <> " " <> ToString[timeused] <> " " <> txt;
       {{projects0[[ip, 1, 1]], ToString@Total[ToExpression[StringSplit[#][[3]]] & /@ (projects0[[ip, 2]][[2 ;;]])]},
        If[AddNewQ && proj == name, Append[projects0[[ip, 2]], entry0], Table[If[WorkingQ && i == h, entry, projects0[[ip, 2]][[i]]], {i, h}]]}, {ip, Length@projects0}]];

  {TabData, ChartData} = Table[name = projects[[ip]][[1, 1]];
                     AccumulateTime = ToExpression[projects[[ip]][[1, 2]]];
                              pdata = StringSplit[#] & /@ (projects[[ip, 2]]);
                              items = {#[[1]] <> " " <> #[[2]], #[[3]], StringRiffle[#[[4 ;;]], " "]} & /@ pdata;
                              {name -> MenuView[#1 -> #2 <> "\n" <> #3 & @@@ items, Length[items]], {name, {#1, ToExpression@#2, #3} & @@@ items}}, {ip, Length@projects}]\[Transpose];

  pltdata = Table[{Total[#2 & @@@ (ChartData[[ip, 2]][[2 ;;]])], Length[ChartData[[ip, 2]]] - 1}, {ip, Length@ChartData}];

  spoly = Max[pltdata\[Transpose][[1]]]/Max[Table[ChartData[[ip, 2, 1]][[2]], {ip, Length[ChartData]}]];

  PriorityPlt = Polygon[Join[Table[{-spoly*ChartData[[ip, 2, 1]][[2]]/2, Length[ChartData] - ip}, {ip, Length[ChartData]}],
                     Reverse@Table[{spoly*ChartData[[ip, 2, 1]][[2]]/2, Length[ChartData] - ip}, {ip, Length[ChartData]}]]];

  chart = Graphics[{Table[{l, h} = pltdata[[ip]]; 
                        priority = Length@ChartData - ip;
                          {EdgeForm[Directive[Cyan, Thickness[0.02]]], Rectangle[#1, #2, RoundingRadius -> 0.1]} & @@ Transpose[{{-l/2, l/2}, {priority - 1/2, priority + 1/2}}], {ip, Length@ChartData}],
                    Table[{l, h} = pltdata[[ip]]; 
                        priority = Length@ChartData - ip;
                          {White, Table[Line[{#1, #2}] & @@ ({{-l/2, n/(h + 1) + priority - 1/2}, {l/2, n/(h + 1) + priority - 1/2}}), {n, h}]}, {ip, Length@ChartData}],
                    Table[{l, h} = pltdata[[ip]]; 
                        priority = Length@ChartData - ip;
                        Text[Column[Style[#, Black, Bold, 18] & /@ {ChartData[[ip, 1]], "Time: " <> ToString[l], "Frequency: " <> ToString[h]}, Left], {Max[pltdata\[Transpose][[1]]], priority}, {Left, Center}, Background -> LightRed], {ip, Length@ChartData}],
                    {FaceForm[], EdgeForm[Red], PriorityPlt}}, ImageSize -> {300, 300}, AspectRatio -> GoldenRatio];

  wedge = Grid@{{TabView[TabData, ImageSize -> {500, 300}], "        ", chart}};
  Print[wedge];
  Export[dir0 <> "/MyTime.dat", Flatten[Table[Transpose[{#}] & /@ Insert[projects[[ip]], {"-----"}, {{1}, {2}}], {ip, Length[projects]}], 2]];
]

ColorPots[v_, n_, tol_ : 10^-6] := Module[{c, test, vc}, 
   c = Association[Append[{#1, #2, #3} -> RGBColor[(#1 + 1)/2, (#2 + 1)/2, (#3 + 1)/2, 0.75] & @@@ Select[Tuples[{0, 1, -1}, 3], Count[Abs@#, 1] == n &], {0, 0, 0} -> RGBColor[0.5, 0.5, 0.5, 0.75]]];
   test = Min[TakeLargest[Abs@Chop[v, tol], n]];
   vc = If[Abs@# >= test, Sign[#], 0] & /@ v;
   Return[c[vc]]
];

Irrationalize[x_, prec_, OptionsPattern[{"NCUT" -> 1000}]] := Module[{pool, eps, n},
  If[ListQ[x], Irrationalize[#, prec, "NCUT" -> OptionValue["NCUT"]] & /@ x,
     If[Head[x] === Real, 
        pool = Flatten@Table[Power[Range[OptionValue["NCUT"]], 1/2]/n, {n,1,10}];
        eps = First@MinimalBy[{Sign[x] #, Abs@N[# - Abs@x]} & /@ pool, #[[2]] &];
        If[eps[[2]] < prec, eps[[1]], False], x]]
]

NumberCommonDivisor[NumList_,prec_:10^-12] := Module[{TempList, DenominatorLCM},
 TempList = Which[Head[#] === Real, Round[#,prec], Head[#] === Integer, #, Head[#] === Times, First@Level[#, {1}], Head[#] === Power, 1, Head[#] === Rational, #, (Head[#] === Complex)&&(Re[#]!=0), Abs[#], (Head[#] === Complex)&&(Re[#]==0), #] &/@ NumList;
 DenominatorLCM = If[MemberQ[TempList, _Rational], LCM @@ (Denominator[#] & /@ Extract[TempList, Position[TempList, _Rational]]), 1];
 DenominatorLCM = If[AllTrue[TempList, Negative], -DenominatorLCM, DenominatorLCM];
 Return[{DenominatorLCM, GCD @@ (TempList DenominatorLCM)}]
]

GetConstantFactor[expr_] := Module[{},
  If[ListQ[expr],  GetConstantFactor[#] & /@ expr, Return[(If[MatchQ[Expand@expr, Plus[_, __]], Level[Expand@expr, {1}], Level[Expand@expr, {0}]] /. Thread[Variables[Expand@expr] -> ConstantArray[1, Length[Variables[Expand@expr]]]])]]
]

SimplifyIrrationalNumber[x_, prec_, OptionsPattern["NCUT" -> 1000]] := Module[{coeff, expr, f, flist},
  If[ListQ[x], SimplifyIrrationalNumber[#, prec, "NCUT" -> OptionValue["NCUT"]] & /@ x, 
     flist = {#, GetConstantFactor[x]/#} &/@ (GetConstantFactor[x]);
     f = First[If[AllTrue[#2, Rationalize[Abs@#] >= 1 &], #1, ## &[]] & @@@ flist];
     {coeff, expr} = DecomposeExpr[Expand[x/f]];
     Total[Irrationalize[If[Head[#]===Rational,1.0 #, #] &/@ coeff, prec, "NCUT" -> OptionValue["NCUT"]] expr]]
]

DecomposeExpr[x_] := Module[{sub, expr, coeff},
  sub = Thread[Variables[x] -> ConstantArray[1, Length[Variables[x]]]];
  expr = If[Head[x] === Plus, Level[x, {1}], {x}];
  coeff = If[Head[x] === Plus, Level[x, {1}], {x}] /. sub;
  Return[{coeff, Expand[expr/coeff]}]
]

SimplifyCommonFactor[expr_, prec_] := Module[{factorLCM, factorGCD},
  Which[ListQ[expr],
        SimplifyCommonFactor[#, prec] & /@ expr,
        expr === 0,
        Return[0],
        True,
        {factorLCM, factorGCD} = NumberCommonDivisor[GetConstantFactor[Expand[expr]],prec prec];
        Return[Expand[expr factorLCM/factorGCD]/. x_Real :> Round[x,prec]]
       ]
]

SimplifyTensorCommonFactor[mat_] := Module[{factorLCM, factorGCD, out, coeff, f},
  coeff = Flatten@GetConstantFactor[mat];
  f = If[DeleteCases[Sign@coeff, 0] ==={}, 1, First@DeleteCases[Sign@coeff, 0]];
  {factorLCM, factorGCD} = NumberCommonDivisor[f coeff];
  If[factorGCD == 0, factorGCD = 1];
  out = Expand[factorLCM mat/factorGCD];
  Return[out]
]

SimplifyTensor[T0_, TT_] := Module[{tsub, s, solutions, out, MatReIm},
  MatReIm = Table[solutions = Transpose[{Flatten[ComplexExpand[f[TT]]], Flatten[ComplexExpand@f[T0]]}];
                  tsub = Flatten[Table[If[s[[1]] === 0, ## &[], Thread[Flatten[s[[1]] {1, -1}] -> Flatten[s[[2]] {1, -1}]]], {s, solutions}]];
                  ComplexExpand@f[TT] /. tsub, {f, {Re, Im}}];
  out = {1, I} . MatReIm;
  Return[out]
]

StrainInvQ[x_, e_] := Module[{}, If[ListQ[x], StrainInvQ[#, e] & /@ x, MemberQ[Level[x, All], e]]]

XInvQ[x_, e_, bool_ : ContainsExactly] := Module[{v},
  If[ListQ[x],
     XInvQ[#, e, bool] & /@ x,
     bool[Cases[Variables[x], Subscript[v_, __] -> v], e]]
]

TSGrid[ts_] := Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ (Sub2List[First[ts][[5]]][[1]])];
VarSubGrid[vs_] := Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ (Sub2List[vs][[1]])];

TSVars[ts_] := Module[{grid, vars, svars, out},
  grid = TSGrid[ts];

  vars  = Sub2List[ts[[1]][[5]]][[1]];
  svars = Sub2List[ts[[1]][[6]]][[1]];
  out   = Join[Partition[vars, Length@vars/(Times@@grid)], {svars}];

  Return[out]
]

T5Vars[t_]  := Sub2List[t[[5]]][[1]]
T6Vars[t_]  := Sub2List[t[[5]]][[1]]
T5Data[t_]  := Sub2List[t[[5]]][[2]]
T6Data[t_]  := Sub2List[t[[6]]][[2]]
T56Data[t_] := Join[T5Data[t], T6Data[t]]

VarHead[var_] := Module[{},
  If[ListQ[var],
     VarHead[#] &/@ var,
     First@Level[var, Infinity]]
]

VarInd[var_] := Module[{info},
  If[ListQ[var],
     VarInd[#] &/@ var,
     info = Level[var, Infinity];
     If[Length@info == 3, info[[2;;3]], info[[2]]]]
]

VarSite[var_] := Module[{info},
  If[ListQ[var],
    VarSite[#] &/@ var,
    info = Level[var, Infinity];
    If[Length@info == 5, info[[3;;5]], Print["Error: var not on grid"];Abort[]]]
]

VarBare[var_] := Module[{info},
   If[ListQ[var],
      VarBare[#] &/@ var,
      info = Level[var, Infinity];
      If[Length@info > 3, Subscript @@ (info[[1;;2]]), var]]
]

OpenExtract[f_, str_, i_, nlines_ : ListQ] := Module[{s, t, num, nskip, nread, out, TotNumLine, StreamPositionList},
  {nskip, nread} = nlines;
  s = OpenRead[f];
  t = ReadString[s];
  t = StringSplit[t, "\n"];
  StreamPositionList = StringLength[#] & /@ t;
  TotNumLine = Length@t;
  SetStreamPosition[s, 0];

  out = Which[StringQ[str], num = Length[FindList[s, str]]; 
                            SetStreamPosition[s, 0]; 
                            Table[Find[s, str]; Do[ReadLine[s], {nskip}]; Table[ReadLine[s], {nread}], {Min[{i, num}]}],
              IntegerQ[str], SetStreamPosition[s, Total[StreamPositionList[[;; str - 1]]]+i]; Table[ReadLine[s], {nread}], 
              True, Abort[]];

  Close[s];
  Return[out]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

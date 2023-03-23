BeginPackage["LINVARIANT`MathematicaPlus`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
NOrderResponse         ::usage "NOrderResponse[eqn, var, n]"
GetVariationVar        ::usage "GetVariationVar[var]"
Complex2Exp            ::usage "Complex2Exp[exp]"
Exp2Complex            ::usage "Exp2Complex[exp]"
GetSubscriptInfo       ::usage "GetSubscriptInfo[atom]"
ReadNonEmptyLine       ::usage "ReadNonEmptyLine[stream]"
PermuteThrough         ::usage "PermuteThrough[tab]"
Var2Var                ::usage "Var2Var[var1, var2, dir]"
ParseFortranNumber     ::usage "ParseFortranNumber[stream]"
MakeMatrixBlock        ::usage "MakeMatrixBlock[mat, dim]"
SimpsonIntegrate       ::usage "SimpsonIntegrate[f, x]"
RectangleIntPath       ::usage "ComplexIntegratePath[spt, ept, npt, ratio]"
partitionBy            ::usage "partitionBy[l, p]"
CloneReshape2D         ::usage "CloneReshape2D[array0, array1]"
Model2P1               ::usage "Model2P1[model]"
FixDoubleCounting      ::usage "FixDoubleCounting[expr]"
OnSiteExprQ            ::usage "OnSiteExprQ[expr]"
MyTime                 ::usage "MyTime[inp, NewProjQ]"

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

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

BeginPackage["LINVARIANT`MathematicaPlus`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
NOrderResponse       ::usage "NOrderResponse[eqn, var, n]"
GetVariationVar      ::usage "GetVariationVar[var]"
Complex2Exp          ::usage "Complex2Exp[exp]"
Exp2Complex          ::usage "Exp2Complex[exp]"
GetSubscripteInfo    ::usage "GetSubscripteInfo[atom]"
ReadNonEmptyLine     ::usage "ReadNonEmptyLine[stream]"
PermuteThrough       ::usage "PermuteThrough[tab]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
Complex2Exp[exp_] := Module[{},
  Expand[exp /. z_?NumericQ :> Abs[z] Exp[I Arg[z]]]
]

Exp2Complex[exp_] := Module[{},
  Expand[ReIm[exp].{1, I}]
]

GetVariationVar[var_] := Module[{\[Delta]var, varnew},
  \[Delta]var = ToString[#] & /@ Level[var, 1];
  varnew = "\!\(\*SubscriptBox[\(" <> "\[Delta]" <> \[Delta]var[[1]] <> "\),\(" <> StringJoin[Riffle[\[Delta]var[[2 ;;]], ","]] <> "\)]\)";
  Return[varnew]
]

NOrderResponse[eqn_, var_, n_] := Module[{\[Epsilon], exp, varnew},
  varnew = ToExpression[GetVariationVar[#]] &/@ var;
  exp = Normal[Series[eqn /. Thread[var -> (var + \[Epsilon] varnew)], {\[Epsilon], 0, n}]];
  Expand[(exp /. {\[Epsilon] -> 1}) - (exp /. {\[Epsilon] -> 0})]
]

GetSubscripteInfo[atom_] := Module[{},
  If[MatchQ[atom, _List], GetSubscripteInfo[#] & /@ atom, Which[AtomQ[atom], ## &[], MatchQ[atom, _Subscript], Level[atom, 1], MatchQ[atom, _Power], If[MatchQ[Level[atom, 1][[2]], _Integer], GetSubscripteInfo[#] & /@ ConstantArray[Level[atom, 1][[1]], Level[atom, 1][[2]]], ##&[]], True, GetSubscripteInfo[#] & /@ Level[atom, 1]]]
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
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

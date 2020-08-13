BeginPackage["LINVARIANT`MathematicaPlus`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
NOrderResponse       ::usage "NOrderResponse[eqn, var, n]"
GetVariationVar      ::usage "GetVariationVar[var]"
Complex2Exp          ::usage "Complex2Exp[exp]"
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

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

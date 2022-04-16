BeginPackage["LINVARIANT`Parser`", {"LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ParseXML                     ::usage "ParseXML[xml, tag, label, level]"
ParseXMLData                 ::usage "ParseXMLData[xml, DataType]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseXML[xml_, tag_, label_, level_: Infinity] := Module[{xmldata, DataType, x},
  xmldata = Cases[xml, XMLElement[tag, Flatten[#1 -> #2 & @@@ If[label==={},label,If[Length@Dimensions[label]==1, {label}, label]]], x_] :> x,   level];
  Return[xmldata]
]

ParseXMLData[xml_, DataType_] := Module[{xmldata},
  Flatten[Cases[#, XMLElement[DataType, {}, x_] :> ParseFortranNumber[x]], 1] & /@ xml
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

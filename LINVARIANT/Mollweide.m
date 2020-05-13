BeginPackage["LINVARIANT`Mollweide`", {"LINVARIANT`Boracite`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetMollweideEne      ::usage "GetImage[Eonsite, vars, init]"
Cart2MollWeide       ::usage "Cart2MollWeide[{theta, phi}]"
MollWeide2Cart       ::usage "MollWeide2Cart[{vx, vy, vz}, vars, image]"
GetPES               ::usage "GetPES[HOnsite, vars, Xaxelabel, texts]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)
rr;
(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetMollweideEne[HOnsite_, mapvars_, \[Theta]_, \[Phi]_, Q3D_] := Module[{sub, vars, X5, init, min, x, op}, 
  X5 = {aa, AA, bb, BB, CC, cc};
  vars = Variables[HOnsite];
  AxesHybrid[s_] := Return[Apply[Plus, (({(1 + #)/2, (1 - #)/2} & /@ s)*{{aa, CC}, {bb, AA}, {cc, BB}})\[Transpose]]];
  Which[Length@mapvars == 6,
        sub = Join[Thread[{aa, bb, cc} -> 1/Sqrt[2] Abs[CoordinateTransformData["Spherical" -> "Cartesian", "Mapping", {rr, \[Theta], \[Phi]}]]], Thread[{AA, BB, CC} -> 1/Sqrt[2] CoordinateTransformData["Spherical" -> "Cartesian", "Mapping", {rr, \[Theta], \[Phi]}]]];
        init = {Select[vars, FreeQ[mapvars, #] &], ConstantArray[0, Length@Select[vars, FreeQ[mapvars, #] &]]}\[Transpose],
        Length@mapvars == 7,
        s = Sign[#] & /@ CoordinateTransformData["Spherical" -> "Cartesian", "Mapping", {1, \[Theta], \[Phi]}];
        sub = If[Length[#[[1]]] == 0, #, Unevaluated[Sequence[]]] & /@ Thread[AxesHybrid[s] -> CoordinateTransformData["Spherical" -> "Cartesian", "Mapping", {rr, \[Theta], \[Phi]}]];
        sub = Table[If[MemberQ[Table[x[[1]], {x, sub}], op], op -> (op /. sub), op -> 0], {op, X5}];
        init = {Select[vars, FreeQ[mapvars, #] &], ConstantArray[0, Length@Select[vars, FreeQ[mapvars, #] &]]}\[Transpose];
        init =  Table[{x[[1]], Which[x[[1]]===Subscript[P,1], -Sign[bb*BB]/.sub, x[[1]]===Subscript[P,2], -Sign[aa*AA]/.sub, x[[1]]===Subscript[P,3], -Sign[cc*CC]/.sub, True, 0]}, {x,init}]/.{rr->0.0001},
        Length@mapvars == 3,
        sub = Thread[mapvars -> CoordinateTransformData["Spherical" -> "Cartesian", "Mapping", {rr, \[Theta], \[Phi]}]];
        init = {Select[vars, FreeQ[mapvars, #] &], ConstantArray[0, Length@Select[vars, FreeQ[mapvars, #] &]]}\[Transpose];

        init = Join[{{Subscript[\[Epsilon], 1, 1], Subscript[\[Epsilon], 2, 2], Subscript[\[Epsilon], 3, 3], Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3], Subscript[\[Epsilon], 1, 2]},ConstantArray[0,6]}\[Transpose], {{aa, 0.1 Sign[Subscript[P,2]]}, {AA,-0.1}, {bb, 0.1 Sign[Subscript[P,1]]}, {BB, -0.1}, {cc, 0.1 Sign[Subscript[P,3]]}, {CC, -0.1}}/.sub/.{rr->1}]
        ];
  min = If[Q3D, 
           FindMinimum[HOnsite /. sub, Join[init, {{rr, 1}}], MaxIterations -> 10000], 
           Minimize[HOnsite /. sub, (Join[init, {{rr, 1}}]\[Transpose])[[1]]]
           ];
  Return[min[[1]]];
]


Cart2MollWeide[{theta_, phi_}]:=Module[{xy, lambda},
  lambda = If[Abs[phi] == Pi/2, phi, lambda /. FindRoot[2 lambda + Sin[2 lambda] == Pi Sin[phi], {lambda, phi}]];
  xy = {2/Pi*theta Cos[lambda], -Sin[lambda]};
  Return[xy];
];

MollWeide2Cart[mepvars_, vars_, image_]:=Module[{rthetaphi},
  rthetaphi = CoordinateTransformData["Cartesian" -> "Spherical", "Mapping", mepvars] /. Thread[vars -> image];
  Return[rthetaphi];
];

GetPES[HOnsite_, vars_, texts_, OptionsPattern[{"3D"->True, "NumDeltaAngle"->18, "NumContour"->50, "NebImages"->{}}]] := Module[{grid, DeltaAngle, theta, phi, pdata, disprange, PES, gridline0, gridline1, LonLat, Xaxelabel, image, ipath, neb, nebpath},
  DeltaAngle = Pi/2/OptionValue["NumDeltaAngle"];
  XX123[var_]:=Which[var===aa, Subscript["a",1], var===AA, Subscript["a",2],
                     var===bb, Subscript["b",1], var===BB, Subscript["b",2],
                     var===cc, Subscript["C",1], var===CC, Subscript["C",2],
                     True, var];
  BarVar[var_]:=Subscript[OverBar[var[[1]]], var[[2]]];
  Which[Length@vars==3||Length@vars==6,
        Xaxelabel = {BarVar[XX123[vars[[1]]]],
                     BarVar[XX123[vars[[2]]]],
                     XX123[vars[[2]]],
                     XX123[vars[[2]]],
                     BarVar[XX123[vars[[1]]]]},
        Length@vars==7,
        Xaxelabel = {Subscript["b",1],
                     Subscript["a",1],
                     Subscript["a",2],
                     Subscript["c",2],
                     Subscript["b",1]}
               ];
  grid = Table[{phi, theta}, {theta, -Pi/2 + 10^-6, Pi/2, DeltaAngle}, {phi, -Pi + 10^-6, Pi + 10^-10, DeltaAngle}];
  nebpath = Table[neb={#3, #2 - Pi/2} &@@@ Table[MollWeide2Cart[vars, Variables[HOnsite], image], {image, OptionValue["NebImages"][[ipath]]}];
                  neb=Pi/2 Cart2MollWeide[#] & /@ neb;
                  {{Black, White}[[ipath]], Line[neb], PointSize[Medium], {Black, White}[[ipath]], Point[neb]}, {ipath, Length@OptionValue["NebImages"]}];
  pdata = {#1, #2, GetMollweideEne[HOnsite, vars, #2 + Pi/2, #1, OptionValue["3D"]]} & @@@ Flatten[grid, 1];
  disprange = MinMax@(pdata\[Transpose][[3]]);
  (* Plot grids *)

  gridline0 = Graphics[{AbsoluteThickness[0.05], Line /@ grid, Line /@ Transpose[grid]}, AspectRatio -> 1/2];
  gridline1 = Flatten[{gridline0[[1, 2]][[Range[OptionValue["NumDeltaAngle"]/2]*4 - 1]], gridline0[[1, 3]][[Range[OptionValue["NumDeltaAngle"]]*4 - 3]]}] // Graphics[{AbsoluteThickness[0.2], #}] &;
  LonLat = gridline1 /. Line[pts_] :> Line[Pi/2 Cart2MollWeide /@ pts] /. Point[pts_] :> Point[Pi/2 Cart2MollWeide[pts]];
  If[Length@vars==7,
     texts = Join[texts,Graphics@Text[Style[Subscript["b",2], 12], {0.05, -0.15}, {-1, 0}]],
     Unevaluated[Sequence[]]
    ];
  txtorigin = Graphics@Text[Style[XX123[vars[[1]]], 12], {0.05, -0.15}, {-1, 0}];
  (* Plot Contour *)

  PES = ListContourPlot[pdata,
                        Epilog -> nebpath,
                        PlotRange->All,
                        Contours -> OptionValue["NumContour"],
                        InterpolationOrder -> 2,
                        ColorFunctionScaling -> False,
                        ColorFunction -> (ColorData["Rainbow"][Rescale[#, disprange]] &),
                        PlotLegends -> Placed[BarLegend[{(ColorData["Rainbow"][Rescale[#, disprange]] &), disprange},
                        LabelStyle -> {Black, Bold, 10}, LegendLayout -> "Column",
                        LegendMargins -> 0,
                        LegendLabel -> Placed[Style["Energy (meV)", 20], Top]], Right],
                        AxesLabel -> {"", XX123[vars[[3]]]},
                        Frame -> None,
                        Axes -> True,
                        AxesStyle -> Directive[Black, 18],
                        Ticks -> {{{-Pi, -(Pi/2), 0, Pi/2, Pi}, Xaxelabel}\[Transpose], {{-(Pi/2), ""}, {Pi/2, ""}}},
                        TicksStyle -> Directive[Black, 12],
                        AspectRatio -> 1/2,
                        ImageSize -> 600];
  mollweide = Pi/2 Cart2MollWeide[#] & /@ PES[[1]][[1]][[1]];
  PES[[1]][[1]][[1]] = mollweide;
  Print[Show[Join[{PES, txtorigin}, texts, {LonLat}]]];
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

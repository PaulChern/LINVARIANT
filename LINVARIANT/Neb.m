BeginPackage["LINVARIANT`Neb`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetImage      ::usage "GetImage[Eonsite, vars, init]"
MPInverse     ::usage "MPInverse[M]"
GetMEP        ::usage "GetMEP[Eonsite, image0, imageN, Nimage, Niter]"
GetMEP0       ::usage "GetMEP[Eonsite, image0, imageN, Nimage, Niter]"
GetMEP1       ::usage "GetMEP[Eonsite, image0, imageN, Nimage, Niter]"
GetMEP3D      ::usage "GetMEP[Eonsite, image0, imageN, Nimage, Niter]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetImage[Eonsite_, vars_, init_] := Module[{image, varsfix, min},
  varsfix = vars;
  Do[varsfix = Delete[varsfix, Position[varsfix, v]], {v, init[[1]]}];
  min = FindMinimum[Eonsite /. Thread[varsfix -> ConstantArray[0, Length@varsfix]], init\[Transpose]];
  image = vars /. Join[min[[2]], Thread[varsfix -> ConstantArray[0, Length@varsfix]]];
  Return[image]
]

GetMEP0[Eonsite_, AllImages_, Nimage_, Niter_, Kratio_, OptionsPattern[{"etol"->0, "Kpara"->10^{4}}]] := Module[{GF, KK, vars, NebVars, \[CapitalDelta]x, \[CapitalDelta]xNormal, \[CapitalDelta]xPara, ineb, i, n, t, Etol, disp, images, Eimages, BondImages, TangentImages, StepImages, error, out=<|"op"->{},"dispara"->{}, "disperp"->{},"dis"->{},"\[CapitalDelta]"->{}|>},
  KK = OptionValue["Kpara"];
  images=AllImages[[-1]];
  vars = Variables@Eonsite;
  Eimages = Eonsite /. (Thread[vars -> #] & /@ images);
  
  Print[ListPlot[{Range[Nimage + 2], Eimages}\[Transpose], 
                         PlotRange -> Full, 
                         LabelingFunction -> (Callout[#1[[2]], Automatic] &), 
                         PlotMarkers -> {Automatic, 8},
                         Joined -> True,
                         Frame -> True, 
                         GridLines -> Automatic, 
                         ImageSize -> 250]];
  Etol=Max@Eimages;
  \[CapitalDelta]x = Table[#[n] & /@ vars, {n, Nimage + 2}] /. Thread[Flatten[Table[#[n] & /@ vars, {n, {1, Nimage + 2}}]] -> ConstantArray[0, 2 Length@vars]];
  NebVars = Flatten[\[CapitalDelta]x[[2 ;; Nimage + 1]]];

  Do[
   BondImages = Differences[images];
   TangentImages = Table[Normalize[BondImages[[n]] + BondImages[[n + 1]]], {n, Range[Nimage]}];
   \[CapitalDelta]xNormal = Table[\[CapitalDelta]x[[n + 1]] - \[CapitalDelta]x[[n + 1]].TangentImages[[n]] TangentImages[[n]], {n, 1, Nimage}];
   \[CapitalDelta]xPara = Table[\[CapitalDelta]x[[n + 1]].TangentImages[[n]] TangentImages[[n]], {n, 1, Nimage}];
   GF = Sum[Eonsite/.(Thread[vars->images[[n+1]]+\[CapitalDelta]xNormal[[n]]]),{n,1,Nimage}] 
      + Sum[KK((BondImages[[n+1]]-BondImages[[n]]+\[CapitalDelta]x[[n]]+\[CapitalDelta]x[[n+2]]-2\[CapitalDelta]x[[n+1]]).TangentImages[[n]])^2,{n,1,Nimage},{n,1,Nimage}]
      + Sum[KK*Kratio[[ineb]]*((BondImages[[n+1]]-BondImages[[n]]+\[CapitalDelta]x[[n]]+\[CapitalDelta]x[[n+2]]-2\[CapitalDelta]x[[n+1]])-(BondImages[[n+1]]-BondImages[[n]]+ \[CapitalDelta]x[[n]]+\[CapitalDelta]x[[n+2]]-2\[CapitalDelta]x[[n+1]]).TangentImages[[n]]*TangentImages[[n]]).((BondImages[[n+1]]-BondImages[[n]]+\[CapitalDelta]x[[n]]+\[CapitalDelta]x[[n+2]]-2\[CapitalDelta]x[[n+1]])-(BondImages[[n+1]]-BondImages[[n]]+\[CapitalDelta]x[[n]]+\[CapitalDelta]x[[n+2]]-2\[CapitalDelta]x[[n+1]]).TangentImages[[n]]*TangentImages[[n]]),{n,1,Nimage}];
   If[Nimage>0,
      {error, StepImages} = FindMinimum[GF, {NebVars, ConstantArray[0, Length@NebVars]}\[Transpose],MaxIterations->100000],
      {error, StepImages} = Minimize[GF, NebVars]];
   images = images + 1/10(\[CapitalDelta]x /. StepImages);
   Eimages = Eonsite /. Thread[vars -> #] & /@ images;
   AppendTo[out["op"],images];
   AppendTo[out["\[CapitalDelta]"],\[CapitalDelta]x /. Join[StepImages, Thread[Flatten[Table[#[n] & /@ vars, {n, {1, Nimage + 2}}]] -> ConstantArray[0, 2 Length@vars]]]];
   AppendTo[out["dis"],Norm[#]&/@BondImages];
   AppendTo[out["dispara"],Norm[#]&/@Table[disp=(BondImages[[n + 1]] - BondImages[[n]]); disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}]];
   AppendTo[out["disperp"],Norm[#]&/@Table[disp=(BondImages[[n + 1]] - BondImages[[n]]); disp - disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}]];
   If[Abs[Max@Eimages-Etol] < OptionValue["etol"], Print["Tolerence "<>ToString[OptionValue["etol"]]<>" met."];Break[], Etol=Max@Eimages];
   , {ineb, Niter}];
  Print[ListPlot[{Range[Nimage + 2], Eimages}\[Transpose], 
                 PlotRange -> Full, 
                 LabelingFunction -> (Callout[#1[[2]], Automatic] &), 
                 PlotMarkers -> {Automatic, 8}, 
                 Joined -> True,
                 Frame -> True, 
                 GridLines -> Automatic, 
                 ImageSize -> 250]];
  EmitSound@Play[Sin[300 t Sin[20 t]], {t, 0, 3}];
  Return[{out, Join[AllImages,out["op"]]}];
]

GetMEP3D[Eonsite_, image0_, imageN_, Nimage_, Niter_] := Module[{vars, images, i, n, Eimages, UVecImages, BondImages, TangentImages, HessianImages, GradImages, ForceSpringImages, ForceKinkImages, ForceNormalImages, StepImages},
  vars = Variables@Eonsite;
  images = Table[image0 + n/(Nimage + 1) (imageN - image0), {n, 0, Nimage + 1}];
  Eimages = Eonsite /. (Thread[vars -> #] & /@ images);
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> Full, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]];
  Do[
   GradImages = Table[D[Eonsite, {vars, 1}] /. Thread[vars -> im], {im, images}];
   HessianImages = Table[Chop[D[Eonsite, {vars, 2}] /. Thread[vars -> im]], {im, images}];
   (*
   Print[LinearAlgebra`Private`MatrixConditionNumber[#]&/@HessianImages];
   Print[HessianImages[[2]]//MatrixForm];
   *)
   UVecImages = Normalize[#] & /@ Differences[images];
   BondImages = # & /@ Differences[images];
   TangentImages = Table[Normalize[UVecImages[[i]] + UVecImages[[i + 1]]], {i, Range[Nimage]}];
   ForceSpringImages = 1.0 Table[(BondImages[[n + 1]] - BondImages[[n]]).TangentImages[[n]] TangentImages[[n]], {n, Range[Nimage]}];
   ForceKinkImages = 0.1 Table[(BondImages[[n + 1]] - BondImages[[n]]) - (BondImages[[n + 1]] - BondImages[[n]]).TangentImages[[n]] TangentImages[[n]], {n, Range[Nimage]}];
   ForceNormalImages = -Table[GradImages[[n + 1]] - GradImages[[n + 1]].TangentImages[[n]] TangentImages[[n]], {n, Nimage}];
   StepImages = Table[LinearSolve[HessianImages[[n+1]], (ForceNormalImages[[n]] + ForceSpringImages[[n]] + ForceKinkImages[[n]])], {n, Nimage}];
   images = Join[{image0}, 0.1 StepImages + images[[2 ;; Nimage + 1]], {imageN}];
   Eimages = Eonsite /. Thread[vars -> #] & /@ images;
   (*Print[images//MatrixForm];
   Print[Norm[#]&/@BondImages];
   Print[ListPlot[{Range[0,Nimage+1],Eimages}\[Transpose],PlotRange\[Rule]Full,LabelingFunction\[Rule](Callout[#1[[2]],Automatic]&),PlotMarkers\[Rule]{Automatic,8},Frame\[Rule]True,GridLines\[Rule]Automatic,ImageSize\[Rule]250]]*), {i, Niter}];
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> {-900,900}, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]]
]

MPInverse[M_] := Inverse[M\[Transpose] M] M\[Transpose]

GetMEP1[Eonsite_, image0_, imageN_, Nimage_, Niter_] := Module[{vars, images, disp, i, n, Eimages, UVecImages, BondImages, TangentImages, HessianImages, GradImages, ForceImages, ForceSpringImages, ForceKinkImages, ForceNormalImages, StepImages, StepSpringImages, StepKinkImages, StepNormalImages, debugarray=<|"op"->{},"dispara"->{},"disperp"->{},"dis"->{},"step"->{}|>},
  vars = Variables@Eonsite;
  images = Table[image0 + n/(Nimage + 1) (imageN - image0), {n, 0, Nimage + 1}];
  Eimages = Eonsite /. (Thread[vars -> #] & /@ images);
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> Full, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]];
  Do[
     GradImages = Table[D[Eonsite, {vars, 1}] /. Thread[vars -> im], {im, images}];
     HessianImages = Table[Chop[D[Eonsite, {vars, 2}] /. Thread[vars -> im]], {im, images}];
     UVecImages = Normalize[#] & /@ Differences[images];
     BondImages = # & /@ Differences[images];
     TangentImages = Table[Normalize[UVecImages[[i]] + UVecImages[[i + 1]]], {i, Range[Nimage]}];
     StepSpringImages = Table[disp=(10^-3)*(BondImages[[n + 1]] - BondImages[[n]]); disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}];
     StepKinkImages = Table[disp=(10^-5)*(BondImages[[n + 1]] - BondImages[[n]]); disp - disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}];
     StepNormalImages = Table[disp=(If[Norm[#]>10,0,#]&/@(-MPInverse[HessianImages[[n+1]]].GradImages[[n + 1]])); disp - disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}];
     StepImages = StepSpringImages + StepKinkImages + StepNormalImages;
     images = Join[{image0}, StepImages + images[[2 ;; Nimage + 1]], {imageN}];
     Eimages = Eonsite /. Thread[vars -> #] & /@ images;
     AppendTo[debugarray["op"],images];
     AppendTo[debugarray["step"],Norm[#]&/@StepNormalImages];
     AppendTo[debugarray["dis"],Norm[#]&/@BondImages];
     AppendTo[debugarray["dispara"],Norm[#]&/@Table[disp=(BondImages[[n + 1]] - BondImages[[n]]); disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}]];
     AppendTo[debugarray["disperp"],Norm[#]&/@Table[disp=(BondImages[[n + 1]] - BondImages[[n]]); disp - disp.TangentImages[[n]] TangentImages[[n]], {n, Nimage}]], 
     {i, Niter}];
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> Full, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]];
  Return[debugarray]
]


GetMEP[Eonsite_, image0_, imageN_, Nimage_, Niter_] := Module[{vars, images, i, n, Eimages, UVecImages, BondImages, TangentImages, HessianImages, GradImages, ForceSpringImages, ForceKinkImages, ForceNormalImages, StepImages}, 
  vars = Variables@Eonsite;
  images = Table[image0 + n/(Nimage + 1) (imageN - image0), {n, 0, Nimage + 1}];
  Eimages = Eonsite /. (Thread[vars -> #] & /@ images);
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> Full, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]];
  Do[
   GradImages = Table[D[Eonsite, {vars, 1}] /. Thread[vars -> im], {im, images}];
   HessianImages = Table[Chop[D[Eonsite, {vars, 2}] /. Thread[vars -> im]], {im, images}];
   UVecImages = Normalize[#] & /@ Differences[images];
   BondImages = # & /@ Differences[images];
   TangentImages = Table[Normalize[UVecImages[[i]] + UVecImages[[i + 1]]], {i, Range[Nimage]}];
   ForceSpringImages = 1.0 Table[(BondImages[[n + 1]] - BondImages[[n]]).TangentImages[[n]] TangentImages[[n]], {n, Range[Nimage]}];
   ForceKinkImages = 0.1 Table[(BondImages[[n + 1]] - BondImages[[n]]) - (BondImages[[n + 1]] - BondImages[[n]]).TangentImages[[n]] TangentImages[[n]], {n, Range[Nimage]}];
   ForceNormalImages = -Table[GradImages[[n + 1]] - GradImages[[n + 1]].TangentImages[[n]] TangentImages[[n]], {n, Nimage}];
   StepImages = Table[LinearSolve[HessianImages[[n+1]], (ForceNormalImages[[n]] + ForceSpringImages[[n]] + ForceKinkImages[[n]])], {n, Nimage}];
   images = Join[{image0}, 0.1 StepImages + images[[2 ;; Nimage + 1]], {imageN}];
   Eimages = Eonsite /. Thread[vars -> #] & /@ images, 
   {i, Niter}];
  Print[ListPlot[{Range[0, Nimage + 1], Eimages}\[Transpose], PlotRange -> {-900,900}, LabelingFunction -> (Callout[#1[[2]], Automatic] &), PlotMarkers -> {Automatic, 8}, Frame -> True, GridLines -> Automatic, ImageSize -> 250]]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

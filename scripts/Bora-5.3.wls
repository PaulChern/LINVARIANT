#!/usr/bin/env wolframscript
Needs["INVARIANT`"];


TempCoeff={Subscript[C, 1111]->2630.720199808822`,Subscript[C, 1122]->787.9532029181885`,Subscript[C, 1212]->1039.2852937577954`,Subscript[\[Beta], 11]->-0.5044768267332222`,Subscript[\[Beta], 1111]->0.17459101764256235`,Subscript[\[Beta], 1122]->0.16707550949168`,Subscript[\[Gamma], p]->0.11360338810072582`,Subscript[Q, 1111]->0.0015617734347307362`,Subscript[Q, 1122]->0.002096435259444331`,Subscript[Q, 1212]->0.0045759294858504366`,Subscript[\[Gamma], pc]->-6.822450975792265`,Subscript[\[Alpha], xx]->-0.14161324961352054`,Subscript[\[Alpha], xxxx]->0.11750744312274372`,Subscript[\[Alpha], xxXX]->-0.19262613956826413`,Subscript[\[Alpha], xxyy]->-0.26144973824015766`,Subscript[\[Alpha], xxYY]->-0.1833140004264515`,Subscript[\[Gamma], x5]->-0.08265877106249529`,Subscript[\[Gamma], xc]->-5.765892474809912`,Subscript[\[Eta], a1b2p3]->7.453047068990745`,Subscript[\[Eta], c1c2p3]->-0.08069700048989199`,Subscript[\[Eta], a2b1p3]->7.408472757634671`,Subscript[\[Gamma], xp]->0.14902633780562038`};


Print["Setup simulations, mesh and discretization."];


param=Join[{\[CapitalDelta]1->2,\[CapitalDelta]2->2,\[CapitalDelta]3->2 ,Subscript[Gx, 1111]->1.0,Subscript[Gx, 1122]->0,Subscript[Gx, 1212]->1.0,Subscript[Gp, 1111]->0.1,Subscript[Gp, 1122]->0,Subscript[Gp, 1212]->0.1},TempCoeff];
Lx=ToExpression[$ScriptCommandLine[[3]]];
Ly=ToExpression[$ScriptCommandLine[[4]]];
Lz=ToExpression[$ScriptCommandLine[[5]]];
Print["Coefficients: ",param];


zbc=Dispatch[{Subscript[u, 3,x_,y_,Lz+1]:>Subscript[u, 3,x,y,1]}];
xbc=Dispatch[{Subscript[u, 3,Lx+1,y_,z_]:>Subscript[u, 3,1,y,z]}];


meshtype="square";
fDisc0=LandauInvariant["Discretized"];
fDisc = Sum[MeshMask[meshtype,{Lx,Ly,Lz},x,y,z]fDisc0, {x, Lx}, {y, Ly}, {z, Lz}]/.param//Expand;
Length[vars=Variables@fDisc];


dir0=Directory[];
case=$ScriptCommandLine[[2]];
NT=ToExpression[$ScriptCommandLine[[6]]];
dir=dir0<>"/Bora-5.3-T"<>case<>"/";
Print["case: "<>case<>"\n dimension: "<>ToString[Lx]<>"x"<>ToString[Ly]<>"x"<>ToString[Lz]];
Print["Start minimizationL ("<>ToString[NT]<>"):"];
If[DirectoryQ[dir],dir<>" exist, will over write!",CreateDirectory[dir]];


Efields=Table[-0.2{Cos[(2 Pi)/NT \[CapitalDelta]t],Cos[(2 Pi)/NT \[CapitalDelta]t],Cos[(2 Pi)/NT \[CapitalDelta]t]}*{0,0,1},{\[CapitalDelta]t,0,NT}];


{pdata,min}=MinimizeBoracite[case,{Lx,Ly,Lz},fDisc,100000,"square",Efields,100000,"TimeSeries"->True];


Export[dir<>"Bora-5.3-T"<>case<>".XPu.gif",Table[Grid[Partition[Flatten[Plotop[pdata[[n]],#,"FrameNum"->n,"RangeType"->"separated"]&/@{"X5","P","u"}],6],Alignment->Left,ItemSize->20],{n,NT+1}],"DisplayDurations"->1,"AnimationRepetitions"->\[Infinity]];

Export[dir<>"Bora-5.3-T"<>case<>".P.gif",Table[PlotPvector[Dispatch[min[[n]]],{Lx,Ly,Lz},"FrameNum"->n],{n,NT+1}],"DisplayDurations"->1,"AnimationRepetitions"->\[Infinity]];


Do[fname=dir<>"Bora-5.3-T"<>case<>"."<>ToString[n]; Export[fname<>".txt",{#1,#2}&@@@min[[n]],"Table"],{n,Length@min}]

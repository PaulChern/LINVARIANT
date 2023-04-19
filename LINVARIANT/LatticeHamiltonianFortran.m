BeginPackage["LINVARIANT`LatticeHamiltonianFortran`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`","LINVARIANT`Fortran`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
WriteLatticeModelF90         ::usage "WriteLatticeModelF90[F90dir, FuncName, FuncType, InvariantList, MMAVars, MMA2F90Vars]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)

WriteLatticeModelF90[F90dir_, FuncName_?StringQ, FuncType_?StringQ, InvariantList_?ListQ, MMAVars_, OptionsPattern[{"LineLimit" -> 500, "AllSites" -> False, "TylorOrder" -> 6}]] := Module[{FortranVar, FortranArrayVar, Invariants, HamList, ham, CoeffType, FortranVarList, ifield, i, inv, MMA, Fortran, f90, body, head, tail, CaseExprTable, var, ii, MMA2F90Vars, nn, ix, iy, iz, FortranVarLength, ExprFunc, ExprFuncBlock, CollectData},

  MMA2F90Vars = FortranVarSub[MMAVars];
  ExprFunc = {{}};
 
  FortranVarList = FuncType <> # & /@ First[Transpose[InvariantList]];
  HamList = Table[{inv[[1]], Expand@If[ListQ[inv[[2]]], inv[[2]].(ToExpression["Coeff"<>inv[[1]]][#] & /@ Range[Length@inv[[2]]]), inv[[2]]]}, {inv, InvariantList}];

  nn = Max@Abs@Append[Flatten[Cases[Variables[HamList], #] & /@ DeleteDuplicates[(Subscript[First[#], _, ix_, iy_, iz_] -> {ix, iy, iz}) & /@ Flatten[MMAVars[[1 ;; -2]]]]], 0];

  Which[
  (*Forces*)
    FuncType == "Forces", 
    Fortran = Flatten[
              Table[MMA = Flatten[Table[Table[FortranArrayVar = If[OptionValue["AllSites"], 
                                              FortranFuncArg[FuncType <> ham[[1]] <> ToString[i] <> ToString[ifield], {"x0", "y0", "z0"}], 
                                              FuncType <> ham[[1]] <> ToString[i] <> ToString[ifield]];
                         {FortranArrayVar, -D[Expand[ham[[2]]], MMAVars[[ifield, i]]]}, {i, Length@MMAVars[[ifield]]}], {ifield, Length@MMAVars}], 1];
                         {#1, Expr2Fortran[#2, MMA2F90Vars]} & @@@ MMA, {ham, HamList}], 1];

    ExprFuncBlock = FortranExprFunc[#1 <> "(nexpr)", FuncType, #2, nn, 2, "AllSites" -> OptionValue["AllSites"], "LineLimit"->OptionValue["LineLimit"]] & @@@ Fortran;
    ExprFunc = Flatten[Append[Flatten[#, 1], {"  "}] & /@ (Transpose[ExprFuncBlock][[2]]), 1];

    ii = 0;
    CollectData = Flatten[Table[Table[ii = ii + 1; 
                      Join[{{"    " <> FuncType <> #1 <> "(" <> ToString[i] <> "," <> ToString[ifield] <> ",:" <> ",:" <> ",:" <> ") = &"}}, 
                      FortranSumRiffle[Table[FuncType <> #1 <> ToString[i] <> ToString[ifield] <> "n" <> ToString[j], 
                      {j, Length[ExprFuncBlock[[ii, 2]]]}], 3, 4]], {i, Length@MMAVars[[ifield]]}], {ifield, Length@MMAVars}] & @@@ HamList, 3];

    body = Join[{#} & /@ Flatten[ExprFuncBlock\[Transpose][[1]]], {{" "}}, CollectData];
    FortranVar = {"Real*8", StringTrim[#], {"cgrid_a%n1+cgrid_b%n1", "cgrid_a%n2+cgrid_b%n2", "cgrid_a%n3+cgrid_b%n3"}} & /@ (StringSplit[Flatten[ExprFuncBlock\[Transpose][[1]]], "="]\[Transpose][[1]]);
    {head, tail} = HeadTailForces[FuncName, If[Length@FortranVarList == 1, First@FortranVarList, {FuncType, FortranVarList}], nn, "ExprFunc" -> False, "AllSites" -> OptionValue["AllSites"], "variables" -> FortranVar],

  (* Variation *)
    FuncType == "Variation",
    Fortran = Table[(*MMA = NOrderResponse[ham[[2]], #, OptionValue["TylorOrder"]] &/@ MMAVars;*)
                    MMA = Table[ParallelMap[NOrderResponse[#, var, OptionValue["TylorOrder"]]&, ham[[2]], 
                                            Method -> Automatic, 
                                            DistributedContexts -> {"LINVARIANT`LatticeHamiltonianFortran`Private`"}], {var, MMAVars}];
                    Expr2Fortran[#, MMA2F90Vars] &/@ MMA, {ham, HamList}];

    CaseExprTable = Transpose[Table[FortranExprFunc[#1 <> ToString[i] <> "(nexpr)", FuncType, #2[[i]], nn, 2, "AllSites" -> OptionValue["AllSites"], "LineLimit"->OptionValue["LineLimit"]], {i, Length@#2}] & @@@Transpose[{FortranVarList, Fortran}]];
    ExprFunc = Flatten[Flatten[Transpose[#][[2]], 2] & /@ CaseExprTable, 1];

    FortranVarLength = Table[Length[#[[1]]] & /@ (CaseExprTable\[Transpose][[i]]), {i, Length@FortranVarList}];
    FortranVar = Flatten[Table[Table[#1 <> ToString[i] <> "n" <> ToString[ii], {ii, #2[[i]]}], {i, Length[#2]}] & @@@ Transpose[{FortranVarList, FortranVarLength}]];
    body = FortranCaseBlock["idelta", {ToString[#] & /@ Range[Length[First@Fortran]], Table[Join[{#} &/@ Flatten[First@Transpose[CaseExprTable[[i]]], 1], {{"  "}, {"    " <> FuncType <> " = &"}}, FortranSumRiffle[StringSplit[Flatten[First@Transpose[CaseExprTable[[i]]], 1], "="]\[Transpose][[1]], 3, 4], {{"  "}}], {i, Length@CaseExprTable}]}\[Transpose], GetCaseDefaults[{{"write(*,*) \"mode out of range!\""}, {"call abort"}}, 1], 1];
    {head, tail} = HeadTailVariation[FuncName, {FuncType, FortranVar}, nn, "AllSites" -> OptionValue["AllSites"], "ExprFunc" -> False], 

    (*SiteEnergy*)
    FuncType == "SiteEnergy", 
    Fortran = Flatten[Table[FortranArrayVar = If[OptionValue["AllSites"], 
                                                 FortranFuncArg[FuncType <> ham[[1]], {"x0", "y0", "z0"}], 
                                                 FuncType <> ham[[1]]];
                            MMA = {{FortranArrayVar, ham[[2]]}};
                            {#1, Expr2Fortran[#2, MMA2F90Vars]} & @@@ MMA, {ham, HamList}], 1];

    ExprFuncBlock = FortranExprFunc[#1 <> "(nexpr)", FuncType, #2, nn, 2, "AllSites" -> OptionValue["AllSites"], "LineLimit"->OptionValue["LineLimit"]] & @@@ Fortran;
    ExprFunc = Flatten[Append[Flatten[#, 1], {"  "}] & /@ (Transpose[ExprFuncBlock][[2]]), 1];

    CollectData = Flatten[Table[
                    Join[{{If[OptionValue["AllSites"], 
                              "    " <> FuncType <> HamList[[i, 1]] <> " = &", 
                              "    " <> FuncType <> HamList[[i, 1]] <> " = &"]}}, 
                         FortranSumRiffle[Table[FuncType <> HamList[[i, 1]] <> "n" <> ToString[j], 
                    {j, Length[ExprFuncBlock[[i, 2]]]}], 3, 4]], {i, Length@HamList}], 1];

    body = Join[{#} & /@ Flatten[ExprFuncBlock\[Transpose][[1]]], {{" "}}, CollectData];
    FortranVar = {"Real*8", StringTrim[#]} & /@ (StringSplit[Flatten[ExprFuncBlock\[Transpose][[1]]], "="]\[Transpose][[1]]);

    {head, tail} = HeadTailSiteEnergy[FuncName, If[Length@FortranVarList == 1, First@FortranVarList, {FuncType, FortranVarList}], nn, "ExprFunc" -> False, "AllSites" -> OptionValue["AllSites"], "variables" -> FortranVar],

  (* HessianOnSite *)
    FuncType == "SiteHessian",
    Fortran = Flatten[Table[MMA = Flatten[Table[FortranArrayVar = FortranFuncArg[FuncType<>ham[[1]], {ToString[i], ToString[j]}];
                                    {FortranArrayVar, Expand@D[ham[[2]], {MMAVars[[i]],1}, {MMAVars[[j]],1}]}, {i, Length@MMAVars}, {j, Length@MMAVars}], 1];
                            {#1, Expr2Fortran[#2, MMA2F90Vars]} & @@@ MMA, {ham, HamList}], 1];
    body = Flatten[FortranExprBlock[#1, #2, 1, "LineLimit"->OptionValue["LineLimit"]] & @@@ Fortran, 1];
    {head, tail} = HeadTailHessianOnSite[FuncName, If[Length@FortranVarList == 1, First@FortranVarList, {FuncType, FortranVarList}], nn]];
  MMA2FORTRAN[F90dir <> "/" <> FuncName, Join[head, body, tail]];
  Return[ExprFunc]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]

BeginPackage["LINVARIANT`Fortran`", {"LINVARIANT`MathematicaPlus`"}]

FortranTools::usage = "The \"FortranTools\" package provides
a set of tools for reading and writing Fortran-formatted files.
It defines two main functions, FWrite and FReadList, that can handle
integer, real, logical and string data.";

FWrite::usage = "FWrite[channel, descriptor, expr1, expr2, ... ]
writes the expressions expri in sequence to the specified output channel
with descriptor specifying the Fortran output format. The descriptor
\"*\" defaults to \"A\" for strings, \"L2\" for booleans, \"I15\" for
integers, and \"E15.7\" for reals. No newline is added to the output.";

MMA2FORTRAN              ::usage = "MMA2FORTRAN[FortranName_String, Fout]"
EncodeExpression         ::usage = "EncodeExpression[expr]"
DecodeExpression         ::usage = "DecodeExpression[code]"
EncodeFortranPolynomial  ::usage = "EncodeFortranPolynomial[mma]"
HopingCode               ::usage = "HopingCode[x, i]"
HopingCodeBlock          ::usage = "HopingCodeBlock[nn]"
FortranExprBlock         ::usage = "FortranExprBlock[expr_, level_] "
FortranExprFunc          ::usage = "FortranExprFunc[expr_, level_] "
FortranSumRiffle         ::usage = "FortranSumRiffle[vars, ind, n]"
HeadTailHJijOnSite       ::usage = "HeadTailHJijOnSite[FunctionName, ReturnVar, FieldDim, nn]"
FortranCaseBlock         ::usage = "FortranCaseBlock[selector, cases, statements, level]"
FortranBodyArray         ::usage = "FortranBodyArray[selector, cases, statements, level]"
GetCaseDefaults          ::usage = "GetCaseDefaults[statement, level]"
GenerateCoefficientsFile ::usage = "GenerateCoefficientsFile[CoeffList, latt, dir]"
FortranParamModule       ::usage = "FortranParamModule[CoeffList]"
FortranReadCoeff         ::usage = "FortranReadCoeff[CoeffList]"
TimesFactor2Real         ::usage = "TimesFactor2Real[expr]"
StrainFromuF90           ::usage = "StrainFromu90[epsilon]"
LINVARIANTModule         ::usage = "LINVARIANTModule[FunctionList]"
FortranMain              ::usage = "FortranMain[]"
FortranConstants         ::usage = "FortranConstants[]"
Makefile                 ::usage = "Makefile[dir]"
Expr2Fortran             ::usage = "Expr2Fortran[expr, FortranVarSub]"
FortranFuncArg           ::usage = "FortranFuncArg[var, arg]"
FortranVarSub            ::usage = "FortranVarSub[field]"
HeadTailVariation        ::usage "HeadTailVariation[FunctionName, ReturnVar, FieldDim, nn]"
HeadTailForces           ::usage "HeadTailForces[FunctionName, ReturnVar, FieldDim, nn]"
HeadTailHessianOnSite    ::usage "HeadTailHessianOnSite[FunctionName, ReturnVar, FieldDim, nn]"
HeadTailSiteEnergy       ::usage "HeadTailSiteEnergy[FunctionName, ReturnVar, FieldDim, nn]"


FWriteArray::usage = "FWriteArray[channel, descriptor, matrix]
writes the two-dimensional array matrix to the specified output channel
with descriptor specifying the Fortran output format. The descriptor
\"*\" defaults to \"A\" for strings, \"L2\" for booleans, \"I15\" for
integers, and \"E15.7\" for reals. Rows are separated by newlines.";

ReadFString::usage = "ReadFString[string] reads one expression from string,
and returns the expression. ReadFString[string, type] reads one object of the
specified type. ReadFString[string, {type1, type2, ... }] reads a sequence of
objects of the specified types.";

ReadListFString::usage = "ReadListFString[string] reads all the remaining
expressions in string, and returns a list of them.
ReadListFString[string, type] reads objects of the specified type from string,
until the end of the string is reached. The list of objects read is returned.
ReadListFString[string, {type1, type2, ... }] reads objects with a sequence
of types, until the end of the string is reached.
ReadListFString[string, types, n] reads only the first n objects	of the
specified types.";

DropNonNumericElements::usage = "DropNonNumericElements[data] drops elements
that contain non-numeric elements in data.";

FToExpression::usage = "FToExpression[string] gives the Mathematica expression
obtained by interpreting string as Fortran-formatted output. Numbers are
converted to Integer or Real expressions, strings representing boolean values
to True or False. All other strings are returned unchanged."

FReadList::usage = "FReadList[\"file\"] reads all the remaining words
in a file interpreting them as Fortran-formatted output, and returns a list
of them. FReadList[\"file\", n] reads only the first n words.";

FReadLine::usage = "FReadLine[\"file\"] reads a single line from a file 
interpreting it as Fortran-formatted output, and returns a list of 
corresponding Mathematica expressions. FReadLine[\"file\", n] reads a total
of n lines.";


(* error messages for the exported objects *)

FWrite::invalid = "Invalid or incompatible format descriptor encountered.";
FWrite::unknown = "Unknown data type encountered.";
FWrite::expbig = 
  "Format error: The exponent `1` is not in the range (-99, 99).";


Begin["`Private`"]    (* begin the private context (implementation part) *)


(* definition of auxiliary functions and local (static) variables *)


(* parse format strings for descriptors like "E12.5". Returns a list of 
   three elements {w, d, e} where w represents the total width of the field, 
   d the number of digits and e the string denoting the exponent ("E" or 
   "D" or "F", i.e. no exponent). If no period is found, or if the width is 
   too small, a warning is issued. *)

parseRealFormat[s_String] := 
    Module[{s1, s2, p, r},
        s1 = StringTake[s, 1];
        s2 = StringDrop[s, 1];
        p = StringPosition[s2, "."];
        r = {Null, Null, "F"};
        If[Dimensions[p] == {1, 2}, (* exactly one . found *)
            p = p[[1, 1]];
            r = {ToExpression[StringTake[s2, p - 1]], 
                 ToExpression[StringDrop[s2, p]], s1} /. Null -> 0;
        ];
        If[(s1 == "F" && r[[1]] - r[[2]] > 2) || (r[[1]] - r[[2]] > 6),
            r,
            Message[FWrite::invalid];
            {Null, Null, s1}
        ]
    ]

(* parse Fortran format descriptors "Aw", "Lw", "Iw", "Fw.d", "Ew.d", "Dw.d", 
   and "*". Returns either the field width w, or a list {w, d, e} where e 
   contains the exponent string, or Null representing default values. 
   If the descriptor does not match the data type of the second argument, 
   a warning is issued. *)

parseFormat[s_String, i_] := 
    Module[{s1, s2, h},
        s1 = StringTake[s, 1];
        s2 = StringDrop[s, 1];
        h = Head[i];
        Which[
            s1 == "A" && h === String, ToExpression[s2],
            s1 == "L" && (i === True || i === False), ToExpression[s2],
            s1 == "I" && h === Integer, ToExpression[s2],
            (s1 == "F" || s1 == "E" || s1 == "D") && h === Real, 
                parseRealFormat[s],
            s1 == "*", Null,
            True, Message[FWrite::invalid]; Null
        ]
    ];

(* pad a string with blanks at the left. If the field width is too small 
   to output the string as a whole, just print *'s. *)

padString = "                                                                                                                                                                                                                                                               ";
errString = "***************************************************************************************************************************************************************************************************************************************************************";

padLeft[s_String, l_Integer] :=
    If[StringLength[s] > l,
        StringTake[errString, -l],
        StringTake[StringJoin[padString, s], -l]
    ];

(* NumberFormat-compatible routine: 
   Output mantissa and exponent in Fortran - like "E" format *)

fortranFormat[ m_, base_, exp_, s_String:"E"] :=
    Which[
        exp === "",
            SequenceForm[m, s, "+00"],
        StringTake[exp, 1] =!= "-" && StringLength[exp] === 1,
            SequenceForm[m, s, "+0", exp],
        StringTake[exp, 1] =!= "-" && StringLength[exp] === 2,
            SequenceForm[m, s, "+", exp],
        StringTake[exp, 1] === "-" && StringLength[exp] === 2,
            SequenceForm[m, s, StringInsert[exp, "0", 2]],
        StringTake[exp, 1] === "-" && StringLength[exp] === 3,
            SequenceForm[m, s, exp],
        True, (* Format error : exponent out of range *)
            Message[FortranTools::expbig, exp];
            $Failed
    ];

(* string conversion routines for data of type string, boolean, 
   integer and real. Every conversion routine expects the variable to be 
   printed and a format specification. If the format specification is Null, 
   it uses appropriate default values. *)

fwrite[s_String, Null] := s;
fwrite[s_String, l_Integer] := padLeft[s, l];

fwrite[False, Null] := padLeft["F", 2];
fwrite[True, Null] := padLeft["T", 2];
fwrite[False, l_Integer] := padLeft["F", l];
fwrite[True, l_Integer] := padLeft["T", l];

fwrite[i_Integer, Null] := fwrite[i, 15];
fwrite[i_Integer, l_Integer] := padLeft[ToString[i], l];

fwrite[r_Real, Null] := fwrite[r, {15, 7, "E"}];
fwrite[r_Real, {Null, _, exp___String}] := fwrite[r, {15, 7, exp}];
fwrite[r_Real, {l_Integer, d_Integer, "F"}] := padLeft[ToString[
    NumberForm[r, {l - 2, d}, NumberPadding -> {" ", "0"}, 
        ExponentFunction -> (Null &)]], 
    l];
fwrite[r_Real, {l_Integer, d_Integer, exp_String}] := padLeft[ToString[
    ScientificForm[r, {l - 6, d}, NumberPadding -> {" ", "0"}, 
        NumberFormat -> (fortranFormat[#1, #2, #3, exp] &)]], 
    l];

fwrite[___] := (Message[FWrite::unknown]; ""); (* default handler *)


(* Conversion routines: Convert strings containing Fortran expressions
   to Mathematica expressions. Numbers are converted to Integer or Real
   expressions, strings representing boolean values to True or False.
   Read::readn is turned off, because otherwise ReadListFString will warn
   about failed conversions. *)

Off[Read::readn];

fToNumber[s_String] :=
    (If[(Length[#2] == 1) && NumberQ[First[#2]], First[#2], #1])&[
        s, ReadListFString[s, Number] ];

(* Code rewritten as pure function for 20% faster execution compared to:
fToNumber[s_String] := 
    Module[{result},
        result = ReadListFString[s, Number];
        If[(Length[result] == 1) && NumberQ[First[result]],
            First[result],
            s
        ]
    ];
*)

(* The Fortran standard allows that real numbers are written without an 
   "E" when the exponent has more than two digits or even as default, 
   e.g. "1.234567-123". To catch this (abnormal) behaviour, fToNumber2 
   checks for a character sequence number - plus/minus - number and 
   tentatively inserts an "E": *)
fToNumber2[s_String] :=
	Module[{result = s, chars = Characters[s]},
		Do[ (* look for patterns like "1-1": *)
			If[ DigitQ[chars[[i-1]]] && DigitQ[chars[[i+1]]] && 
				MemberQ[{"+", "-"}, chars[[i]]],
				result = fToNumber[StringInsert[s, "E", i]];
				Break[]
			],
			{i, 2, Length[chars] - 1}
		];
		result
	];

fToBoolean[s_String] := 
    Module[{u},
        u = ToUpperCase[s];
        Switch[If[StringLength[u] > 0, StringTake[u, {1}], ""],
            "T", True,
            "F", False,
            ".", Switch[If[StringLength[u] > 1, StringTake[u, {2}], ""],
                    "T", True,
                    "F", False,
                    _, s
                ],
            _, s
        ]
    ];



(* definition of the exported functions *)

(* FWrite routine : 
   Write to a file as specified by Fortran format descriptors. FWrite is 
   recursively defined and terminates with a Null statement as soon as 
   all atomic or list arguments are processed. 
   Attention: All arguments must be compatible with the format descriptor 
   used. FWrite does not append any newlines! *)

FWrite[file_String, format_String, rest___] := 
    Module[{stream},
        stream = OpenWrite[file];
        FWrite[stream, format, rest];
        Close[stream];
    ];

FWrite[stream_OutputStream, format_String] := Null;
FWrite[stream_OutputStream, format_String, i_ /; AtomQ[i], rest___] :=
    (WriteString[stream, fwrite[i, parseFormat[format, i]]];
     FWrite[stream, format, rest]);
FWrite[stream_OutputStream, format_String, {list__}, rest___] := 
    FWrite[stream, format, list, rest];


(* FWriteArray routine : 
   Write an array to a file. Rows are separated by newlines. *)

FWriteArray[file_String, format_String, rest___] := 
    Module[{stream},
        stream = OpenWrite[file];
        FWriteArray[stream, format, rest];
        Close[stream];
    ];

FWriteArray[stream_OutputStream, format_String, m_ /; MatrixQ[m]] :=
    Scan[(FWrite[stream, format, #]; Write[stream]) &, m];


(* ReadFString/ReadListFString: Input data reading from string variable. *)

ReadFString[string_String, rest___] := 
    Module[{stream, result},
        stream = StringToStream[string];
        result = Read[stream, rest];
        Close[stream];
        result
    ];
ReadListFString[string_String, rest___] := 
    Module[{stream, result},
        stream = StringToStream[string];
        result = ReadList[stream, rest];
        Close[stream];
        result
    ];


(* DropNonNumericElements: Remove non-numeric elements from a list *)

Attributes[DropNonNumericElements] = {Listable};
DropNonNumericElements[e_ /; NumberQ[e]] := e;
(* If e is not numeric, drop it completely by inserting an empty sequence 
   instead *)
DropNonNumericElements[e_] := Sequence[];


(* FToExpression: Convert strings to Mathematica expressions by interpreting
   them as Fortran output. *)

Attributes[FToExpression] = {Listable};
FToExpression[s_String] := 
    Module[{r},
		r = fToNumber[s];
		If[Head[r] === String,
			r = fToNumber2[s];
			If[Head[r] === String,
				r = fToBoolean[s]
			]
		];
		r
    ];


(* FReadList: Arguments can either be number of words to be 
   read followed by options or only options. Caveat: Multiple spaces are
   counted as multiple word separators, so specifying n might result in
   fewer items than expected! *)

FReadList[file_, n___Integer, opts___] :=
    FToExpression[ReadList[file, Word, n, opts, RecordLists -> True]];

FReadLine[file_, opts___] := First[FReadLine[file, 1, opts]];
FReadLine[file_, n_Integer, opts___] := 
    FToExpression[ReadListFString[#, Word, opts]]& /@ 
        ReadList[file, Record, n, RecordLists -> False, opts];

MMA2FORTRAN[FortranName_String, Fout_, OptionsPattern["PageWidth"->72]] := Module[{LeftPad, FunctionBegin, FunctionEnd, stream, ss},
  LeftPad = "      ";
  stream = OpenWrite[FortranName <> ".f90", FormatType -> OutputForm, PageWidth -> OptionValue["PageWidth"]];
  FWriteArray[stream, "*", Fout];
  Close[stream];
]

EncodeExpression[expr_] := Module[{exprlist, atom, FindUnsubdividable},
  exprlist = If[MatchQ[expr, Plus[_, __]], Level[expr, {1}], {Level[expr, {0}]}];
  FindUnsubdividable[atom_] := If[MatchQ[atom, _List], FindUnsubdividable[#] & /@ atom, If[AtomQ[atom], atom, <|Head[atom] -> FindUnsubdividable[Level[atom, {1}]]|>]];
  Return[<|Plus -> FindUnsubdividable[exprlist]|>]
]

DecodeExpression[code_] := Module[{},
  Which[MatchQ[code, _Association], 
          First[Keys[code]] @@ DecodeExpression[First[Values[code]]], 
        MatchQ[code, _List], 
          DecodeExpression[#] & /@ code, 
        AtomQ[code], 
          code
        ]
]

EncodeFortranPolynomial[mma_] := Module[{},
  Which[
   MatchQ[mma, _Symbol], mma,
   MatchQ[mma, _Integer], mma,
   MatchQ[mma, _Real], {mma},
   MatchQ[mma, _List], EncodeFortranPolynomial[#] & /@ mma,
   MatchQ[mma, _Association],
   Which[
    Keys[mma] === {Plus}, 
    EncodeFortranPolynomial[#] & /@ First[Values[mma]],
    Keys[mma] === {Subscript}, First[Values[mma]],
    Keys[mma] === {Power}, 
    ConstantArray[EncodeFortranPolynomial[First[Values[mma]][[1]]], 
     First[Values[mma]][[2]]],
    Keys[mma] === {Times}, 
    EncodeFortranPolynomial[#] & /@ First[Values[mma]]]]
]

HopingCode[xx_, dd_] := Module[{d},
  d = If[NumberQ[dd], dd, If[AtomQ[dd], 0, First@Level[dd, {1}]]];
  ToExpression[StringReplace[xx <> ToString[d], {"-" -> "i"}]]
]

HopingCodeBlock[nn_, OptionsPattern["int"->1]] := Module[{Fint, n, ij, xyz, varxyz, s, vars, dvars, expr, cgrid},
 cgrid=<|"x"->"1","y"->"2","z"->"3"|>;
 Fint[n_] := StringRepeat["  ", n];
 varxyz = Table[StringJoin[Riffle[Flatten[Table[xyz <> i <> ToString[j], {i, {"", "i"}}, {xyz, {"x", "y", "z"}}]], ","]], {j, nn}];
 dvars = Table[{Fint[1] <> "Integer             :: " <> xyz}, {xyz, varxyz}];
 expr = Flatten[Table[s = If[i == "i", "-", "+"]; Table[{Fint[OptionValue["int"]] <> xyz <> i <> ToString[j] <> " = " <> "(" <> xyz <> "0" <> s <> ToString[j] <> ")" <> "-floor(real(" <> xyz <> "0" <> s <> ToString[j] <> "-1" <> ")/real(cgrid%n" <> cgrid[xyz] <> "))*cgrid%n" <> cgrid[xyz]}, {j, nn}, {xyz, {"x", "y", "z"}}], {i, {"", "i"}}], 2];
 Return[{varxyz, dvars, expr}]
]

FortranSumRiffle[vars_, ind_, n_] := Module[{Fint, x, out},
  Fint = StringRepeat["  ", ind];
  out = StringRiffle[#, {Fint <> "+", "+", " &"}] & /@ Partition[StringTrim[#] & /@ vars, UpTo[n]];
  out[[-1]] = StringDelete[out[[-1]], "&"];
  Return[Transpose@{out}]
]

FortranExprFunc[f_?StringQ, FuncType_, expr0_, nn_, level_, OptionsPattern[{"LineLimit" -> 2000, "AllSites" -> False}]] := Module[{Fint, i, n, xyz, term, expr, FoutList, Fout, head, tail, subfunc, test, HeadTail, fvar, euij}, 
  HeadTail = Which[FuncType === "Variation",     HeadTailVariation, 
                   FuncType === "Forces",        HeadTailForces, 
                   FuncType === "SiteEnergy",    HeadTailSiteEnergy,
                   FuncType === "HessianOnSite", HeadTailHessianOnSite];

  euij = If[FuncType === "Variation", True, False];

  Fint[n_] := StringRepeat["  ", n];
  expr = Partition[expr0, UpTo@OptionValue["LineLimit"]];
  FoutList = Table[term = ToString[FortranForm[atom]];
                   {If[StringTake[term, 1] != "-", 
                       Fint[level + 1] <> "+" <> StringDelete[term, Whitespace] <> " &", 
                       Fint[level + 1] <> StringDelete[term, Whitespace] <> " &"]}, {atom, #}] & /@ expr;

  Do[FoutList[[i]][[-1]][[1]] = StringDrop[FoutList[[i]][[-1]][[1]], -2], {i, Length@FoutList}];

  Fout = Table[subfunc = StringDelete[StringReplace[f, "(nexpr)" -> "n_" <> ToString[i]], RegularExpression["\\(.+\\)"]];
               fvar = StringDelete[f, RegularExpression["\\(\\w+\\)"]];
               {head, tail} = HeadTail[subfunc, StringDelete[fvar, RegularExpression["\\(.+\\)"]], nn, "euij" -> euij, "array" -> False, "ExprFunc" -> True, "AllSites" -> OptionValue["AllSites"]]; 
               {Fint[level] <> StringDelete[StringReplace[f, "(nexpr)" -> "n" <> ToString[i]], RegularExpression["\\(.+\\)"]] <> " = " <> StringJoin[StringSplit[head[[1, 1]]][[2 ;; -2]]], 
                Join[head, {{Fint[level] <> fvar <> " = &"}, {Fint[level + 1] <> fvar <> " &"}}, FoutList[[i]], {{Fint[level]}}, tail]}, {i, Length[FoutList]}];

  Return[Transpose@Fout]
]

FortranExprBlock[f_?StringQ, expr0_, nn_, level_, OptionsPattern[{"LineLimit" -> 100}]] := Module[{Fint, i, n, xyz, term, expr, FoutList, Fout},
  Fint[n_] := StringRepeat["  ", n];
  expr = DeleteCases[Join[Partition[expr0, OptionValue["LineLimit"]], {Complement[expr0, Flatten[Partition[expr0, OptionValue["LineLimit"]]]]}], {}];
  FoutList = Table[term = ToString[FortranForm[atom]];
                   {If[StringTake[term, 1] != "-", 
                       Fint[level + 1] <> "+" <> StringDelete[term, Whitespace] <> " &", 
                       Fint[level + 1] <> StringDelete[term, Whitespace] <> " &"]}, {atom, #}] & /@ expr;

  Do[FoutList[[i]][[-1]][[1]] = StringDrop[FoutList[[i]][[-1]][[1]], -2], {i, Length@FoutList}];

  Fout = Flatten[Join[{{Fint[level] <> f <> " = &"}, {Fint[level+1] <> f <> " &"}}, #, {{Fint[level]}}] & /@ FoutList, 1];
  Return[Fout]
]

HeadTailHJijOnSite[FunctionName_?StringQ, ReturnVar_?StringQ, nn_, OptionsPattern[{"AllSites" -> False, "ExprFunc" -> False, "array" -> True}]] := Module[{Fint, n, head, tail}, 
  Fint[n_] := StringRepeat["  ", n];
  head = Join[{{"Function " <> FunctionName <> "(x0, y0, z0, Fields) Result(" <> ReturnVar <> ")"}, 
               {Fint[1]}, 
               {Fint[1] <> "Implicit none"}, 
               {Fint[1] <> "Integer, Intent(in) :: x0, y0, z0"}, 
               {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"}, 
               {Fint[1] <> "Real*8              :: " <> ReturnVar}}, 
              HopingCodeBlock[nn][[2]], 
              {{Fint[1]}},
              HopingCodeBlock[nn][[3]]];
  tail = {{Fint[1]}, {"End Function " <> FunctionName}};
  Return[{head, tail}]
]

StrainFromuF90[epsilon_, nn_] := Module[{head, tail, Fint, n},
  Fint[n_] := StringRepeat["  ", n];

  head = Join[{{"Function GetHeterostructureStrain(Fields, xx, yy, zz, ncut) Result(euij)"},
               {Fint[1]},
               {Fint[1] <> "Use Parameters"},
               {Fint[1] <> "Implicit none"},
               {Fint[1] <> "Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"},
               {Fint[1] <> "Integer, Intent(in)    :: xx, yy, zz, ncut"},
               {Fint[1] <> "Integer                :: dx, dy, dz"},
               {Fint[1] <> "Integer                :: x0, y0, z0"},
               {Fint[1] <> "Real*8                 :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"}},
               HopingCodeBlock[1][[2]],
              {{Fint[1]},
               {Fint[1] <> "euij = 0.0d0"},
               {Fint[1]},
               {Fint[1] <> "do dx = -ncut, ncut"},
               {Fint[2] <> "do dy = -ncut, ncut"},
               {Fint[3] <> "do dz = -ncut, ncut"},
               {Fint[4] <> "x0 = (xx+dx)-floor(real(xx+dx-1)/real(cgrid%n1))*cgrid%n1"},
               {Fint[4] <> "y0 = (yy+dy)-floor(real(yy+dy-1)/real(cgrid%n2))*cgrid%n2"},
               {Fint[4] <> "z0 = (zz+dz)-floor(real(zz+dz-1)/real(cgrid%n3))*cgrid%n3"},
               {Fint[4]}},
               {Fint[3] <> First[#]} &/@ (HopingCodeBlock[1][[3]]),
              {{Fint[1]}}];

  tail = {{Fint[1]}, 
          {Fint[3] <> "end do"},
          {Fint[2] <> "end do"},
          {Fint[1] <> "end do"},
          {"End Function GetHeterostructureStrain"}};

  Return[Join[head, {StringRepeat["  ", 3] <> First[#]} &/@ epsilon, tail]]
]

GetCaseDefaults[statement_, level_] := Module[{Fint, n, default},
  Fint[n_] := StringRepeat["  ", n];
  default = {Fint[level+2] <> #} &/@ statement;
  Return[{"DEFAULT", default}]
]

FortranCaseBlock[selector_?StringQ, CasesStates_, default_, level_] := Module[{Fint, n, head, tail, body, cases, statements},
  Fint[n_] := StringRepeat["  ", n];
  {cases, statements} = Join[CasesStates, {default}]\[Transpose];
  head = {{Fint[level] <> "SELECT CASE (" <> selector <> ")"}};
  body = Flatten[Join[{{Fint[level] <> "CASE " <> Which[NumberQ[ToExpression[ToLowerCase[#1]]], "("<>#1<>")", #1==="DEFAULT", #1, True, "(\"" <> ToLowerCase[#1] <> "\")"]}}, #2] & @@@ ({cases, statements}\[Transpose]), 1];
  tail = {{Fint[level] <> "END SELECT"}};
  Return[Join[head, body, tail]]
]

FortranBodyArray[var_?StringQ, ExprArray_, level_] := Module[{len, Fint, n, i, body, tail},
  Fint[n_] := StringRepeat["  ", n];
  len = Length[ExprArray];
  body = Flatten[Table[FortranExprBlock[var <> ToString[i], ExprArray[[i]], level], {i, len}], 1];
  tail = {{Fint[level] <> var <> " = (/" <> StringJoin[Riffle[var <> ToString[#] &/@ Range[len], ","]] <> "/)"}};
  Return[Join[body, tail]]
]

GenerateCoefficientsFile[dir_, latt_, CoeffList_] := Module[{CoeffOut, c},
  CoeffOut = Flatten[Table[Join[{{c[[1]] <> "    :    ", Length[c[[2]]]}}, Partition[Flatten[c[[2]]], 1], {{""}}], {c, CoeffList}], 1];
  Export[dir<>"/database/"<>"Coefficients.dat", Flatten[#] & /@ Join[Join[{{"lattice    :    "}}, Map[ToString[NumberForm[#, {3, 16}]] &, #] & /@ latt, {{}}], CoeffOut]];
]

FortranParamModule[CoeffList_, vars_] := Module[{Fint, n, head, tail, body1, body2, c, Fout, FieldDim, FieldBinary, i, NumField},
  Fint[n_] := StringRepeat["  ", n];
  FieldDim = Max[Append[Length[#] &/@ vars, 3]];
  NumField = Length[vars] + 1;
  FieldBinary = Flatten@Append[Table[If[MemberQ[Cases[#, Subscript[__, n_] -> n], i], 1, 0], {i, FieldDim}] & /@ vars, {1, 1, 1}];
  head ={{"Module Parameters"},
         {Fint[1] <> "use Constants"}, 
         {Fint[1] <> "use common"}, 
         {Fint[1] <> "implicit none"}, 
         {Fint[1] <> "integer                 :: NODE_ME, NCPU, IERROR, IONODE"}, 
         {Fint[1] <> "Integer, parameter      :: stdin = 5"}, 
         {Fint[1] <> "Integer, parameter      :: stdout = 6"}, 
         {Fint[1] <> "Integer, parameter      :: ifileno = 55                          !< File handle number for input files"}, 
         {Fint[1] <> "Integer, parameter      :: ofileno = 66                          !< File handle number for output files"},
         {Fint[1] <> "Integer                 :: NumField = " <> ToString[NumField] <> ", FieldDim = " <> ToString[FieldDim]},
         {Fint[1] <> "Integer, dimension("<>ToString[FieldDim]<>","<>ToString[NumField]<>") :: FieldsBinary = " <> "Reshape((/" <> StringRiffle[ToString[#] & /@ FieldBinary, ", "] <> "/), (/"<>ToString[FieldDim]<>","<>ToString[NumField]<>"/))"},
         {Fint[1] <> "Integer                 :: OnSiteDim = " <> ToString[Total[FieldBinary]]},
         {Fint[1] <> "Real*8                  :: latt(3,3)"},
         {Fint[1] <> "Real*8, dimension (:,:,:,:,:), allocatable :: EwaldMat"},
         {Fint[1] <> "Real*8, dimension (:,:,:),     allocatable :: EwaldHessian"},
         {Fint[1]}};
  tail = {{Fint[1]}, { "End Module Parameters"}};
  body1 = Table[{Fint[1] <> "Real*8                  :: " <> c[[1]] <> "(" <> ToString[Length[c[[2]]]] <> ")"}, {c, CoeffList}];
  body2 = {{Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! system"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "type(s_unitcell) :: unitcell"},
           {Fint[1] <> "Real(dp)         :: symprec"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! parallelization"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "integer        :: nproc_k, nproc_ob, nproc_rgrid"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! DFT"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "integer        :: iperiodic"},
           {Fint[1] <> "integer        :: layout_multipole"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Geometry and composition"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "type(s_grid) :: cgrid"},
           {Fint[1] <> "type(s_grid) :: rgrid"},
           {Fint[1] <> "type(s_grid) :: kgrid"},
           {Fint[1] <> "Real(dp), allocatable  ::  kbz(:,:)"},
           {Fint[1] <> "Real(dp)     :: DWq(3,3), alat"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! External prob"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Real*8                :: EAmp(4)"},
           {Fint[1] <> "Real*8                :: EPhi(4)"},
           {Fint[1] <> "Real*8                :: GateField(4)"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! TB and Jij data"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Integer                      :: SpinDim, NumWann, OrbMul"},
           {Fint[1] <> "Integer                      :: NumWannSites, ContourNPoints(3)"},
           {Fint[1] <> "Integer                      :: supercell(3), k_mesh(3), fft_grid(3), Jij_R(3)"},
           {Fint[1] <> "Integer,allocatable          :: NumRpts(:), WannSiteInd(:), WannBlockInd(:,:), WannDist(:,:,:)"},
           {Fint[1] <> "real(dp)                     :: ContourMin, ContourMax, ContourHeight"},
           {Fint[1] <> "Integer                      :: potim"},
           {Fint[1] <> "complex(dp), allocatable     :: WannFunc(:,:,:,:,:,:)"},
           {Fint[1] <> "complex(dp), allocatable     :: WannBasis(:,:,:,:,:)"},
           {Fint[1] <> "real(dp), allocatable        :: redcoord(:,:,:)       !< Coordinates for Heisenberg exchange couplings"},
           {Fint[1] <> "real(dp), allocatable        :: jc(:,:,:,:,:)         !< Exchange couplings"},
           {Fint[1] <> "Integer,  allocatable        :: SiteOrbInfo(:,:)"},
           {Fint[1] <> "Integer                      :: JijSites(2)"},
           {Fint[1] <> "complex(dp), allocatable     :: HRmn(:,:,:,:), ContourPath(:)"},
           {Fint[1] <> "integer, allocatable         :: Ti0(:,:,:)"},
           {Fint[1] <> "Real(dp),allocatable         :: Efermi(:)"},
           {Fint[1] <> "logical                      :: WriteG"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Simulation paramters"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "character(len=40) :: RestartFields    !< restart file"},
           {Fint[1] <> "character(len=40) :: RestartVelocity  !< restart file"},
           {Fint[1] <> "character(len=10) :: NameSim          !< Name of simulation"},
           {Fint[1] <> "character(len=20) :: CoeffFile        !< Name of coefficient data file"},
           {Fint[1] <> "character(len=20) :: TrajectoryFile   !< Name of Config file"},
           {Fint[1] <> "character(len=1)  :: aunits           !< Atomic units to simulate model Hamiltonians (Y/N)"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Solvers"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "character(len=20) :: Solver                 !< Model solver"},
           {Fint[1] <> "character(len=20) :: EfieldType             !< Electric field type"},
           {Fint[1] <> "integer           :: ThermoSteps            !< Thermo up"},
           {Fint[1] <> "integer           :: CoolingSteps           !< Thermo up"},
           {Fint[1] <> "integer           :: NumSteps               !< Number of Monte Carlo steps"},
           {Fint[1] <> "integer           :: TapeRate               !< Number of steps in measurement phase"},
           {Fint[1] <> "integer           :: TrainRate              !< Number of steps in measurement phase"},
           {Fint[1] <> "Logical           :: DipoleQ                !< dipole dipole interaction"},
           {Fint[1] <> "Logical           :: EfieldQ                !< dipole dipole interaction"},
           {Fint[1] <> "Logical           :: TrainQ                 !< dipole dipole interaction"},
           {Fint[1] <> "Logical           :: CLAMPQ(6)              !< strain frozen"},
           {Fint[1] <> "Logical           :: FrozenQ(" <> ToString[NumField] <> ")             !< parameters frozen"},
           {Fint[1] <> "real*8            :: Temp                   !< Temperature"},
           {Fint[1] <> "real*8            :: Pressure               !< Pressure"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Markov Chain Monte Carlo"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Real*8  :: damp                   !< damp size, move scale"},
           {Fint[1] <> "Real*8  :: DampRatio              !< damp size, move scale"},
           {Fint[1] <> "Real*8  :: AcceptRatio            !< accept ratio"},
           {Fint[1] <> "integer :: seed                   !< random seed"},
           {Fint[1] <> "integer :: AvrgInterval           !< Sampling interval for averages"},
           {Fint[1] <> "integer :: BuffMcAvrg             !< Buffer size for averages"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Wang Landau Monte Carlo"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Real*8  :: wl_a                 !< reducing factor"},
           {Fint[1] <> "integer :: wl_rate              !< spiking rate"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Parallel Tempering Monte Carlo"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Real*8  :: ReplicaT0"},
           {Fint[1] <> "Real*8  :: ReplicaTN"},
           {Fint[1] <> "integer :: SwapRate               !< Number of sweeps for one swap"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Molecular Dynamics"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Real*8            :: deltaT"},
           {Fint[1] <> "character(len=20) :: ThermoState"},
           {Fint[1] <> "Logical           :: ReservoirQ(2)"},
           {Fint[1] <> "Integer           :: ReservoirRate"},
           {Fint[1] <> "Real*8            :: ReservoirRatio"},
           {Fint[1] <> "Real*8            :: ReservoirTau"},
           {Fint[1] <> "Real*8            :: ReservoirMass"},
           {Fint[1]},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "! Minimization"},
           {Fint[1] <> "!---------------------------------------------------------------------------------"},
           {Fint[1] <> "Character(len=20) :: OptAlgo"},
           {Fint[1] <> "Real*8            :: opt_tol"},
           {Fint[1] <> "Integer           :: optdim"},
           {Fint[1] <> "Integer           :: opt_iter=0"}};

  Fout = Join[head, body1, body2, tail];
  Return[Fout]
]

FortranReadCoeff[CoeffList_] := Module[{Fint, n, c, CaseBlock, default, head, tail, body, Fout},
  Fint[n_] := StringRepeat["  ", n];
  head = {{"Subroutine ReadCoefficients"}, 
          {Fint[1] <> "use Constants"}, 
          {Fint[1] <> "use Parameters"}, 
          {Fint[1] <> "use FileParser"}, 
          {Fint[1] <> "implicit none"}, 
          {Fint[1] <> "character(len=50)   :: keyword, cache "}, 
          {Fint[1] <> "integer             :: rd_len,i_err,i,j,i_errb, pos, Ndim"}, 
          {Fint[1] <> "Real(dp)            :: tmpdp"}, 
          {Fint[1] <> "logical             :: comment"}, 
          {Fint[1]}, 
          {Fint[1] <> "open(ifileno,file=CoeffFile)"}, 
          {Fint[1]}, 
          {Fint[1] <> "do"}, 
          {Fint[2] <> "10 continue"}, 
          {Fint[2] <> "keyword=\"\""}, 
          {Fint[2] <> "call bytereader(keyword,rd_len,ifileno,i_errb)"}, 
          {Fint[2] <> "call caps2small(keyword)"}, 
          {Fint[2] <> "comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&"}, 
          {Fint[2] <> "(scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&"}, {Fint[2] <> "(scan(trim(keyword),'!')==1))"}, 
          {Fint[2]}, 
          {Fint[2] <> "if (comment) then"}, 
          {Fint[3] <> "read(ifileno,*) cache"}, 
          {Fint[2] <> "else"}, 
          {Fint[2]}, 
          {Fint[3] <> "keyword=trim(keyword)"}, 
          {Fint[3]}};
  tail = {{Fint[2] <> "end if"}, 
          {Fint[2]}, 
          {Fint[2] <> "if (i_errb==20) goto 20"}, 
          {Fint[2] <> "if (i_errb==10) goto 10"}, 
          {Fint[2]}, 
          {Fint[1] <> "end do"}, 
          {Fint[1]}, 
          {Fint[1] <> "20 continue"}, 
          {Fint[1]}, 
          {Fint[1] <> "mass = mass*mpme"}, 
          {Fint[1]}, 
          {"End Subroutine ReadCoefficients"}};
  CaseBlock = Join[
       {{"lattice", 
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "read(ifileno, *, iostat=i_err)((latt(i,j),i=1,3),j=1,3)"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3] <> "alat = (norm2(latt(:,1)) + norm2(latt(:,2)) + norm2(latt(:,3)))/3"},
          {Fint[3]}}},
        {"dwq",
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "read(ifileno,*,iostat=i_err) ((DWq(j,i), j=1,3), i=1,3)"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3]}}},
        {"efield",
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "do i = 1, 4"},
          {Fint[3] <> "read(ifileno,*,iostat=i_err) EAmp(i), EPhi(i), GateField(i)"},
          {Fint[3] <> "end do"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3]}}},
        {"contour",
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "pos = scan(cache, ':')"},
          {Fint[3] <> "cache = trim(cache(pos+1:))"},
          {Fint[3] <> "read(cache, *, iostat=i_err) Ndim"},
          {Fint[3] <> "read(ifileno,*,iostat=i_err) ContourMin, ContourMax, ContourHeight"},
          {Fint[3] <> "read(ifileno,*,iostat=i_err) ContourNPoints(:)"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3]}}},
        {"wanniersite",
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "pos = scan(cache, ':')"},
          {Fint[3] <> "cache = trim(cache(pos+1:))"},
          {Fint[3] <> "read(cache, *, iostat=i_err) NumWannSites"},
          {Fint[3] <> "Allocate(SiteOrbInfo(3,NumWannSites))"},
          {Fint[3] <> "do i = 1, NumWannSites"},
          {Fint[3] <> "Read(ifileno,*,iostat=i_err) SiteOrbInfo(1,i), SiteOrbInfo(2,i), SiteOrbInfo(3,i)"},
          {Fint[3] <> "end do"},
          {Fint[3] <> "do i = 1, 2"},
          {Fint[3] <> "  do j = 1, NumWannSites"},
          {Fint[3] <> "    If(SiteOrbInfo(3,j).gt.0) then"},
          {Fint[3] <> "      JijSites(i) = j"},
          {Fint[3] <> "      SiteOrbInfo(3,j) = SiteOrbInfo(3,j) - 1"},
          {Fint[3] <> "      exit"},
          {Fint[3] <> "    end if"},
          {Fint[3] <> "  end do"},
          {Fint[3] <> "end do"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3]}}},
        {"poscar",
         {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
          {Fint[3] <> "pos = scan(cache, ':')"},
          {Fint[3] <> "cache = trim(cache(pos+1:))"},
          {Fint[3] <> "Read(cache, *, iostat=i_err) unitcell%nion"},
          {Fint[3] <> "Allocate(unitcell%site(3,unitcell%nion))"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err)"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err) tmpdp"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err) ((unitcell%latt_a(i,j), j=1,3), i=1,3)"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err)"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err)"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err)"},
          {Fint[3] <> "unitcell%latt_a = tmpdp*unitcell%latt_a"},
          {Fint[3] <> "Read(ifileno, *, iostat=i_err) ((unitcell%site(j,i), j=1,3), i=1,unitcell%nion)"},
          {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[3]}}}},
     Table[{c[[1]],
     {{Fint[3] <> "read(ifileno, '(A)', iostat=i_err) cache"},
      {Fint[3] <> "pos = scan(cache, ':')"},
      {Fint[3] <> "cache = trim(cache(pos+1:))"},
      {Fint[3] <> "read(cache, *, iostat=i_err) Ndim"},
      {Fint[3] <> "do i = 1, Ndim"},
      {Fint[3] <> "read(ifileno,*,iostat=i_err) " <> c[[1]] <> "(i)"},
      {Fint[3] <> "end do"},
      {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
      {Fint[3]}}}, {c, CoeffList}]];
  default = GetCaseDefaults[{{"if(len(trim(keyword))>0) then"}, Fint[1] <> {"read(ifileno,*)"}, {"end if"}}, 1];
  body = FortranCaseBlock["keyword", CaseBlock, default, 1];
  Fout = Join[head, body, tail];
  Return[Fout]
]

TimesFactor2Real[expr_] := Which[MatchQ[expr, Times[_, __]], Construct @@ Flatten[{Times, Map[(If[NumberQ[#] && #+1 != 0, N[#], #] &), Level[expr, {1}]]}], IntegerQ[expr], N[expr], True, expr]

FortranMain[] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"Program main"},
          {Fint[1] <> "Use LINVARIANT"},
          {Fint[1] <> "Use Parameters"},
          {Fint[1] <> "Use Inputs"},
          {Fint[1] <> "Use Outputs"},
          {Fint[1] <> "Use mpi"},
          {Fint[1]},
          {Fint[1] <> "Implicit none"},
          {Fint[1]},
          {Fint[1] <> "integer :: PROCESS_RANK, SIZE_CLUSTER, IERROR, MESSAGE_ITEM"},
          {Fint[1]},
          {Fint[1] <> "Integer           :: i, j, istep, ireplica, ix, iy, iz"},
          {Fint[1] <> "Integer           :: N1, N2, N3"},
          {Fint[1] <> "Real*8, dimension (:,:,:,:,:), allocatable :: Fields"},
          {Fint[1] <> "Real*8, dimension (:,:,:,:,:), allocatable :: dFieldsdt"},
          {Fint[1] <> "Real*8, dimension (:,:,:,:,:), allocatable :: EwaldField"},
          {Fint[1] <> "Real*8, dimension (:), allocatable         :: udamp"},
          {Fint[1] <> "Real*8, dimension (:), allocatable         :: TempList"},
          {Fint[1] <> "Integer,dimension (:,:), allocatable       :: Replicas"},
          {Fint[1] <> "Real*8            :: gm, etadamp"},
          {Fint[1] <> "Real*4            :: tstart, tend"},
          {Fint[1] <> "Real*8            :: e0ij(3,3), de0ijdt(3,3), TK, dene"},
          {Fint[1] <> "Real*8            :: Etot, Epot, Ekin"},
          {Fint[1] <> "character(len=10) :: NameEwald"},
          {Fint[1] <> "character(len=10) :: FileCode"},
          {Fint[1] <> "logical           :: file_exists"},
          {Fint[1]},
          {Fint[1]},
          {Fint[1] <> "call MPI_INIT(IERROR)"},
          {Fint[1] <> "call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE_CLUSTER, IERROR)"},
          {Fint[1] <> "call MPI_COMM_RANK(MPI_COMM_WORLD, PROCESS_RANK, IERROR)"},
          {Fint[1]},
          {Fint[1] <> "if(SIZE_CLUSTER.gt.1) then"},
          {Fint[2] <> "FileCode = int2str(PROCESS_RANK)"},
          {Fint[1] <> "else"},
          {Fint[2] <> "FileCode = int2str(0)"},
          {Fint[1] <> "end if"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Read control parameters !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1]},
          {Fint[1] <> "INQUIRE(FILE=\"LINVARIANT.inp\", EXIST=file_exists)"},
          {Fint[1] <> "if (file_exists) then"},
          {Fint[2] <> "open(ifileno,file='LINVARIANT.inp')"},
          {Fint[2] <> "call ReadParameters(ifileno)"},
          {Fint[2] <> "close(ifileno)"},
          {Fint[1] <> "else"},
          {Fint[2] <> "call set_input_defaults()"},
          {Fint[1] <> "end if"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Read Coefficients !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> ""},
          {Fint[1] <> "open(ifileno,file=CoeffFile)"},
          {Fint[1] <> "Call ReadCoefficients(ifileno)"},
          {Fint[1] <> "close(ifileno)"},
          {Fint[1] <> ""},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Read or Write Ewald Matrix !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1]},
          {Fint[1] <> "if(PROCESS_RANK == 0) then"},
          {Fint[2] <> "INQUIRE(FILE=\"EwaldMat.dat\", EXIST=file_exists)"},
          {Fint[2] <> "if (file_exists) then"},
          {Fint[3] <> "write(*,*) 'EwaldMat.dat exist, checking ...'"},
          {Fint[3] <> "open(ifileno,file='EwaldMat.dat',form='unformatted',status='old')"},
          {Fint[3] <> "read(ifileno) NameEwald, N1, N2, N3"},
          {Fint[3] <> "close(ifileno)"},
          {Fint[3] <> "if(N1.ne.cgrid%n1 .or. N2.ne.cgrid%n2 .or. N3.ne.cgrid%n3 .or. NameEwald .ne. NameSim) then"},
          {Fint[4] <> "write(*,*) 'Dipole file mismatch, current simulation: ', NameSim, cgrid%n1, cgrid%n2, cgrid%n3"},
          {Fint[4] <> "write(*,*) 'Ewald file: ', NameEwald, N1, N2, N3"},
          {Fint[4] <> "write(*,*) 'Now, regenerate Ewald Matrix!!!'"},
          {Fint[4] <> "call EwaldMatrix(latt)"},
          {Fint[3] <> "else"},
          {Fint[4] <> "write(*,*) 'Dipole file matches and will be used.'"},
          {Fint[3] <> "end if"},
          {Fint[3] <> "close(ifileno)"},
          {Fint[2] <> "else"},
          {Fint[3] <> "write(*,*) 'No EwaldMat.dat is found, generating ...'"},
          {Fint[3] <> "call EwaldMatrix(latt)"},
          {Fint[2] <> "end if"},
          {Fint[1] <> "end if"},
          {Fint[1]},
          {Fint[1] <> "allocate(EwaldMat(3,3,cgrid%n1,cgrid%n2,cgrid%n3))"},
          {Fint[1] <> "open(ifileno,file='EwaldMat.dat',form='unformatted',status='old')"},
          {Fint[1] <> "read(ifileno) NameEwald, N1, N2, N3"},
          {Fint[1] <> "read(ifileno) (((((EwaldMat(j,i,ix,iy,iz),j=1,3),i=1,3),ix=1,cgrid%n1),iy=1,cgrid%n2),iz=1,cgrid%n3)"},
          {Fint[1] <> "close(ifileno)"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Initialize Fields !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1]},
          {Fint[1] <> "allocate(Fields(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3))"},
          {Fint[1] <> "Call GetInitConfig(Fields, e0ij, RestartFields, FileCode)"},
          {Fint[1]},
          {Fint[1] <> "allocate(dFieldsdt(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3))"},
          {Fint[1] <> "Call GetInitConfig(dFieldsdt, de0ijdt, RestartVelocity, FileCode)"},
          {Fint[1]},
          {Fint[1] <> "allocate(EwaldField(3,NumField,cgrid%n1,cgrid%n2,cgrid%n3))"},
          {Fint[1] <> "Call GetEwaldField(Fields, EwaldField)"},
          {Fint[1]},
          {Fint[1] <> "Call RemoveGridDrifts(Fields)"},
          {Fint[1] <> "Call RemoveGridDrifts(dFieldsdt)"},
          {Fint[1]},
          {Fint[1] <> "gm = 0.0D0"},
          {Fint[1] <> "Tk = Thermometer(dFieldsdt, de0ijdt)"},
          {Fint[1] <> "dFieldsdt = Sqrt(Temp/Tk)*dFieldsdt"},
          {Fint[1] <> "de0ijdt = Sqrt(Temp/Tk)*de0ijdt"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Start Simulations !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "if(PROCESS_RANK == 0) then"},
          {Fint[2] <> "Call fmkdir(trim(Solver)//\".out\")"},
          {Fint[1] <> "end if"},
          {Fint[1] <> "Call get_walltime(tend)"},
          {Fint[1] <> "allocate(udamp(NumField))"},
          {Fint[1] <> "allocate(TempList(SIZE_CLUSTER))"},
          {Fint[1] <> "allocate(Replicas(SIZE_CLUSTER,5))"},
          {Fint[1] <> "if(trim(Solver).eq.\"MD\") then"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "! Initialize MD !"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "open(ifileno, &"},
          {Fint[2] <> "file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileCode)//'.dat', &"},
          {Fint[2] <> "form='unformatted',status='unknown')"},
          {Fint[2] <> "!!!!!!!!!!!"},
          {Fint[2] <> "! MD LOOP !"},
          {Fint[2] <> "!!!!!!!!!!!"},
          {Fint[2] <> "do istep = 1, NumSteps"},
          {Fint[3] <> "Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, Temp, gm)"},
          {Fint[3] <> "if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then"},
          {Fint[4] <> "tstart = tend"},
          {Fint[4] <> "Call get_walltime(tend)"},
          {Fint[4] <> "Epot = GetEtot(Fields, e0ij)"},
          {Fint[4] <> "Ekin = GetEkin(dFieldsdt, de0ijdt)"},
          {Fint[4] <> "Etot = Epot + Ekin"},
          {Fint[4] <> "Tk = Thermometer(dFieldsdt, de0ijdt)"},
          {Fint[4] <> "write(*,*) trim(Solver)//" step: ", istep, tend - tstart, \"(second)\", Tk, gm"},
          {Fint[4] <> "write(*,*) \"Eenergy (tot, pot, kin): \", Etot, Epot, Ekin"},
          {Fint[4] <> "write(*,*) e0ij(1,1), e0ij(2,2), e0ij(3,3)"},
          {Fint[4] <> "write(*,*) e0ij(2,3), e0ij(1,3), e0ij(1,2)"},
          {Fint[4] <> "Call WriteBinary(ifileno, Fields, dFieldsdt, e0ij, de0ijdt)"},
          {Fint[3] <> "end if"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "close(ifileno)"},
          {Fint[2] <> "Call WriteFinal('FinalConfig-'//trim(FileCode)//'.dat', Fields, e0ij)"},
          {Fint[2] <> "Call WriteFinal('FinalVelocity-'//trim(FileCode)//'.dat', dFieldsdt, de0ijdt)"},
          {Fint[2] <> "Call BinaryToData(TrajectoryFile, 0)"},
          {Fint[2] <> "Call GetObservables(\"Observables\", 1, 0)"},
          {Fint[1] <> "else if(trim(Solver).eq.\"MCMC\") then"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "! Initialize MCMC !"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "open(ifileno, &"},
          {Fint[2] <> "file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileCode)//'.dat', &"},
          {Fint[2] <> "form='unformatted',status='unknown')"},
          {Fint[2] <> "dene = 0.0D0"},
          {Fint[2] <> "do i = 1, NumField"},
          {Fint[3] <> "udamp(i) = damp"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "etadamp = damp"},
          {Fint[2] <> "!!!!!!!!!!!!!"},
          {Fint[2] <> "! MCMC LOOP !"},
          {Fint[2] <> "!!!!!!!!!!!!!"},
          {Fint[2] <> "Call random_seed()"},
          {Fint[2] <> "do istep = 1, NumSteps"},
          {Fint[3] <> "Call MCMCStep(istep, Fields, EwaldField, e0ij, Temp, udamp, etadamp, dene)"},
          {Fint[3] <> "if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then"},
          {Fint[4] <> "tstart = tend"},
          {Fint[4] <> "Call get_walltime(tend)"},
          {Fint[4] <> "Etot = GetEtot(Fields, e0ij)"},
          {Fint[4] <> "write(*,*) trim(Solver)//\" step: \", istep, tend - tstart, \"(second)\", dene"},
          {Fint[4] <> "write(*,*) \"Eenergy (pot): \", Etot"},
          {Fint[4] <> "write(*,*) e0ij(1,1), e0ij(2,2), e0ij(3,3)"},
          {Fint[4] <> "write(*,*) e0ij(2,3), e0ij(1,3), e0ij(1,2)"},
          {Fint[4] <> "Call WriteBinary(ifileno, Fields, dFieldsdt, e0ij, de0ijdt)"},
          {Fint[3] <> "end if"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "close(ifileno)"},
          {Fint[2] <> "Call WriteFinal('FinalConfig-'//trim(FileCode)//'.dat', Fields, e0ij)"},
          {Fint[2] <> "Call WriteFinal('FinalVelocity-'//trim(FileCode)//'.dat', dFieldsdt, de0ijdt)"},
          {Fint[2] <> "Call BinaryToData(TrajectoryFile, 0)"},
          {Fint[2] <> "Call GetObservables(\"Observables\", 1, 0)"},
          {Fint[1] <> "else if((trim(Solver).eq.\"PTMC\").or.(trim(Solver).eq.\"PTMD\")) then"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "! Initialize PT !"},
          {Fint[2] <> "!!!!!!!!!!!!!!!!!"},
          {Fint[2] <> "do ireplica = 1, SIZE_CLUSTER"},
          {Fint[3] <> "TempList(ireplica) = ReplicaT0 + (ireplica-1)*(ReplicaTN-ReplicaT0)/(SIZE_CLUSTER-1)"},
          {Fint[3] <> "Replicas(ireplica,:) = (/ireplica, ireplica, 0, 0, 0/)"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "de0ijdt = Sqrt(TempList(PROCESS_RANK+1)/Temp)*de0ijdt"},
          {Fint[2] <> "dFieldsdt = Sqrt(TempList(PROCESS_RANK+1)/Temp)*dFieldsdt"},
          {Fint[2] <> "if(PROCESS_RANK == 0) then"},
          {Fint[3] <> "open(11,file=trim(Solver)//'.out/'//'REPLICAS.dat',status='unknown')"},
          {Fint[3] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"E15.6)\") (TempList(i), i=1,SIZE_CLUSTER)"},
          {Fint[2] <> "end if"},
          {Fint[2] <> "open(ifileno+100*PROCESS_RANK,&"},
          {Fint[2] <> "file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileCode)//'.dat',&"},
          {Fint[2] <> "form='unformatted',status='unknown')"},
          {Fint[2] <> "Replicas(1,3) = 1"},
          {Fint[2] <> "Replicas(SIZE_CLUSTER,3) = -1"},
          {Fint[2] <> "dene = 0.0D0"},
          {Fint[2] <> "do i = 1, NumField"},
          {Fint[3] <> "udamp(i) = damp"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "etadamp = damp"},
          {Fint[2] <> "!!!!!!!!!!"},
          {Fint[2] <> "! PT LOOP !"},
          {Fint[2] <> "!!!!!!!!!!!"},
          {Fint[2] <> "Call random_seed()"},
          {Fint[2] <> "if(PROCESS_RANK == 0) then"},
          {Fint[3] <> "write(*,*) trim(Solver)//\" solver runing in \", SIZE_CLUSTER, \" mpi processors.\""},
          {Fint[3] <> "write(*,*) \"Temperauture set: \", TempList"},
          {Fint[2] <> "end if"},
          {Fint[2]},
          {Fint[2] <> "if(CoolingSteps.gt.0) then"},
          {Fint[3] <> "do ireplica = SIZE_CLUSTER, 1, -1"},
          {Fint[4] <> "if(trim(Solver).eq.\"PTMC\") then"},
          {Fint[5] <> "do istep = 1, CoolingSteps"},
          {Fint[6] <> "Call MCMCStep(istep, Fields, EwaldField, e0ij, TempList(PROCESS_RANK+1), udamp, etadamp, dene)"},
          {Fint[5] <> "end do"},
          {Fint[4] <> "else"},
          {Fint[5] <> "do istep = 1, CoolingSteps"},
          {Fint[6] <> "Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, TempList(PROCESS_RANK+1), gm)"},
          {Fint[5] <> "end do"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "if(ireplica.eq.PROCESS_RANK+1) then"},
          {Fint[5] <> "Call WriteFinal('FinalConfig-'//trim(FileCode)//'.dat', Fields, e0ij)"},
          {Fint[5] <> "Call WriteFinal('FinalVelocity-'//trim(FileCode)//'.dat', dFieldsdt, de0ijdt)"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, SIZE_CLUSTER, Replicas, PROCESS_RANK+1, ireplica)"},
          {Fint[4] <> "Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)"},
          {Fint[3] <> "end do"},
          {Fint[3] <> "RestartFields = 'FinalConfig'"},
          {Fint[3] <> "Call GetInitConfig(Fields, e0ij, RestartFields, FileCode)"},
          {Fint[2] <> "end if"},
          {Fint[2] <> "CoolingSteps = 0"},
          {Fint[2] <> "gm = 0.0D0"},
          {Fint[2]},
          {Fint[2] <> "if(trim(Solver).eq.\"PTMC\") then"},
          {Fint[3] <> "do istep = 1, NumSteps"},
          {Fint[4] <> "Call MCMCStep(istep, Fields, EwaldField, e0ij, TempList(PROCESS_RANK+1), udamp, etadamp, dene)"},
          {Fint[4] <> "if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then"},
          {Fint[5] <> "if(PROCESS_RANK == 0) then"},
          {Fint[6] <> "tstart = tend"},
          {Fint[6] <> "Call get_walltime(tend)"},
          {Fint[6] <> "write(*,*) trim(Solver)//\" step: \", istep, tend - tstart, \"(second)\""},
          {Fint[5] <> "end if"},
          {Fint[5] <> "Call WriteBinary(ifileno+100*PROCESS_RANK, Fields, dFieldsdt, e0ij, de0ijdt)"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "if((mod(istep,SwapRate).eq.0).and.(istep.lt.ThermoSteps)) then"},
          {Fint[5] <> "Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, SIZE_CLUSTER, Replicas, PROCESS_RANK+1, istep)"},
          {Fint[5] <> "if(PROCESS_RANK == 0) then"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"E15.6)\") (TempList(i), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,1), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,2), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,3), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,4), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,5), i=1,SIZE_CLUSTER)"},
          {Fint[5] <> "end if"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "else"},
          {Fint[3] <> "do istep = 1, NumSteps"},
          {Fint[4] <> "Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, TempList(PROCESS_RANK+1), gm)"},
          {Fint[4] <> "Epot = GetEtot(Fields, e0ij)"},
          {Fint[4] <> "if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then"},
          {Fint[5] <> "if(PROCESS_RANK == 0) then"},
          {Fint[6] <> "tstart = tend"},
          {Fint[6] <> "Call get_walltime(tend)"},
          {Fint[6] <> "write(*,*) trim(Solver)//\" step: \", istep, tend - tstart, \"(second)\""},
          {Fint[5] <> "end if"},
          {Fint[5] <> "Tk = Thermometer(dFieldsdt, de0ijdt)"},
          {Fint[5] <> "do i = 1, SIZE_CLUSTER"},
          {Fint[6] <> "if(i == PROCESS_RANK+1) then"},
          {Fint[7] <> "write(*,'(I5,A1,I10,F10.4,A1,F10.4,A12,F10.6)') i, \"@\", istep, Tk, \"/\", TempList(i), \"(K)  Epot: \", Epot"},
          {Fint[6] <> "end if"},
          {Fint[6] <> "call MPI_BARRIER(MPI_COMM_WORLD, IERROR)"},
          {Fint[5] <> "end do"},
          {Fint[5] <> "Call WriteBinary(ifileno+100*PROCESS_RANK, Fields, dFieldsdt, e0ij, de0ijdt)"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "if((mod(istep,SwapRate).eq.0).and.(istep.lt.ThermoSteps)) then"},
          {Fint[5] <> "Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, SIZE_CLUSTER, Replicas, PROCESS_RANK+1, istep)"},
          {Fint[5] <> "if(PROCESS_RANK == 0) then"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"E15.6)\") (TempList(i), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,1), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,2), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,3), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,4), i=1,SIZE_CLUSTER)"},
          {Fint[6] <> "write(11, \"(\"//trim(int2str(SIZE_CLUSTER))//\"I8)\") (Replicas(i,5), i=1,SIZE_CLUSTER)"},
          {Fint[5] <> "end if"},
          {Fint[4] <> "end if"},
          {Fint[4] <> "Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "end if"},
          {Fint[2] <> "write(*,*) \"PT process: \", PROCESS_RANK, \"of\", SIZE_CLUSTER, \"done!\""},
          {Fint[2] <> "close(ifileno+100*PROCESS_RANK)"},
          {Fint[2] <> "Call WriteFinal('FinalConfig-'//trim(FileCode)//'.dat', Fields, e0ij)"},
          {Fint[2] <> "Call WriteFinal('FinalVelocity-'//trim(FileCode)//'.dat', dFieldsdt, de0ijdt)"},
          {Fint[2] <> "Call BinaryToData(TrajectoryFile, PROCESS_RANK)"},
          {Fint[2] <> "Call GetObservables(\"Observables\", 1, PROCESS_RANK)"},
          {Fint[2] <> "if(PROCESS_RANK == 0) then"},
          {Fint[3] <> "close(11)"},
          {Fint[2] <> "end if"},
          {Fint[2] <> "call MPI_FINALIZE(IERROR)"},
          {Fint[1] <> "else"},
          {Fint[2] <> "write(*,*) trim(Solver)//\" is not implemented yet!\""},
          {Fint[2] <> "call abort"},
          {Fint[1] <> "end if"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!"},
          {Fint[1] <> "! End !"},
          {Fint[1] <> "!!!!!!!"},
          {Fint[1] <> "deallocate(udamp)"},
          {Fint[1] <> "deallocate(TempList)"},
          {Fint[1] <> "deallocate(Replicas)"},
          {Fint[1] <> "deallocate(Fields)"},
          {Fint[1] <> "deallocate(EwaldMat)"},
          {Fint[1] <> "deallocate(EwaldField)"},
          {Fint[1]},
          {Fint[1] <> "End Program main"}};
  Return[Fout]
]

FortranConstants[] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"!-------------------------------------------------------------------------------"},
          {"! MODULE: Constants"},
          {"!> Physical constants"},
          {"!> The magnitude of the parameters has been chosen so that it mtaches"},
          {"!> the values given in the NIST database with the CODATA recommended 2014 values"},
          {"!> https://physics.nist.gov/cuu/Constants/"},
          {"!-------------------------------------------------------------------------------"},
          {"module Constants"},
          {Fint[1] <> "implicit none"},
          {Fint[1] <> "Integer, parameter    :: snglprec = selected_real_kind(6, 37)  !< define precision for single reals"},
          {Fint[1] <> "Integer, parameter    :: dblprec = selected_real_kind(15, 307) !< define precision for double reals"},
          {Fint[1] <> "!.. Scalar parameters"},
          {Fint[1] <> "real(dblprec),parameter :: pi = 3.141592653589793_dblprec"},
          {Fint[1] <> "real(dblprec) :: gama         = 1.760859644d11     ! s^(-1)*T^(-1)"},
          {Fint[1] <> "real(dblprec) :: k_bolt       = 1.38064852d-23     ! J/K"},
          {Fint[1] <> "real(dblprec) :: k_bolt_ev    = 8.6173303d-5       ! eV/K"},
          {Fint[1] <> "real(dblprec) :: mub          = 5.7883818012d-5    ! eV/T"},
          {Fint[1] <> "real(dblprec) :: mu0          = 1.2566370614d-6    ! N/A^2"},
          {Fint[1] <> "real(dblprec) :: mry          = 2.179872325d-21    ! J"},
          {Fint[1] <> "real(dblprec) :: hbar_mev     = 6.582119514d-13    ! meV*s"},
          {Fint[1] <> "real(dblprec) :: hbar_ev      = 6.582119514d-16    ! eV*s"},
          {Fint[1] <> "real(dblprec) :: Joule_ev     = 6.241509126d18     ! eV"},
          {Fint[1] <> "real(dblprec) :: Hartree      = 27.211386245       ! eV"},
          {Fint[1] <> "real(dblprec) :: ry_ev        = 13.605693009_dblprec     ! eV"},
          {Fint[1] <> "real(dblprec) :: hbar         = 1.054571800e-34    ! J*s"},
          {Fint[1] <> "real(dblprec) :: ev           = 1.6021766208d-19   ! J"},
          {Fint[1] <> "real(dblprec) :: amu          = 1.660539040d-27    ! kg"},
          {Fint[1] <> "real(dblprec) :: ame          = 9.1093837015d-31   ! kg"},
          {Fint[1] <> "real(dblprec) :: mpme         = 1.83615267343d3    ! mp/me ratio"},
          {Fint[1] <> "real(dblprec) :: angstrom     = 1.0d-10            ! m"},
          {Fint[1] <> "real(dblprec) :: bohr         = 0.52917721067d-10  ! m"},
          {Fint[1] <> "real(dblprec) :: time_fs      = 2.4188843265857d-2 ! fs"},
          {Fint[1] <> "real(dblprec) :: g_e_abs      = 2.00231930436182"},
          {"end module Constants"}};
  Return[Fout]
]

Makefile[dir_, prefix_] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{".SUFFIXES: .f .f90 .F90 .SUFFIXES .prj"},
          {" "},
          {"FC = mpif90 -fopenmp -O3"},
          {"FC = mpiifort -fopenmp -O3"},
          {" "},
          {"PREFIX = " <> prefix},
          {"exe = LINVARIANT.x"},
          {" "},
          {"SRC_DIR   = $(PREFIX)/src"},
          {"BUILD_DIR = $(PREFIX)/build"},
          {"BIN_DIR   = $(PREFIX)/bin"},
          {" "},
          {"#SRCS := $(shell find $(SRC_DIR) -name \"*.f90\")"},
          {"SRCS = Constants.f90 Parameters.f90 FileParser.f90 \\"},
          {"       Inputs.f90 \\"},
          {"       LINVARIANT.f90 \\"},
          {"       EwaldMatrix.f90 \\"},
          {"       Outputs.f90 \\"},
          {"       mc.f90 md.f90 pt.f90 \\"},
          {"       main.f90"},
          {" "},
          {"OBJS = $(foreach obj, $(SRCS:.f90=.o), $(BUILD_DIR)/$(obj))"},
          {" "},
          {"#"},
          {"$(exe):   $(OBJS)"},
          {"	$(FC) $(LDFLAGS) -o $(exe) $(OBJS) $(LIBS) -O2"},
          {"	chmod +x $(exe); cp $(exe) $(BIN_DIR); rm -f *.f90 *.mod $(OBJS)"},
          {"#"},
          {"#"},
          {"%.o: %.f90"},
          {"	$(FC) -c $(FFLAGS) $<"},
          {"#"},
          {"%.f90:"},
          {"	cp $(shell find $(SRC_DIR) -name \"*.f90\") $(BUILD_DIR)"},
          {" "},
          {"install:     $(exe)"},
          {"	cp $(exe) $(BIN_DIR)"},
          {"	cd $(BIN_DIR); chmod +x $(exe)"},
          {"	    echo \"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\""},
          {"	    echo \"To decouple the dipole-dipole between different type of fields, modify GetSiteEnergyEwald.f90, GetEwaldForces.f90, GetEwaldField.f90, UpdateEwaldField.f90\""},
          {"	    echo \"To decouple microscopic and macroscopic strain, modify GetDeltaHEps.f90, GetSiteENergyGeps.f90, GetForcesEps.f90\""},
          {"	    echo \"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\""},
          {" "},
          {"clean: "},
          {"	rm -f *.f90 *.mod $(OBJS) $(exe)"},
          {"#"}};
  Export[dir, Fout, "Table", "TextDelimiters"->None];
]

Expr2Fortran[expr_, FortranVarSub_] := Module[{},
  TimesFactor2Real[#] & /@ If[MatchQ[expr, Plus[_, __]], Level[expr, {1}], Level[expr, {0}]] /. FortranVarSub
]

FortranFuncArg[var_, arg_] := Module[{},
  ToString@FortranForm@ToExpression[var <> "[" <> StringRiffle[arg, ","] <> "]"]
]

FortranVarSub[MMAVars_] := Module[{var, ix, iy, iz},
  var = DeleteCases[DeleteDuplicates[GetSubscriptInfo[#][[1]] & /@ Flatten[MMAVars]], ToExpression["\[Epsilon]0"]];
  Join[Table[(Subscript @@ Join[{var[[i]]}, ToExpression["{a_,ix_,iy_,iz_}"]] :> Fields[a, ii, HopingCode["x", ix], HopingCode["y", iy], HopingCode["z", iz]]) /. ii -> i, {i, Length@var}],
       Table[(Subscript @@ Join[{var[[i]]}, ToExpression["{a_,ix_,iy_,iz_,0}"]] :> Fields[a, ii, HopingCode["x", ix], HopingCode["y", iy], HopingCode["z", iz]]) /. ii -> i, {i, Length@var}],
       Subscript @@ Join[{ToExpression["\[Delta]" <> ToString[#]]}, ToExpression["{a_,0,0,0}"]] :> delta[a] & /@ var,
       {Subscript[ToExpression["\[Epsilon]"], ToExpression["i_"], ToExpression["j_"]] :> eij[i, j]},
       {Subscript[ToExpression["\[Epsilon]0"], ToExpression["i_"], ToExpression["j_"]] :> e0ij[i, j]},
       {Subscript[ToExpression["\[Epsilon]0"], ToExpression["i_"], ToExpression["j_"], 0] :> e0ij[i, j]},
       {Subscript[ToExpression["\[Epsilon]u"], ToExpression["i_"], ToExpression["j_"], ToExpression["ix_"], ToExpression["iy_"], ToExpression["iz_"]] :> euij[i, j, HopingCode["x", ix], HopingCode["y", iy], HopingCode["z", iz]]},
       {Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 1] :> delta[1],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 2, 2] :> delta[2],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 3, 3] :> delta[3],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 2, 3] :> delta[4],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 3] :> delta[5],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 2] :> delta[6]}] /. {a -> ToExpression["a"], ix -> ToExpression["ix"], iy -> ToExpression["iy"], iz -> ToExpression["iz"], x -> ToExpression["x"], y -> ToExpression["y"], z -> ToExpression["z"], delta -> ToExpression["delta"],  Fields -> ToExpression["Fields"], eij -> ToExpression["eij"], e0ij -> ToExpression["e0ij"], euij -> ToExpression["euij"], i -> ToExpression["i"], j -> ToExpression["j"]}
]

HeadTailVariation[FunctionName_?StringQ, ReturnVar_, nn_, OptionsPattern[{"AllSites" -> False, "ExprFunc" -> False, "euij" -> False, "array" -> True}]] := Module[{Fint, n, head, tail, res, reslist},

  Fint[n_] := StringRepeat["  ", n];

  {res, reslist} = Which[StringQ[ReturnVar], {ReturnVar, {ReturnVar}},
                         ListQ[ReturnVar] && Length[ReturnVar] == 2, {First[ReturnVar], ReturnVar[[2]]},
                         True, Print["ReturnVar either a string or a list with 2 elements"; Abort[]]];

  head = Join[{{"Function " <> FunctionName <> "(x0, y0, z0, Fields, e0ij,"<>If[OptionValue["euij"], " euij, ", " "]<>"idelta, delta) Result(" <> res <> ")"},
               {Fint[1]},
               {Fint[1] <> "Implicit none"},
               {Fint[1] <> "Integer, Intent(in) :: x0, y0, z0, idelta"},
               {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"},
               {Fint[1] <> "Real*8,  Intent(in) :: delta(FieldDim)"},
               {Fint[1] <> "Real*8,  Intent(in) :: e0ij(3,3)"},
               If[OptionValue["euij"],
               {Fint[1] <> "Real*8,  Intent(in) :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"},
               {Fint[1] <> "Real*8              :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"}],
               {Fint[1] <> "Real*8              :: eij(3,3)"},
               {Fint[1] <> "Integer             :: ncut = " <> ToString[nn]}},
               {Fint[1] <> "Real*8              :: " <> #} & /@ DeleteDuplicates[Join[{res}, reslist]],
              If[OptionValue["euij"], Join[HopingCodeBlock[nn][[2]], {{Fint[1]}}], {{Fint[1]}}],
              If[OptionValue["euij"], {{Fint[1]}}, {{Fint[1] <> "euij = GetHeterostructureStrain(Fields, x0, y0, z0, ncut)"}, {Fint[1]}}],
              If[OptionValue["euij"], Join[HopingCodeBlock[nn][[3]], {{Fint[1]}}], {{Fint[1]}}],
               {Fint[1] <> # <> " = 0.0D0"} & /@ DeleteDuplicates[Join[{res}, reslist]],
               {{Fint[1]}}];

  tail = {{Fint[1]}, {"End Function " <> FunctionName}, {Fint[1]}};
  Return[{head, tail}]
]

HeadTailForces[FunctionName_?StringQ, ReturnVar_, nn_, OptionsPattern[{"AllSites" -> True, "ExprFunc" -> False, "array" -> True, "euij" -> False, "variables"->{}}]] := Module[{Fint, n, head, tail, res, reslist, func},

  Fint[n_] := StringRepeat["  ", n];
  {res, reslist} = Which[StringQ[ReturnVar], {ReturnVar, {ReturnVar}}, ListQ[ReturnVar] && Length[ReturnVar] == 2, {First[ReturnVar], ReturnVar[[2]]}, True, Print["ReturnVar either a string or a list with 2 elements"; Abort[]]];

  func = If[OptionValue["AllSites"], {{"Function " <> FunctionName <> "(Fields, e0ij) Result(" <> res <> ")"}}, {{"Function " <> FunctionName <> "(x0, y0, z0, Fields, e0ij) Result(" <> res <> ")"}}];

  head = Join[
      func,
      {{Fint[1]},
       {Fint[1] <> "Implicit none"},
       {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8,  Intent(in) :: e0ij(3,3)"},
       If[OptionValue["euij"],
       {Fint[1] <> "Real*8,  Intent(in) :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8              :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"}],
       {Fint[1] <> If[OptionValue["AllSites"], "Integer             :: ", "Integer, Intent(in) :: "] <> "x0, y0, z0"},
       {Fint[1] <> "Real*8              :: eij(3,3)"},
       {Fint[1] <> "Integer             :: ncut = " <> ToString[nn]}},
       {Fint[1] <> "Real*8              :: " <> # <> If[OptionValue["array"], If[OptionValue["AllSites"], "(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)", "(Max(FieldDim, 6), NumField+1)"], If[OptionValue["AllSites"], "(cgrid%n1, cgrid%n2, cgrid%n3)", " "]]} & /@ DeleteDuplicates[Join[{res}, reslist]],
       {Fint[1] <> #1 <> "              :: " <> #2 <> "(" <> StringRiffle[#3, ","] <> ")"} &@@@ OptionValue["variables"],
      {{Fint[1]}},
      HopingCodeBlock[nn][[2]],
      {{Fint[1]}},
      {Fint[1] <> # <> " = 0.0D0"} & /@ DeleteDuplicates[Join[{res}, reslist]],
      If[OptionValue["ExprFunc"],
         If[OptionValue["AllSites"],
           {{Fint[1] <> "!$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(" <> StringJoin[Riffle[Join[{"eij,euij,x0","y0","z0"},HopingCodeBlock[nn][[1]]], ","]] <> ")"},
            {Fint[1] <> "!$OMP    DO COLLAPSE(3)"},
            {Fint[1] <> "do z0 = 1, cgrid%n3"},
            {Fint[1] <> "do y0 = 1, cgrid%n2"},
            {Fint[1] <> "do x0 = 1, cgrid%n1"},
            {Fint[2]}, 
            {Fint[1] <> "euij = GetHeterostructureStrain(Fields, x0, y0, z0, ncut)"},
            {Fint[2]}}, 
           {{Fint[2]},
             {Fint[1] <> "euij = GetHeterostructureStrain(Fields, x0, y0, z0, ncut)"},
             {Fint[2]}}], {{Fint[1]}}],
      {{Fint[1]}}];
  tail = Join[
    If[OptionValue["ExprFunc"],
      If[OptionValue["AllSites"],
        {{Fint[1] <> "end do"},
         {Fint[1] <> "end do"},
         {Fint[1] <> "end do"},
         {Fint[1] <> "!$OMP    END DO"},
         {Fint[1] <> "!$OMP    END PARALLEL"}}, {}],
        {{Fint[1]}}],
    If[ListQ[ReturnVar], {{Fint[1] <> res <> "=" <> StringJoin[Riffle[reslist, "+"]]}}, {}],
    {{Fint[1]},
     {"End Function " <> FunctionName}}];
  Return[{head, tail}]]

HeadTailHessianOnSite[FunctionName_?StringQ, ReturnVar_, nn_, OptionsPattern[{"euij" -> False, "array" -> True, "ExprFunc" -> False, "AllSites" -> False}]] :=
 Module[{Fint, n, head, tail, res, reslist, func},
  Fint[n_] := StringRepeat["  ", n];
  {res, reslist} = Which[StringQ[ReturnVar], {ReturnVar, {ReturnVar}}, ListQ[ReturnVar] && Length[ReturnVar] == 2, {First[ReturnVar], ReturnVar[[2]]}, True, Print["ReturnVar either a string or a list with 2 elements"; Abort[]]];
  func = {{"Function " <> FunctionName <> "(x0, y0, z0, Fields, e0ij) Result(" <> res <> ")"}};
  head = Join[
      func,
      {{Fint[1]},
       {Fint[1] <> "use Parameters"},
       {Fint[1] <> "use Inputs"}},
      {{Fint[1]},
       {Fint[1] <> "Implicit none"},
       {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8,  Intent(in) :: e0ij(3,3)"},
       {Fint[1] <> "Integer, Intent(in) :: x0, y0, z0"},
       If[OptionValue["euij"],
       {Fint[1] <> "Real*8,  Intent(in) :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8              :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"}],
       {Fint[1] <> "Real*8              :: eij(3,3)"},
       {Fint[1] <> "Integer             :: ncut = " <> ToString[nn]}},
      {Fint[1] <> "Real*8              :: " <> # <> "(OnSiteDim, OnSiteDim)"} & /@ DeleteDuplicates[Join[{res}, reslist]],
      {{Fint[1]}},
      HopingCodeBlock[nn][[2]],
      {{Fint[1]}},
      {Fint[1] <> # <> " = 0.0D0"} & /@ DeleteDuplicates[Join[{res}, reslist]],
      {{Fint[1] <> "euij = GetHeterostructureStrain(Fields, x0, y0, z0, ncut)"},
       {Fint[1]}},
      HopingCodeBlock[nn, "int" -> 1][[3]],
      {{Fint[1]}}];
  tail = Join[
    If[ListQ[ReturnVar], {{Fint[1] <> res <> "=" <> StringJoin[Riffle[reslist, "+"]]}}, {}],
    {{Fint[1]},
     {"End Function " <> FunctionName}}];
  Return[{head, tail}]
]

HeadTailSiteEnergy[FunctionName_?StringQ, ReturnVar_, nn_, OptionsPattern[{"AllSites" -> False, "ExprFunc" -> False, "array" -> True, "euij" -> False, "variables"->{}}]] := Module[{Fint, n, head, tail, res, reslist, func},
  Fint[n_] := StringRepeat["  ", n];
  {res, reslist} = Which[StringQ[ReturnVar], {ReturnVar, {ReturnVar}}, ListQ[ReturnVar] && Length[ReturnVar] == 2, {First[ReturnVar], ReturnVar[[2]]}, True, Print["ReturnVar either a string or a list with 2 elements"; Abort[]]];
  func = If[OptionValue["AllSites"], {{"Function " <> FunctionName <> "(Fields, e0ij) Result(" <> res <> ")"}}, {{"Function " <> FunctionName <> "(x0, y0, z0, Fields, e0ij) Result(" <> res <> ")"}}];
  head = Join[
      func,
      {{Fint[1]},
       {Fint[1] <> "Implicit none"},
       {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8,  Intent(in) :: e0ij(3,3)"},
       If[OptionValue["euij"],
       {Fint[1] <> "Real*8,  Intent(in) :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"},
       {Fint[1] <> "Real*8              :: euij(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)"}],
       {Fint[1] <> If[OptionValue["AllSites"], "Integer             :: ", "Integer, Intent(in) :: "] <> "x0, y0, z0"},
       {Fint[1] <> "Real*8              :: eij(3,3)"},
       {Fint[1] <> "Integer             :: ncut = " <> ToString[nn]}},
       {Fint[1] <> "Real*8              :: " <> # <> If[OptionValue["AllSites"], "(cgrid%n1, cgrid%n2, cgrid%n3)", ""]} & /@ DeleteDuplicates[Join[{res}, reslist]],
       {Fint[1] <> #1 <> "              :: " <> #2} &@@@ OptionValue["variables"],
      {{Fint[1]}},
      HopingCodeBlock[nn][[2]],
      {{Fint[1]}},
      {Fint[1] <> # <> " = 0.0D0"} & /@ DeleteDuplicates[Join[{res}, reslist]],
      {{Fint[1]}},
      If[OptionValue["AllSites"],
         {{Fint[1] <> "!$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(" <> StringJoin[Riffle[Join[{"eij,euij,x0","y0","z0"},HopingCodeBlock[nn][[1]]], ","]] <> ")"},
          {Fint[1] <> "!$OMP    DO COLLAPSE(3)"},
          {Fint[1] <> "do z0 = 1, cgrid%n3"},
          {Fint[1] <> "do y0 = 1, cgrid%n2"},
          {Fint[1] <> "do x0 = 1, cgrid%n1"}}, {}],
      {{Fint[1]},
       {Fint[1] <> "euij = GetHeterostructureStrain(Fields, x0, y0, z0, ncut)"},
       {Fint[1]}},
      HopingCodeBlock[nn, "int" -> If[OptionValue["AllSites"], 2, 1]][[3]],
      {{Fint[1]}}];
  tail = Join[
    If[OptionValue["AllSites"],
      {{Fint[1] <> "end do"},
       {Fint[1] <> "end do"},
       {Fint[1] <> "end do"},
       {Fint[1] <> "!$OMP    END DO"},
       {Fint[1] <> "!$OMP    END PARALLEL"}}, {}],
    If[ListQ[ReturnVar], {{Fint[1] <> res <> "=" <> StringJoin[Riffle[reslist, "+"]]}}, {}],
    {{Fint[1]},
     {"End Function " <> FunctionName}}];
  Return[{head, tail}]
]

End[ ]         (* end the private context *)

(* protect exported symbols *)

Protect[
    FWrite, FWriteArray,
    ReadFString, ReadListFString,
    DropNonNumericElements,
    FToExpression,
    FReadList, FReadLine
    ];

EndPackage[ ]  (* end the package context *)

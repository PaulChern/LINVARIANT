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
FortranParsers           ::usage = "FortranParsers[]"
FortranInputs            ::usage = "FortranInputs[]"
FortranOutputs           ::usage = "FortranOutputs[]"
Makefile                 ::usage = "Makefile[dir]"
Expr2Fortran             ::usage = "Expr2Fortran[expr, FortranVarSub]"
FortranVarStr            ::usage = "FortranVarStr[var, arg]"
FortranVarSub            ::usage = "FortranVarSub[field]"

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

HopingCode[x_, i_] := ToExpression[StringReplace[ToString[x] <> ToString[i], "-" -> "i"]]

HopingCodeBlock[nn_, OptionsPattern["int"->1]] := Module[{Fint, n, ij, xyz, varxyz, s, vars, dvars, expr},
 Fint[n_] := StringRepeat["  ", n];
 varxyz = Table[StringJoin[Riffle[Flatten[Table[xyz <> i <> ToString[j], {i, {"", "i"}}, {xyz, {"x", "y", "z"}}]], ","]], {j, nn}];
 dvars = Table[{Fint[1] <> "Integer             :: " <> xyz}, {xyz, varxyz}];
 expr = Flatten[Table[s = If[i == "i", "-", "+"]; Table[{Fint[OptionValue["int"]] <> xyz <> i <> ToString[j] <> " = " <> "(" <> xyz <> "0" <> s <> ToString[j] <> ")" <> "-floor(real(" <> xyz <> "0" <> s <> ToString[j] <> "-1" <> ")/real(NGrid" <> xyz <> "))*NGrid" <> xyz}, {j, nn}, {xyz, {"x", "y", "z"}}], {i, {"", "i"}}], 2];
 Return[{varxyz, dvars, expr}]
]

FortranExprBlock[f_?StringQ, expr0_, level_, OptionsPattern[{"LineLimit" -> 300}]] := Module[{Fint, i, n, xyz, term, expr, FoutList, Fout},
  Fint[n_] := StringRepeat["  ", n];
  expr = DeleteCases[Join[Partition[expr0, OptionValue["LineLimit"]], {Complement[expr0, Flatten[Partition[expr0, OptionValue["LineLimit"]]]]}], {}];
  FoutList = Table[term = ToString[FortranForm[atom]];
                   {If[StringTake[term, 1] != "-", Fint[level + 1] <> "+" <> StringDelete[term, Whitespace] <> " &", Fint[level + 1] <> StringDelete[term, Whitespace] <> " &"]}, {atom, #}] & /@ expr;
  Do[FoutList[[i]][[-1]][[1]] = StringDrop[FoutList[[i]][[-1]][[1]], -2], {i, Length@FoutList}];
  Fout = Flatten[Join[{{Fint[level] <> f <> " = &"}, {Fint[level+1] <> f <> " &"}}, #, {{Fint[level]}}] & /@ FoutList, 1];
  Return[Fout]
]

HeadTailHJijOnSite[FunctionName_?StringQ, ReturnVar_?StringQ, nn_] := Module[{Fint, n, head, tail}, 
  Fint[n_] := StringRepeat["  ", n];
  head = Join[{{"Function " <> FunctionName <> "(x0, y0, z0, Fields) Result(" <> ReturnVar <> ")"}, 
               {Fint[1]}, 
               {Fint[1] <> "Implicit none"}, 
               {Fint[1] <> "Integer, Intent(in) :: x0, y0, z0"}, 
               {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
               {Fint[1] <> "Real*8              :: " <> ReturnVar}}, 
              HopingCodeBlock[nn][[2]], 
              {{Fint[1]}},
              HopingCodeBlock[nn][[3]]];
  tail = {{Fint[1]}, {"End Function " <> FunctionName}};
  Return[{head, tail}]
]

StrainFromuF90[epsilon_] := Module[{head, tail, Fint, n},
  Fint[n_] := StringRepeat["  ", n];
  head = Join[{{"Function StrainFromu(x0, y0, z0, Fields) Result(euij)"},
               {Fint[1]},
               {Fint[1] <> "Use Parameters"},
               {Fint[1] <> "Implicit none"},
               {Fint[1] <> "Integer, Intent(in) :: x0, y0, z0"}, 
               {Fint[1] <> "Real*8,  Intent(in) :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
               {Fint[1] <> "Real*8              :: euij(3,3)"}},
              HopingCodeBlock[1][[2]],
              {{Fint[1]}},
              HopingCodeBlock[1][[3]],
              {{Fint[1]},
               {Fint[1] <> "euij = 0.0d0"},
               {Fint[1]}}];
  tail = {{Fint[1]}, {"End Function StrainFromu"}};
  Return[Join[head, epsilon, tail]]
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

GenerateCoefficientsFile[CoeffList_, latt_, dir_] := Module[{CoeffOut, c},
  CoeffOut = Flatten[Table[Join[{{c[[1]] <> "    :    ", Length[c[[2]]]}}, Partition[Flatten[c[[2]]], 1], {{""}}], {c, CoeffList}], 1];
  Export[dir<>"/"<>"Coefficients.dat", Flatten[#] & /@ Join[Join[{{"alat    :    "}}, Map[ToString[NumberForm[#, {3, 16}]] &, #] & /@ latt, {{}}], CoeffOut]];
]

FortranParamModule[CoeffList_, FieldsDef_] := Module[{Fint, n, head, tail, body, c, Fout, FieldDim, Charges, mass},
  Fint[n_] := StringRepeat["  ", n];
  FieldDim = FieldsDef\[Transpose][[2]];
  Charges = FieldsDef\[Transpose][[3]]; 
  mass = FieldsDef\[Transpose][[4]];
  head ={{"Module Parameters"},
         {Fint[1] <> "use Constants"}, 
         {Fint[1] <> "implicit none"}, 
         {Fint[1] <> "Integer, parameter    :: ifileno = 55                          !< File handle number for input files"}, 
         {Fint[1] <> "Integer, parameter    :: ofileno = 66                          !< File handle number for output files"},
         {Fint[1] <> "Integer               :: NumField = " <> ToString[Length[FieldDim]] <> ", FieldDim = " <> ToString[Max[FieldDim]]},
         {Fint[1] <> "Real*8, dimension("<>ToString[Length[FieldDim]]<>")  :: FieldCharge = (/" <> StringJoin[Riffle[ToString[#] & /@ Charges, ", "]] <> "/)"},
         {Fint[1] <> "Integer, dimension("<>ToString[Length[FieldDim]+1]<>") :: FieldDimList = " <> "(/" <> StringJoin[Riffle[ToString[#] & /@ Join[FieldDim, {6}], ", "]] <> "/)"},
         {Fint[1] <> "Integer               :: OnSiteDim = " <> ToString[Total[FieldDim]]},
         {Fint[1] <> "Real*8                :: alat(3,3)"},
         {Fint[1] <> "Real*8, dimension("<> ToString[Length[FieldDim]] <>")  :: mass = " <> "(/" <> StringJoin[Riffle[mass, ", "]] <> "/)"},
         {Fint[1] <> "Real*8, dimension (:,:,:,:,:), allocatable :: EwaldMat"},
         {Fint[1]}};
  tail = {{Fint[1]}, { "End Module Parameters"}};
  body = Table[{Fint[1] <> "Real*8                :: " <> c[[1]] <> "(" <> ToString[Length[c[[2]]]] <> ")"}, {c, CoeffList}];
  Fout = Join[head, body, tail];
  Return[Fout]
]

FortranReadCoeff[CoeffList_] := Module[{Fint, n, c, CaseBlock, default, head, tail, body, Fout},
  Fint[n_] := StringRepeat["  ", n];
  head = {{"Subroutine ReadCoefficients(ifile)"}, 
          {Fint[1] <> "use Parameters"}, 
          {Fint[1] <> "use FileParser"}, 
          {Fint[1] <> "implicit none"}, 
          {Fint[1] <> "integer, intent(in) :: ifile"}, 
          {Fint[1] <> "character(len=50)   :: keyword, cache "}, 
          {Fint[1] <> "integer             :: rd_len,i_err,i,j,i_errb, pos, Ndim"}, 
          {Fint[1] <> "logical             :: comment"}, 
          {Fint[1]}, 
          {Fint[1] <> "do"}, 
          {Fint[2] <> "10 continue"}, 
          {Fint[2] <> "keyword=\"\""}, 
          {Fint[2] <> "call bytereader(keyword,rd_len,ifile,i_errb)"}, 
          {Fint[2] <> "call caps2small(keyword)"}, 
          {Fint[2] <> "comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&"}, 
          {Fint[2] <> "(scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&"}, {Fint[2] <> "(scan(trim(keyword),'!')==1))"}, 
          {Fint[2]}, 
          {Fint[2] <> "if (comment) then"}, 
          {Fint[3] <> "read(ifile,*) cache"}, 
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
          {"End Subroutine ReadCoefficients"}};
  CaseBlock = Join[{{"alat", {{Fint[3] <> "read(ifile, '(A)', iostat=i_err) cache"},
        {Fint[3] <> "read(ifile, *, iostat=i_err)((alat(i,j),i=1,3),j=1,3)"},
        {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},{Fint[3]}}}},
     Table[{c[[1]],
     {{Fint[3] <> "read(ifile, '(A)', iostat=i_err) cache"},
      {Fint[3] <> "pos = scan(cache, ':')"},
      {Fint[3] <> "cache = trim(cache(pos+1:))"},
      {Fint[3] <> "read(cache, *, iostat=i_err) Ndim"},
      {Fint[3] <> "do i = 1, Ndim"},
      {Fint[3] <> "read(ifile,*,iostat=i_err) " <> c[[1]] <> "(i)"},
      {Fint[3] <> "end do"},
      {Fint[3] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
      {Fint[3]}}}, {c, CoeffList}]];
  default = GetCaseDefaults[{{"if(len(trim(keyword))>0) then"}, Fint[1] <> {"read(ifile,*)"}, {"end if"}}, 1];
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
          {Fint[3] <> "if(N1.ne.NGridx .or. N2.ne.NGridy .or. N3.ne.NGridz .or. NameEwald .ne. NameSim) then"},
          {Fint[4] <> "write(*,*) 'Dipole file mismatch, current simulation: ', NameSim, NGridx, NGridy, NGridz"},
          {Fint[4] <> "write(*,*) 'Ewald file: ', NameEwald, N1, N2, N3"},
          {Fint[4] <> "write(*,*) 'Now, regenerate Ewald Matrix!!!'"},
          {Fint[4] <> "call EwaldMatrix(alat)"},
          {Fint[3] <> "else"},
          {Fint[4] <> "write(*,*) 'Dipole file matches and will be used.'"},
          {Fint[3] <> "end if"},
          {Fint[3] <> "close(ifileno)"},
          {Fint[2] <> "else"},
          {Fint[3] <> "write(*,*) 'No EwaldMat.dat is found, generating ...'"},
          {Fint[3] <> "call EwaldMatrix(alat)"},
          {Fint[2] <> "end if"},
          {Fint[1] <> "end if"},
          {Fint[1]},
          {Fint[1] <> "allocate(EwaldMat(3,3,NGridx,NGridy,NGridz))"},
          {Fint[1] <> "open(ifileno,file='EwaldMat.dat',form='unformatted',status='old')"},
          {Fint[1] <> "read(ifileno) NameEwald, N1, N2, N3"},
          {Fint[1] <> "read(ifileno) (((((EwaldMat(j,i,ix,iy,iz),j=1,3),i=1,3),ix=1,NGridx),iy=1,NGridy),iz=1,NGridz)"},
          {Fint[1] <> "close(ifileno)"},
          {Fint[1]},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1] <> "! Initialize Fields !"},
          {Fint[1] <> "!!!!!!!!!!!!!!!!!!!!!"},
          {Fint[1]},
          {Fint[1] <> "allocate(Fields(FieldDim,NumField,NGridx,NGridy,NGridz))"},
          {Fint[1] <> "Call GetInitConfig(Fields, e0ij, RestartFields, FileCode)"},
          {Fint[1]},
          {Fint[1] <> "allocate(dFieldsdt(FieldDim,NumField,NGridx,NGridy,NGridz))"},
          {Fint[1] <> "Call GetInitConfig(dFieldsdt, de0ijdt, RestartVelocity, FileCode)"},
          {Fint[1]},
          {Fint[1] <> "allocate(EwaldField(FieldDim,NumField,NGridx,NGridy,NGridz))"},
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

FortranParsers[] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"Module Fileparser"},
          {Fint[1] <> "implicit none"},
          {Fint[1] <> "Contains"},
          {Fint[1] <> ""},
          {Fint[1] <> "Subroutine bytereader(keyword,rd_len,ifile,i_err)"},
          {Fint[2] <> "implicit none"},
          {Fint[2] <> "character(len=*), intent(out) :: keyword  !< Parsed keyword"},
          {Fint[2] <> "integer, intent(out) :: rd_len !< Length of parsed keyword"},
          {Fint[2] <> "integer, intent(in) :: ifile  !< File to read from"},
          {Fint[2] <> "integer, intent(out) :: i_err  !< Error status of reading"},
          {Fint[2] <> "logical :: rd_done,rd_start"},
          {Fint[2]},
          {Fint[2] <> "rd_done=.false."},
          {Fint[2] <> "rd_start=.false."},
          {Fint[2] <> "rd_len=0"},
          {Fint[2] <> "do while(.not.rd_done.and.rd_len<len(keyword))"},
          {Fint[3] <> "rd_len=rd_len+1"},
          {Fint[3] <> "read(ifile,'(a1)',advance='no',end=20,eor=10) keyword(rd_len:rd_len)"},
          {Fint[3] <> "rd_start=rd_start.or.keyword(rd_len:rd_len)/=\" \".or.keyword(rd_len:rd_len)/=\":\""},
          {Fint[3] <> "rd_done=rd_start.and.(keyword(rd_len:rd_len)==\" \".or.keyword(rd_len:rd_len)==\":\")"},
          {Fint[3] <> "if(keyword(rd_len:rd_len)==\":\") keyword(rd_len:rd_len)=\"\""},
          {Fint[2] <> "end do"},
          {Fint[2]},
          {Fint[2] <> "i_err=0"},
          {Fint[2] <> "keyword=adjustl(keyword(1:rd_len)//'')"},
          {Fint[2] <> "return"},
          {Fint[2] <> "! final word"},
          {Fint[2] <> "10  continue"},
          {Fint[2] <> "i_err=10"},
          {Fint[2] <> "keyword=adjustl(keyword(1:rd_len)//'')"},
          {Fint[2] <> "return"},
          {Fint[2] <> "! end of file"},
          {Fint[2] <> "20  continue"},
          {Fint[2] <> "i_err=20"},
          {Fint[2] <> "keyword=adjustl(keyword(1:rd_len)//'')"},
          {Fint[2] <> "return"},
          {Fint[2]},
          {Fint[1] <> "End Subroutine bytereader"},
          {Fint[1]},
          {Fint[1]},
          {Fint[1] <> "!> Convert lower case characters to upper case"},
          {Fint[1] <> "Subroutine small2caps(str)"},
          {Fint[2] <> "implicit none"},
          {Fint[2] <> "character(len=*),intent(inout):: str  !< string to convert"},
          {Fint[2]},
          {Fint[2] <> "integer i"},
          {Fint[2]},
          {Fint[2] <> "do i=1,len(str)"},
          {Fint[3] <> "if(str(i:i)>=\"a\" .and. str(i:i)<= \"z\") str(i:i)=achar(iachar(str(i:i))-32)"},
          {Fint[2] <> "end do"},
          {Fint[2]},
          {Fint[1] <> "End Subroutine small2caps"},
          {Fint[1]},
          {Fint[1]},
          {Fint[1] <> "!> Convert upper case characters to lower case"},
          {Fint[1] <> "Subroutine caps2small(str)"},
          {Fint[2] <> "implicit none"},
          {Fint[2] <> "character(len=*),intent(inout):: str  !< string to convert"},
          {Fint[2]},
          {Fint[2] <> "integer i"},
          {Fint[2]},
          {Fint[2] <> "do i=1,len(str)"},
          {Fint[3] <> "if(str(i:i)>=\"A\" .and. str(i:i)<= \"Z\") str(i:i)=achar(iachar(str(i:i))+32)"},
          {Fint[2] <> "end do"},
          {Fint[2]},
          {Fint[1] <> "End Subroutine caps2small"},
          {Fint[1]},
          {Fint[1]},
          {"End Module Fileparser"}};
  Return[Fout]
]

FortranInputs[] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"Module Inputs"},
          {Fint[1] <> "Use Parameters"},
          {Fint[1] <> "implicit none"},
          {Fint[1]},
          {Fint[1] <> "logical :: sane_input = .true."},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "! Geometry and composition"},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "integer :: NGridx              !< Number of cell repetitions in x direction"},
          {Fint[1] <> "integer :: NGridy              !< Number of cell repetitions in y direction"},
          {Fint[1] <> "integer :: NGridz              !< Number of cell repetitions in z direction"},
          {Fint[1]},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "! Exchange data"},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "real(dblprec), dimension(:,:,:), allocatable :: redcoord       !< Coordinates for Heisenberg exchange couplings"},
          {Fint[1] <> "real(dblprec), dimension(:,:,:,:,:), allocatable :: jc         !< Exchange couplings"},
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
          {Fint[1] <> "character(len=10) :: Solver                 !< Model solver"},
          {Fint[1] <> "integer           :: ThermoSteps            !< Thermo up"},
          {Fint[1] <> "integer           :: CoolingSteps           !< Thermo up"},
          {Fint[1] <> "integer           :: NumSteps               !< Number of Monte Carlo steps"},
          {Fint[1] <> "integer           :: TapeRate               !< Number of steps in measurement phase"},
          {Fint[1] <> "Logical           :: DipoleQ                !< dipole dipole interaction"},
          {Fint[1] <> "real*8            :: Temp                   !< Temperature"},
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
          {Fint[1] <> "! Parallel Tempering Monte Carlo"},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "Real*8  :: ReplicaT0"},
          {Fint[1] <> "Real*8  :: ReplicaTN"},
          {Fint[1] <> "integer :: SwapRate               !< Number of sweeps for one swap"},
          {Fint[1]},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "! Molecular Dynamics"},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "Real*8  :: deltaT"},
          {Fint[1] <> "Real*8  :: NoseMass"},
          {Fint[1]},
          {Fint[1] <> "Contains"},
          {Fint[1]},
          {Fint[1] <> "include \"ReadCoefficients.f90\""},
          {Fint[1]},
          {Fint[1] <> "Subroutine ReadParameters(ifile)"},
          {Fint[2] <> "Use FileParser"},
          {Fint[2] <> "Use Parameters"},
          {Fint[2] <> "Use Constants"},
          {Fint[2]},
          {Fint[2] <> "implicit none"},
          {Fint[2]},
          {Fint[2] <> "integer, intent(in) :: ifile   !< File to read from"},
          {Fint[2] <> "character(len=50)   :: keyword, cache "},
          {Fint[2] <> "integer             :: rd_len,i_err,i,i_stat,i_errb,ii, i_all, pos"},
          {Fint[2] <> "logical             :: comment"},
          {Fint[2]},
          {Fint[2] <> "do"},
          {Fint[3] <> "10     continue"},
          {Fint[3] <> "! Read file character for character until first whitespace"},
          {Fint[3] <> "keyword="""},
          {Fint[3] <> "call bytereader(keyword,rd_len,ifile,i_errb)"},
          {Fint[3]},
          {Fint[3] <> "! converting Capital letters"},
          {Fint[3] <> "call caps2small(keyword)"},
          {Fint[3] <> "! check for comment markers (currently % and #)"},
          {Fint[3] <> "comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&"},
          {Fint[3] <> "(scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&"},
          {Fint[3] <> "(scan(trim(keyword),'!')==1))"},
          {Fint[3] <> "if (comment) then"},
          {Fint[4] <> "read(ifile,*)"},
          {Fint[3] <> "else"},
          {Fint[4] <> "! Parse keyword"},
          {Fint[4] <> "keyword=trim(keyword)"},
          {Fint[4] <> "select case(keyword)"},
          {Fint[5] <> "case('restartfields')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) RestartFields"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('restartvelocity')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) RestartVelocity"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('namesim')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) NameSim"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('solver')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) Solver"},
          {Fint[5] <> "TrajectoryFile = trim(Solver)"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('numsteps')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) NumSteps"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('thermosteps')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) ThermoSteps"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('coolingsteps')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) CoolingSteps"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('dipoleq')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) DipoleQ"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('deltat')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) DeltaT"},
          {Fint[5] <> "DeltaT = DeltaT/time_fs"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('replicat0')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) ReplicaT0"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('replicatn')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) ReplicaTN"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('nosemass')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) NoseMass"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('temp')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) Temp"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('taperate')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) TapeRate"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err             "},
          {Fint[4] <> "case('swaprate')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) SwapRate"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('seed')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) seed"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('damp')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) damp"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('dampratio')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) DampRatio"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('acceptratio')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) AcceptRatio"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('coefffile')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) CoeffFile"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case('mesh')"},
          {Fint[5] <> "read(ifile, '(A)', iostat=i_err) cache"},
          {Fint[5] <> "pos = scan(cache, '=')"},
          {Fint[5] <> "cache = trim(cache(pos+1:))"},
          {Fint[5] <> "read(cache, *, iostat=i_err) NGridx, NGridy, NGridz"},
          {Fint[5] <> "if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err"},
          {Fint[4] <> "case default"},
          {Fint[5] <> "if(len(trim(keyword))>0) then"},
          {Fint[6] <> "read(ifile,*)"},
          {Fint[5] <> "end if"},
          {Fint[5]},
          {Fint[4] <> "end select"},
          {Fint[3] <> "end if"},
          {Fint[3]},
          {Fint[3] <> "! End of file"},
          {Fint[3] <> "if (i_errb==20) goto 20"},
          {Fint[3] <> "! End of row"},
          {Fint[3] <> "if (i_errb==10) goto 10"},
          {Fint[2] <> "end do"},
          {Fint[2]},
          {Fint[2] <> "20  continue"},
          {Fint[2]},
          {Fint[2] <> "return"},
          {Fint[1] <> "End Subroutine ReadParameters"},
          {Fint[1]},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "! SUBROUTINE: set_input_defaults"},
          {Fint[1] <> "!> Sets default values to input variables"},
          {Fint[1] <> "!---------------------------------------------------------------------------------"},
          {Fint[1] <> "Subroutine set_input_defaults()"},
          {Fint[2] <> "Use FileParser"},
          {Fint[2]},
          {Fint[2] <> "implicit none"},
          {Fint[2]},
          {Fint[2] <> "real(dblprec) :: one=1.0_dblprec"},
          {Fint[2] <> "real(dblprec) :: zero=0.0_dblprec"},
          {Fint[2]},
          {Fint[2] <> "!Geometry and composition"},
          {Fint[2] <> "NGridx                = 1"},
          {Fint[2] <> "NGridy                = 1"},
          {Fint[2] <> "NGridz                = 1"},
          {Fint[2]},
          {Fint[2] <> "!Solvers"},
          {Fint[2] <> "Solver            = \"MCMC\""},
          {Fint[2] <> "ThermoSteps       = 10000"},
          {Fint[2] <> "CoolingSteps      = 0"},
          {Fint[2] <> "NumSteps          = 10000"},
          {Fint[2] <> "TapeRate          = 1000"},
          {Fint[2] <> "DipoleQ           = .true."},
          {Fint[2] <> "call caps2small(Solver)"},
          {Fint[2]},
          {Fint[2] <> "!Molecular Dynamics"},
          {Fint[2] <> "DeltaT            =1.0D-6"},
          {Fint[2] <> "NoseMass          =1.0D1"},
          {Fint[2]},
          {Fint[2] <> "!Markov Chain Monte Carlo"},
          {Fint[2] <> "seed              = 0"},
          {Fint[2] <> "damp              = 0.001"},
          {Fint[2] <> "dampRatio         = 1.5"},
          {Fint[2] <> "AcceptRatio       = 0.3"},
          {Fint[2] <> "AvrgInterval      = 100"},
          {Fint[2] <> "BuffMcAvrg        = 0"},
          {Fint[2] <> "BuffMcAvrg        = 0.001"},
          {Fint[2]},
          {Fint[2] <> "!Parallel Tempering Monte Carlo"},
          {Fint[2] <> "ReplicaT0         = 0.0D0"},
          {Fint[2] <> "ReplicaTN         = 500.0D0"},
          {Fint[2] <> "SwapRate          = 1"},
          {Fint[2]},
          {Fint[2] <> "!Simulation"},
          {Fint[2] <> "RestartFields     = \"random\""},
          {Fint[2] <> "RestartVelocity   = \"random\""},
          {Fint[2] <> "NameSim           = \"LINVARIANT\""},
          {Fint[2] <> "CoeffFile         = \"Coefficients.dat\""},
          {Fint[2] <> "TrajectoryFile    = trim(Solver)"},
          {Fint[2] <> "aunits            = \"N\""},
          {Fint[1] <> "End Subroutine set_input_defaults"},
          {Fint[1]},
          {Fint[1] <> "Subroutine InitFromFile(filename, Fields, e0ij)"},
          {Fint[2] <> "Implicit None"},
          {Fint[2]},
          {Fint[2] <> "character(*), Intent(in) :: filename"},
          {Fint[2] <> "Integer                  :: FileHandle = 1111"},
          {Fint[2] <> "Integer                  :: ix, iy, iz, i, ifield"},
          {Fint[2] <> "Real*8                   :: e"},
          {Fint[2] <> "Real*8, Intent(inout)    :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
          {Fint[2] <> "Real*8, Intent(inout)    :: e0ij(3,3)"},
          {Fint[2]},
          {Fint[2] <> "e0ij = 0.0D0"},
          {Fint[2]},
          {Fint[2] <> "open(FileHandle,file=filename,form='formatted',status='old')"},
          {Fint[2] <> "Read(FileHandle, \"(3E25.15)\") e0ij(1,1), e0ij(2,2), e0ij(3,3)"},
          {Fint[2] <> "Read(FileHandle, \"(3E25.15)\") e0ij(2,3), e0ij(1,3), e0ij(1,2)"},
          {Fint[2] <> "do While (.True.)"},
          {Fint[3] <> "Read(FileHandle, \"(3I10)\", END=999)      ix, iy, iz"},
          {Fint[3] <> "do ifield = 1, NumField"},
          {Fint[4] <> "Read(FileHandle, \"(3E25.15)\", END=999) (Fields(i, ifield, ix, iy, iz), i=1,FieldDimList(ifield))"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "999  close(FileHandle)"},
          {Fint[2]},
          {Fint[1] <> "End Subroutine InitFromFile"},
          {"End Module Inputs"}};
  Return[Fout]
]

FortranOutputs[] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"Module Outputs "},
          {Fint[1] <> "Use LINVARIANT"},
          {Fint[1] <> "Use Parameters"},
          {Fint[1] <> "Use Inputs"},
          {Fint[1] <> "implicit none"},
          {Fint[1]},
          {Fint[1] <> "Contains"},
          {Fint[1]},
          {Fint[1] <> "subroutine fmkdir(newDirPath)"},
          {Fint[2] <> "implicit none"},
          {Fint[2]},
          {Fint[2] <> "character(len=*), intent(in) :: newDirPath"},
          {Fint[2] <> "character(len=256)           :: mkdirCmd"},
          {Fint[2] <> "logical                      :: dirExists"},
          {Fint[2]},
          {Fint[2] <> "! Check if the directory exists first"},
          {Fint[2] <> "! inquire(file=trim(newDirPath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort"},
          {Fint[2] <> "inquire(directory=newDirPath, exist=dirExists)         ! Works with ifort, but not gfortran"},
          {Fint[2]},
          {Fint[2] <> "if (dirExists) then"},
          {Fint[3] <> "write (*,*) \"Directory already exists: '\"//trim(newDirPath)//\"'\""},
          {Fint[2] <> "else"},
          {Fint[3] <> "mkdirCmd = 'mkdir -p '//trim(newDirPath)"},
          {Fint[3] <> "write(*,'(a)') \"Creating new directory: '\"//trim(mkdirCmd)//\"'\""},
          {Fint[3] <> "call system(mkdirCmd)"},
          {Fint[2] <> "endif"},
          {Fint[1] <> "end subroutine fmkdir"},
          {Fint[1]},
          {Fint[1] <> "Subroutine WriteBinary(FileHandle, Fields, dFieldsdt, e0ij, de0ijdt)"},
          {Fint[2] <> "Implicit None"},
          {Fint[2]},
          {Fint[2] <> "Integer        :: FileHandle"},
          {Fint[2] <> "Integer        :: ix, iy, iz, i, ifield"},
          {Fint[2] <> "Real*8         :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
          {Fint[2] <> "Real*8         :: dFieldsdt(FieldDim, NumField, NGridx, NGridy, NGridz)"},
          {Fint[2] <> "Real*8         :: Etot, Epot, Ekin, e0ij(3,3), de0ijdt(3,3)"},
          {Fint[2]},
          {Fint[2] <> "Epot = GetEtot(Fields, e0ij)"},
          {Fint[2] <> "Ekin = GetEkin(dFieldsdt, de0ijdt)"},
          {Fint[2] <> "Etot = Epot + Ekin"},
          {Fint[2]},
          {Fint[2] <> "write(FileHandle) Etot, Epot, Ekin"},
          {Fint[2] <> "write(FileHandle) e0ij(1,1), e0ij(2,2), e0ij(3,3)"},
          {Fint[2] <> "write(FileHandle) e0ij(2,3), e0ij(1,3), e0ij(1,2)"},
          {Fint[2] <> "do iz = 1, NGridz"},
          {Fint[3] <> "do iy = 1, NGridy"},
          {Fint[4] <> "do ix = 1, NGridx"},
          {Fint[5] <> "write(FileHandle) ix, iy, iz"},
          {Fint[5] <> "!          write(FileHandle) EOnSite(ix,iy,iz,Fields)"},
          {Fint[5] <> "do ifield = 1, NumField"},
          {Fint[6] <> "write(FileHandle) (Fields(i,ifield,ix,iy,iz), i=1,FieldDimList(ifield))"},
          {Fint[5] <> "end do"},
          {Fint[4] <> "end do"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "end do"},
          {Fint[1] <> "End Subroutine"},
          {Fint[1]},
          {Fint[1] <> "Subroutine WriteFinal(filename, Fields, e0ij)"},
          {Fint[2] <> "Implicit None"},
          {Fint[2]},
          {Fint[2] <> "character(*)   :: filename"},
          {Fint[2] <> "Integer        :: FileHandle = 1111"},
          {Fint[2] <> "Integer        :: ix, iy, iz, i, ifield"},
          {Fint[2] <> "Real*8         :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
          {Fint[2] <> "Real*8         :: e0ij(3,3)"},
          {Fint[2]},
          {Fint[2] <> "open(FileHandle,file=trim(Solver)//'.out/'//filename,form='formatted',status='unknown')"},
          {Fint[2] <> "write(FileHandle, \"(3E25.15)\") e0ij(1,1), e0ij(2,2), e0ij(3,3)"},
          {Fint[2] <> "write(FileHandle, \"(3E25.15)\") e0ij(2,3), e0ij(1,3), e0ij(1,2)"},
          {Fint[2] <> "do iz = 1, NGridz"},
          {Fint[3] <> "do iy = 1, NGridy"},
          {Fint[4] <> "do ix = 1, NGridx"},
          {Fint[5] <> "write(FileHandle, \"(3I10)\")    ix, iy, iz"},
          {Fint[5] <> "!          write(FileHandle, \"(E25.15)\")  EOnSite(ix,iy,iz,Fields)"},
          {Fint[5] <> "do ifield = 1, NumField"},
          {Fint[6] <> "write(FileHandle, \"(3E25.15)\") (Fields(i,ifield,ix,iy,iz), i=1,FieldDimList(ifield))"},
          {Fint[5] <> "end do"},
          {Fint[4] <> "end do"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "end do"},
          {Fint[2] <> "close(FileHandle)"},
          {Fint[1] <> "End Subroutine"},
          {Fint[1]},
          {Fint[1] <> "Subroutine BinaryToData(filename, processor)"},
          {Fint[2] <> "Implicit None"},
          {Fint[2]},
          {Fint[2] <> "character(*)    :: filename"},
          {Fint[2] <> "Integer         :: io_mode, io_out, processor"},
          {Fint[2] <> "Integer         :: ix, iy, iz, imc=0"},
          {Fint[2] <> "Integer         :: ifield, i, igrid"},
          {Fint[2] <> "Real*8          :: Etot, Epot, Ekin"},
          {Fint[2] <> "Real*8          :: e11, e22, e33, e23, e13, e12"},
          {Fint[2] <> "Real*8          :: field(FieldDim, NumField)"},
          {Fint[2]},
          {Fint[2] <> "io_mode = 11"},
          {Fint[2] <> "io_out = 12"},
          {Fint[2]},
          {Fint[2] <> "Open(io_mode,file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(int2str(processor))//'.dat',&"},
          {Fint[2] <> "form='unformatted',status='old')"},
          {Fint[2] <> "Open(io_out,file=trim(Solver)//'.out/'//trim(filename)//'-'//trim(int2str(processor))//'.dat',&"},
          {Fint[2] <> "status='unknown')"},
          {Fint[2]},
          {Fint[2] <> "do While (.True.)"},
          {Fint[3] <> "imc = imc + 1"},
          {Fint[3]},
          {Fint[3] <> "Read(io_mode, END=999) Etot, Epot, Ekin"},
          {Fint[3] <> "Read(io_mode, END=999) e11, e22, e33"},
          {Fint[3] <> "Read(io_mode, END=999) e23, e13, e12"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") Etot, Epot, Ekin"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") e11, e22, e33"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") e23, e13, e12"},
          {Fint[3]},
          {Fint[3] <> "do igrid = 1, NGridz*NGridy*NGridx"},
          {Fint[4] <> "Read(io_mode, END=999) ix, iy, iz"},
          {Fint[4] <> "do ifield = 1, NumField"},
          {Fint[5] <> "Read(io_mode, END=999) (field(i, ifield), i=1,FieldDimList(ifield))"},
          {Fint[4] <> "end do"},
          {Fint[4] <> "Write(io_out, \"(3I10)\") ix, iy, iz"},
          {Fint[4] <> "do ifield = 1, NumField"},
          {Fint[5] <> "Write(io_out, \"(3E25.15)\") (field(i, ifield), i=1,FieldDimList(ifield))"},
          {Fint[4] <> "end do"},
          {Fint[3] <> "end do"},
          {Fint[2] <> "End do"},
          {Fint[2]},
          {Fint[2] <> "999 Close(io_mode)"},
          {Fint[2] <> "Close(io_out)"},
          {Fint[2]},
          {Fint[1] <> "END Subroutine BinaryToData"},
          {Fint[1]},
          {Fint[1] <> "Subroutine GetObservables(filename, order, processor)"},
          {Fint[2] <> "Implicit None"},
          {Fint[2]},
          {Fint[2] <> "character(*)    :: filename"},
          {Fint[2] <> "Integer         :: order, io_mode, io_out, processor"},
          {Fint[2] <> "Integer         :: ix, iy, iz, imc=0"},
          {Fint[2] <> "Integer         :: ifield, i, j, igrid"},
          {Fint[2] <> "Real*8          :: Etot, Epot, Ekin"},
          {Fint[2] <> "Real*8          :: e11, e22, e33, e23, e13, e12"},
          {Fint[2] <> "Real*8          :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
          {Fint[2] <> "Real*8          :: Observables(FieldDim, NumField)"},
          {Fint[2]},
          {Fint[2] <> "io_mode = 11"},
          {Fint[2] <> "io_out = 12"},
          {Fint[2]},
          {Fint[2] <> "Open(io_mode,file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(int2str(processor))//'.dat',&"},
          {Fint[2] <> "form='unformatted',status='old')"},
          {Fint[2] <> "Open(io_out,file=trim(Solver)//'.out/'//trim(filename)//'-'//trim(int2str(processor))//'.dat',&"},
          {Fint[2] <> "status='unknown')"},
          {Fint[2]},
          {Fint[2] <> "do While (.True.)"},
          {Fint[3] <> "imc = imc + 1"},
          {Fint[3]},
          {Fint[3] <> "Read(io_mode, END=999) Etot, Epot, Ekin"},
          {Fint[3] <> "Read(io_mode, END=999) e11, e22, e33"},
          {Fint[3] <> "Read(io_mode, END=999) e23, e13, e12"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") Etot, Epot, Ekin"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") e11, e22, e33"},
          {Fint[3] <> "Write(io_out, \"(3E25.15)\") e23, e13, e12"},
          {Fint[3]},
          {Fint[3] <> "do igrid = 1, NGridz*NGridy*NGridx"},
          {Fint[4] <> "Read(io_mode, END=999) ix, iy, iz"},
          {Fint[4] <> "do ifield = 1, NumField"},
          {Fint[5] <> "Read(io_mode, END=999) (Fields(i, ifield, ix, iy, iz), i=1,FieldDimList(ifield))"},
          {Fint[4] <> "end do"},
          {Fint[3] <> "end do"},
          {Fint[3]},
          {Fint[3] <> "do ifield = 1, NumField"},
          {Fint[4] <> "do j = 1, FieldDimList(ifield)"},
          {Fint[5] <> "Observables(j, ifield) = Sum(Fields(j,ifield,:,:,:)**order)/(NGridx*NGridy*NGridz)"},
          {Fint[4] <> "end do"},
          {Fint[4] <> "Write(io_out, \"(3E25.15)\") (Observables(i, ifield), i=1,FieldDimList(ifield))"},
          {Fint[3] <> "end do"},
          {Fint[3]},
          {Fint[2] <> "End do"},
          {Fint[2]},
          {Fint[2] <> "999 Close(io_mode)"},
          {Fint[2] <> "Close(io_out)"},
          {Fint[2]},
          {Fint[1] <> "END Subroutine GetObservables"},
          {Fint[1]},
          {"End Module Outputs"}};
  Return[Fout]
]

Makefile[dir_] := Module[{Fint, Fout},
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{".SUFFIXES: .f .f90 .F90 .SUFFIXES .prj"},
          {" "},
          {"FC = mpif90 -fopenmp -O3"},
          {"FC = mpiifort -fopenmp -O3"},
          {" "},
          {"PREFIX = /home/xxx/directory/to/LINVARIANT"},
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
          {"	    echo \"To decouple the dipole-dipole between different type of fields, modify GetEOnSiteEwald.f90, GetEwaldForces.f90, GetEwaldField.f90, UpdateEwaldField.f90\""},
          {"	    echo \"To decouple microscopic and macroscopic strain, modify GetDeltaHEps.f90, GetEOnSiteEps.f90, GetForcesEps.f90\""},
          {"	    echo \"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\""},
          {" "},
          {"clean: "},
          {"	rm -f *.f90 *.mod $(OBJS) $(exe)"},
          {"#"}};
  Export[dir, Fout, "Table", "TextDelimiters"->None];
]

LINVARIANTModule[DeltaFunctionList_?ListQ, OtherFunctionList_:{}, OptionsPattern[{"EwaldField"->True}]] := Module[{Fint, n, Fout, IncFuncBlock, DeltaHBlock, EwaldFieldBlock, EwaldFieldArg},
  Fint[n_] := StringRepeat["  ", n];
  EwaldFieldArg = If[OptionValue["EwaldField"], "EwaldField,", ""];
  EwaldFieldBlock = If[OptionValue["EwaldField"], {{Fint[2] <> " Real*8,  Intent(in) :: EwaldField(NumField, NGridx, NGridy, NGridz, FieldDim)"}}, {}];
  IncFuncBlock = Join[Flatten[{{Fint[1] <> " Include \"GetDeltaH" <> # <> ".f90\""}, 
                               {Fint[1] <> " Include \"GetForces" <> # <> ".f90\""},
                               {Fint[1] <> " Include \"GetEOnSite" <> # <> ".f90\""}} & /@ Join[DeltaFunctionList, {""}], 1],
                      {Fint[1] <> " Include \"" <> # <> ".f90\""} &/@ OtherFunctionList
                     ];
  DeltaHBlock = {Fint[3] <> "+ DeltaH" <> # <> "(ix, iy, iz, Fields, e0ij, euij, idelta, delta) &"} & /@ DeltaFunctionList;
  DeltaHBlock[[1]] = StringReplace[DeltaHBlock[[1]], "+ " -> ""];
  DeltaHBlock[[-1]] = StringReplace[DeltaHBlock[[-1]], " &" -> ""];
  Fout = Join[{{"Module LINVARIANT"},
     {Fint[1]},
     {Fint[1] <> "Use Parameters"},
     {Fint[1] <> "Use Constants"},
     {Fint[1] <> "Use Inputs"},
     {Fint[1]},
     {Fint[1] <> "Implicit none"},
     {Fint[1] <> "Contains"},
     {Fint[1]}},
    IncFuncBlock,
    {{Fint[1]},
     {Fint[1] <> "Function GetEtot(Fields, e0ij) Result(Etot)"}, 
     {Fint[2] <> "Implicit none"}, 
     {Fint[2] <> "Real*8,  Intent(in)    :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8,  Intent(in)    :: e0ij(3,3)"}, 
     {Fint[2] <> "Real*8                 :: Etot"}, 
     {Fint[2] <> "Integer                :: ix, iy, iz"}, 
     {Fint[2]}, 
     {Fint[2] <> "Etot = 0.0D0"}, 
     {Fint[2]}, 
     {Fint[2] <> "if (DipoleQ) then"}, 
     {Fint[3] <> "do iz = 1, NGridz"}, 
     {Fint[4] <> "do iy = 1, NGridy"}, 
     {Fint[5] <> "do ix = 1, NGridx"}, 
     {Fint[6] <> "Etot = Etot + GetEOnSite(ix, iy, iz, Fields, e0ij) + GetEOnSiteEwald(ix, iy, iz, Fields)"}, 
     {Fint[5] <> "end do"}, 
     {Fint[4] <> "end do"}, 
     {Fint[3] <> "end do"}, 
     {Fint[2] <> "else"}, 
     {Fint[3] <> "do iz = 1, NGridz"},
     {Fint[4] <> "do iy = 1, NGridy"},
     {Fint[5] <> "do ix = 1, NGridx"},
     {Fint[6] <> "Etot = Etot + GetEOnSite(ix, iy, iz, Fields, e0ij)"},
     {Fint[5] <> "end do"},
     {Fint[4] <> "end do"},
     {Fint[3] <> "end do"},
     {Fint[2] <> "end if"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Function GetEtot"},
     {Fint[1]},
     {Fint[1] <> "Function GridPbc(a,p) Result(b)"},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Integer::a,p,b"},
     {Fint[2]},
     {Fint[2] <> "b=a-floor((real(a)-1.0)/real(p))*p"},
     {Fint[1] <> "End Function GridPbc"},
     {Fint[1]},
     {Fint[1] <> "Subroutine RemoveGridDrifts(Fields)"}, 
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(inout) :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8                 :: GridDriftx, GridDrifty, GridDriftz"}, 
     {Fint[2] <> "Integer                :: ix, iy, iz"}, 
     {Fint[2]}, 
     {Fint[2] <> "GridDriftx = Sum(Fields(1,NumField,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2] <> "GridDrifty = Sum(Fields(2,NumField,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2] <> "GridDriftz = Sum(Fields(3,NumField,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2]}, 
     {Fint[2] <> "do iz = 1, NGridz"}, 
     {Fint[3] <> "do iy = 1, NGridy"}, 
     {Fint[4] <> "do ix = 1, NGridx"}, 
     {Fint[5] <> "Fields(1,NumField,ix,iy,iz) = Fields(1,NumField,ix,iy,iz) - GridDriftx"}, 
     {Fint[5] <> "Fields(2,NumField,ix,iy,iz) = Fields(2,NumField,ix,iy,iz) - GridDrifty"}, 
     {Fint[5] <> "Fields(3,NumField,ix,iy,iz) = Fields(3,NumField,ix,iy,iz) - GridDriftz"}, 
     {Fint[4] <> "end do"}, 
     {Fint[3] <> "end do"}, 
     {Fint[2] <> "end do"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Subroutine RemoveGridDrifts"},
     {Fint[1]},
     {Fint[1] <> "Function FieldsTo1D(Fields) Result(Fields1D)"},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in)    :: Fields(FieldDim, NumField, NGridx, NGridy, NGridz)"},
     {Fint[2] <> "Real*8                 :: Fields1D(FieldDim*NumField*NGridx*NGridy*NGridz)"},
     {Fint[2] <> "Integer                :: ix, iy, iz, ifield, i, id"},
     {Fint[2]},
     {Fint[2] <> "do iz = 1, NGridz"},
     {Fint[3] <> "do iy = 1, NGridy"},
     {Fint[4] <> "do ix = 1, NGridx"},
     {Fint[5] <> "do ifield = 1, NumField"},
     {Fint[6] <> "do i = 1, FieldDim"},
     {Fint[7] <> "id = (iz-1)*NGridy*NGridx*NumField*FieldDim &"},
     {Fint[7] <> "+ (iy-1)*NGridx*NumField*FieldDim &"},
     {Fint[7] <> "+ (ix-1)*NumField*FieldDim &"},
     {Fint[7] <> "+ (ifield-1)*FieldDim &"},
     {Fint[7] <> "+ i"},
     {Fint[7] <> "Fields1D(id) = Fields(i,ifield,ix,iy,iz)"},
     {Fint[6] <> "end do"},
     {Fint[5] <> "end do"},
     {Fint[4] <> "end do"},
     {Fint[3] <> "end do"},
     {Fint[2] <> "end do"},
     {Fint[2]},
     {Fint[1] <> "End Function FieldsTo1D"},
     {Fint[1]},
     {Fint[1] <> "Function FieldsToND(Fields1D) Result(Fields)"},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in)    :: Fields1D(FieldDim*NumField*NGridx*NGridy*NGridz)"},
     {Fint[2] <> "Real*8                 :: Fields(FieldDim,NumField,NGridx,NGridy,NGridz)"},
     {Fint[2]},
     {Fint[2] <> "Fields = Reshape(Fields1D, (/FieldDim, NumField, NGridx, NGridy, NGridz/))"},
     {Fint[2]},
     {Fint[1] <> "End Function FieldsToND"},
     {Fint[1]},
     {Fint[1] <> "Function eij2eta(eij) Result(eta)"},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in)    :: eij(3,3)"},
     {Fint[2] <> "Real*8                 :: eta(6)"},
     {Fint[2]},
     {Fint[2] <> "eta(1) = eij(1,1)"},
     {Fint[2] <> "eta(2) = eij(2,2)"},
     {Fint[2] <> "eta(3) = eij(3,3)"},
     {Fint[2] <> "eta(4) = eij(2,3)"},
     {Fint[2] <> "eta(5) = eij(1,3)"},
     {Fint[2] <> "eta(6) = eij(1,2)"},
     {Fint[2]},
     {Fint[1] <> "End Function eij2eta"},
     {Fint[1]},
     {Fint[1] <> "Function eta2eij(eta) Result(eij)"},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in)    :: eta(6)"},
     {Fint[2] <> "Real*8                 :: eij(3,3)"},
     {Fint[2]},
     {Fint[2] <> "eij(1,1) = eta(1)"},
     {Fint[2] <> "eij(2,2) = eta(2)"},
     {Fint[2] <> "eij(3,3) = eta(3)"},
     {Fint[2] <> "eij(2,3) = eta(4) "},
     {Fint[2] <> "eij(1,3) = eta(5) "},
     {Fint[2] <> "eij(1,2) = eta(6)"},
     {Fint[2]},
     {Fint[2] <> "eij(3,2) = eta(4)"},
     {Fint[2] <> "eij(3,1) = eta(5)"},
     {Fint[2] <> "eij(2,1) = eta(6)"},
     {Fint[2]},
     {Fint[1] <> "End Function eta2eij"},
     {Fint[1]},
     {Fint[1] <> "Function int2str(i) Result(str)"},
     {Fint[2] <> "integer, intent(in) :: i"},
     {Fint[2] <> "character(len=20)   :: str"},
     {Fint[2]},
     {Fint[2] <> "write (str, '(I5.5)') i"},
     {Fint[2] <> "str = trim(str)"},
     {Fint[1] <> "end function int2str"},
     {Fint[1]},
     {Fint[1] <> "Function Thermometer(dFieldsdt, de0ijdt) Result(TK)"}, 
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in) :: dFieldsdt(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8,  Intent(in) :: de0ijdt(3,3)"}, 
     {Fint[2] <> "Real*8              :: TK, Ekin"}, 
     {Fint[2] <> "Integer             :: i, ifield, ix, iy, iz, TotalDim"}, 
     {Fint[2]}, 
     {Fint[2] <> "TotalDim = NGridx*NGridy*NGridz*OnSiteDim + 6"}, 
     {Fint[2] <> "Ekin = GetEkin(dFieldsdt, de0ijdt)"}, 
     {Fint[2] <> "TK = 2*Ekin*Hartree/k_bolt_ev/TotalDim"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Function Thermometer"},
     {Fint[1]}, 
     {Fint[1] <> "Function GetEkin(dFieldsdt, de0ijdt) Result(Ekin)"}, 
     {Fint[2] <> "Implicit none"}, 
     {Fint[2] <> "Real*8,  Intent(in)  :: dFieldsdt(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8,  Intent(in)  :: de0ijdt(3,3)"}, 
     {Fint[2] <> "Real*8               :: Ekin, detadt(6)"}, 
     {Fint[2] <> "Integer              :: i, ifield, ix, iy, iz"}, 
     {Fint[2]}, 
     {Fint[2] <> "detadt = eij2eta(de0ijdt)"}, 
     {Fint[2] <> "Ekin = 0.0D0"}, 
     {Fint[2]}, 
     {Fint[2] <> "do iz = 1, NGridz"}, 
     {Fint[3] <> "do iy = 1, NGridy"}, 
     {Fint[4] <> "do ix = 1, NGridx"}, 
     {Fint[5] <> "do ifield = 1, NumField"}, 
     {Fint[6] <> "do i = 1, FieldDimList(ifield)"}, 
     {Fint[7] <> "Ekin = Ekin + 0.5d0*mass(ifield)*mpme*dFieldsdt(i, ifield, ix, iy, iz)**2"}, 
     {Fint[6] <> "end do"}, 
     {Fint[5] <> "end do"}, 
     {Fint[4] <> "end do"}, 
     {Fint[3] <> "end do"}, 
     {Fint[2] <> "end do"}, 
     {Fint[2] <> "do i = 1, 6"}, 
     {Fint[3] <> "Ekin = Ekin + 0.5d0*NGridx*NGridy*NGridz*mass(NumField+1)*mpme*detadt(i)**2"}, 
     {Fint[2] <> "end do"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Function GetEkin"},
     {Fint[1]}, 
     {Fint[1] <> "Subroutine NoseHooverUpdate(gm, T0, dFieldsdt, de0ijdt)"}, 
     {Fint[2] <> "Implicit none"}, 
     {Fint[2] <> "Real*8,  Intent(in)     :: dFieldsdt(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8,  Intent(in)     :: de0ijdt(3,3), T0"}, 
     {Fint[2] <> "Real*8,  Intent(inout)  :: gm"}, 
     {Fint[2] <> "Real*8                  :: Ekin"}, 
     {Fint[2] <> "Integer                 :: TotalDim"}, 
     {Fint[2]}, 
     {Fint[2] <> "TotalDim = NGridx*NGridy*NGridz*OnSiteDim"}, 
     {Fint[2] <> "Ekin = GetEkin(dFieldsdt, de0ijdt)"}, 
     {Fint[2] <> "gm = gm + (Ekin - 0.5D0*TotalDim*k_bolt_ev*T0/Hartree)*DeltaT/7.464D0/NoseMass"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Subroutine NoseHooverUpdate"},
     {Fint[1]}, 
     {Fint[1] <> "Subroutine get_walltime(wctime)"}, 
     {Fint[2] <> "Implicit none"}, 
     {Fint[2] <> "Real*4, Intent(out) :: wctime"}, 
     {Fint[2] <> "Integer             :: r, c"}, 
     {Fint[2]}, 
     {Fint[2] <> "call system_clock(c, r)"}, 
     {Fint[2]}, 
     {Fint[2] <> "wctime = Real(c) / r"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Subroutine get_walltime"},
     {Fint[1]},
     {Fint[1] <> " Function Cell2Volume(A) Result(VOLCELL)"},
     {Fint[1]},
     {Fint[2] <> "! CALCULATES THE VOLUME OF THE UNIT CELL,POSSIBLY WITH A MINUS SIGN"},
     {Fint[1]},
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8::A(3,3)"},
     {Fint[2] <> "Real*8::VOLCELL"},
     {Fint[2]},
     {Fint[2] <> "VOLCELL=(A(2,1)*A(3,2)-A(3,1)*A(2,2))*A(1,3)+&"},
     {Fint[2] <> "(A(3,1)*A(1,2)-A(1,1)*A(3,2))*A(2,3)+&"},
     {Fint[2] <> "(A(1,1)*A(2,2)-A(2,1)*A(1,2))*A(3,3)"},
     {Fint[2]},
     {Fint[1] <> "End Function Cell2Volume"},
     {Fint[1]},
     {Fint[1] <> "Subroutine RECLAT(A,B,IOPT)"},
     {Fint[2]},
     {Fint[2] <> "! CALCULATES RECIPROCAL LATTICE VECTORS.THEIR PRODUCT WITH DIRECT"},
     {Fint[2] <> "! LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1"},
     {Fint[2]},
     {Fint[2] <> "IMPLICIT NONE"},
     {Fint[2] <> "Real*8::A(3,3),B(3,3)"},
     {Fint[2] <> "integer::iopt,i"},
     {Fint[2] <> "Real*8::pi,c,ci"},
     {Fint[2]},
     {Fint[2] <> "PI=ACOS(-1.D0)"},
     {Fint[2] <> "B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)"},
     {Fint[2] <> "B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)"},
     {Fint[2] <> "B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)"},
     {Fint[2] <> "B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)"},
     {Fint[2] <> "B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)"},
     {Fint[2] <> "B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)"},
     {Fint[2] <> "B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)"},
     {Fint[2] <> "B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)"},
     {Fint[2] <> "B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)"},
     {Fint[2] <> "C=1.D0"},
     {Fint[2]},
     {Fint[2] <> "if(IOPT.eq.1) then"},
     {Fint[3] <> "C=2.D0*PI"},
     {Fint[2] <> "endif"},
     {Fint[2]},
     {Fint[2] <> "do i=1,3"},
     {Fint[3] <> "CI=C/(A(1,i)*B(1,i)+A(2,i)*B(2,i)+A(3,i)*B(3,i))"},
     {Fint[3] <> "B(1,i)=B(1,i)*CI"},
     {Fint[3] <> "B(2,i)=B(2,i)*CI"},
     {Fint[3] <> "B(3,i)=B(3,i)*CI"},
     {Fint[2] <> "end do"},
     {Fint[2]},
     {Fint[1] <> "End Subroutine"},
     {Fint[1]},
     {"End Module LINVARIANT"}}];
  Return[Fout]
]

Expr2Fortran[expr_, FortranVarSub_] := Module[{},
  TimesFactor2Real[#] & /@ If[MatchQ[expr, Plus[_, __]], Level[expr, {1}], Level[expr, {0}]] /. FortranVarSub
]

FortranVarStr[var_, arg_] := Module[{},
  ToString@FortranForm@ToExpression[var <> "[" <> StringRiffle[arg, ","] <> "]"]
]

FortranVarSub[MMAVars_] := Module[{var},
  var = DeleteCases[DeleteDuplicates[GetSubscriptInfo[#][[1]] & /@ Flatten[MMAVars]], ToExpression["\[Epsilon]0"]];
  Join[Table[(Subscript @@ Join[{var[[i]]}, ToExpression["{a_,ix_,iy_,iz_}"]] :> Fields[a, ii, HopingCode[x, ix], HopingCode[y, iy],         HopingCode[z, iz]]) /. ii -> i, {i, Length@var}],
       Table[(Subscript @@ Join[{var[[i]]}, ToExpression["{a_,ix_,iy_,iz_,0}"]] :> Fields[a, ii, HopingCode[x, ix], HopingCode[y, iy],       HopingCode[z, iz]]) /. ii -> i, {i, Length@var}],
       Subscript @@ Join[{ToExpression["\[Delta]" <> ToString[#]]}, ToExpression["{a_,0,0,0}"]] :> delta[a] & /@ var,
       {Subscript[ToExpression["\[Epsilon]"], ToExpression["i_"], ToExpression["j_"]] :> eij[i, j]},
       {Subscript[ToExpression["\[Epsilon]0"], ToExpression["i_"], ToExpression["j_"]] :> e0ij[i, j]},
       {Subscript[ToExpression["\[Epsilon]0"], ToExpression["i_"], ToExpression["j_"], 0] :> e0ij[i, j]},
       {Subscript[ToExpression["\[Epsilon]u"], ToExpression["i_"], ToExpression["j_"]] :> euij[i, j]},
       {Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 1] :> delta[1],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 2, 2] :> delta[2],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 3, 3] :> delta[3],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 2, 3] :> delta[4],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 3] :> delta[5],
        Subscript[ToExpression["\[Delta]\[Epsilon]0"], 1, 2] :> delta[6]}] /. {a -> ToExpression["a"], ix -> ToExpression["ix"], iy ->        ToExpression["iy"], iz -> ToExpression["iz"], x -> ToExpression["x"], y -> ToExpression["y"], z -> ToExpression["z"], delta -> ToExpression["delta"],  Fields -> ToExpression["Fields"], eij -> ToExpression["eij"], e0ij -> ToExpression["e0ij"], euij -> ToExpression["euij"], i -> ToExpression["i"], j -> ToExpression["j"]}
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

BeginPackage["LINVARIANT`Fortran`"]

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
Expr2Fortran             ::usage = "Expr2Fortran[expr, FortranVarSub]"
FortranVarStr            ::usage = "FortranVarStr[var, arg]"

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
         {Fint[1] <> "implicit none"}, 
         {Fint[1] <> "Integer, parameter    :: snglprec = selected_real_kind(6, 37)  !< define precision for single reals"}, 
         {Fint[1] <> "Integer, parameter    :: dblprec = selected_real_kind(15, 307) !< define precision for double reals"}, 
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
     {Fint[2] <> "GridDriftx = Sum(Fields(1,3,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2] <> "GridDrifty = Sum(Fields(2,3,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2] <> "GridDriftz = Sum(Fields(3,3,:,:,:))/(NGridx*NGridy*NGridz)"}, 
     {Fint[2]}, 
     {Fint[2] <> "do iz = 1, NGridz"}, 
     {Fint[3] <> "do iy = 1, NGridy"}, 
     {Fint[4] <> "do ix = 1, NGridx"}, 
     {Fint[5] <> "Fields(1,3,ix,iy,iz) = Fields(1,3,ix,iy,iz) - GridDriftx"}, 
     {Fint[5] <> "Fields(2,3,ix,iy,iz) = Fields(2,3,ix,iy,iz) - GridDrifty"}, 
     {Fint[5] <> "Fields(3,3,ix,iy,iz) = Fields(3,3,ix,iy,iz) - GridDriftz"}, 
     {Fint[4] <> "end do"}, 
     {Fint[3] <> "end do"}, 
     {Fint[2] <> "end do"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Subroutine RemoveGridDrifts"},
     {Fint[1]},
     {Fint[1] <> "Function Thermometer(dFieldsdt, e0ij) Result(TK)"}, 
     {Fint[2] <> "Implicit none"},
     {Fint[2] <> "Real*8,  Intent(in) :: dFieldsdt(FieldDim, NumField, NGridx, NGridy, NGridz)"}, 
     {Fint[2] <> "Real*8,  Intent(in) :: e0ij(3,3)"}, 
     {Fint[2] <> "Real*8              :: TK, Ekin"}, 
     {Fint[2] <> "Integer             :: i, ifield, ix, iy, iz, TotalDim"}, 
     {Fint[2]}, 
     {Fint[2] <> "Ekin = 0.0D0"}, 
     {Fint[2] <> "TotalDim = NGridx*NGridy*NGridz*OnSiteDim"}, 
     {Fint[2]}, 
     {Fint[2] <> "do iz = 1, NGridz"}, 
     {Fint[3] <> "do iy = 1, NGridy"}, 
     {Fint[4] <> "do ix = 1, NGridx"}, 
     {Fint[5] <> "do ifield = 1, NumField"}, 
     {Fint[6] <> "do i = 1, FieldDimList(ifield)"}, 
     {Fint[7] <> "Ekin = Ekin + 0.5d0*mass(ifield)*dFieldsdt(i, ifield, ix, iy, iz)**2"}, 
     {Fint[6] <> "end do"}, 
     {Fint[5] <> "end do"}, 
     {Fint[4] <> "end do"}, 
     {Fint[3] <> "end do"}, 
     {Fint[2] <> "end do"}, 
     {Fint[2]}, 
     {Fint[2] <> "TK = 0.5D0*Ekin/TotalDim"}, 
     {Fint[2]}, 
     {Fint[1] <> "End Function Thermometer"},
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

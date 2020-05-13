
(* :Title: FortranTools.m -- A collection of tools for reading and
           writing Fortran-formatted files *)

(* :Context: Utilities`FortranTools` *)

(* :Author: Markus Lischka (mlischka@physik.tu-muenchen.de) *)

(* :Summary:
   The "FortranTools" package provides a set of tools for reading
   and writing Fortran-formatted files. It defines two main
   functions, FWrite and FReadList, that can handle integer, real, logical
   and string data.
*)

(* :Copyright: (c) 2000 by Markus Lischka

   Permission is granted to distribute this file for any purpose except 
   for inclusion in commercial software or program collections. 
   This copyright notice must remain intact.
*)

(* :Package Version: 1.2 *)

(* :Mathematica Version: 4.0 *)

(* :History:
   1.0 Initial version (2000-07-21).
   1.1 Update: Handle special Fortran numbers with missing "E" (2001-12-04).
   1.2 Added some documentation. Initial submission to MathSource (2003-03-30).
*)

(* :Keywords: Fortran, file input/output *)

(* :Sources:
   Thomas B. Bahder. Mathematica for Scientists and Engineers. 
       Addison-Wesley, 1995 (ch. 11).
   Roman E. Maeder. Programming in Mathematica, 3rd ed. Addison-Wesley, 1996.
*)

(* :Warnings:
   Package turns off the error message Read::readn.
*)

(* :Limitations:
   The field width used in FWrite is limited to 255 characters 
   (padString, errString).
*)

(* :Discussion:
   See some usage examples at the end of the package.
*)

(* :Requirements:
   None
*)


(* set up the package context, including public imports *)

BeginPackage["LINVARIANT`Fortran`"]

(* usage messages for the exported functions and the context itself *)

FortranTools::usage = "The \"FortranTools\" package provides
a set of tools for reading and writing Fortran-formatted files.
It defines two main functions, FWrite and FReadList, that can handle
integer, real, logical and string data.";

FWrite::usage = "FWrite[channel, descriptor, expr1, expr2, ... ]
writes the expressions expri in sequence to the specified output channel
with descriptor specifying the Fortran output format. The descriptor
\"*\" defaults to \"A\" for strings, \"L2\" for booleans, \"I15\" for
integers, and \"E15.7\" for reals. No newline is added to the output.";

FWriteArray::usage = "FWriteArray[channel, descriptor, matrix]
writes the two-dimensional array matrix to the specified output channel
with descriptor specifying the Fortran output format. The descriptor
\"*\" defaults to \"A\" for strings, \"L2\" for booleans, \"I15\" for
integers, and \"E15.7\" for reals. Rows are separated by newlines.";

ReadString::usage = "ReadString[string] reads one expression from string,
and returns the expression. ReadString[string, type] reads one object of the
specified type. ReadString[string, {type1, type2, ... }] reads a sequence of
objects of the specified types.";

ReadListString::usage = "ReadListString[string] reads all the remaining
expressions in string, and returns a list of them.
ReadListString[string, type] reads objects of the specified type from string,
until the end of the string is reached. The list of objects read is returned.
ReadListString[string, {type1, type2, ... }] reads objects with a sequence
of types, until the end of the string is reached.
ReadListString[string, types, n] reads only the first n objects	of the
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
   Read::readn is turned off, because otherwise ReadListString will warn
   about failed conversions. *)

Off[Read::readn];

fToNumber[s_String] :=
    (If[(Length[#2] == 1) && NumberQ[First[#2]], First[#2], #1])&[
        s, ReadListString[s, Number] ];

(* Code rewritten as pure function for 20% faster execution compared to:
fToNumber[s_String] := 
    Module[{result},
        result = ReadListString[s, Number];
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


(* ReadString/ReadListString: Input data reading from string variable. *)

ReadString[string_String, rest___] := 
    Module[{stream, result},
        stream = StringToStream[string];
        result = Read[stream, rest];
        Close[stream];
        result
    ];
ReadListString[string_String, rest___] := 
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
    FToExpression[ReadListString[#, Word, opts]]& /@ 
        ReadList[file, Record, n, RecordLists -> False, opts];


End[ ]         (* end the private context *)


(* protect exported symbols *)

Protect[
    FWrite, FWriteArray,
    ReadString, ReadListString,
    DropNonNumericElements,
    FToExpression,
    FReadList, FReadLine
    ];

EndPackage[ ]  (* end the package context *)


(* :Examples:

Needs["Utilities`FortranTools`"]

?FortranTools

?Utilities`FortranTools`*

?FWrite

(* unsupported data type gives an error message *)
FWrite["demo.txt", "*", Null]

FWrite["demo.txt", "*", N[Pi, 8], 13!, True, " String"]
!! "demo.txt"

(* overflow in last column *)
FWrite["demo.txt", "I4", Table[i!, {i, 0, 8}]]
!! "demo.txt"

FWrite["demo.txt", "*", Table[Random[], {3}, {3}]]
!! "demo.txt"

?FWriteArray

FWriteArray["demo.txt", "F10.5", Table[Random[], {3}, {3}]];
!! "demo.txt"

stream = OpenWrite["demo.txt"];
FWriteArray[stream, "*", Table[Random[], {3}, {3}]];
FWriteArray[stream, "L14", Table[EvenQ[i + j], {i, 3}, {j, 3}]];
Close[stream];
!! "demo.txt"

?FToExpression

FToExpression[{{"1", "1.0", "1.0E+00"}, {".FALSE.", "F", "String"}, 
{".TRUE.", "T", "String"}}]

?DropNonNumericElements

DropNonNumericElements[%%]

?FReadList

FReadList["demo.txt"]

?FReadLine

FReadLine["demo.txt"]

FReadLine["demo.txt", 3]

*)

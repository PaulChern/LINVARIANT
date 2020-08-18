#!/usr/bin/env wolframscript

tmp = URLSave["https://github.com/PaulChern/LINVARIANT/archive/master.zip"];
dest = FileNameJoin[{$UserBaseDirectory, "Peng", "LINVARIANT"}];
tmpExpanded = CreateDirectory[];
ExtractArchive[tmp, tmpExpanded];
CopyDirectory[FileNameJoin[{tmpExpanded, "LINVARIANT-master", "LINVARIANT"}], dest];
DeleteFile[tmp];
DeleteDirectory[tmpExpanded, DeleteContents -> True];
Print["LINVARIANT at current stage preserve all the rights to Peng Chen!"]
Print["Installed LINVARIANT successfully!"]
Print["Please restart Mathematica or the kernel."]

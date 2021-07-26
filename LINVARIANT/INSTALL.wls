#!/usr/bin/env wolframscript
tmp = URLSave["https://github.com/PaulChern/LINVARIANT/archive/refs/tags/0.9-beta.zip"];
dest = FileNameJoin[{$UserBaseDirectory, "Applications", "LINVARIANT"}];
tmpExpanded = CreateDirectory[];
ExtractArchive[tmp, tmpExpanded];
CopyDirectory[FileNameJoin[{tmpExpanded, "LINVARIANT-0.9-beta", "LINVARIANT"}], dest];
DeleteFile[tmp];
DeleteDirectory[tmpExpanded, DeleteContents -> True];
Print["LINVARIANT at current stage preserve all the rights to Peng Chen!"]
Print["Use as your own risk!"]
Print["Installed LINVARIANT successfully!"]
Print["Please restart Mathematica or the kernel."]

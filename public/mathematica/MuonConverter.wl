(* ::Package:: *)

(* Initialization file for the package MuonConverter` *)

(* SetEnvironment["MU2E_ELASTIC"-> NotebookDirectory[]<>"/Elastic"]  - Not needed, Mu2eNRET finds automatically *)
(* root is above the 4 repository directories *)
root=ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
mu2eNRETDir = FileNameJoin[{root, "Mu2e_NRET", "v2", "mathematica"}]

Get["Mu2eNRET`", Path->mu2eNRETDir];
Get["hadronization`"];
Get["parameters`"];
Get["formfactors`"];
Get["interface`"];

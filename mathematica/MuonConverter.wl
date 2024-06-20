(* 
  ::  MuonConverter.wl is a part of the MuonConverter package.                      
  ::  Copyright (C) 2024 MuonConverter authors (see AUTHORS for details).
  ::  MuonConverter is licenced under the BSD 3-Clause, see LICENSE for details.
*)

(* Initialization file for MuonConverter` *)

(* Fix the relative paths for Mu2e_NRET  *)
MuonBridgePATH = ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]];
Mu2eNRETPATH = FileNameJoin[{MuonBridgePATH, "Mu2e_NRET", "v2", "mathematica"}]

Get["Mu2eNRET`", Path -> Mu2eNRETPATH];
Get["hadronization`"];
Get["parameters`"];
Get["formfactors`"];
Get["interface`"];
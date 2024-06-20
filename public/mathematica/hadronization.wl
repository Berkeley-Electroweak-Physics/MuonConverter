BeginPackage["MuonConverter`"];

(* HadronizeWET *)
HadronizeWET::usage = "Hadronize matches or 'hadronizes' or the EFT of relativistic quarks to an EFT with relativistic nucleons";

Begin["`Private`"];

HadronizeWET[WCs_, q_] := Module[{WCkeys, CmueqqAssoc, RETnucleonAssoc},
(* Supported keys - see __insert arXiv no.__ *)
WCkeys = {"C51", "C52", 
          "C61u", "C61d", "C61s", "C62u", "C62d", "C62s",
          "C63u", "C63d", "C63s", "C64u", "C64d", "C64s",
          "C65u", "C65d", "C65s", "C66u", "C66d", "C66s",
          "C67u", "C67d", "C67s", "C68u", "C68d", "C68s",
          "C69u", "C69d", "C69s", "C610u", "C610d", "C610s",
          "C71", "C72", "C73", "C74", "C75", "C76", "C77", "C78",
          "C79u", "C79d", "C79s", "C710u", "C710d", "C710s",
          "C711u", "C711d", "C711s", "C712u", "C712d", "C712s",
          "C713u", "C713d", "C713s", "C714u", "C714d", "C714s",
          "C715u", "C715d", "C715s", "C716u", "C716d", "C716s"};

(* Initialize WC association *)
Zerofication[x_]:= 0.;
CmueqqAssoc = AssociationMap[Zerofication, WCkeys];

(* Initialize an association with all possible keys 
   and fill specified keys with assigned values *)
CmueqqAssoc = KeyDrop[CmueqqAssoc, Keys[WCs]];
CmueqqAssoc = AssociateTo[CmueqqAssoc, WCs];

(* Initialize form factors *)
Do[F1[qq,nn]  = FormFactor[q, "1", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[F2[qq,nn]  = FormFactor[q, "2", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[F3[qq,nn]  = FormFactor[q, "3", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FA[qq,nn]  = FormFactor[q, "A", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FPp[qq,nn] = FormFactorPprimed[q, qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FS[qq,nn]  = FormFactor[q, "S", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FP[qq,nn]  = FormFactorP[q, qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FA3[qq,nn] = FormFactor[q, "A3", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT0[qq,nn] = FormFactor[q, "T0", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT1[qq,nn] = FormFactor[q, "T1", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT2[qq,nn] = FormFactor[q, "T2", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT3[qq,nn] = FormFactor[q, "T3", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FG[nn]     = FormFactorNucleon[q,"FG",nn], {nn,{"p","n"}}];
Do[FGt[nn]    = FormFactorNucleon[q,"FGt",nn], {nn,{"p","n"}}];
Do[FF[nn]     = FormFactorNucleon[q,"FF",nn], {nn,{"p","n"}}];
Do[FFt[nn]    = FormFactorNucleon[q,"FFt",nn], {nn,{"p","n"}}];

(* Override defaults if provided in parameters.wl *)
If[F1up != None, F1["u","p"] = F1up];
If[F1dp != None, F1["d","p"] = F1dp];
If[F1sp != None, F1["s","p"] = F1sp];
If[F1un != None, F1["u","n"] = F1un];
If[F1dn != None, F1["d","n"] = F1dn];
If[F1sn != None, F1["s","n"] = F1sn];

If[F2up != None, F2["u","p"] = F2up];
If[F2dp != None, F2["d","p"] = F2dp];
If[F2sp != None, F2["s","p"] = F2sp];
If[F2un != None, F2["u","n"] = F2un];
If[F2dn != None, F2["d","n"] = F2dn];
If[F2sn != None, F2["s","n"] = F2sn];

If[F3up != None, F3["u","p"] = F3up];
If[F3dp != None, F3["d","p"] = F3dp];
If[F3sp != None, F3["s","p"] = F3sp];
If[F3un != None, F3["u","n"] = F3un];
If[F3dn != None, F3["d","n"] = F3dn];
If[F3sn != None, F3["s","n"] = F3sn];

If[FAup != None, FA["u","p"] = FAup];
If[FAdp != None, FA["d","p"] = FAdp];
If[FAsp != None, FA["s","p"] = FAsp];
If[FAun != None, FA["u","n"] = FAun];
If[FAdn != None, FA["d","n"] = FAdn];
If[FAsn != None, FA["s","n"] = FAsn];

If[FPpup != None, FPp["u","p"] = FPpup];
If[FPpdp != None, FPp["d","p"] = FPpdp];
If[FPpsp != None, FPp["s","p"] = FPpsp];
If[FPpun != None, FPp["u","n"] = FPpun];
If[FPpdn != None, FPp["d","n"] = FPpdn];
If[FPpsn != None, FPp["s","n"] = FPpsn];

If[FSup != None, FS["u","p"] = FSup];
If[FSdp != None, FS["d","p"] = FSdp];
If[FSsp != None, FS["s","p"] = FSsp];
If[FSun != None, FS["u","n"] = FSun];
If[FSdn != None, FS["d","n"] = FSdn];
If[FSsn != None, FS["s","n"] = FSsn];

If[FPup != None, FP["u","p"] = FPup];
If[FPdp != None, FP["d","p"] = FPdp];
If[FPsp != None, FP["s","p"] = FPsp];
If[FPun != None, FP["u","n"] = FPun];
If[FPdn != None, FP["d","n"] = FPdn];
If[FPsn != None, FP["s","n"] = FPsn];

If[FT0up != None, FT0["u","p"] = FT0up];
If[FT0dp != None, FT0["d","p"] = FT0dp];
If[FT0sp != None, FT0["s","p"] = FT0sp];
If[FT0un != None, FT0["u","n"] = FT0un];
If[FT0dn != None, FT0["d","n"] = FT0dn];
If[FT0sn != None, FT0["s","n"] = FT0sn];

If[FT1up != None, FT1["u","p"] = FT1up];
If[FT1dp != None, FT1["d","p"] = FT1dp];
If[FT1sp != None, FT1["s","p"] = FT1sp];
If[FT1un != None, FT1["u","n"] = FT1un];
If[FT1dn != None, FT1["d","n"] = FT1dn];
If[FT1sn != None, FT1["s","n"] = FT1sn];

If[FT2up != None, FT2["u","p"] = FT2up];
If[FT2dp != None, FT2["d","p"] = FT2dp];
If[FT2sp != None, FT2["s","p"] = FT2sp];
If[FT2un != None, FT2["u","n"] = FT2un];
If[FT2dn != None, FT2["d","n"] = FT2dn];
If[FT2sn != None, FT2["s","n"] = FT2sn];

If[FT3up != None, FT3["u","p"] = FT3up];
If[FT3dp != None, FT3["d","p"] = FT3dp];
If[FT3sp != None, FT3["s","p"] = FT3sp];
If[FT3un != None, FT3["u","n"] = FT3un];
If[FT3dn != None, FT3["d","n"] = FT3dn];
If[FT3sn != None, FT3["s","n"] = FT3sn];

If[FGp != None, FG["p"] = FGp];
If[FGn != None, FG["n"] = FGn];

If[FGtp != None, FGt["p"] = FGtp];
If[FGtn != None, FGt["n"] = FGtn];

If[FFp != None, FF["p"] = FFp];
If[FFn != None, FF["n"] = FFn];

If[FFtp != None, FFt["p"] = FFtp];
If[FFtn != None, FFt["n"] = FFtn];


(* Define the momentum transfer q^2 for select form factors *)
qSq = q[[1]] * q[[1]] - q[[2]] * q[[2]] - q[[3]] * q[[3]] - q[[4]] * q[[4]];

(* Hadronize to the nucleon RET basis *)
RETnucleonAssoc = <| "d1p" -> + (1. / mu) * CmueqqAssoc["C65u"] * FS["u","p"] + (1. / md) * CmueqqAssoc["C65d"] * FS["d","p"] + (1. / ms) * CmueqqAssoc["C65s"] * FS["s","p"]
							  + CmueqqAssoc["C75"] * FF["p"] + CmueqqAssoc["C71"] * FG["p"]
							  (* T-odd contribution (F3 is set to zero by default) *)
							  - I * (mMinus / mp) * (CmueqqAssoc["C61u"] * F3["u","p"] + CmueqqAssoc["C61d"] * F3["d","p"] + CmueqqAssoc["C61s"] * F3["s","p"])\
							  - I * (mmu**2 - me**2) * (CmueqqAssoc["C69u"] * F3["u","p"] + CmueqqAssoc["C69d"] * F3["d","p"] + CmueqqAssoc["C69s"] * F3["s","p"]),

					 "d1n" -> + (1. / mu) * CmueqqAssoc["C65u"] * FS["u","n"] + (1. / md) * CmueqqAssoc["C65d"] * FS["d","n"] + (1. / ms) * CmueqqAssoc["C65s"] * FS["s","n"]
							  + CmueqqAssoc["C75"] * FF["n"] + CmueqqAssoc["C71"] * FG["n"]
							  (* T-odd contribution (F3 is set to zero by default) *)
							  - I * (mMinus / mn) * (CmueqqAssoc["C61u"] * F3["u","n"] + CmueqqAssoc["C61d"] * F3["d","n"] + CmueqqAssoc["C61s"] * F3["s","n"])\
							  - I * (mmu**2 - me**2) * (CmueqqAssoc["C69u"] * F3["u","n"] + CmueqqAssoc["C69d"] * F3["d","n"] + CmueqqAssoc["C69s"] * F3["s","n"]),
				
					 "d2p" -> + (1. / mu) * CmueqqAssoc["C67u"] * FP["u","p"] + (1. / md) * CmueqqAssoc["C67d"] * FP["d","p"] + (1. / ms) * CmueqqAssoc["C67s"] * FP["s","p"]
							  + CmueqqAssoc["C73"] * FGt["p"] + CmueqqAssoc["C77"] * FFt["p"]
							  - I * (mMinus / (2 * mp)) * (CmueqqAssoc["C63u"] * FPp["u","p"] + CmueqqAssoc["C63d"] * FPp["d","p"] + CmueqqAssoc["C63s"] * FPp["s","p"])
							  - I * ((mPlus * mMinus) / (2 * mp)) * (CmueqqAssoc["C711u"] * FPp["u","p"] + CmueqqAssoc["C711d"] * FPp["d","p"] + CmueqqAssoc["C711s"] * FPp["s","p"])
							  (* T-odd contribution (F3 and FT3 are set to zero by default) *)
							  - 4 * I * mMinus * ((1 / mu) * CmueqqAssoc["C715u"] * FT3["u","p"] + (1 / md) * CmueqqAssoc["C715d"] * FT3["d","p"] + (1 / ms) * CmueqqAssoc["C715s"] * FT3["s","p"])
							  + (mPlus / mp) * (CmueqqAssoc["C62u"] * F3["u", "p"] + CmueqqAssoc["C62d"] * F3["d", "p"] + CmueqqAssoc["C62s"] * F3["s", "p"]),
				
					 "d2n" -> + (1. / mu) * CmueqqAssoc["C67u"] * FP["u","n"] + (1. / md) * CmueqqAssoc["C67d"] * FP["d","n"] + (1. / ms) * CmueqqAssoc["C67s"] * FP["s","n"]
							  + CmueqqAssoc["C73"] * FGt["n"] + CmueqqAssoc["C77"] * FFt["n"]
							  - I * (mMinus / (2 * mp)) * (CmueqqAssoc["C63u"] * FPp["u","n"] + CmueqqAssoc["C63d"] * FPp["d","n"] + CmueqqAssoc["C63s"] * FPp["s","n"])
							  - I * ((mPlus * mMinus) / (2 * mp)) * (CmueqqAssoc["C711u"] * FPp["u","n"] + CmueqqAssoc["C711d"] * FPp["d","n"] + CmueqqAssoc["C711s"] * FPp["s","n"])
							  (* T-odd contribution (F3 and FT3 are set to zero by default) *)
							  - 4 * I * mMinus * ((1 / mu) * CmueqqAssoc["C715u"] * FT3["u","n"] + (1 / md) * CmueqqAssoc["C715d"] * FT3["d","n"] + (1 / ms) * CmueqqAssoc["C715s"] * FT3["s","n"])
							  + (mPlus / mn) * (CmueqqAssoc["C62u"] * F3["u", "n"] + CmueqqAssoc["C62d"] * F3["d", "n"] + CmueqqAssoc["C62s"] * F3["s", "n"]),
				
					"d3p" -> + (1. / mu) * CmueqqAssoc["C66u"] * FS["u","p"] + (1. / md) * CmueqqAssoc["C66d"] * FS["d","p"] + (1. / ms) * CmueqqAssoc["C66s"] * FS["s","p"]
							 + CmueqqAssoc["C72"] * FG["p"] + CmueqqAssoc["C76"] * FF["p"]
							 (* T-odd contribution (F3 is set to zero by default) *)
							 - I * (mmu^2 - me^2) * (CmueqqAssoc["C710u"] * F3["u","p"] + CmueqqAssoc["C710d"] * F3["d","p"] + CmueqqAssoc["C710s"] * F3["s","p"]),
				
					"d3n" -> + (1. / mu) * CmueqqAssoc["C66u"] * FS["u","n"] + (1. / md) * CmueqqAssoc["C66d"] * FS["d","n"] + (1. / ms) * CmueqqAssoc["C66s"] * FS["s","n"]
							 + CmueqqAssoc["C72"] * FG["n"] + CmueqqAssoc["C76"] * FF["n"]
							 (* T-odd contribution (F3 is set to zero by default) *)
							 - I * (mmu^2 - me^2) * (CmueqqAssoc["C710u"] * F3["u","n"] + CmueqqAssoc["C710d"] * F3["d","n"] + CmueqqAssoc["C710s"] * F3["s","n"]),
				
					"d4p" -> + (1. / mu) * CmueqqAssoc["C68u"] * FP["u","p"] + (1. / md) * CmueqqAssoc["C68d"] * FP["d","p"] + (1. / ms) * CmueqqAssoc["C68s"] * FP["s","p"]
							 + CmueqqAssoc["C74"] * FGt["p"] + CmueqqAssoc["C78"] * FFt["p"]
							 + (mPlus / (2 * mp)) * (CmueqqAssoc["C64u"] * FPp["u","p"] + CmueqqAssoc["C64d"] * FPp["d","p"] + CmueqqAssoc["C64s"] * FPp["s","p"])
							 - I * ((mPlus * mMinus) / (2 * mp)) * (CmueqqAssoc["C712u"] * FPp["u","p"] + CmueqqAssoc["C712d"] * FPp["d","p"] + CmueqqAssoc["C712s"] * FPp["s","p"])
							 (* T-odd contribution (FT3 is set to zero by default) *)
							 + 4 * mPlus * (CmueqqAssoc["C716u"] * FT3["u","p"] + CmueqqAssoc["C716d"] * FT3["d","p"] + CmueqqAssoc["C716s"] * FT3["s","p"]),
				
					"d4n" -> + (1. / mu) * CmueqqAssoc["C68u"] * FP["u","n"] + (1. / md) * CmueqqAssoc["C68d"] * FP["d","n"] + (1. / ms) * CmueqqAssoc["C68s"] * FP["s","n"]
							 + CmueqqAssoc["C74"] * FGt["n"] + CmueqqAssoc["C78"] * FFt["n"]
							 + (mPlus / (2 * mn)) * (CmueqqAssoc["C64u"] * FPp["u","n"] + CmueqqAssoc["C64d"] * FPp["d","n"] + CmueqqAssoc["C64s"] * FPp["s","n"])
							 - I * ((mPlus * mMinus) / (2 * mn)) * (CmueqqAssoc["C712u"] * FPp["u","n"] + CmueqqAssoc["C712d"] * FPp["d","n"] + CmueqqAssoc["C712s"] * FPp["s","n"])
							 (* T-odd contribution (FT3 is set to zero by default) *)
							 + 4 * mPlus * (CmueqqAssoc["C716u"] * FT3["u","n"] + CmueqqAssoc["C716d"] * FT3["d","n"] + CmueqqAssoc["C716s"] * FT3["s","n"]),
				
					"d5p" -> + CmueqqAssoc["C61u"] * F1["u","p"] + CmueqqAssoc["C61d"] * F1["d","p"] + CmueqqAssoc["C61s"] * F1["s","p"]
							 + mPlus * (CmueqqAssoc["C79u"] * F1["u","p"] + CmueqqAssoc["C79d"] * F1["d","p"] + CmueqqAssoc["C79s"] * F1["s","p"])
							 - (qSq / (2 * mp)) * (CmueqqAssoc["C713u"] * (FT1["u","p"] - 4 * FT2["u","p"]) + CmueqqAssoc["C713d"] * (FT1["d","p"] - 4 * FT2["d","p"]) + CmueqqAssoc["C713s"] * (FT1["s","p"] - 4 * FT2["s","p"])),
		
					"d5n" -> + CmueqqAssoc["C61u"] * F1["u","n"] + CmueqqAssoc["C61d"] * F1["d","n"] + CmueqqAssoc["C61s"] * F1["s","n"]
							 + mPlus * (CmueqqAssoc["C79u"] * F1["u","n"] + CmueqqAssoc["C79d"] * F1["d","n"] + CmueqqAssoc["C79s"] * F1["s","n"])
							 - (qSq / (2 * mn)) * (CmueqqAssoc["C713u"] * (FT1["u","n"] - 4 * FT2["u","n"]) + CmueqqAssoc["C713d"] * (FT1["d","n"] - 4 * FT2["d","n"]) + CmueqqAssoc["C713s"] * (FT1["s","n"] - 4 * FT2["s","n"])),
				
					"d6p" -> - (1./2.) * (CmueqqAssoc["C61u"] * F2["u","p"] + CmueqqAssoc["C61d"] * F2["d","p"] + CmueqqAssoc["C61s"] * F2["s","p"])
							 - (1./2.) * mPlus * (CmueqqAssoc["C79u"] * F2["u","p"] + CmueqqAssoc["C79d"] * F2["d","p"] + CmueqqAssoc["C79s"] * F2["s","p"]),
		
					"d6n" -> - (1./2.) * (CmueqqAssoc["C61u"] * F2["u","n"] + CmueqqAssoc["C61d"] * F2["d","n"] + CmueqqAssoc["C61s"] * F2["s","n"])
							 - (1./2.) * mPlus * (CmueqqAssoc["C79u"] * F2["u","n"] + CmueqqAssoc["C79d"] * F2["d","n"] + CmueqqAssoc["C79s"] * F2["s","n"]),
		
					"d7p" -> + CmueqqAssoc["C63u"] * FA["u","p"] + CmueqqAssoc["C63d"] * FA["d","p"] + CmueqqAssoc["C63s"] * FA["s","p"]
							 + mPlus * (CmueqqAssoc["C711u"] * FA["u","p"] + CmueqqAssoc["C711d"] * FA["d","p"] + CmueqqAssoc["C711s"] * FA["s","p"])
							 (* T-odd contribution (FT3 is set to zero by default) *)
							 - 2 * I * (qSq / mp) * (CmueqqAssoc["C715u"] * FT3["u","p"] + CmueqqAssoc["C715d"] * FT3["d","p"] + CmueqqAssoc["C715s"] * FT3["s","p"]),
		
					"d7n" -> + CmueqqAssoc["C63u"] * FA["u","n"] + CmueqqAssoc["C63d"] * FA["d","n"] + CmueqqAssoc["C63s"] * FA["s","n"]
							 + mPlus * (CmueqqAssoc["C711u"] * FA["u","n"] + CmueqqAssoc["C711d"] * FA["d","n"] + CmueqqAssoc["C711s"] * FA["s","n"])
							 (* T-odd contribution (FT3 is set to zero by default) *)
							 - 2 * I * (qSq / mn) * (CmueqqAssoc["C715u"] * FT3["u","n"] + CmueqqAssoc["C715d"] * FT3["d","n"] + CmueqqAssoc["C715s"] * FT3["s","n"]),
					
							 (* T-odd contribution (FA3 is set to zero by default) *)
					"d8p" -> + CmueqqAssoc["C63u"] * FA3["u","p"] + CmueqqAssoc["C63d"] * FA3["d","p"] + CmueqqAssoc["C63s"] * FA3["s","p"]
							 + mPlus * (CmueqqAssoc["C711u"] * FA3["u","p"] + CmueqqAssoc["C711d"] * FA3["d","p"] + CmueqqAssoc["C711s"] * FA3["s","p"]),

							 (* T-odd contribution (FA3 is set to zero by default) *)
					"d8n" -> + CmueqqAssoc["C63u"] * FA3["u","n"] + CmueqqAssoc["C63d"] * FA3["d","n"] + CmueqqAssoc["C63s"] * FA3["s","n"]
							 + mPlus * (CmueqqAssoc["C711u"] * FA3["u","n"] + CmueqqAssoc["C711d"] * FA3["d","n"] + CmueqqAssoc["C711s"] * FA3["s","n"]),
				
					"d9p" -> - (\[Alpha] / \[Pi]) * CmueqqAssoc["C51"] * (mL / qSq) * (qu * F1["u","p"] + qd * F1["d","p"] + qs * F1["s","p"])
							 - mL * (CmueqqAssoc["C79u"] * F1["u","p"] + CmueqqAssoc["C79d"] * F1["d","p"] + CmueqqAssoc["C79s"] * F1["s","p"]),
		
					"d9n" -> - (\[Alpha] / \[Pi]) * CmueqqAssoc["C51"] * (mL / qSq) * (qu * F1["u","n"] + qd * F1["d","n"] + qs * F1["s","n"])
							 - mL * (CmueqqAssoc["C79u"] * F1["u","n"] + CmueqqAssoc["C79d"] * F1["d","n"] + CmueqqAssoc["C79s"] * F1["s","n"]),
		
					"d10p" -> + (\[Alpha] / (2 * \[Pi])) * CmueqqAssoc["C51"] * (mL / qSq) * (qu * F2["u","p"] + qd * F2["d","p"] + qs * F2["s","p"])
							  + (mL / 2.) * (CmueqqAssoc["C79u"] * F2["u","p"] + CmueqqAssoc["C79d"] * F2["d","p"] + CmueqqAssoc["C79s"] * F2["s","p"]),
		
					"d10n" -> + (\[Alpha] / (2 * \[Pi])) * CmueqqAssoc["C51"] * (mL / qSq) * (qu * F2["u","n"] + qd * F2["d","n"] + qs * F2["s","n"])
							  + (mL / 2.) * (CmueqqAssoc["C79u"] * F2["u","n"] + CmueqqAssoc["C79d"] * F2["d","n"] + CmueqqAssoc["C79s"] * F2["s","n"]),
		
					"d11p" -> - mL * (CmueqqAssoc["C711u"] * FA["u","p"] + CmueqqAssoc["C711d"] * FA["d","p"] + CmueqqAssoc["C711s"] * FA["s","p"]),
		
					"d11n" -> - mL * (CmueqqAssoc["C711u"] * FA["u","n"] + CmueqqAssoc["C711d"] * FA["d","n"] + CmueqqAssoc["C711s"] * FA["s","n"]),
		
							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d12p" -> - mL * (CmueqqAssoc["C711u"] * FA3["u","p"] + CmueqqAssoc["C711d"] * FA3["d","p"] + CmueqqAssoc["C711s"] * FA3["s","p"]),

							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d12n" -> - mL * (CmueqqAssoc["C711u"] * FA3["u","n"] + CmueqqAssoc["C711d"] * FA3["d","n"] + CmueqqAssoc["C711s"] * FA3["s","n"]),
		
					"d13p" -> + CmueqqAssoc["C62u"] * F1["u","p"] + CmueqqAssoc["C62d"] * F1["d","p"] + CmueqqAssoc["C62s"] * F1["s","p"]
							  - I * Mminus * (CmueqqAssoc["C710u"] * F1["u","p"] + CmueqqAssoc["C710d"] * F1["d","p"] + CmueqqAssoc["C710s"] * F1["s","p"])
							  - (qSq / (2 * mp)) * (CmueqqAssoc["C714u"] * (FT1["u","p"] - 4 * FT2["u","p"]) + CmueqqAssoc["C714d"] * (FT1["d","p"] - 4 * FT2["d","p"]) + CmueqqAssoc["C714s"] * (FT1["s","p"] - 4 * FT2["s","p"])),
		
					"d13n" -> + CmueqqAssoc["C62u"] * F1["u","n"] + CmueqqAssoc["C62d"] * F1["d","n"] + CmueqqAssoc["C62s"] * F1["s","n"]
							  - I * Mminus * (CmueqqAssoc["C710u"] * F1["u","n"] + CmueqqAssoc["C710d"] * F1["d","n"] + CmueqqAssoc["C710s"] * F1["s","n"])
							  - (qSq / (2 * mn)) * (CmueqqAssoc["C714u"] * (FT1["u","n"] - 4 * FT2["u","n"]) + CmueqqAssoc["C714d"] * (FT1["d","n"] - 4 * FT2["d","n"]) + CmueqqAssoc["C714s"] * (FT1["s","n"] - 4 * FT2["s","n"])),
				
					"d14p" -> - (1./2.) * (CmueqqAssoc["C62u"] * F2["u","p"] + CmueqqAssoc["C62d"] * F2["d","p"] + CmueqqAssoc["C62s"] * F2["s","p"])
							  + (I / 2.) * mMinus * (CmueqqAssoc["C710u"] * F2["u","p"] + CmueqqAssoc["C710d"] * F2["d","p"] + CmueqqAssoc["C710s"] * F2["s","p"]),
		
					"d14n" -> - (1./2.) * (CmueqqAssoc["C62u"] * F2["u","n"] + CmueqqAssoc["C62d"] * F2["d","n"] + CmueqqAssoc["C62s"] * F2["s","n"])
							  + (I / 2.) * mMinus * (CmueqqAssoc["C710u"] * F2["u","n"] + CmueqqAssoc["C710d"] * F2["d","n"] + CmueqqAssoc["C710s"] * F2["s","n"]),
		
					"d15p" -> + CmueqqAssoc["C64u"] * FA["u","p"] + CmueqqAssoc["C64d"] * FA["d","p"] + CmueqqAssoc["C64s"] * FA["s","p"]
							  - I * mMinus * (CmueqqAssoc["C712u"] * FA["u","p"] + CmueqqAssoc["C712d"] * FA["d","p"] + CmueqqAssoc["C712s"] * FA["s","p"])
							  (* T-odd contribution (FT3 is set to zero by default) *)
							  - 2 * (qSq / mp) * (CmueqqAssoc["C716u"] * FT3["u","p"] + CmueqqAssoc["C716d"] * FT3["d","p"] + CmueqqAssoc["C716s"] * FT3["s","p"]),
		
					"d15n" -> + CmueqqAssoc["C64u"] * FA["u","n"] + CmueqqAssoc["C64d"] * FA["d","n"] + CmueqqAssoc["C64s"] * FA["s","n"]
							  - I * mMinus * (CmueqqAssoc["C712u"] * FA["u","n"] + CmueqqAssoc["C712d"] * FA["d","n"] + CmueqqAssoc["C712s"] * FA["s","n"])
							  (* T-odd contribution (FT3 is set to zero by default) *)
							  - 2 * (qSq / mn) * (CmueqqAssoc["C716u"] * FT3["u","n"] + CmueqqAssoc["C716d"] * FT3["d","n"] + CmueqqAssoc["C716s"] * FT3["s","n"]),
		
							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d16p" -> + CmueqqAssoc["C64u"] * FA3["u","p"] + CmueqqAssoc["C64d"] * FA3["d","p"] + CmueqqAssoc["C64s"] * FA3["s","p"]
							  - I * mMinus * (CmueqqAssoc["C712u"] * FA3["u","p"] + CmueqqAssoc["C712d"] * FA3["d","p"] + CmueqqAssoc["C712s"] * FA3["s","p"]),

							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d16n" -> + CmueqqAssoc["C64u"] * FA3["u","n"] + CmueqqAssoc["C64d"] * FA3["d","n"] + CmueqqAssoc["C64s"] * FA3["s","n"]
							  - I * mMinus * (CmueqqAssoc["C712u"] * FA3["u","n"] + CmueqqAssoc["C712d"] * FA3["d","n"] + CmueqqAssoc["C712s"] * FA3["s","n"]),
		
					"d17p" -> + (\[Alpha] / \[Pi]) * CmueqqAssoc["C52"] * (mL / qSq) * (qu * F1["u","p"] + qd * F1["d","p"] + qs * F1["s","p"])
							  + mL * (CmueqqAssoc["C710u"] * F1["u","p"] + CmueqqAssoc["C710d"] * F1["d","p"] + CmueqqAssoc["C710s"] * F1["s","p"]),
		
					"d17n" -> + (\[Alpha] / \[Pi]) * CmueqqAssoc["C52"] * (mL / qSq) * (qu * F1["u","n"] + qd * F1["d","n"] + qs * F1["s","n"])
							  + mL * (CmueqqAssoc["C710u"] * F1["u","n"] + CmueqqAssoc["C710d"] * F1["d","n"] + CmueqqAssoc["C710s"] * F1["s","n"]),
		
					"d18p" -> - (\[Alpha] / (2 * \[Pi])) * CmueqqAssoc["C52"] * (mL / qSq) * (qu * F2["u","p"] + qd * F2["d","p"] + qs * F2["s","p"])
							  - (mL / 2.) * (CmueqqAssoc["C710u"] * F2["u","p"] + CmueqqAssoc["C710d"] * F2["d","p"] + CmueqqAssoc["C710s"] * F2["s","p"]),
		
					"d18n" -> - (\[Alpha] / (2 * \[Pi])) * CmueqqAssoc["C52"] * (mL / qSq) * (qu * F2["u","n"] + qd * F2["d","n"] + qs * F2["s","n"])
							  - (mL / 2.) * (CmueqqAssoc["C710u"] * F2["u","n"] + CmueqqAssoc["C710d"] * F2["d","n"] + CmueqqAssoc["C710s"] * F2["s","n"]),
		
					"d19p" -> + mL * (CmueqqAssoc["C712u"] * FA["u","p"] + CmueqqAssoc["C712d"] * FA["d","p"] + CmueqqAssoc["C712s"] * FA["s","p"]),
		
					"d19n" -> + mL * (CmueqqAssoc["C712u"] * FA["u","n"] + CmueqqAssoc["C712d"] * FA["d","n"] + CmueqqAssoc["C712s"] * FA["s","n"]),
				
							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d20p" -> + mL * (CmueqqAssoc["C712u"] * FA3["u","p"] + CmueqqAssoc["C712d"] * FA3["d","p"] + CmueqqAssoc["C712s"] * FA3["s","p"]),

							  (* T-odd contribution (FA3 is set to zero by default) *)
					"d20n" -> + mL * (CmueqqAssoc["C712u"] * FA3["u","n"] + CmueqqAssoc["C712d"] * FA3["d","n"] + CmueqqAssoc["C712s"] * FA3["s","n"]),
		
					"d21p" -> + CmueqqAssoc["C69u"] * FT0["u","p"] + CmueqqAssoc["C69d"] * FT0["d","p"] + CmueqqAssoc["C69s"] * FT0["s","p"],
		
					"d21n" -> + CmueqqAssoc["C69u"] * FT0["u","n"] + CmueqqAssoc["C69d"] * FT0["d","n"] + CmueqqAssoc["C69s"] * FT0["s","n"],
		
					"d22p" -> - (CmueqqAssoc["C69u"] * FT1["u","p"] + CmueqqAssoc["C69d"] * FT1["d","p"] + CmueqqAssoc["C69s"] * FT1["s","p"]),
		
					"d22n" -> - (CmueqqAssoc["C69u"] * FT1["u","n"] + CmueqqAssoc["C69d"] * FT1["d","n"] + CmueqqAssoc["C69s"] * FT1["s","n"]),
		
					"d23p" -> - (CmueqqAssoc["C69u"] * FT2["u","p"] + CmueqqAssoc["C69d"] * FT2["d","p"] + CmueqqAssoc["C69s"] * FT2["s","p"]),
		
					"d23n" -> - (CmueqqAssoc["C69u"] * FT2["u","n"] + CmueqqAssoc["C69d"] * FT2["d","n"] + CmueqqAssoc["C69s"] * FT2["s","n"]),
		
							  (* T-odd contribution (FT3 set to zero by default) *)
					"d24p" -> - (CmueqqAssoc["C69u"] * FT3["u","p"] + CmueqqAssoc["C69d"] * FT3["d","p"] + CmueqqAssoc["C69s"] * FT3["s","p"]),
					
							  (* T-odd contribution (FT3 set to zero by default) *)
					"d24n" -> - (CmueqqAssoc["C69u"] * FT3["u","n"] + CmueqqAssoc["C69d"] * FT3["d","n"] + CmueqqAssoc["C69s"] * FT3["s","n"]),
		
					"d25p" -> + CmueqqAssoc["C610u"] * FT0["u","p"] + CmueqqAssoc["C610d"] * FT0["d","p"] + CmueqqAssoc["C610s"] * FT0["s","p"],
		
					"d25n" -> + CmueqqAssoc["C610u"] * FT0["u","n"] + CmueqqAssoc["C610d"] * FT0["d","n"] + CmueqqAssoc["C610s"] * FT0["s","n"],
		
					"d26p" -> - (CmueqqAssoc["C610u"] * FT1["u","p"] + CmueqqAssoc["C610d"] * FT1["d","p"] + CmueqqAssoc["C610s"] * FT1["s","p"]),
		
					"d26n" -> - (CmueqqAssoc["C610u"] * FT1["u","n"] + CmueqqAssoc["C610d"] * FT1["d","n"] + CmueqqAssoc["C610s"] * FT1["s","n"]),
		
					"d27p" -> - (CmueqqAssoc["C610u"] * FT2["u","p"] + CmueqqAssoc["C610d"] * FT2["d","p"] + CmueqqAssoc["C610s"] * FT2["s","p"]),
		
					"d27n" -> - (CmueqqAssoc["C610u"] * FT2["u","n"] + CmueqqAssoc["C610d"] * FT2["d","n"] + CmueqqAssoc["C610s"] * FT2["s","n"]),
		
							  (* T-odd contribution (FT3 set to zero by default) *)
					"d28p" -> - (CmueqqAssoc["C610u"] * FT3["u","p"] + CmueqqAssoc["C610d"] * FT3["d","p"] + CmueqqAssoc["C610s"] * FT3["s","p"]),
		
							  (* T-odd contribution (FT3 set to zero by default) *)
					"d28n" -> - (CmueqqAssoc["C610u"] * FT3["u","n"] + CmueqqAssoc["C610d"] * FT3["d","n"] + CmueqqAssoc["C610s"] * FT3["s","n"]),
		
					"d29p" -> - mL * (CmueqqAssoc["C713u"] * (FT0["u","p"] - (qSq / mp^2) * FT2["u","p"]) + CmueqqAssoc["C713d"] * (FT0["d","p"] - (qSq / mp^2) * FT2["d","p"]) + CmueqqAssoc["C713s"] * (FT0["s","p"] - (qSq / mp^2) * FT2["s","p"])),
		
					"d29n" -> - mL * (CmueqqAssoc["C713u"] * (FT0["u","n"] - (qSq / mn^2) * FT2["u","n"]) + CmueqqAssoc["C713d"] * (FT0["d","n"] - (qSq / mn^2) * FT2["d","n"]) + CmueqqAssoc["C713s"] * (FT0["s","n"] - (qSq / mn^2) * FT2["s","n"])),
		
					"d30p" -> - mL * (CmueqqAssoc["C714u"] * (FT0["u","p"] - (qSq / mp^2) * FT2["u","p"]) + CmueqqAssoc["C714d"] * (FT0["d","p"] - (qSq / mp^2) * FT2["d","p"]) + CmueqqAssoc["C714s"] * (FT0["s","p"] - (qSq / mp^2) * FT2["s","p"])),
		
					"d30n" -> - mL * (CmueqqAssoc["C714u"] * (FT0["u","n"] - (qSq / mn^2) * FT2["u","n"]) + CmueqqAssoc["C714d"] * (FT0["d","n"] - (qSq / mn^2) * FT2["d","n"]) + CmueqqAssoc["C714s"] * (FT0["s","n"] - (qSq / mn^2) * FT2["s","n"])),
		
					"d31p" -> + (1 / 4.) * mL * (CmueqqAssoc["C716u"] * FT0["u","p"] + CmueqqAssoc["C716d"] * FT0["d","p"] + CmueqqAssoc["C716s"] * FT0["s","p"]),
		
					"d31n" -> + (1 / 4.) * mL * (CmueqqAssoc["C716u"] * FT0["u","n"] + CmueqqAssoc["C716d"] * FT0["d","n"] + CmueqqAssoc["C716s"] * FT0["s","n"]),
		
					"d32p" -> + (1 / 4.) * mL * (CmueqqAssoc["C715u"] * FT0["u","p"] + CmueqqAssoc["C715d"] * FT0["d","p"] + CmueqqAssoc["C715s"] * FT0["s","p"]),
		
					"d32n" -> + (1 / 4.) * mL * (CmueqqAssoc["C715u"] * FT0["u","n"] + CmueqqAssoc["C715d"] * FT0["d","n"] + CmueqqAssoc["C715s"] * FT0["s","n"]) |>;

(* Convert to the isospin basis (written in two-component form {isoscalar, isovector}) *)
RETisospinAssoc = <| "d1"  -> {(RETnucleonAssoc["d1p"] + RETnucleonAssoc["d1n"])/2., (RETnucleonAssoc["d1p"] - RETnucleonAssoc["d1n"])/2.},
					 "d2"  -> {(RETnucleonAssoc["d2p"] + RETnucleonAssoc["d2n"])/2., (RETnucleonAssoc["d2p"] - RETnucleonAssoc["d2n"])/2.},
                  	 "d3"  -> {(RETnucleonAssoc["d3p"] + RETnucleonAssoc["d3n"])/2., (RETnucleonAssoc["d3p"] - RETnucleonAssoc["d3n"])/2.},
					 "d4"  -> {(RETnucleonAssoc["d4p"] + RETnucleonAssoc["d4n"])/2., (RETnucleonAssoc["d4p"] - RETnucleonAssoc["d4n"])/2.},
					 "d5"  -> {(RETnucleonAssoc["d5p"] + RETnucleonAssoc["d5n"])/2., (RETnucleonAssoc["d5p"] - RETnucleonAssoc["d5n"])/2.},
					 "d6"  -> {(RETnucleonAssoc["d6p"] + RETnucleonAssoc["d6n"])/2., (RETnucleonAssoc["d6p"] - RETnucleonAssoc["d6n"])/2.},
					 "d7"  -> {(RETnucleonAssoc["d7p"] + RETnucleonAssoc["d7n"])/2., (RETnucleonAssoc["d7p"] - RETnucleonAssoc["d7n"])/2.},
					 "d8"  -> {(RETnucleonAssoc["d8p"] + RETnucleonAssoc["d8n"])/2., (RETnucleonAssoc["d8p"] - RETnucleonAssoc["d8n"])/2.},
					 "d9"  -> {(RETnucleonAssoc["d9p"] + RETnucleonAssoc["d9n"])/2., (RETnucleonAssoc["d9p"] - RETnucleonAssoc["d9n"])/2.},
					 "d10" -> {(RETnucleonAssoc["d10p"] + RETnucleonAssoc["d10n"])/2., (RETnucleonAssoc["d10p"] - RETnucleonAssoc["d10n"])/2.},
					 "d11" -> {(RETnucleonAssoc["d11p"] + RETnucleonAssoc["d11n"])/2., (RETnucleonAssoc["d11p"] - RETnucleonAssoc["d11n"])/2.},
					 "d12" -> {(RETnucleonAssoc["d12p"] + RETnucleonAssoc["d12n"])/2., (RETnucleonAssoc["d12p"] - RETnucleonAssoc["d12n"])/2.},
					 "d13" -> {(RETnucleonAssoc["d13p"] + RETnucleonAssoc["d13n"])/2., (RETnucleonAssoc["d13p"] - RETnucleonAssoc["d13n"])/2.},
					 "d14" -> {(RETnucleonAssoc["d14p"] + RETnucleonAssoc["d14n"])/2., (RETnucleonAssoc["d14p"] - RETnucleonAssoc["d14n"])/2.},
					 "d15" -> {(RETnucleonAssoc["d15p"] + RETnucleonAssoc["d15n"])/2., (RETnucleonAssoc["d15p"] - RETnucleonAssoc["d15n"])/2.},
					 "d16" -> {(RETnucleonAssoc["d16p"] + RETnucleonAssoc["d16n"])/2., (RETnucleonAssoc["d16p"] - RETnucleonAssoc["d16n"])/2.},
					 "d17" -> {(RETnucleonAssoc["d17p"] + RETnucleonAssoc["d17n"])/2., (RETnucleonAssoc["d17p"] - RETnucleonAssoc["d17n"])/2.},
					 "d18" -> {(RETnucleonAssoc["d18p"] + RETnucleonAssoc["d18n"])/2., (RETnucleonAssoc["d18p"] - RETnucleonAssoc["d18n"])/2.},
					 "d19" -> {(RETnucleonAssoc["d19p"] + RETnucleonAssoc["d19n"])/2., (RETnucleonAssoc["d19p"] - RETnucleonAssoc["d19n"])/2.},
					 "d20" -> {(RETnucleonAssoc["d20p"] + RETnucleonAssoc["d20n"])/2., (RETnucleonAssoc["d20p"] - RETnucleonAssoc["d20n"])/2.},
					 "d21" -> {(RETnucleonAssoc["d21p"] + RETnucleonAssoc["d21n"])/2., (RETnucleonAssoc["d21p"] - RETnucleonAssoc["d21n"])/2.},
					 "d22" -> {(RETnucleonAssoc["d22p"] + RETnucleonAssoc["d22n"])/2., (RETnucleonAssoc["d22p"] - RETnucleonAssoc["d22n"])/2.},
					 "d23" -> {(RETnucleonAssoc["d23p"] + RETnucleonAssoc["d23n"])/2., (RETnucleonAssoc["d23p"] - RETnucleonAssoc["d23n"])/2.},
					 "d24" -> {(RETnucleonAssoc["d24p"] + RETnucleonAssoc["d24n"])/2., (RETnucleonAssoc["d24p"] - RETnucleonAssoc["d24n"])/2.},
					 "d25" -> {(RETnucleonAssoc["d25p"] + RETnucleonAssoc["d25n"])/2., (RETnucleonAssoc["d25p"] - RETnucleonAssoc["d25n"])/2.},
					 "d26" -> {(RETnucleonAssoc["d26p"] + RETnucleonAssoc["d26n"])/2., (RETnucleonAssoc["d26p"] - RETnucleonAssoc["d26n"])/2.},
					 "d27" -> {(RETnucleonAssoc["d27p"] + RETnucleonAssoc["d27n"])/2., (RETnucleonAssoc["d27p"] - RETnucleonAssoc["d27n"])/2.},
					 "d28" -> {(RETnucleonAssoc["d28p"] + RETnucleonAssoc["d28n"])/2., (RETnucleonAssoc["d28p"] - RETnucleonAssoc["d28n"])/2.},
					 "d29" -> {(RETnucleonAssoc["d29p"] + RETnucleonAssoc["d29n"])/2., (RETnucleonAssoc["d29p"] - RETnucleonAssoc["d29n"])/2.},
					 "d30" -> {(RETnucleonAssoc["d30p"] + RETnucleonAssoc["d30n"])/2., (RETnucleonAssoc["d30p"] - RETnucleonAssoc["d30n"])/2.},
					 "d31" -> {(RETnucleonAssoc["d31p"] + RETnucleonAssoc["d31n"])/2., (RETnucleonAssoc["d31p"] - RETnucleonAssoc["d31n"])/2.},
					 "d32" -> {(RETnucleonAssoc["d32p"] + RETnucleonAssoc["d32n"])/2., (RETnucleonAssoc["d32p"] - RETnucleonAssoc["d32n"])/2.} |>;

(* Return in list format compatible with Mu2e_NRET *)
RETisospinList = {{1,  RETisospinAssoc["d1"][[1]], RETisospinAssoc["d1"][[2]]},
				  {2,  RETisospinAssoc["d2"][[1]], RETisospinAssoc["d2"][[2]]},
				  {3,  RETisospinAssoc["d3"][[1]], RETisospinAssoc["d3"][[2]]},
				  {4,  RETisospinAssoc["d4"][[1]], RETisospinAssoc["d4"][[2]]},
				  {5,  RETisospinAssoc["d5"][[1]], RETisospinAssoc["d5"][[2]]},
				  {6,  RETisospinAssoc["d6"][[1]], RETisospinAssoc["d6"][[2]]},
				  {7,  RETisospinAssoc["d7"][[1]], RETisospinAssoc["d7"][[2]]},
				  {8,  RETisospinAssoc["d8"][[1]], RETisospinAssoc["d8"][[2]]},
				  {9,  RETisospinAssoc["d9"][[1]], RETisospinAssoc["d9"][[2]]},
				  {10, RETisospinAssoc["d10"][[1]], RETisospinAssoc["d10"][[2]]},
				  {11, RETisospinAssoc["d11"][[1]], RETisospinAssoc["d11"][[2]]},
				  {12, RETisospinAssoc["d12"][[1]], RETisospinAssoc["d12"][[2]]},
				  {13, RETisospinAssoc["d13"][[1]], RETisospinAssoc["d13"][[2]]},
				  {14, RETisospinAssoc["d14"][[1]], RETisospinAssoc["d14"][[2]]},
				  {15, RETisospinAssoc["d15"][[1]], RETisospinAssoc["d15"][[2]]},
				  {16, RETisospinAssoc["d16"][[1]], RETisospinAssoc["d16"][[2]]},
				  {17, RETisospinAssoc["d17"][[1]], RETisospinAssoc["d17"][[2]]},
				  {18, RETisospinAssoc["d18"][[1]], RETisospinAssoc["d18"][[2]]},
				  {19, RETisospinAssoc["d19"][[1]], RETisospinAssoc["d19"][[2]]},
				  {20, RETisospinAssoc["d20"][[1]], RETisospinAssoc["d20"][[2]]},
				  {21, RETisospinAssoc["d21"][[1]], RETisospinAssoc["d21"][[2]]},
				  {22, RETisospinAssoc["d22"][[1]], RETisospinAssoc["d22"][[2]]},
				  {23, RETisospinAssoc["d23"][[1]], RETisospinAssoc["d23"][[2]]},
				  {24, RETisospinAssoc["d24"][[1]], RETisospinAssoc["d24"][[2]]},
				  {25, RETisospinAssoc["d25"][[1]], RETisospinAssoc["d25"][[2]]},
				  {26, RETisospinAssoc["d26"][[1]], RETisospinAssoc["d26"][[2]]},
				  {27, RETisospinAssoc["d27"][[1]], RETisospinAssoc["d27"][[2]]},
				  {28, RETisospinAssoc["d28"][[1]], RETisospinAssoc["d28"][[2]]},
				  {29, RETisospinAssoc["d29"][[1]], RETisospinAssoc["d29"][[2]]},
				  {30, RETisospinAssoc["d30"][[1]], RETisospinAssoc["d30"][[2]]},
				  {31, RETisospinAssoc["d31"][[1]], RETisospinAssoc["d31"][[2]]},
				  {32, RETisospinAssoc["d32"][[1]], RETisospinAssoc["d32"][[2]]}};

RETisospinList]
End[]
EndPackage[]
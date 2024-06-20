(* 
  ::  formfactors.wl is a part of the MuonConverter package.                      
  ::  Copyright (C) 2024 MuonConverter authors (see AUTHORS for details).
  ::  MuonConverter is licenced under the BSD 3-Clause, see LICENSE for details.
*)

BeginPackage["MuonConverter`"]

(* HadronizeWET *)
PrintFormFactors::usage = "Check the numerical values of the form-factors at q^2";

Begin["`Private`"]

(* -------------------------------------------------------------------------- * 
 *  Form factors
 * -------------------------------------------------------------------------- *)

(* Define functions to initialize full q-dependent form factor *)

FormFactor[q_, op_, parton_, nucleon_] := Module[{qSq},
qSq = q[[1]] * q[[1]] - q[[2]] * q[[2]] - q[[3]] * q[[3]] - q[[4]] * q[[4]];
FF0[op, parton, nucleon] + qSq * FFPrime[op, parton, nucleon]
]

FormFactorNucleon[q_, op_, nucleon_] := Module[{qSq},
qSq = q[[1]] * q[[1]] - q[[2]] * q[[2]] - q[[3]] * q[[3]] - q[[4]] * q[[4]];
FF0[op, nucleon] + qSq * FFPrime[op, nucleon]
]

FormFactorPprimed[q_, parton_, nucleon_] := Module[{qSq, PP, EP, b},
qSq = q[[1]] * q[[1]] - q[[2]] * q[[2]] - q[[3]] * q[[3]] - q[[4]] * q[[4]];
PP = (mN^2 / (m\[Pi]^2 - qSq)) * FFPion["Pp", parton, nucleon];
EP = (mN^2 / (m\[Eta]^2 - qSq)) * FFEta["Pp", parton, nucleon];
b = FF0["Pp", parton, nucleon] - mN^2 * ((FFPion["Pp", parton, nucleon] / m\[Pi]^2) + (FFEta["Pp", parton, nucleon] / m\[Eta]^2));
PP + EP + b
]

FormFactorP[q_, parton_, nucleon_] := Module[{qSq, mtilde, Deltau, Deltad, PP, EP, bP},
qSq = q[[1]] * q[[1]] - q[[2]] * q[[2]] - q[[3]] * q[[3]] - q[[4]] * q[[4]];
mtilde = ((1 / mu) + (1 / md) + (1 / ms))^-1;
Deltau = (1 / 2) * (gA + \[CapitalDelta]\[CapitalSigma]ud);
Deltad = (1 / 2) * (- gA + \[CapitalDelta]\[CapitalSigma]ud);
PP = (mN^2 / (m\[Pi]^2 - qSq)) * FFPion["P", parton, nucleon];
EP = (mN^2 / (m\[Eta]^2 - qSq)) * FFEta["P", parton, nucleon];
bP["u"] = bP["d"] = (mN / 3) * (\[CapitalDelta]\[CapitalSigma]ud - 2 * \[CapitalDelta]s) * ((6 * rs - 3 * \[CapitalDelta]GT8) / (4 * rs + 2.)) + mN * (-(1 / 3) * gA + \[CapitalDelta]s) - mN * mtilde * ((Deltau / mu) + (Deltad / md) + \[CapitalDelta]s / ms);
bP["s"] = bP["u"] + (1 / 2) * mN * \[CapitalDelta]GT8 * (Deltau + Deltad - 2 *  \[CapitalDelta]s); 
PP + EP + bP[parton]
]

(* Define separate function for form factors with non-zero pion and eta poles *)

(* ======================= Vector ======================= *)

(* F1 *)

(* Set the value of the F1 form factor at zero momentum transfer *)
If[F1up0 == None, FF0["1","u","p"] = 2, FF0["1","u","p"] = F1up0];
If[F1dp0 == None, FF0["1","d","p"] = 1, FF0["1","d","p"] = F1dp0];
If[F1sp0 == None, FF0["1","s","p"] = 0, FF0["1","s","p"] = F1sp0];

If[F1un0 == None, FF0["1","u","n"] = 1, FF0["1","u","n"] = F1un0];
If[F1dn0 == None, FF0["1","d","n"] = 2, FF0["1","d","n"] = F1dn0];
If[F1sn0 == None, FF0["1","s","n"] = 0, FF0["1","s","n"] = F1sn0];

(* Set first derivative of F1 at zero momentum transfer *)
If[F1upPrime == None, FFPrime["1","u","p"] = (1. / 6.) * (2. * rEp2 + rEn2 + rEs2) - (1. / (4. * mN^2)) * (2. * \[Mu]p + \[Mu]n + \[Mu]s - 2.), FFPrime["1","u","p"] = F1upPrime];
If[F1dpPrime == None, FFPrime["1","d","p"] = (1. / 6.) * (rEp2 + 2. * rEn2 + rEs2) - (1. / (4. * mN^2)) * (\[Mu]p + 2 * \[Mu]n + \[Mu]s - 1.), FFPrime["1","d","p"] = F1dpPrime];
If[F1spPrime == None, FFPrime["1","s","p"] = (1. / 6.) * rEs2 - (\[Mu]s / (4 * mN^2)), FFPrime["1","s","p"] = F1spPrime];

If[F1unPrime == None, FFPrime["1","u","n"] = (1. / 6.) * (rEp2 + 2. * rEn2 + rEs2) - (1. / (4. * mN^2)) * (\[Mu]p + 2 * \[Mu]n + \[Mu]s - 1.), FFPrime["1","u","n"] = F1unPrime];
If[F1dnPrime == None, FFPrime["1","d","n"] = (1. / 6.) * (2. * rEp2 + rEn2 + rEs2) - (1. / (4. * mN^2)) * (2. * \[Mu]p + \[Mu]n + \[Mu]s - 2.), FFPrime["1","d","n"] = F1dnPrime];
If[F1snPrime == None, FFPrime["1","s","n"] = (1. / 6.) * rEs2 - (\[Mu]s / (4 * mN^2)), FFPrime["1","s","n"] = F1snPrime];

(* F2 *)

(* Set the value of the F2 form factor at zero momentum transfer *)
If[F2up0 == None, FF0["2","u","p"] = 2. * (\[Mu]p - 1.) + \[Mu]n + \[Mu]s, FF0["2","u","p"] = F2up0];
If[F2dp0 == None, FF0["2","d","p"] = 2. * \[Mu]n  + (\[Mu]p - 1.) + \[Mu]s, FF0["2","d","p"] = F2dp0];
If[F2sp0 == None, FF0["2","s","p"] = \[Mu]s, FF0["2","s","p"] = F2sp0];

If[F2un0 == None, FF0["2","u","n"] = 2. * \[Mu]n  + (\[Mu]p - 1.) + \[Mu]s, FF0["2","u","n"] = F2un0];
If[F2dn0 == None, FF0["2","d","n"] = 2. * (\[Mu]p - 1.) + \[Mu]n + \[Mu]s, FF0["2","d","n"] = F2dn0];
If[F2sn0 == None, FF0["2","s","n"] = \[Mu]s, FF0["2","s","n"] = F2sn0];

(* Set first derivative of F2 at zero momentum transfer *)
If[F2upPrime == None, FFPrime["2","u","p"] = (1. / 6.) * (2 * (rMp2 - rEp2) + rMn2 - rEn2 + rMs2 - rEs2) + (1 / (4 * mN^2)) * (2 * \[Mu]p + \[Mu]n + \[Mu]s - 2.), FFPrime["2","u","p"] = F2upPrime];
If[F2dpPrime == None, FFPrime["2","d","p"] = (1. / 6.) * (rMp2 - rEp2 + 2 * (rMn2 - rEn2) + rMs2 - rEs2) + (1 / (4 * mN^2)) * (\[Mu]p + 2 * \[Mu]n + \[Mu]s - 1.), FFPrime["2","d","p"] = F2dpPrime];
If[F2spPrime == None, FFPrime["2","s","p"] = (1. / 6.) * (rMs2 - rEs2) + (\[Mu]s / (4 * mN^2)), FFPrime["2","s","p"] = F2spPrime];

If[F2unPrime == None, FFPrime["2","u","n"] = (1. / 6.) * (rMp2 - rEp2 + 2 * (rMn2 - rEn2) + rMs2 - rEs2) + (1 / (4 * mN^2)) * (\[Mu]p + 2 * \[Mu]n + \[Mu]s - 1.), FFPrime["2","u","n"] = F2unPrime];
If[F2dnPrime == None, FFPrime["2","d","n"] = (1. / 6.) * (2 * (rMp2 - rEp2) + rMn2 - rEn2 + rMs2 - rEs2) + (1 / (4 * mN^2)) * (2 * \[Mu]p + \[Mu]n + \[Mu]s - 2.), FFPrime["2","d","n"] = F2dnPrime];
If[F2snPrime == None, FFPrime["2","s","n"] = (1. / 6.) * (rMs2 - rEs2) + (\[Mu]s / (4 * mN^2)), FFPrime["2","s","n"] = F2snPrime];

(* F3 (set to zero by default) *)

(* Set the value of the F3 form factor at zero momentum transfer *)
If[F3up0 == None, FF0["3","u","p"] = 0, FF0["3","u","p"] = F3up0];
If[F3dp0 == None, FF0["3","d","p"] = 0, FF0["3","d","p"] = F3dp0];
If[F3sp0 == None, FF0["3","s","p"] = 0, FF0["3","s","p"] = F3sp0];

If[F3un0 == None, FF0["3","u","n"] = 0, FF0["3","u","n"] = F3un0];
If[F3dn0 == None, FF0["3","d","n"] = 0, FF0["3","d","n"] = F3dn0];
If[F3sn0 == None, FF0["3","s","n"] = 0, FF0["3","s","n"] = F3sn0];

(* Set first derivative of F3 at zero momentum transfer *)
If[F3upPrime == None, FFPrime["3","u","p"] = 0, FFPrime["3","u","p"] = F3upPrime];
If[F3dpPrime == None, FFPrime["3","d","p"] = 0, FFPrime["3","d","p"] = F3dpPrime];
If[F3spPrime == None, FFPrime["3","s","p"] = 0, FFPrime["3","s","p"] = F3spPrime];

If[F3unPrime == None, FFPrime["3","u","n"] = 0, FFPrime["3","u","n"] = F3unPrime];
If[F3dnPrime == None, FFPrime["3","d","n"] = 0, FFPrime["3","d","n"] = F3dnPrime];
If[F3snPrime == None, FFPrime["3","s","n"] = 0, FFPrime["3","s","n"] = F3snPrime];


(* ======================= Axial-vector ======================= *)

(* FA *)

(* Set the value of the F3 form factor at zero momentum transfer *)
If[FAup0 == None, FF0["A","u","p"] = (1 / 2) * (gA + \[CapitalDelta]\[CapitalSigma]ud), FF0["A","u","p"] = FAup0];
If[FAdp0 == None, FF0["A","d","p"] = (1 / 2) * (- gA + \[CapitalDelta]\[CapitalSigma]ud), FF0["A","d","p"] = FAdp0];
If[FAsp0 == None, FF0["A","s","p"] = \[CapitalDelta]s, FF0["A","s","p"] = FAsp0];

If[FAun0 == None, FF0["A","u","n"] = (1 / 2) * (- gA + \[CapitalDelta]\[CapitalSigma]ud), FF0["A","u","n"] = FAun0];
If[FAdn0 == None, FF0["A","d","n"] = (1 / 2) * (gA + \[CapitalDelta]\[CapitalSigma]ud), FF0["A","d","n"] = FAdn0];
If[FAsn0 == None, FF0["A","s","n"] = \[CapitalDelta]s, FF0["A","s","n"] = FAsn0];

(* Set first derivative of F3 at zero momentum transfer *)
If[FAupPrime == None, FFPrime["A","u","p"] = (1 / 12) * (gA * rA2uMinusd + \[CapitalDelta]\[CapitalSigma]ud * rA2uPlusd), FFPrime["A","u","p"] = FAupPrime];
If[FAdpPrime == None, FFPrime["A","d","p"] = (1 / 12) * (- gA * rA2uMinusd + \[CapitalDelta]\[CapitalSigma]ud * rA2uPlusd), FFPrime["A","d","p"] = FAdpPrime];
If[FAspPrime == None, FFPrime["A","s","p"] = (1 / 6) * \[CapitalDelta]s * rA2s, FFPrime["A","s","p"] = FAspPrime];

If[FAunPrime == None, FFPrime["A","u","n"] = (1 / 12) * (- gA * rA2uMinusd + \[CapitalDelta]\[CapitalSigma]ud * rA2uPlusd), FFPrime["A","u","n"] = FAunPrime];
If[FAdnPrime == None, FFPrime["A","d","n"] = (1 / 12) * (gA * rA2uMinusd + \[CapitalDelta]\[CapitalSigma]ud * rA2uPlusd), FFPrime["A","d","n"] = FAdnPrime];
If[FAsnPrime == None, FFPrime["A","s","n"] = (1 / 6) * \[CapitalDelta]s * rA2s, FFPrime["A","s","n"] = FAsnPrime];

(* FPPrimed *)

(* Set the value of the F3 form factor at zero momentum transfer *)
FF0["Pp","u","p"] = FPprimedup0;
FF0["Pp","d","p"] = FPprimeddp0;
FF0["Pp","s","p"] = FPprimedsp0;

FF0["Pp","u","n"] = FPprimedun0;
FF0["Pp","d","n"] = FPprimeddn0;
FF0["Pp","s","n"] = FPprimedsn0;

(* Set the value of the form factor at the pion-pole *)
aPprimepiu = 2 * gA;
aPprimepid = -2 * gA;

If[FPprimedupPion == None, FFPion["Pp", "u", "p"] = aPprimepiu, FFPion["Pp", "u", "p"] = FPprimedupPion];
If[FPprimeddpPion == None, FFPion["Pp", "d", "p"] = aPprimepid, FFPion["Pp", "d", "p"] = FPprimeddpPion];
If[FPprimedspPion == None, FFPion["Pp", "s", "p"] = 0,          FFPion["Pp", "s", "p"] = FPprimedspPion];

If[FPprimedunPion == None, FFPion["Pp", "u", "n"] = aPprimepid, FFPion["Pp", "u", "n"] = FPprimedunPion];
If[FPprimeddnPion == None, FFPion["Pp", "d", "n"] = aPprimepiu, FFPion["Pp", "d", "n"] = FPprimeddnPion];
If[FPprimedsnPion == None, FFPion["Pp", "s", "n"] = 0,          FFPion["Pp", "s", "n"] = FPprimedsnPion];

(* Set the value of the form factor at the eta-pole *)
aPprimeetau = (2. / 3.) * (\[CapitalDelta]\[CapitalSigma]ud + 2. * \[CapitalDelta]s * (1. + \[CapitalDelta]GT8));
aPprimeetad = (2. / 3.) * (\[CapitalDelta]\[CapitalSigma]ud + 2. * \[CapitalDelta]s * (1. + \[CapitalDelta]GT8));
aPprimeetas = - (1. / 3.) * (\[CapitalDelta]\[CapitalSigma]ud + 2. * \[CapitalDelta]s * (1. + \[CapitalDelta]GT8));

If[FPprimedupEta == None, FFEta["Pp", "u", "p"] = aPprimeetau, FFEta["Pp", "u", "p"] = FPprimedupEta];
If[FPprimeddpEta == None, FFEta["Pp", "d", "p"] = aPprimeetad, FFEta["Pp", "d", "p"] = FPprimeddpEta];
If[FPprimedspEta == None, FFEta["Pp", "s", "p"] = aPprimeetas, FFEta["Pp", "s", "p"] = FPprimedspEta];

If[FPprimedunEta == None, FFEta["Pp", "u", "n"] = aPprimeetad, FFEta["Pp", "u", "n"] = FPprimedunEta];
If[FPprimeddnEta == None, FFEta["Pp", "d", "n"] = aPprimeetau, FFEta["Pp", "d", "n"] = FPprimeddnEta];
If[FPprimedsnEta == None, FFEta["Pp", "s", "n"] = aPprimeetas, FFEta["Pp", "s", "n"] = FPprimedsnEta];

(* FA3 (set to zero by default) *)

If[FA3up0 == None, FF0["A3","u","p"] = 0, FF0["A3","u","p"] = FA3up0];
If[FA3dp0 == None, FF0["A3","d","p"] = 0, FF0["A3","d","p"] = FA3dp0];
If[FA3sp0 == None, FF0["A3","s","p"] = 0, FF0["A3","s","p"] = FA3sp0];

If[FA3un0 == None, FF0["A3","u","n"] = 0, FF0["A3","u","n"] = FA3un0];
If[FA3dn0 == None, FF0["A3","d","n"] = 0, FF0["A3","d","n"] = FA3dn0];
If[FA3sn0 == None, FF0["A3","s","n"] = 0, FF0["A3","s","n"] = FA3sn0];

(* Set first derivative of F3 at zero momentum transfer *)
If[FA3upPrime == None, FFPrime["A3","u","p"] = 0, FFPrime["A3","u","p"] = FA3upPrime];
If[FA3dpPrime == None, FFPrime["A3","d","p"] = 0, FFPrime["A3","d","p"] = FA3dpPrime];
If[FA3spPrime == None, FFPrime["A3","s","p"] = 0, FFPrime["A3","s","p"] = FA3spPrime];

If[FA3unPrime == None, FFPrime["A3","u","n"] = 0, FFPrime["A3","u","n"] = FA3unPrime];
If[FA3dnPrime == None, FFPrime["A3","d","n"] = 0, FFPrime["A3","d","n"] = FA3dnPrime];
If[FA3snPrime == None, FFPrime["A3","s","n"] = 0, FFPrime["A3","s","n"] = FA3snPrime];

(* ======================= Scalar ======================= *)

(* FS *)
xi = (1 - rud) / (1 + rud);
sigmaup = (1 / 2) * \[Sigma]\[Pi]Ntilde * (1 - xi) + c5hat * (1 - (1 / xi));
sigmadp = (1 / 2) * \[Sigma]\[Pi]Ntilde * (1 + xi) + c5hat * (1 + (1 / xi));
sigmaun = (1 / 2) * \[Sigma]\[Pi]Ntilde * (1 - xi) - c5hat * (1 - (1 / xi));
sigmadn = (1 / 2) * \[Sigma]\[Pi]Ntilde * (1 + xi) - c5hat * (1 + (1 / xi));

(* Set the value of the F3 form factor at zero momentum transfer *)
If[FSup0 == None, FF0["S","u","p"] = sigmaup, FF0["S","u","p"] = FSup0];
If[FSdp0 == None, FF0["S","d","p"] = sigmadp, FF0["S","d","p"] = FSdp0];
If[FSsp0 == None, FF0["S","s","p"] = \[Sigma]s, FF0["S","s","p"] = FSsp0];

If[FSun0 == None, FF0["S","u","n"] = sigmaun, FF0["S","u","n"] = FSun0];
If[FSdn0 == None, FF0["S","d","n"] = sigmadn, FF0["S","d","n"] = FSdn0];
If[FSsn0 == None, FF0["S","s","n"] = \[Sigma]s, FF0["S","s","n"] = FSsn0];

(* Set first derivative of F3 at zero momentum transfer *)
FFPrime["S","u","p"] = FSupPrime;
FFPrime["S","d","p"] = FSdpPrime;
FFPrime["S","s","p"] = FSspPrime;

FFPrime["S","u","n"] = FSunPrime;
FFPrime["S","d","n"] = FSdnPrime;
FFPrime["S","s","n"] = FSsnPrime;

(* ======================= Psuedo-scalar ======================= *)

(* FP *)

(* Set the value of the form factor at the pion-pole *)
aPpiu = (m\[Pi]^2 / mN) * (1 / (1 + (1 / rud))) * gA;
aPpid = (m\[Pi]^2 / mN) * (1 / (1 + rud)) * gA;

If[FPupPion == None, FFPion["P","u","p"] = aPpiu, FFPion["P","u","p"] = FPupPion];
If[FPdpPion == None, FFPion["P","d","p"] = aPpid, FFPion["P","d","p"] = FPdpPion];
If[FPspPion == None, FFPion["P","s","p"] = 0.,    FFPion["P","s","p"] = FPspPion];

If[FPunPion == None, FFPion["P","u","n"] = aPpid, FFPion["P","u","n"] = FPunPion];
If[FPdnPion == None, FFPion["P","d","n"] = aPpiu, FFPion["P","d","n"] = FPdnPion];
If[FPsnPion == None, FFPion["P","s","n"] = 0.,    FFPion["P","s","n"] = FPsnPion];

(* Set the value of the form factor at the eta-pole *)
aPetau = (m\[Eta]^2 / mN) * (1 / (1 + (1 / rud))) * (1 / (1 + 2 * rs)) * (\[CapitalDelta]\[CapitalSigma]ud - 2 * \[CapitalDelta]s) * (1 + \[CapitalDelta]GT8);
aPetad = (m\[Eta]^2 / mN) * (1 / (1 + rud)) * (1 / (1 + 2 * rs)) * (\[CapitalDelta]\[CapitalSigma]ud - 2 * \[CapitalDelta]s) * (1 + \[CapitalDelta]GT8);
aPetas = - (m\[Eta]^2 / mN) * (1 / (2 + (1 / rs))) * (\[CapitalDelta]\[CapitalSigma]ud - 2 * \[CapitalDelta]s) * (1 + \[CapitalDelta]GT8);

If[FPupEta == None, FFEta["P","u","p"] = aPetau, FFEta["P","u","p"] = FPupEta];
If[FPdpEta == None, FFEta["P","d","p"] = aPetad, FFEta["P","d","p"] = FPdpEta];
If[FPspEta == None, FFEta["P","s","p"] = aPetas, FFEta["P","s","p"] = FPspEta];

If[FPunEta == None, FFEta["P","u","n"] = aPetad, FFEta["P","u","n"] = FPunEta];
If[FPdnEta == None, FFEta["P","d","n"] = aPetau, FFEta["P","d","n"] = FPdnEta];
If[FPsnEta == None, FFEta["P","s","n"] = aPetas, FFEta["P","s","n"] = FPsnEta];

(* ======================= Tensor ======================= *)

(* FT0 *)

(* Set the value of the T0 form factor at zero momentum transfer *)
FF0["T0","u","p"] = FT0up0;
FF0["T0","d","p"] = FT0dp0;
FF0["T0","s","p"] = FT0sp0;

FF0["T0","u","n"] = FT0un0;
FF0["T0","d","n"] = FT0dn0;
FF0["T0","s","n"] = FT0sn0;

(* Set first derivative of T0 at zero momentum transfer *)
FFPrime["T0","u","p"] = FT0upPrime;
FFPrime["T0","d","p"] = FT0dpPrime;
FFPrime["T0","s","p"] = FT0spPrime;

FFPrime["T0","u","n"] = FT0unPrime;
FFPrime["T0","d","n"] = FT0dnPrime;
FFPrime["T0","s","n"] = FT0snPrime;

(* FT1 *)

(* Set the value of the T0 form factor at zero momentum transfer *)
FF0["T1","u","p"] = FT1up0;
FF0["T1","d","p"] = FT1dp0;
FF0["T1","s","p"] = FT1sp0;

FF0["T1","u","n"] = FT1un0;
FF0["T1","d","n"] = FT1dn0;
FF0["T1","s","n"] = FT1sn0;

(* Set first derivative of T0 at zero momentum transfer *)
FFPrime["T1","u","p"] = FT1upPrime;
FFPrime["T1","d","p"] = FT1dpPrime;
FFPrime["T1","s","p"] = FT1spPrime;

FFPrime["T1","u","n"] = FT1unPrime;
FFPrime["T1","d","n"] = FT1dnPrime;
FFPrime["T1","s","n"] = FT1snPrime;

(* FT2 *)

(* Set the value of the T0 form factor at zero momentum transfer *)
FF0["T2","u","p"] = FT2up0;
FF0["T2","d","p"] = FT2dp0;
FF0["T2","s","p"] = FT2sp0;

FF0["T2","u","n"] = FT2un0;
FF0["T2","d","n"] = FT2dn0;
FF0["T2","s","n"] = FT2sn0;

(* Set first derivative of T0 at zero momentum transfer *)
FFPrime["T2","u","p"] = FT2upPrime;
FFPrime["T2","d","p"] = FT2dpPrime;
FFPrime["T2","s","p"] = FT2spPrime;

FFPrime["T2","u","n"] = FT2unPrime;
FFPrime["T2","d","n"] = FT2dnPrime;
FFPrime["T2","s","n"] = FT2snPrime;

(* FT3 *)

(* Set the value of the T0 form factor at zero momentum transfer *)
FF0["T3","u","p"] = FT3up0;
FF0["T3","d","p"] = FT3dp0;
FF0["T3","s","p"] = FT3sp0;

FF0["T3","u","n"] = FT3un0;
FF0["T3","d","n"] = FT3dn0;
FF0["T3","s","n"] = FT3sn0;

(* Set first derivative of T0 at zero momentum transfer *)
FFPrime["T3","u","p"] = FT3upPrime;
FFPrime["T3","d","p"] = FT3dpPrime;
FFPrime["T3","s","p"] = FT3spPrime;

FFPrime["T3","u","n"] = FT3unPrime;
FFPrime["T3","d","n"] = FT3dnPrime;
FFPrime["T3","s","n"] = FT3snPrime;

(* ======================= Gluonic ======================= *)

(* FG *)

(* Set the value of the FG form factor at zero momentum transfer *)
FF0["FG","p"] = FGp0;
FF0["FG","n"] = FGn0;

(* Set first derivative of FG at zero momentum transfer *)
FFPrime["FG","p"] = FGpPrime;
FFPrime["FG","n"] = FGnPrime;

(* FGtilde *)

(* Set the value of the FG form factor at zero momentum transfer *)
FF0["FGt","p"] = FGtildep0; (* Need to implement expression *)
FF0["FGt","n"] = FGtilden0;

(* Set first derivative of FG at zero momentum transfer *)
FFPrime["FGt","p"] = FGtildepPrime;
FFPrime["FGt","n"] = FGtildenPrime;

(* ======================= Rayleigh ======================= *)

(* FF *)

(* Set the value of the FG form factor at zero momentum transfer *)
FF0["FF","p"] = FFp0;
FF0["FF","n"] = FFn0;

(* Set first derivative of FG at zero momentum transfer *)
FFPrime["FF","p"] = FFpPrime;
FFPrime["FF","n"] = FFnPrime;

(* FFtilde *)

(* Set the value of the FG form factor at zero momentum transfer *)
FF0["FFt","p"] = FFtildep0;
FF0["FFt","n"] = FFtilden0;

(* Set first derivative of FG at zero momentum transfer *)
FFPrime["FFt","p"] = FFtildepPrime;
FFPrime["FFt","n"] = FFtildenPrime;

(* Check the form factors *)
PrintFormFactors[q_] := Module[{},
Print["q = ", q];
Do[F1[qq,nn] = FormFactor[q, "1", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[F2[qq,nn] = FormFactor[q, "2", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[F3[qq,nn] = FormFactor[q, "3", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FA[qq,nn] = FormFactor[q, "A", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FPp[qq,nn] = FormFactorPprimed[q, qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FS[qq,nn] = FormFactor[q, "S", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FP[qq,nn] = FormFactorP[q, qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT0[qq,nn] = FormFactor[q, "T0", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT1[qq,nn] = FormFactor[q, "T1", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT2[qq,nn] = FormFactor[q, "T2", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FT3[qq,nn] = FormFactor[q, "T3", qq, nn], {qq,{"u","d","s"}},{nn,{"p","n"}}];
Do[FG[nn] = FormFactorNucleon[q,"FG",nn], {nn,{"p","n"}}];
Do[FGt[nn] = FormFactorNucleon[q,"FGt",nn], {nn,{"p","n"}}];
Do[FF[nn] = FormFactorNucleon[q,"FF",nn], {nn,{"p","n"}}];
Do[FFt[nn] = FormFactorNucleon[q,"FFt",nn], {nn,{"p","n"}}];

(* Vector F1 (proton) *)
Print["F1up(qSq) = ", F1["u", "p"]];
Print["F1dp(qSq) = ", F1["d", "p"]];
Print["F1sp(qSq) = ", F1["s", "p"]];
Print[];
(* Vector F1 (neutron) *)
Print["F1un(qSq) = ", F1["u", "n"]];
Print["F1dn(qSq) = ", F1["d", "n"]];
Print["F1sn(qSq) = ", F1["s", "n"]];
Print[];
(* Vector F2 (proton) *)
Print["F2up(qSq) = ", F2["u", "p"]];
Print["F2dp(qSq) = ", F2["d", "p"]];
Print["F2sp(qSq) = ", F2["s", "p"]];
Print[];
(* Vector F2 (neutron) *)
Print["F2un(qSq) = ", F2["u", "n"]];
Print["F2dn(qSq) = ", F2["d", "n"]];
Print["F2sn(qSq) = ", F2["s", "n"]];
Print[];
(* Vector F3 (proton) *)
Print["F3up(qSq) = ", F3["u", "p"]];
Print["F3dp(qSq) = ", F3["d", "p"]];
Print["F3sp(qSq) = ", F3["s", "p"]];
Print[];
(* Vector F3 (neutron) *)
Print["F3un(qSq) = ", F3["u", "n"]];
Print["F3dn(qSq) = ", F3["d", "n"]];
Print["F3sn(qSq) = ", F3["s", "n"]];
Print[];
(* Axial-vector FA (proton) *)
Print["FAup(qSq) = ", FA["u", "p"]];
Print["FAdp(qSq) = ", FA["d", "p"]];
Print["FAsp(qSq) = ", FA["s", "p"]];
Print[];
(* Axial-vector FA (neutron) *)
Print["FAun(qSq) = ", FA["u", "n"]];
Print["FAdn(qSq) = ", FA["d", "n"]];
Print["FAsn(qSq) = ", FA["s", "n"]];
Print[];
(* Axial-vector FPprimed (proton) *)
Print["FPprimedup(qSq) = ", FPp["u", "p"]];
Print["FPprimeddp(qSq) = ", FPp["d", "p"]];
Print["FPprimedsp(qSq) = ", FPp["s", "p"]];
Print[];
(* Axial-vector FPprimed (neutron) *)
Print["FPprimedun(qSq) = ", FPp["u", "n"]];
Print["FPprimeddn(qSq) = ", FPp["d", "n"]];
Print["FPprimedsn(qSq) = ", FPp["s", "n"]];
Print[];
(* Scalar FS (proton) *)
Print["FSup(qSq) = ", FS["u", "p"]];
Print["FSdp(qSq) = ", FS["d", "p"]];
Print["FSsp(qSq) = ", FS["s", "p"]];
Print[];
(* Scalar FS (neutron) *)
Print["FSun(qSq) = ", FS["u", "n"]];
Print["FSdn(qSq) = ", FS["d", "n"]];
Print["FSsn(qSq) = ", FS["s", "n"]];
Print[];
(* Psuedo-scalar FP (proton) *)
Print["FPup(qSq) = ", FP["u", "p"]];
Print["FPdp(qSq) = ", FP["d", "p"]];
Print["FPsp(qSq) = ", FP["s", "p"]];
Print[];
(* Psuedo-scalar FP (neutron) *)
Print["FPun(qSq) = ", FP["u", "n"]];
Print["FPdn(qSq) = ", FP["d", "n"]];
Print["FPsn(qSq) = ", FP["s", "n"]];
Print[];
(* Tensor FT0 (proton) *)
Print["FT0up(qSq) = ", FT0["u", "p"]];
Print["FT0dp(qSq) = ", FT0["d", "p"]];
Print["FT0sp(qSq) = ", FT0["s", "p"]];
Print[];
(* Tensor FT0 (neutron) *)
Print["FT0un(qSq) = ", FT0["u", "n"]];
Print["FT0dn(qSq) = ", FT0["d", "n"]];
Print["FT0sn(qSq) = ", FT0["s", "n"]];
Print[];
(* Tensor FT1 (proton) *)
Print["FT1up(qSq) = ", FT1["u", "p"]];
Print["FT1dp(qSq) = ", FT1["d", "p"]];
Print["FT1sp(qSq) = ", FT1["s", "p"]];
Print[];
(* Tensor FT1 (neutron) *)
Print["FT1un(qSq) = ", FT1["u", "n"]];
Print["FT1dn(qSq) = ", FT1["d", "n"]];
Print["FT1sn(qSq) = ", FT1["s", "n"]];
Print[];
(* Tensor FT2 (proton) *)
Print["FT2up(qSq) = ", FT2["u", "p"]];
Print["FT2dp(qSq) = ", FT2["d", "p"]];
Print["FT2sp(qSq) = ", FT2["s", "p"]];
Print[];
(* Tensor FT2 (neutron) *)
Print["FT2un(qSq) = ", FT2["u", "n"]];
Print["FT2dn(qSq) = ", FT2["d", "n"]];
Print["FT2sn(qSq) = ", FT2["s", "n"]];
Print[];
(* Tensor FT3 (proton) *)
Print["FT3up(qSq) = ", FT3["u", "p"]];
Print["FT3dp(qSq) = ", FT3["d", "p"]];
Print["FT3sp(qSq) = ", FT3["s", "p"]];
Print[];
(* Tensor FT3 (neutron) *)
Print["FT3un(qSq) = ", FT3["u", "n"]];
Print["FT3dn(qSq) = ", FT3["d", "n"]];
Print["FT3sn(qSq) = ", FT3["s", "n"]];
Print[];
(* Gluonic-even (proton) *)
Print["FGup(qSq) = ", FG["p"]];
Print[];
(* Gluonic-even (neutron) *)
Print["FGun(qSq) = ", FG["n"]];
Print[];
(* Gluonic-odd (proton) *)
Print["FGtup(qSq) = ", FGt["p"]];
Print[];
(* Gluonic-odd (neutron) *)
Print["FGtun(qSq) = ", FGt["n"]];
Print[];
(* Rayleigh-even (proton) *)
Print["FFup(qSq) = ", FF["p"]];
Print[];
(* Rayleigh-even (neutron) *)
Print["FFun(qSq) = ", FF["n"]];
Print[];
(* Rayleigh-odd (proton) *)
Print["FFtup(qSq) = ", FFt["p"]];
Print[];
(* Rayleigh-odd (neutron) *)
Print["FFtun(qSq) = ", FFt["n"]];
Print[];

];

End[]
EndPackage[]
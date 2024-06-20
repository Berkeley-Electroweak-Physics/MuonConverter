BeginPackage["MuonConverter`"]
Begin["`Private`"]

(* -------------------------------------------------------------------------- * 
 *  The u, d, and s masses are MSbar at 2 GeV from 2016 PDG
 * -------------------------------------------------------------------------- *)
mu = 2.16*^-3 (* GeV *); qu = 2./3.;
md = 4.67*^-3 (* GeV *); qd = -1./3.;
ms = 93*^-3  (* GeV *); qs = -1./3.;
mtilde = 1/(1/mu + 1/md + 1/ms) (* Eq. B17 in [1611.00368] *);

(* -------------------------------------------------------------------------- * 
 *  Quark masses at matching scales at all other scales are obtained by       *
 *  one-loop running with the appropriate number of active flavors            *
 * -------------------------------------------------------------------------- *)
MQ["u","MZ"] = 0.0015; MQ["u","MB"] = 0.0020; MQ["u","2GeV"] = mu;
MQ["d","MZ"] = 0.0031; MQ["d","MB"] = 0.0042; MQ["d","2GeV"] = md; 
MQ["s","MZ"] = 0.063;  MQ["s","MB"] = 0.086;  MQ["s","2GeV"] = ms;
MQ["c","MZ"] = 0.78;   MQ["c","MB"] = 1.06;
MQ["b","MZ"] = 3.08;

(* -------------------------------------------------------------------------- * 
 *  The average nucleon mass and the meson masses from PDG 2016
 * -------------------------------------------------------------------------- *)
mp = 0.938272; (* GeV *)
mn = 0.939565; (* GeV *)
(* Average of m_p and m_n *)
mN = (0.938272 + 0.939565) / 2.; (* GeV *)
m\[Pi]  = 134.9766*^-3; (* GeV *)
m\[Eta] = 547.862*^-3; (* GeV *)

(* -------------------------------------------------------------------------- * 
 *  The tensor charges and form factor coefficients @ 2 GeV
 * -------------------------------------------------------------------------- *)
gT["u"] = + 0.794;
gT["d"] = - 0.204;
gT["s"] = + 3.2*^-4;

BT10["u","p"] = BT10["d","n"] = 3.0;
BT10["d","p"] = BT10["u","n"] = 0.24;
BT10["s","p"] = BT10["s","n"] = 0.0;

AT10t["u","p"] = AT10t["d","n"] = -0.50;
AT10t["d","p"] = AT10t["u","n"] = 0.46;
AT10t["s","p"] = AT10t["s","n"] = 0.0;

(* -------------------------------------------------------------------------- * 
 *  Electroweak inputs
 * -------------------------------------------------------------------------- *)
\[Alpha]     = 1./137.035999084;  (* PDG 2018 -- EW SM review *)
MH           = 125.10;      (* [GeV] *)
MZ           = 91.1876;     (* [GeV] *)
MW           = 80.379;      (* [GeV] *)
MT           = 173.;        (* [GeV] *)
GF           = 1.1663787*^-5;
vev          = 1/Sqrt[Sqrt[2.]*GF] (* [GeV] *);
SW2MZ        = 0.23122 (* From PDG 2018; MSbar at MZ *);
CW2MZ        = 1-SW2MZ;
SW           = Sqrt[SW2MZ];
sw           = SW;
\[Lambda]    = 2*MH^2/vev^2  (* Dimensionless Higgs trilinear *);

(* -------------------------------------------------------------------------- * 
 *  Lepton masses
 * -------------------------------------------------------------------------- *)
me = 0.000510998928 (* GeV *);
mmu = 0.1056583715 (* GeV *);
mPlus = mmu + me;
mMinus = mmu - me;

(* Arbitrary scale *)
mL = mN;

(* -------------------------------------------------------------------------- * 
 *  Form factor parameters
 * -------------------------------------------------------------------------- *)

(* ======================= Global ======================= *)

rud = 0.474; (* Ratio of up and down quark masses (rud = mu / md) *)
rs = 27.33; (* Ratio of strange quark mass with sum of up and down masses (rs = (2ms) / (mu + md)) *)

(* ======================= Vector ======================= *)

(* F1(q^2) *)
F1up = None;
F1dp = None;
F1sp = None;

F1un = None;
F1dn = None;
F1sn = None;

(* F1 evaluated at zero momentum transfer -- F1 @ q^2 = 0 *)
F1up0 = None;
F1dp0 = None;
F1sp0 = None;

F1un0 = None;
F1dn0 = None;
F1sn0 = None;

(* F1 first derivative evaluated at zero momentum transfer -- (d [F1] / d q^2) @ q^2 = 0 *)
F1upPrime = None;
F1dpPrime = None;
F1spPrime = None;

F1unPrime = None;
F1dnPrime = None;
F1snPrime = None;

(* F2(q^2) *)
F2up = None;
F2dp = None;
F2sp = None;

F2un = None;
F2dn = None;
F2sn = None;

(* F2 evaluated at zero momentum transfer -- F2 @ q^2 = 0 *)
F2up0 = None;
F2dp0 = None;
F2sp0 = None;

F2un0 = None;
F2dn0 = None;
F2sn0 = None;

(* F2 first derivative evaluated at zero momentum transfer -- (d [F2] / d q^2) @ q^2 = 0 *)
F2upPrime = None;
F2dpPrime = None;
F2spPrime = None;

F2unPrime = None;
F2dnPrime = None;
F2snPrime = None;

(* F3 evaluated at zero momentum transfer -- F3 @ q^2 = 0 *)
F3up0 = None;
F3dp0 = None;
F3sp0 = None;

F3un0 = None;
F3dn0 = None;
F3sn0 = None;

(* F3 first derivative evaluated at zero momentum transfer -- (d [F3] / d q^2) @ q^2 = 0 *)
F3upPrime = None;
F3dpPrime = None;
F3spPrime = None;

F3unPrime = None;
F3dnPrime = None;
F3snPrime = None;

(*Proton nuclear dipole moment*)
\[Mu]p = 2.792847;

(* Neutron nuclear dipole moment *)
\[Mu]n = -1.91304;

(* Strange nuclear dipole moment *) 
\[Mu]s = -0.036;

(* Proton electric charge radius squared [1/GeV^2] *) 
rEp2 = 18.16;

(* Neutron electric charge radius squared [1/GeV^2] *) 
rEn2 = -2.966;

(* Strange electric charge radius squared [1/GeV^2] *) 
rEs2 = -0.115;

(* Proton magnetic charge radius squared [1/GeV^2] *) 
rMp2 = 18.6;

(* Neutron magnetic charge radius squared [1/GeV^2] *) 
rMn2 = 19.1;

(* Strange magnetic charge radius squared [1/GeV^2] *) 
rMs2 = -0.26;

(* ======================= Axial vector ======================= *)

(* FA(q^2) *)
FAup = None;
FAdp = None;
FAsp = None;

FAun = None;
FAdn = None;
FAsn = None;

(* FA evaluated at zero momentum transfer -- FA @ q^2 = 0 *)
FAup0 = None;
FAdp0 = None;
FAsp0 = None;

FAun0 = None;
FAdn0 = None;
FAsn0 = None;

(* FA first derivative evaluated at zero momentum transfer -- (d [FA] / d q^2) @ q^2 = 0 *)
FAupPrime = None;
FAdpPrime = None;
FAspPrime = None;

FAunPrime = None;
FAdnPrime = None;
FAsnPrime = None;

(* FPprimed(q^2) *)
FPprimedup = None;
FPprimeddp = None;
FPprimedsp = None;

FPprimedun = None;
FPprimeddn = None;
FPprimedsn = None;

(* FPprimed residua of pion pole contribution *)
FPprimedupPion = None;
FPprimeddpPion = None;
FPprimedspPion = None;

FPprimedunPion = None;
FPprimeddnPion = None;
FPprimedsnPion = None;

(* FPprimed residua of eta pole contribution *)
FPprimedupEta = None;
FPprimeddpEta = None;
FPprimedspEta = None;

FPprimedunEta = None;
FPprimeddnEta = None;
FPprimedsnEta = None;

(* FPprimed evaluated at zero momentum transfer -- FPprimed @ q^2 = 0 *)
FPprimedup0 = 119;
FPprimeddp0 = -130;
FPprimedsp0 = -1.6;

FPprimedun0 = -130;
FPprimeddn0 = 119;
FPprimedsn0 = -1.6;

(* FA3 evaluated at zero momentum transfer -- FA3 @ q^2 = 0 *)
FA3up0 = None;
FA3dp0 = None;
FA3sp0 = None;

FA3un0 = None;
FA3dn0 = None;
FA3sn0 = None;

(* FA3 first derivative evaluated at zero momentum transfer -- (d [FA3] / d q^2) @ q^2 = 0 *)
FA3upPrime = None;
FA3dpPrime = None;
FA3spPrime = None;

FA3unPrime = None;
FA3dnPrime = None;
FA3snPrime = None;

(* gA *)
gA = 1.2754;

(* Delta u + Delta d *)
\[CapitalDelta]\[CapitalSigma]ud = 0.397;

(* Deltas *)
\[CapitalDelta]s = -0.045;

(* Axial charge radius (u - d) [1/GeV^2] *)
rA2uMinusd = 10.1;

(* Axial charge radius (u + d) [1/GeV^2] *)
rA2uPlusd = 12.58;

(* Axial charge radius (s) [1/GeV^2] *)
rA2s = 12.;

(* Axial FPprime correction *)
\[CapitalDelta]GT8 = 0.50;

(* ======================= Scalar ======================= *)

(* FS(q^2) *)
FSup = None;
FSdp = None;
FSsp = None;

FSun = None;
FSdn = None;
FSsn = None;

(* FS evaluated at zero momentum transfer -- FS @ q^2 = 0 *)
FSup0 = None;
FSdp0 = None;
FSsp0 = None;

FSun0 = None;
FSdn0 = None;
FSsn0 = None;

(* FS first derivative evaluated at zero momentum transfer -- (d [FS] / d q^2) @ q^2 = 0 *)
FSupPrime = 0.72;
FSdpPrime = 0.59;
FSspPrime = 0.17;

FSunPrime = 0.59;
FSdnPrime = 0.72;
FSsnPrime = 0.17;

(* FP(q^2) *)
FPup = None;
FPdp = None;
FPsp = None;

FPun = None;
FPdn = None;
FPsn = None;

(* FP residua of pion pole contribution *)
FPupPion = None;
FPdpPion = None;
FPspPion = None;

FPunPion = None;
FPdnPion = None;
FPsnPion = None;

(* FP residua of eta pole contribution *)
FPupEta = None;
FPdpEta = None;
FPspEta = None;

FPunEta = None;
FPdnEta = None;
FPsnEta = None;

\[Sigma]\[Pi]Ntilde = 48*^-3;

c5hat = -0.51*^-3;

\[Sigma]s = 43.3*^-3;

(* ======================= Gluonic ======================= *)

(* FG(q^2) *)
FGp = None;
FGn = None;

(* FG evaluated at zero momentum transfer -- FG @ q^2 = 0 *)
FGp0 = -50.4*^-3;
FGn0 = -50.4*^-3;

(* FG first derivative -- (d [FG] / d q^2) @ q^2 = 0 *)
FGpPrime = -0.14;
FGnPrime = -0.14;

(* FGtilde(q^2) *)
FGtildep = None;
FGtilden = None;

(* FGtilde evaluated at zero momentum transfer -- FGtilde @ q^2 = 0 *)
FGtildep0 = None; (* Need to implement expression *)
FGtilden0 = None;

FGtildepPrime = 0.;
FGtildenPrime = 0.;

(* FGtilde residua of pion pole contribution *)
FGtildepPion = None;
FGtildenPion = None;

(* FGtilde residua of eta pole contribution *)
FGtildepEta = None;
FGtildenEta = None;

(* Strong coupling at 2 GeV *)
alphaSAt2GeV = 0.297;

(* mG *)
mG = 0.823;

(* ======================= Tensor ======================= *)

(* FT0(q^2) *)
FT0up = None;
FT0dp = None;
FT0sp = None;

FT0un = None;
FT0dn = None;
FT0sn = None;

(* FT0 evaluated at zero momentum transfer -- FT0 @ q^2 = 0 *)
FT0up0 = 0.784; (* gTu *)
FT0dp0 = -0.204; (* gTd *)
FT0sp0 = -2.7*^-3; (* gTs *)

FT0un0 = -0.204; (* gTd *)
FT0dn0 = 0.784; (* gTu *)
FT0sn0 = -2.7*^-3; (* gTs *)

(* FT0 first derivative evaluated at zero momentum transfer -- (d [FT0] / d q^2) @ q^2 = 0 *)
FT0upPrime = 0.54;
FT0dpPrime = -0.11;
FT0spPrime = -0.0014;

FT0unPrime = -0.11;
FT0dnPrime = 0.54;
FT0snPrime = -0.0014;

(* FT1(q^2) *)
FT1up = None;
FT1dp = None;
FT1sp = None;

FT1un = None;
FT1dn = None;
FT1sn = None;

(* FT1 evaluated at zero momentum transfer -- FT1 @ q^2 = 0 *)
FT1up0 = -3.0;
FT1dp0 = 1.0;
FT1sp0 = 0.018;

FT1un0 = 1.0;
FT1dn0 = -3.0;
FT1sn0 = 0.018;

(* FT1 first derivative evaluated at zero momentum transfer -- (d [FT1] / d q^2) @ q^2 = 0 *)
FT1upPrime = -14.0;
FT1dpPrime = 5.0;
FT1spPrime = 0.082;

FT1unPrime = 5.0;
FT1dnPrime = -14.0;
FT1snPrime = 0.082;

(* FT2(q^2) *)
FT2up = None;
FT2dp = None;
FT2sp = None;

FT2un = None;
FT2dn = None;
FT2sn = None;

(* FT2 evaluated at zero momentum transfer -- FT2 @ q^2 = 0 *)
FT2up0 = -0.1;
FT2dp0 = 0.6;
FT2sp0 = 0.004;

FT2un0 = 0.6;
FT2dn0 = -0.1;
FT2sn0 = 0.004;

(* FT2 first derivative evaluated at zero momentum transfer -- (d [FT2] / d q^2) @ q^2 = 0 *)
FT2upPrime = -1.8;
FT2dpPrime = 2.1;
FT2spPrime = 0.015;

FT2unPrime = 2.1;
FT2dnPrime = -1.8;
FT2snPrime = 0.015;

(* FT3 evaluated at zero momentum transfer -- FT3 @ q^2 = 0 *)
FT3up0 = 0.;
FT3dp0 = 0.;
FT3sp0 = 0.;

FT3un0 = 0.;
FT3dn0 = 0.;
FT3sn0 = 0.;

(* FT3 first derivative evaluated at zero momentum transfer -- (d [FT3] / d q^2) @ q^2 = 0 *)
FT3upPrime = 0.;
FT3dpPrime = 0.;
FT3spPrime = 0.;

FT3unPrime = 0.;
FT3dnPrime = 0.;
FT3snPrime = 0.;

(* ======================= Rayleigh ======================= *)

(* FF(q^2) *)
FFp = None;
FFn = None;

(* FF evaluated at zero momentum transfer -- FF @ q^2 = 0 *)
FFp0 = 4.7*^-7;
FFn0 = 1.5*^-6;

(* FF first derivative -- (d [FF] / d q^2) @ q^2 = 0 *)
FFpPrime = 0.0;
FFnPrime = 0.0;

(* FFtilde(q^2) *)
FFtildep = None;
FFtilden = None;

(* FFtilde evaluated at zero momentum transfer -- FFtilde @ q^2 = 0 *)
FFtildep0 = -3.83*^-6;
FFtilden0 = 3.9*^-7;

(* FFtilde first derivative -- (d [FFtilde] / d q^2) @ q^2 = 0 *)
FFtildepPrime = 0.0;
FFtildenPrime = 0.0;

End[]
EndPackage[]
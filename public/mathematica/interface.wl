(* ::Package:: *)

BeginPackage["MuonConverter`"];

(* JMStoMuonConverter *)
JMStoMuonConverterWCxf::usage = "Translates WCs from the JMS -> MuonConverter basis";
MuonConverterWCxfToMuonConverterCemuqq::usage = "Convert from WCxf naming conventions to C^{(d)}_{e,mu,q,q} used internally see App. E of __insert arXiV no.__";
NucleonToIsospinBasis::usage = "Convert from the proton/neutron basis to the isospin basis used by Mu2e_NRET";
JMStoRET::usage = "Compose a series of translations and return RET d-coefficients for input into Mu2e_NRET";
HMMRZtoWET::usage = "Compose a series of translations and return WET d-coefficients for input into Mu2e_NRET";
ComputeRate::usage = "Compute the rate for a given set of WCs";

Begin["`Private`"];

JMStoMuonConverterWCxf[WCsJMS_] := Module[{WCsMC},
(* Initialize a MuonConverter Association *)
WCsMC = <| "Tegamma_12"->0., "ATegamma_12"->0., "VVeu_1211"->0., "VVed_1211"->0., "VVed_1222"->0., "AVVeu_1211"->0., "AVVed_1211"->0., "AVVed_1222"->0.,
           "VAVeu_1211"->0., "VAVed_1211"->0., "VAVed_1222"->0., "AVAVeu_1211"->0., "AVAVed_1211"->0., "AVAVed_1222"->0., "SSeu_1211"->0., "SSed_1211"->0.,
           "SSed_1222"->0., "ASeu_1211"->0., "ASed_1211"->0., "ASed_1222"->0., "SAeu_1211"->0., "SAed_1211"->0., "SAed_1222"->0., "AAeu_1211"->0.,
           "AAed_1211"->0., "AAed_1222"->0., "TTeu_1211"->0., "TTed_1211"->0., "TTed_1222"->0., "ATTeu_1211"->0., "ATTed_1211"->0., "ATTed_1222"->0. |>;

(* Translate WCsJMS to WCsMC *)
WCsMC["Tegamma_12"] = 4. * \[Pi]^2 * WCsJMS["egamma_12"];
WCsMC["ATegamma_12"] = - 4. * \[Pi]^2 * I * WCsJMS["egamma_12"];

WCsMC["VVeu_1211"] = (WCsJMS["VeuLL_1211"] + WCsJMS["VeuRR_1211"] + WCsJMS["VeuLR_1211"] + WCsJMS["VueLR_1112"]) / 4.;
WCsMC["VVed_1211"] = (WCsJMS["VedLL_1211"] + WCsJMS["VedRR_1211"] + WCsJMS["VedLR_1211"] + WCsJMS["VdeLR_1112"]) / 4.;
WCsMC["VVed_1222"] = (WCsJMS["VedLL_1222"] + WCsJMS["VedRR_1222"] + WCsJMS["VedLR_1222"] + WCsJMS["VdeLR_2212"]) / 4.;

WCsMC["AVVeu_1211"] = (-WCsJMS["VeuLL_1211"] + WCsJMS["VeuRR_1211"] - WCsJMS["VeuLR_1211"] + WCsJMS["VueLR_1112"]) / 4.;
WCsMC["AVVed_1211"] = (-WCsJMS["VedLL_1211"] + WCsJMS["VedRR_1211"] - WCsJMS["VedLR_1211"] + WCsJMS["VdeLR_1112"]) / 4.;
WCsMC["AVVed_1222"] = (-WCsJMS["VedLL_1222"] + WCsJMS["VedRR_1222"] - WCsJMS["VedLR_1222"] + WCsJMS["VdeLR_2212"]) / 4.;

WCsMC["VAVeu_1211"] = (-WCsJMS["VeuLL_1211"] + WCsJMS["VeuRR_1211"] + WCsJMS["VeuLR_1211"] - WCsJMS["VueLR_1112"]) / 4.;
WCsMC["VAVed_1211"] = (-WCsJMS["VedLL_1211"] + WCsJMS["VedRR_1211"] + WCsJMS["VedLR_1211"] - WCsJMS["VdeLR_1112"]) / 4.;
WCsMC["VAVed_1222"] = (-WCsJMS["VedLL_1222"] + WCsJMS["VedRR_1222"] + WCsJMS["VedLR_1222"] - WCsJMS["VdeLR_2212"]) / 4.;

WCsMC["AVAVeu_1211"] = (WCsJMS["VeuLL_1211"] + WCsJMS["VeuRR_1211"] - WCsJMS["VeuLR_1211"] - WCsJMS["VueLR_1112"]) / 4.;
WCsMC["AVAVed_1211"] = (WCsJMS["VedLL_1211"] + WCsJMS["VedRR_1211"] - WCsJMS["VedLR_1211"] - WCsJMS["VdeLR_1112"]) / 4.;
WCsMC["AVAVed_1222"] = (WCsJMS["VedLL_1222"] + WCsJMS["VedRR_1222"] - WCsJMS["VedLR_1222"] - WCsJMS["VdeLR_2212"]) / 4.;

WCsMC["SSeu_1211"] = (WCsJMS["SeuRL_1211"] + WCsJMS["SeuRR_1211"]) / 4.;
WCsMC["SSed_1211"] = (WCsJMS["SedRL_1211"] + WCsJMS["SedRR_1211"]) / 4.;
WCsMC["SSed_1222"] = (WCsJMS["SedRL_1211"] + WCsJMS["SedRR_1222"]) / 4.;

WCsMC["ASeu_1211"] = -I * (WCsJMS["SeuRL_1211"] + WCsJMS["SeuRR_1211"]) / 4.;
WCsMC["ASed_1211"] = -I * (WCsJMS["SedRL_1211"] + WCsJMS["SedRR_1211"]) / 4.;
WCsMC["ASed_1222"] = -I * (WCsJMS["SedRL_1222"] + WCsJMS["SedRR_1222"]) / 4.;

WCsMC["SAeu_1211"] = I * (WCsJMS["SeuRL_1211"] - WCsJMS["SeuRR_1211"]) / 4.;
WCsMC["SAed_1211"] = I * (WCsJMS["SedRL_1211"] - WCsJMS["SedRR_1211"]) / 4.;
WCsMC["SAed_1222"] = I * (WCsJMS["SedRL_1222"] - WCsJMS["SedRR_1222"]) / 4.;

WCsMC["AAeu_1211"] = (WCsJMS["SeuRL_1211"] - WCsJMS["SeuRR_1211"]) / 4.;
WCsMC["AAed_1211"] = (WCsJMS["SedRL_1211"] - WCsJMS["SedRR_1211"]) / 4.;
WCsMC["AAed_1222"] = (WCsJMS["SedRL_1222"] - WCsJMS["SedRR_1222"]) / 4.;

WCsMC["TTeu_1211"] = WCsJMS["TeuRR_1211"] / 4.;
WCsMC["TTed_1211"] = WCsJMS["TedRR_1211"] / 4.;
WCsMC["TTed_1222"] = WCsJMS["TedRR_1222"] / 4.;

WCsMC["ATTeu_1211"] = -I * WCsJMS["TeuRR_1211"] / 4.;
WCsMC["ATTed_1211"] = -I * WCsJMS["TedRR_1211"] / 4.;
WCsMC["ATTed_1222"] = -I * WCsJMS["TedRR_1222"] / 4.;

WCsMC];

MuonConverterWCxfToMuonConverterCemuqq[WCsMC_] := Module[{Cemuqq},
(* Convert from WCxf naming conventions to C^6_{e,mu,q,q} used internally see App. E of __insert arXiV no.__ *)
Cemuqq = <| "C51"  -> WCsMC["Tegamma_12"],  "C52"  -> WCsMC["ATegamma_12"], 
            "C61u" -> WCsMC["VVeu_1211"],   "C61d" -> WCsMC["VVed_1211"],   "C61s" -> WCsMC["VVed_1222"], 
            "C62u" -> WCsMC["AVVeu_1211"],  "C62d" -> WCsMC["AVVed_1211"],  "C62s" -> WCsMC["AVVed_1222"], 
            "C63u" -> WCsMC["VAVeu_1211"],  "C63d" -> WCsMC["VAVed_1211"],  "C63s" -> WCsMC["VAVed_1222"], 
            "C64u" -> WCsMC["AVAVeu_1211"], "C64d" -> WCsMC["AVAVed_1211"], "C64s" -> WCsMC["AVAVed_1222"],
            "C65u" -> WCsMC["SSeu_1211"],   "C65d" -> WCsMC["SSed_1211"],   "C65s" -> WCsMC["SSed_1222"],
            "C67u" -> WCsMC["ASeu_1211"],   "C67d" -> WCsMC["ASed_1211"],   "C67s" -> WCsMC["ASed_1222"],
            "C66u" -> WCsMC["SAeu_1211"],   "C66d" -> WCsMC["SAed_1211"],   "C66s" -> WCsMC["SAed_1222"],
            "C68u" -> WCsMC["AAeu_1211"],   "C68d" -> WCsMC["AAed_1211"],   "C68s" -> WCsMC["AAed_1222"],
            "C69u" -> WCsMC["TTeu_1211"],   "C69d" -> WCsMC["TTed_1211"],   "C69s" -> WCsMC["TTed_1222"],
            "C610u"-> WCsMC["ATTeu_1211"],  "C610d"-> WCsMC["ATTed_1211"], "C610s" -> WCsMC["ATTed_1222"] |>;

Cemuqq];

JMStoRET[WETJMS_, q_] := Module[{JMSkeys, WETJMSi, WETWCxf, WETinternal, RETisospinList},
JMSkeys = {"egamma_12", "VeuLL_1211", "VedLL_1211", "VedLL_1222", 
           "VeuRR_1211", "VedRR_1211", "VedRR_1222", "VeuLR_1211",
           "VedLR_1211", "VedLR_1222", "VueLR_1112", "VdeLR_1112",
           "VdeLR_2212", "SeuRL_1211", "SedRL_1211", "SedRL_1222",
           "SeuRR_1211", "SedRR_1211", "SedRR_1222", "TeuRR_1211",
           "TedRR_1211", "TedRR_1222"};
WETJMSi = <|Table[If[Not[KeyExistsQ[WETJMS, JMSkeys[[i]]]], JMSkeys[[i]] -> 0., JMSkeys[[i]] -> WETJMS[JMSkeys[[i]]]], {i, 1, Length[JMSkeys]}]|>;
(* Translate to MuonConverter WCxf WET basis *)
WETWCxf = JMStoMuonConverterWCxf[WETJMSi];
(* Translate to MuonConverter internal Cmueqq WET basis *)
WETinternal = MuonConverterWCxfToMuonConverterCemuqq[WETWCxf];
(* Hadronize to RET and return list for input into Mu2e_NRET *)
RETisospinList = HadronizeWET[WETinternal, q];

RETisospinList]

HMMRZtoWET[WETHMMRZ_, q_] := Module[{WETinternal, RETisospinList},
(* Translate to MuonConverter internal Cmueqq WET basis *)
WETinternal = MuonConverterWCxfToMuonConverterCemuqq[WETWCxf];
(* Hadronize to RET and return list for input into Mu2e_NRET *)
RETisospinList = HadronizeWET[WETinternal, q];

RETisospinList]

ComputeRate[WET_, q_, isotope_, interaction_, oscb_, isochar_, basis_:"JMS", plots_:"none", printDetails_:False] := Block[{data, ds},
(* Hadronize the WET association *)
If[basis == "JMS", ds = JMStoRET[WET, q], ds = HMMRZtoWET[WET, q]];

(* Initialize data Association *)
data = <|"element" -> element, "isotope" -> isotope, "interaction" -> interaction, "oscb" -> oscb,
         "plots" -> plots, "mL" -> mL * 10^3 (* MeV *), "muonlower" -> 1, "ds" -> ds, 
         "dummy" -> 0 (*avoid trailing comma*)|>;

(* Initialize Mu2e_NRET *)
Global`batch[data, printDetails];
(* Compute the rate *)
(*DecayRate = Global`calcDecayRate[False];*)
(* Account for the normalization used in the Mu2eNRET code *)
v = 246.2; (* GeV *)
{v^4 * Global`Decayrate (* 1/s *), v^4 * Global`BranchingRatio}
];

End[]
EndPackage[]
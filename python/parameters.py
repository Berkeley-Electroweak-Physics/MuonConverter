import sys
import numpy as np
import warnings

class NumericalInput:
    def __init__(self, my_input_dict=None):
        """ 
        The numerical input for MuonConverter. 
        All masses are in units of GeV and most recent lattice results (as of March 2024).

        Default values can be overriden by providing the optional 
        dictionary 'my_input_dict' which specifies the numerical 
        values of (a subset of) input parameters.
        """
        self.input_parameters = {}

        # ----------------------- Couplings ----------------------- #

        # The strong coupling constant @ MZ
        self.input_parameters['asMZ'] = 0.1179

        # The Fermi constant
        self.input_parameters['GF'] = 1.1663787e-5

        # The inverse of QED alpha @ MZ
        self.input_parameters['aMZinv'] = 127.952

        # The inverse of QED alpha @ mtau
        self.input_parameters['amtauinv'] = 133.472

        # The inverse of QED alpha at low scales
        self.input_parameters['alowinv'] = 137.035999084

        # Sine-squared of the weak mixing angle (MS-bar)
        self.input_parameters['sw2_MSbar'] = 0.23121

        # ----------------------- Boson masses ----------------------- #

        # Z boson mass 
        self.input_parameters['MZ'] = 91.1876 # GeV

        # Higgs boson mass
        self.input_parameters['Mh'] = 125.10 # GeV

        # W boson mass
        self.input_parameters['MW'] = 80.379 # GeV

        # ----------------------- Lepton masses ----------------------- #

        # Tau mass
        self.input_parameters['mtau'] = 1.77686 # GeV

        # Muon mass
        self.input_parameters['mmu'] = 105.6583715e-3 # GeV

        # Electron mass
        self.input_parameters['me'] = 0.000510998928 # GeV

        # ----------------------- Baryon masses ----------------------- #

        # Proton mass
        self.input_parameters['mproton'] = 938.272081e-3 # GeV

        # Neutron mass
        self.input_parameters['mneutron'] = 939.565413e-3 # GeV

        # ----------------------- Meson masses ----------------------- #

        # Pion mass
        self.input_parameters['mpi0'] =134.98e-3 # GeV

        # Eta mass
        self.input_parameters['meta'] = 547.862e-3 # GeV

        # ----------------------- Quark masses ----------------------- #

        # Strange quark mass, MS-bar at 2 GeV
        self.input_parameters['ms_at_2GeV'] = 0.093 # GeV

        # Down quark mass, MS-bar at 2 GeV
        self.input_parameters['md_at_2GeV'] = 0.00467 # GeV

        # Up quark mass, MS-bar at 2 GeV
        self.input_parameters['mu_at_2GeV'] = 0.00216 # GeV

        # ----------------------- Quark charges ----------------------- #

        # Up quark charge
        self.input_parameters['qu'] = 2./3.

        # Down quark charge
        self.input_parameters['qd'] = -1./3.

        # Strange quark charge
        self.input_parameters['qs'] = -1./3.


        # ----------------------- Leptonic scale ----------------------- #

        # Leptonic scale (set to the nucleon mass by default, final capture and decay rate computations are independent of this value)
        self.input_parameters['mL'] = (1 / 2) * (self.input_parameters['mproton'] + self.input_parameters['mneutron'])

        # ----------------------- Form factor parameters ----------------------- #

        # ======================= Global ======================= #
        
        # Ratio of up and down quark masses (rud = mu / md)
        self.input_parameters['rud'] = 0.474

        # Ratio of strange quark mass with sum of up and down masses (rs = (2 ms) / (mu + md))
        self.input_parameters['rs'] = 27.33

        # ======================= Vector ======================= #

        # F1(q^2)
        self.input_parameters['F1up'] = None
        self.input_parameters['F1dp'] = None
        self.input_parameters['F1sp'] = None

        self.input_parameters['F1un'] = None
        self.input_parameters['F1dn'] = None
        self.input_parameters['F1sn'] = None

        # F1 evaluated at zero momentum transfer -- F1 @ q^2 = 0
        self.input_parameters['F1up_0'] = None
        self.input_parameters['F1dp_0'] = None
        self.input_parameters['F1sp_0'] = None

        self.input_parameters['F1un_0'] = None
        self.input_parameters['F1dn_0'] = None
        self.input_parameters['F1sn_0'] = None

        # F1 first derivative evaluated at zero momentum transfer -- (d [F1] / d q^2) @ q^2 = 0
        self.input_parameters['F1up_prime'] = None
        self.input_parameters['F1dp_prime'] = None
        self.input_parameters['F1sp_prime'] = None

        self.input_parameters['F1un_prime'] = None
        self.input_parameters['F1dn_prime'] = None
        self.input_parameters['F1sn_prime'] = None

        # F2(q^2)
        self.input_parameters['F2up'] = None
        self.input_parameters['F2dp'] = None
        self.input_parameters['F2sp'] = None

        self.input_parameters['F2un'] = None
        self.input_parameters['F2dn'] = None
        self.input_parameters['F2sn'] = None

        # F2 evaluated at zero momentum transfer -- F2 @ q^2 = 0
        self.input_parameters['F2up_0'] = None
        self.input_parameters['F2dp_0'] = None
        self.input_parameters['F2sp_0'] = None

        self.input_parameters['F2un_0'] = None
        self.input_parameters['F2dn_0'] = None
        self.input_parameters['F2sn_0'] = None

        # F2 first derivative evaluated at zero momentum transfer -- (d [F2] / d q^2) @ q^2 = 0
        self.input_parameters['F2up_prime'] = None
        self.input_parameters['F2dp_prime'] = None
        self.input_parameters['F2sp_prime'] = None

        self.input_parameters['F2un_prime'] = None
        self.input_parameters['F2dn_prime'] = None
        self.input_parameters['F2sn_prime'] = None

        # F3(q^2) - set to zero by default
        self.input_parameters['F3up'] = None
        self.input_parameters['F3dp'] = None
        self.input_parameters['F3sp'] = None

        self.input_parameters['F3un'] = None
        self.input_parameters['F3dn'] = None
        self.input_parameters['F3sn'] = None

        # Proton nuclear dipole moment
        self.input_parameters['mup'] = 2.792847

        # Neutron nuclear dipole moment
        self.input_parameters['mun'] = -1.91304

        # Strange nuclear dipole moment
        self.input_parameters['mus'] = -0.036

        # Proton electric charge radius squared [1/GeV^2]
        self.input_parameters['rEp2'] = 18.16

        # Neutron electric charge radius squared [1/GeV^2]
        self.input_parameters['rEn2'] = -2.966

        # Strange electric charge radius squared [1/GeV^2]
        self.input_parameters['rEs2'] = -0.115

        # Proton magnetic charge radius squared [1/GeV^2]
        self.input_parameters['rMp2'] = 18.6

        # Neutron magnetic charge radius squared [1/GeV^2]
        self.input_parameters['rMn2'] = 19.1

        # Strange magnetic charge radius squared [1/GeV^2]
        self.input_parameters['rMs2'] = -0.26

        # ======================= Axial vector ======================= #

        # FA(q^2)
        self.input_parameters['FAup'] = None
        self.input_parameters['FAdp'] = None
        self.input_parameters['FAsp'] = None

        self.input_parameters['FAun'] = None
        self.input_parameters['FAdn'] = None
        self.input_parameters['FAsn'] = None

        # FA evaluated at zero momentum transfer -- FA @ q^2 = 0
        self.input_parameters['FAup_0'] = None
        self.input_parameters['FAdp_0'] = None
        self.input_parameters['FAsp_0'] = None

        self.input_parameters['FAun_0'] = None
        self.input_parameters['FAdn_0'] = None
        self.input_parameters['FAsn_0'] = None

        # FA first derivative evaluated at zero momentum transfer -- (d [FA] / d q^2) @ q^2 = 0
        self.input_parameters['FAup_prime'] = None
        self.input_parameters['FAdp_prime'] = None
        self.input_parameters['FAsp_prime'] = None

        self.input_parameters['FAun_prime'] = None
        self.input_parameters['FAdn_prime'] = None
        self.input_parameters['FAsn_prime'] = None

        # FPprimed(q^2)
        self.input_parameters['FPprimedup'] = None
        self.input_parameters['FPprimeddp'] = None
        self.input_parameters['FPprimedsp'] = None

        self.input_parameters['FPprimedun'] = None
        self.input_parameters['FPprimeddn'] = None
        self.input_parameters['FPprimedsn'] = None

        # FPprimed residua of pion pole contribution
        self.input_parameters['FPprimedup_pion'] = None
        self.input_parameters['FPprimeddp_pion'] = None
        self.input_parameters['FPprimedsp_pion'] = None

        self.input_parameters['FPprimedun_pion'] = None
        self.input_parameters['FPprimeddn_pion'] = None
        self.input_parameters['FPprimedsn_pion'] = None

        # FPprimed residua of eta pole contribution
        self.input_parameters['FPprimedup_eta'] = None
        self.input_parameters['FPprimeddp_eta'] = None
        self.input_parameters['FPprimedsp_eta'] = None

        self.input_parameters['FPprimedun_eta'] = None
        self.input_parameters['FPprimeddn_eta'] = None
        self.input_parameters['FPprimedsn_eta'] = None

        # FPprimed evaluated at zero momentum transfer -- FPprimed @ q^2 = 0
        self.input_parameters['FPprimedup_0'] = 119
        self.input_parameters['FPprimeddp_0'] = -130
        self.input_parameters['FPprimedsp_0'] = -1.6

        self.input_parameters['FPprimedun_0'] = -130
        self.input_parameters['FPprimeddn_0'] = 119
        self.input_parameters['FPprimedsn_0'] = -1.6

        self.input_parameters['FA3up'] = None
        self.input_parameters['FA3dp'] = None
        self.input_parameters['FA3sp'] = None

        self.input_parameters['FA3un'] = None
        self.input_parameters['FA3dn'] = None
        self.input_parameters['FA3sn'] = None

        # gA
        self.input_parameters['gA'] = 1.2754

        # Delta u + Delta d
        self.input_parameters['DeltaSigmaud'] = 0.397

        # Deltas
        self.input_parameters['Deltas'] = -0.045

        # Axial charge radius (u - d) [1/GeV^2]
        self.input_parameters['rA2uMinusd'] = 10.1

        # Axial charge radius (u + d) [1/GeV^2]
        self.input_parameters['rA2uPlusd'] = 12.58

        # Axial charge radius (s) [1/GeV^2]
        self.input_parameters['rA2s'] = 12.

        # Axial FPprime correction
        self.input_parameters['DeltaGT8'] = 0.50

        # ======================= Scalar ======================= #

        # FS(q^2)
        self.input_parameters['FSup'] = None
        self.input_parameters['FSdp'] = None
        self.input_parameters['FSsp'] = None

        self.input_parameters['FSun'] = None
        self.input_parameters['FSdn'] = None
        self.input_parameters['FSsn'] = None

        # FS evaluated at zero momentum transfer -- FS @ q^2 = 0
        self.input_parameters['FSup_0'] = None
        self.input_parameters['FSdp_0'] = None
        self.input_parameters['FSsp_0'] = None

        self.input_parameters['FSun_0'] = None
        self.input_parameters['FSdn_0'] = None
        self.input_parameters['FSsn_0'] = None

        # FS first derivative evaluated at zero momentum transfer -- (d [FS] / d q^2) @ q^2 = 0
        self.input_parameters['FSup_prime'] = 0.72
        self.input_parameters['FSdp_prime'] = 0.59
        self.input_parameters['FSsp_prime'] = 0.17

        self.input_parameters['FSun_prime'] = 0.59
        self.input_parameters['FSdn_prime'] = 0.72
        self.input_parameters['FSsn_prime'] = 0.17

        # FP(q^2)
        self.input_parameters['FPup'] = None
        self.input_parameters['FPdp'] = None
        self.input_parameters['FPsp'] = None

        self.input_parameters['FPun'] = None
        self.input_parameters['FPdn'] = None
        self.input_parameters['FPsn'] = None

        # FP residua of pion pole contribution
        self.input_parameters['FPup_pion'] = None
        self.input_parameters['FPdp_pion'] = None
        self.input_parameters['FPsp_pion'] = None

        self.input_parameters['FPun_pion'] = None
        self.input_parameters['FPdn_pion'] = None
        self.input_parameters['FPsn_pion'] = None

        # FP residua of eta pole contribution
        self.input_parameters['FPup_eta'] = None
        self.input_parameters['FPdp_eta'] = None
        self.input_parameters['FPsp_eta'] = None

        self.input_parameters['FPun_eta'] = None
        self.input_parameters['FPdn_eta'] = None
        self.input_parameters['FPsn_eta'] = None

        self.input_parameters['sigmapiNtilde'] = 48e-3

        self.input_parameters['c5hat'] = -0.51e-3

        # sigmas
        self.input_parameters['sigmas'] = 43.3e-3

        # ======================= Gluonic ======================= #

        # FG(q^2)
        self.input_parameters['FGp'] = None
        self.input_parameters['FGn'] = None

        # FG evaluated at zero momentum transfer -- FG @ q^2 = 0
        self.input_parameters['FGp_0'] = -50.4e-3
        self.input_parameters['FGn_0'] = -50.4e-3

        # FG first derivative -- (d [FG] / d q^2) @ q^2 = 0
        self.input_parameters['FGp_prime'] = -0.14
        self.input_parameters['FGn_prime'] = -0.14

        # FGtilde(q^2)
        self.input_parameters['FGtildep'] = None
        self.input_parameters['FGtilden'] = None

        # FGtilde evaluated at zero momentum transfer -- FGtilde @ q^2 = 0
        self.input_parameters['FGtildep_0'] = None
        self.input_parameters['FGtilden_0'] = None

        # FGtilde residua of pion pole contribution
        self.input_parameters['FGtildep_pion'] = None
        self.input_parameters['FGtilden_pion'] = None

        # FGtilde residua of eta pole contribution
        self.input_parameters['FGtildep_eta'] = None
        self.input_parameters['FGtilden_eta'] = None

        # Strong coupling at 2 GeV
        self.input_parameters['alphas_at_2GeV'] = 0.297

        # mG
        self.input_parameters['mG'] = 0.823

        # ======================= Tensor ======================= #

        # FT0(q^2)
        self.input_parameters['FT0up'] = None
        self.input_parameters['FT0dp'] = None
        self.input_parameters['FT0sp'] = None

        self.input_parameters['FT0un'] = None
        self.input_parameters['FT0dn'] = None
        self.input_parameters['FT0sn'] = None

        # FT0 evaluated at zero momentum transfer -- FT0 @ q^2 = 0
        self.input_parameters['FT0up_0'] = 0.784 # gTu
        self.input_parameters['FT0dp_0'] = -0.204 # gTd
        self.input_parameters['FT0sp_0'] = -2.7e-3 # gTs

        self.input_parameters['FT0un_0'] = -0.204 # gTd
        self.input_parameters['FT0dn_0'] = 0.784 # gTu
        self.input_parameters['FT0sn_0'] = -2.7e-3 # gTs

        # FT0 first derivative evaluated at zero momentum transfer -- (d [FT0] / d q^2) @ q^2 = 0
        self.input_parameters['FT0up_prime'] = 0.54
        self.input_parameters['FT0dp_prime'] = -0.11
        self.input_parameters['FT0sp_prime'] = -0.0014

        self.input_parameters['FT0un_prime'] = -0.11
        self.input_parameters['FT0dn_prime'] = 0.54
        self.input_parameters['FT0sn_prime'] = -0.0014

        # FT1(q^2)
        self.input_parameters['FT1up'] = None
        self.input_parameters['FT1dp'] = None
        self.input_parameters['FT1sp'] = None

        self.input_parameters['FT1un'] = None
        self.input_parameters['FT1dn'] = None
        self.input_parameters['FT1sn'] = None

        # FT1 evaluated at zero momentum transfer -- FT1 @ q^2 = 0
        self.input_parameters['FT1up_0'] = -3.0
        self.input_parameters['FT1dp_0'] = 1.0
        self.input_parameters['FT1sp_0'] = 0.018

        self.input_parameters['FT1un_0'] = 1.0
        self.input_parameters['FT1dn_0'] = -3.0
        self.input_parameters['FT1sn_0'] = 0.018

        # FT1 first derivative evaluated at zero momentum transfer -- (d [FT1] / d q^2) @ q^2 = 0
        self.input_parameters['FT1up_prime'] = -14.0
        self.input_parameters['FT1dp_prime'] = 5.0
        self.input_parameters['FT1sp_prime'] = 0.082

        self.input_parameters['FT1un_prime'] = 5.0
        self.input_parameters['FT1dn_prime'] = -14.0
        self.input_parameters['FT1sn_prime'] = 0.082

        # FT2(q^2)
        self.input_parameters['FT2up'] = None
        self.input_parameters['FT2dp'] = None
        self.input_parameters['FT2sp'] = None

        self.input_parameters['FT2un'] = None
        self.input_parameters['FT2dn'] = None
        self.input_parameters['FT2sn'] = None

        # FT2 evaluated at zero momentum transfer -- FT2 @ q^2 = 0
        self.input_parameters['FT2up_0'] = -0.1
        self.input_parameters['FT2dp_0'] = 0.6
        self.input_parameters['FT2sp_0'] = 0.004

        self.input_parameters['FT2un_0'] = 0.6
        self.input_parameters['FT2dn_0'] = -0.1
        self.input_parameters['FT2sn_0'] = 0.004

        # FT2 first derivative evaluated at zero momentum transfer -- (d [FT2] / d q^2) @ q^2 = 0
        self.input_parameters['FT2up_prime'] = -1.8
        self.input_parameters['FT2dp_prime'] = 2.1
        self.input_parameters['FT2sp_prime'] = 0.015

        self.input_parameters['FT2un_prime'] = 2.1
        self.input_parameters['FT2dn_prime'] = -1.8
        self.input_parameters['FT2sn_prime'] = 0.015

        self.input_parameters['FT3up'] = None
        self.input_parameters['FT3dp'] = None
        self.input_parameters['FT3sp'] = None

        self.input_parameters['FT3un'] = None
        self.input_parameters['FT3dn'] = None
        self.input_parameters['FT3sn'] = None

        # ======================= Rayleigh ======================= #

        # FF(q^2)
        self.input_parameters['FFp'] = None
        self.input_parameters['FFn'] = None 

        # FF evaluated at zero momentum transfer -- FF @ q^2 = 0
        self.input_parameters['FFp_0'] = 4.7e-7
        self.input_parameters['FFn_0'] = 1.5e-6

        # FF first derivative -- (d [FF] / d q^2) @ q^2 = 0
        self.input_parameters['FFp_prime'] = None
        self.input_parameters['FFn_prime'] = None

        # FFtilde(q^2)
        self.input_parameters['FFtildep'] = None
        self.input_parameters['FFtilden'] = None

        # FFtilde evaluated at zero momentum transfer -- FFtilde @ q^2 = 0
        self.input_parameters['FFtildep_0'] = -3.83e-6
        self.input_parameters['FFtilden_0'] = 3.9e-7

        # FFtilde first derivative -- (d [FFtilde] / d q^2) @ q^2 = 0
        self.input_parameters['FFtildep_prime'] = None
        self.input_parameters['FFtilden_prime'] = None


        # Update primary input parameters with user-specified values (optional)
        if my_input_dict is None:
            pass
        else:
            # Issue a user warning if a key is not defined:
            for input_key in my_input_dict.keys():
                if input_key in self.input_parameters.keys():
                    pass
                else:
                    raise Exception(input_key + ' is not a valid key for an input parameter.')
            # Create the dictionary.
            self.input_parameters.update(my_input_dict)
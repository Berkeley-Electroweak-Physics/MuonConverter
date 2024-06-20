import sys
import numpy as np
import warnings
import os.path
from parameters import NumericalInput
from formfactors import *

class Hadronization:
    def __init__(self, Cmueqq_dict, input_dict):
        """
        Class for Wilson coefficients in 3 flavor QCD + QED.

        Cmueqq_dict (dictionary): dictionary containing the initial conditions
        of the dimension-five to dimension-seven three-flavor-QCD Wilson coefficients 
        of the form {'C51' : value, 'C52' : value, ...}. By default all values are
        set to zero. The Wilson coefficients should be specified in the MS-bar scheme 
        at 2 GeV.

        Cmueqq_dict (dictionary): dictionary containing desired input parameters.
        """
        # The dictionary of input parameters
        if input_dict == None:
            self.ip = NumericalInput().input_parameters
        else:
            self.ip = NumericalInput(my_input_dict = input_dict).input_parameters

        # Supported Wilson coefficients (following Cmueqq notation from __insert arXiv no.__)
        self.wc_keys = ['C51', 'C52', 
                        'C61u', 'C61d', 'C61s', 'C62u', 'C62d', 'C62s',
                        'C63u', 'C63d', 'C63s', 'C64u', 'C64d', 'C64s',
                        'C65u', 'C65d', 'C65s', 'C66u', 'C66d', 'C66s',
                        'C67u', 'C67d', 'C67s', 'C68u', 'C68d', 'C68s',
                        'C69u', 'C69d', 'C69s', 'C610u', 'C610d', 'C610s',
                        'C71', 'C72', 'C73', 'C74', 'C75', 'C76', 'C77', 'C78',
                        'C79u', 'C79d', 'C79s', 'C710u', 'C710d', 'C710s',
                        'C711u', 'C711d', 'C711s', 'C712u', 'C712d', 'C712s',
                        'C713u', 'C713d', 'C713s', 'C714u', 'C714d', 'C714s',
                        'C715u', 'C715d', 'C715s', 'C716u', 'C716d', 'C716s']

        # Initialize coefficient dictionary
        self.Cmueqq_dict = {}

        # Issue a user warning if a key is not defined:
        for wc_name in Cmueqq_dict.keys():
            if wc_name in self.wc_keys:
                pass
            else:
                warnings.warn('The key ' + wc_name + ' is not a valid key.')

        # Create the dictionary.
        for wc_name in self.wc_keys:
            if wc_name in Cmueqq_dict.keys():
                self.Cmueqq_dict[wc_name] = Cmueqq_dict[wc_name]
            else:
                self.Cmueqq_dict[wc_name] = 0.

    def hadronize_WET(self, q):
        """
        Match the relativistic 3-flavor parton basis (WET) to the relativistic (RET)
        isospin basis derived in Eq. ___ of __insert arXiv no.__.

        Input:
        q (numpy array): Outgoing 4-momentum of the electron
        mL (float): Arbitrary leptonic scale - see 2208.07945

        Ouput:
        Returns a dictionary contianing the relativistic Wilson ("d") coefficients in the isospin basis
        """
        #####################################################
        # ---------- Define the input parameters ---------- #
        #####################################################
        # Nucleon masses
        mp = self.ip['mproton']
        mn = self.ip['mneutron']

        # Electromagnetic coupling
        alpha = 1/self.ip['alowinv']

        # Quark masses at 2 GeV
        mu = self.ip['mu_at_2GeV']
        md = self.ip['md_at_2GeV']
        ms = self.ip['ms_at_2GeV']
        
        # Quark charged
        qu = self.ip['qu']
        qd = self.ip['qd']
        qs = self.ip['qs']

        # Lepton masses
        me = self.ip['me']
        mmu = self.ip['mmu']

        # Arbitrary leptonic scale
        mL = self.ip['mL']

        # Convinient notation (without E_bind)
        m_plus = mmu + me
        m_minus = mmu - me

        #####################################################
        # ------------ Initialize form factors ------------ #
        #####################################################
        
        # Vector
        F1up = F1('u', 'p', self.ip).form_factor(q)
        F1dp = F1('d', 'p', self.ip).form_factor(q)
        F1sp = F1('s', 'p', self.ip).form_factor(q)

        F1un = F1('u', 'n', self.ip).form_factor(q)
        F1dn = F1('d', 'n', self.ip).form_factor(q)
        F1sn = F1('s', 'n', self.ip).form_factor(q)

        F2up = F2('u', 'p', self.ip).form_factor(q)
        F2dp = F2('d', 'p', self.ip).form_factor(q)
        F2sp = F2('s', 'p', self.ip).form_factor(q)

        F2un = F2('u', 'n', self.ip).form_factor(q)
        F2dn = F2('d', 'n', self.ip).form_factor(q)
        F2sn = F2('s', 'n', self.ip).form_factor(q)

        F3up = F3('u', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        F3dp = F3('d', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        F3sp = F3('s', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        F3un = F3('u', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        F3dn = F3('d', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        F3sn = F3('s', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        # Axial-vector
        FAup = FA('u', 'p', self.ip).form_factor(q)
        FAdp = FA('d', 'p', self.ip).form_factor(q)
        FAsp = FA('s', 'p', self.ip).form_factor(q)

        FAun = FA('u', 'n', self.ip).form_factor(q)
        FAdn = FA('d', 'n', self.ip).form_factor(q)
        FAsn = FA('s', 'n', self.ip).form_factor(q)

        FPpup = FPprimed('u', 'p', self.ip).form_factor(q)
        FPpdp = FPprimed('d', 'p', self.ip).form_factor(q)
        FPpsp = FPprimed('s', 'p', self.ip).form_factor(q)

        FPpun = FPprimed('u', 'n', self.ip).form_factor(q)
        FPpdn = FPprimed('d', 'n', self.ip).form_factor(q)
        FPpsn = FPprimed('s', 'n', self.ip).form_factor(q)

        FA3up = FA3('u', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FA3dp = FA3('d', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FA3sp = FA3('s', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        FA3un = FA3('u', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FA3dn = FA3('d', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FA3sn = FA3('s', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        # Scalar
        FSup = FS('u', 'p', self.ip).form_factor(q)
        FSdp = FS('d', 'p', self.ip).form_factor(q)
        FSsp = FS('s', 'p', self.ip).form_factor(q)

        FSun = FS('u', 'n', self.ip).form_factor(q)
        FSdn = FS('d', 'n', self.ip).form_factor(q)
        FSsn = FS('s', 'n', self.ip).form_factor(q)

        # Psuedoscalar
        FPup = FP('u', 'p', self.ip).form_factor(q)
        FPdp = FP('d', 'p', self.ip).form_factor(q)
        FPsp = FP('s', 'p', self.ip).form_factor(q)

        FPun = FP('u', 'n', self.ip).form_factor(q)
        FPdn = FP('d', 'n', self.ip).form_factor(q)
        FPsn = FP('s', 'n', self.ip).form_factor(q)

        # GG
        FGp = FG('p', self.ip).form_factor(q)
        FGn = FG('n', self.ip).form_factor(q)

        # GGtilde
        FGtildep = FGtilde('p', self.ip).form_factor(q)
        FGtilden = FGtilde('n', self.ip).form_factor(q)

        # Tensor
        FT0up = FT0('u', 'p', self.ip).form_factor(q)
        FT0dp = FT0('d', 'p', self.ip).form_factor(q)
        FT0sp = FT0('s', 'p', self.ip).form_factor(q)

        FT0un = FT0('u', 'n', self.ip).form_factor(q)
        FT0dn = FT0('d', 'n', self.ip).form_factor(q)
        FT0sn = FT0('s', 'n', self.ip).form_factor(q)

        FT1up = FT1('u', 'p', self.ip).form_factor(q)
        FT1dp = FT1('d', 'p', self.ip).form_factor(q)
        FT1sp = FT1('s', 'p', self.ip).form_factor(q)

        FT1un = FT1('u', 'n', self.ip).form_factor(q)
        FT1dn = FT1('d', 'n', self.ip).form_factor(q)
        FT1sn = FT1('s', 'n', self.ip).form_factor(q)
        
        FT2up = FT2('u', 'p', self.ip).form_factor(q)
        FT2dp = FT2('d', 'p', self.ip).form_factor(q)
        FT2sp = FT2('s', 'p', self.ip).form_factor(q)

        FT2un = FT2('u', 'n', self.ip).form_factor(q)
        FT2dn = FT2('d', 'n', self.ip).form_factor(q)
        FT2sn = FT2('s', 'n', self.ip).form_factor(q)

        FT3up = FT3('u', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FT3dp = FT3('d', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FT3sp = FT3('s', 'p', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        FT3un = FT3('u', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FT3dn = FT3('d', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)
        FT3sn = FT3('s', 'n', self.ip).form_factor(q) # T-odd (set to zero by default in formfactor.py)

        # FF
        Fgammap = Fgamma('p', self.ip).form_factor(q) 
        Fgamman = Fgamma('n', self.ip).form_factor(q) 

        # FFtilde
        Fgammatildep = Fgammatilde('p', self.ip).form_factor(q) 
        Fgammatilden = Fgammatilde('n', self.ip).form_factor(q)

        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]

        # Hadronize to the nucleon RET basis
        RET_nucleon_dict = {
            'd1p':  + (1. / mu) * self.Cmueqq_dict['C65u'] * FSup + (1. / md) * self.Cmueqq_dict['C65d'] * FSdp + (1. / ms) * self.Cmueqq_dict['C65s'] * FSsp\
                    + self.Cmueqq_dict['C75'] * Fgammap + self.Cmueqq_dict['C71'] * FGp\
                    # T-odd contributions (zero by default)
                    - 1j * (m_minus / mp) * (self.Cmueqq_dict['C61u'] * F3up + self.Cmueqq_dict['C61d'] * F3dp + self.Cmueqq_dict['C61s'] * F3sp)\
                    - 1j * (mmu**2 - me**2) * (self.Cmueqq_dict['C69u'] * F3up + self.Cmueqq_dict['C69d'] * F3dp + self.Cmueqq_dict['C69s'] * F3sp),

            'd1n':  + (1. / mu) * self.Cmueqq_dict['C65u'] * FSun + (1. / md) * self.Cmueqq_dict['C65d'] * FSdn + (1. / ms) * self.Cmueqq_dict['C65s'] * FSsn\
                    + self.Cmueqq_dict['C75'] * Fgamman + self.Cmueqq_dict['C71'] * FGn\
                    # T-odd contributions (zero by default)
                    - 1j * (m_minus / mn) * (self.Cmueqq_dict['C61u'] * F3un + self.Cmueqq_dict['C61d'] * F3dn + self.Cmueqq_dict['C61s'] * F3sn)\
                    - 1j * (mmu**2 - me**2) * (self.Cmueqq_dict['C69u'] * F3un + self.Cmueqq_dict['C69d'] * F3dn + self.Cmueqq_dict['C69s'] * F3sn),

            'd2p':  + (1. / mu) * self.Cmueqq_dict['C67u'] * FPup + (1. / md) * self.Cmueqq_dict['C67d'] * FPdp + (1. / ms) * self.Cmueqq_dict['C67s'] * FPsp\
                    + self.Cmueqq_dict['C73'] * FGtildep + self.Cmueqq_dict['C77'] * Fgammatildep\
                    - 1j * (m_minus / (2 * mp)) * (self.Cmueqq_dict['C63u'] * FPpup + self.Cmueqq_dict['C63d'] * FPpdp + self.Cmueqq_dict['C63s'] * FPpsp)\
                    - 1j * ((m_plus * m_minus) / (2 * mp)) * (self.Cmueqq_dict['C711u'] * FPpup + self.Cmueqq_dict['C711d'] * FPpdp + self.Cmueqq_dict['C711s'] * FPpsp)\
                    # T-odd contributions (zero by default)
                    + (m_plus / mp) * (self.Cmueqq_dict['C62u'] * F3up + self.Cmueqq_dict['C62d'] * F3dp + self.Cmueqq_dict['C62s'] * F3sp)\
                    - 4j * m_minus * (self.Cmueqq_dict['C715u'] * FT3up + self.Cmueqq_dict['C715d'] * FT3dp + self.Cmueqq_dict['C715s'] * FT3sp),

            'd2n':  + (1. / mu) * self.Cmueqq_dict['C67u'] * FPun + (1. / md) * self.Cmueqq_dict['C67d'] * FPdn + (1. / ms) * self.Cmueqq_dict['C67s'] * FPsn\
                    + self.Cmueqq_dict['C73'] * FGtilden + self.Cmueqq_dict['C77'] * Fgammatilden\
                    - 1j * (m_minus / (2 * mn)) * (self.Cmueqq_dict['C63u'] * FPpun + self.Cmueqq_dict['C63d'] * FPpdn + self.Cmueqq_dict['C63s'] * FPpsn)\
                    - 1j * ((m_plus * m_minus) / (2 * mn)) * (self.Cmueqq_dict['C711u'] * FPpun + self.Cmueqq_dict['C711d'] * FPpdn + self.Cmueqq_dict['C711s'] * FPpsn)
                    # T-odd contributions (zero by default)
                    + (m_plus / mn) * (self.Cmueqq_dict['C62u'] * F3un + self.Cmueqq_dict['C62d'] * F3dn + self.Cmueqq_dict['C62s'] * F3sn)\
                    - 4j * m_minus * (self.Cmueqq_dict['C715u'] * FT3un + self.Cmueqq_dict['C715d'] * FT3dn + self.Cmueqq_dict['C715s'] * FT3sn),

            'd3p':  + (1. / mu) * self.Cmueqq_dict['C66u'] * FSup + (1. / md) * self.Cmueqq_dict['C66d'] * FSdp + (1. / ms) * self.Cmueqq_dict['C66s'] * FSsp\
                    + self.Cmueqq_dict['C72'] * FGp + self.Cmueqq_dict['C76'] * Fgammap\
                    # T-odd contributions (zero by default)
                    -1j * (mmu**2 - me**2) * (self.Cmueqq_dict['C710u'] * F3up + self.Cmueqq_dict['C710d'] * F3dp + self.Cmueqq_dict['C710s'] * F3sp),

            'd3n':  + (1. / mu) * self.Cmueqq_dict['C66u'] * FSun + (1. / md) * self.Cmueqq_dict['C66d'] * FSdn + (1. / ms) * self.Cmueqq_dict['C66s'] * FSsn\
                    + self.Cmueqq_dict['C72'] * FGn + self.Cmueqq_dict['C76'] * Fgamman\
                    # T-odd contributions (zero by default)
                    -1j * (mmu**2 - me**2) * (self.Cmueqq_dict['C710u'] * F3un + self.Cmueqq_dict['C710d'] * F3dn + self.Cmueqq_dict['C710s'] * F3sn),

            'd4p':  + (1. / mu) * self.Cmueqq_dict['C68u'] * FPup + (1. / md) * self.Cmueqq_dict['C68d'] * FPdp + (1. / ms) * self.Cmueqq_dict['C68s'] * FPsp\
                    + self.Cmueqq_dict['C74'] * FGtildep + self.Cmueqq_dict['C78'] * Fgammatildep\
                    + (m_plus / (2 * mp)) * (self.Cmueqq_dict['C64u'] * FPpup + self.Cmueqq_dict['C64d'] * FPpdp + self.Cmueqq_dict['C64s'] * FPpsp)\
                    - 1j * ((m_plus * m_minus) / (2 * mp)) * (self.Cmueqq_dict['C712u'] * FPpup + self.Cmueqq_dict['C712d'] * FPpdp + self.Cmueqq_dict['C712s'] * FPpsp)\
                    # T-odd contributions (zero by default)
                    + 4 * m_plus * (self.Cmueqq_dict['C716u'] * FT3up + self.Cmueqq_dict['C716d'] * FT3dp + self.Cmueqq_dict['C716s'] * FT3sp),

            'd4n':  + (1. / mu) * self.Cmueqq_dict['C68u'] * FPun + (1. / md) * self.Cmueqq_dict['C68d'] * FPdn + (1. / ms) * self.Cmueqq_dict['C68s'] * FPsn\
                    + self.Cmueqq_dict['C74'] * FGtilden + self.Cmueqq_dict['C78'] * Fgammatilden\
                    + (m_plus / (2 * mn)) * (self.Cmueqq_dict['C64u'] * FPpun + self.Cmueqq_dict['C64d'] * FPpdn + self.Cmueqq_dict['C64s'] * FPpsn)\
                    - 1j * ((m_plus * m_minus) / (2 * mn)) * (self.Cmueqq_dict['C712u'] * FPpun + self.Cmueqq_dict['C712d'] * FPpdn + self.Cmueqq_dict['C712s'] * FPpsn)\
                    # T-odd contributions (zero by default)
                    + 4 * m_plus * (self.Cmueqq_dict['C716u'] * FT3un + self.Cmueqq_dict['C716d'] * FT3dn + self.Cmueqq_dict['C716s'] * FT3sn),

            'd5p':  + self.Cmueqq_dict['C61u'] * F1up + self.Cmueqq_dict['C61d'] * F1dp + self.Cmueqq_dict['C61s'] * F1sp\
                    + (m_plus) * (self.Cmueqq_dict['C79u'] * F1up + self.Cmueqq_dict['C79d'] * F1dp + self.Cmueqq_dict['C79s'] * F1sp)\
                    - (qSq / (2 * mp)) * (self.Cmueqq_dict['C713u'] * (FT1up - 4 * FT2up) + self.Cmueqq_dict['C713d'] * (FT1dp - 4 * FT2dp) + self.Cmueqq_dict['C713s'] * (FT1sp - 4 * FT2sp)),

            'd5n':  + self.Cmueqq_dict['C61u'] * F1un + self.Cmueqq_dict['C61d'] * F1dn + self.Cmueqq_dict['C61s'] * F1sn\
                    + (m_plus) * (self.Cmueqq_dict['C79u'] * F1un + self.Cmueqq_dict['C79d'] * F1dn + self.Cmueqq_dict['C79s'] * F1sn)\
                    - (qSq / (2 * mn)) * (self.Cmueqq_dict['C713u'] * (FT1un - 4 * FT2un) + self.Cmueqq_dict['C713d'] * (FT1dn - 4 * FT2dn) + self.Cmueqq_dict['C713s'] * (FT1sn - 4 * FT2sn)),

            'd6p':  - (1./2.) * (self.Cmueqq_dict['C61u'] * F2up + self.Cmueqq_dict['C61d'] * F2dp + self.Cmueqq_dict['C61s'] * F2sp)\
                    - (1./2.) * m_plus * (self.Cmueqq_dict['C79u'] * F2up + self.Cmueqq_dict['C79d'] * F2dp + self.Cmueqq_dict['C79s'] * F2sp),

            'd6n':  - (1./2.) * (self.Cmueqq_dict['C61u'] * F2un + self.Cmueqq_dict['C61d'] * F2dn + self.Cmueqq_dict['C61s'] * F2sn)\
                    - (1./2.) * m_plus * (self.Cmueqq_dict['C79u'] * F2un + self.Cmueqq_dict['C79d'] * F2dn + self.Cmueqq_dict['C79s'] * F2sn),

            'd7p':  + self.Cmueqq_dict['C63u'] * FAup + self.Cmueqq_dict['C63d'] * FAdp + self.Cmueqq_dict['C63s'] * FAsp\
                    + m_plus * (self.Cmueqq_dict['C711u'] * FAup + self.Cmueqq_dict['C711d'] * FAdp + self.Cmueqq_dict['C711s'] * FAsp)\
                    # T-odd contribution (set to zero by default)
                    - 2 * (qSq / mp) * (self.Cmueqq_dict['C715u'] * FT3up + self.Cmueqq_dict['C715d'] * FT3dp + self.Cmueqq_dict['C715s'] * FT3sp),

            'd7n':  + self.Cmueqq_dict['C63u'] * FAun + self.Cmueqq_dict['C63d'] * FAdn + self.Cmueqq_dict['C63s'] * FAsn
                    + m_plus * (self.Cmueqq_dict['C711u'] * FAun + self.Cmueqq_dict['C711d'] * FAdn + self.Cmueqq_dict['C711s'] * FAsn)\
                    # T-odd contribution (set to zero by default)
                    - 2 * (qSq / mn) * (self.Cmueqq_dict['C715u'] * FT3un + self.Cmueqq_dict['C715d'] * FT3dn + self.Cmueqq_dict['C715s'] * FT3sn),

                    # T-odd contribution (set to zero by default)
            'd8p':  + m_plus * (self.Cmueqq_dict['C711u'] * FA3up + self.Cmueqq_dict['C711d'] * FA3dp + self.Cmueqq_dict['C711s'] * FA3sp)\
                    + self.Cmueqq_dict['C63u'] * FA3up + self.Cmueqq_dict['C63d'] * FA3dp + self.Cmueqq_dict['C63s'] * FA3sp,

                    # T-odd contribution (set to zero by default)
            'd8n':  + m_plus * (self.Cmueqq_dict['C711u'] * FA3un + self.Cmueqq_dict['C711d'] * FA3dn + self.Cmueqq_dict['C711s'] * FA3sn)\
                    + self.Cmueqq_dict['C63u'] * FA3un + self.Cmueqq_dict['C63d'] * FA3dn + self.Cmueqq_dict['C63s'] * FA3sn,

            'd9p':  - (alpha / np.pi) * self.Cmueqq_dict['C51'] * (mL / qSq) * (qu * F1up + qd * F1dp + qs * F1sp)\
                    - mL * (self.Cmueqq_dict['C79u'] * F1up + self.Cmueqq_dict['C79d'] * F1dp + self.Cmueqq_dict['C79s'] * F1sp),

            'd9n':  - (alpha / np.pi) * self.Cmueqq_dict['C51'] * (mL / qSq) * (qu * F1un + qd * F1dn + qs * F1sn)\
                    - mL * (self.Cmueqq_dict['C79u'] * F1un + self.Cmueqq_dict['C79d'] * F1dn + self.Cmueqq_dict['C79s'] * F1sn),

            'd10p': + (alpha / (2 * np.pi)) * self.Cmueqq_dict['C51'] * (mL / qSq) * (qu * F2up + qd * F2dp + qs * F2sp)\
                    + (mL / 2.) * (self.Cmueqq_dict['C79u'] * F2up + self.Cmueqq_dict['C79d'] * F2dp + self.Cmueqq_dict['C79s'] * F2sp),

            'd10n': + (alpha / (2 * np.pi)) * self.Cmueqq_dict['C51'] * (mL / qSq) * (qu * F2un + qd * F2dn + qs * F2sn)\
                    + (mL / 2.) * (self.Cmueqq_dict['C79u'] * F2un + self.Cmueqq_dict['C79d'] * F2dn + self.Cmueqq_dict['C79s'] * F2sn),

            'd11p': - mL * (self.Cmueqq_dict['C711u'] * FAup + self.Cmueqq_dict['C711d'] * FAdp + self.Cmueqq_dict['C711s'] * FAsp),

            'd11n': - mL * (self.Cmueqq_dict['C711u'] * FAun + self.Cmueqq_dict['C711d'] * FAdn + self.Cmueqq_dict['C711s'] * FAsn),

                    # T-odd contribution (set to zero by default)
            'd12p': - mL * (self.Cmueqq_dict['C711u'] * FA3up + self.Cmueqq_dict['C711d'] * FA3dp + self.Cmueqq_dict['C711s'] * FA3sp),
            
                    # T-odd contribution (set to zero by default)
            'd12n': - mL * (self.Cmueqq_dict['C711u'] * FA3un + self.Cmueqq_dict['C711d'] * FA3dn + self.Cmueqq_dict['C711s'] * FA3sn),

            'd13p': + self.Cmueqq_dict['C62u'] * F1up + self.Cmueqq_dict['C62d'] * F1dp + self.Cmueqq_dict['C62s'] * F1sp\
                    - 1j * m_minus * (self.Cmueqq_dict['C710u'] * F1up + self.Cmueqq_dict['C710d'] * F1dp + self.Cmueqq_dict['C710s'] * F1sp)\
                    - (qSq / (2 * mp)) * (self.Cmueqq_dict['C714u'] * (FT1up - 4 * FT2up) + self.Cmueqq_dict['C714d'] * (FT1dp - 4 * FT2dp) + self.Cmueqq_dict['C714s'] * (FT1sp - 4 * FT2sp)),

            'd13n': + self.Cmueqq_dict['C62u'] * F1un + self.Cmueqq_dict['C62d'] * F1dn + self.Cmueqq_dict['C62s'] * F1sn\
                    - 1j * m_minus * (self.Cmueqq_dict['C710u'] * F1un + self.Cmueqq_dict['C710d'] * F1dn + self.Cmueqq_dict['C710s'] * F1sn)\
                    - (qSq / (2 * mn)) * (self.Cmueqq_dict['C714u'] * (FT1un - 4 * FT2un) + self.Cmueqq_dict['C714d'] * (FT1dn - 4 * FT2dn) + self.Cmueqq_dict['C714s'] * (FT1sn - 4 * FT2sn)),

            'd14p': - (1./2.) * (self.Cmueqq_dict['C62u'] * F2up + self.Cmueqq_dict['C62d'] * F2dp + self.Cmueqq_dict['C62s'] * F2sp)\
                    + (1j / 2.) * m_minus * (self.Cmueqq_dict['C710u'] * F2up + self.Cmueqq_dict['C710d'] * F2dp + self.Cmueqq_dict['C710s'] * F2sp),

            'd14n': - (1./2.) * (self.Cmueqq_dict['C62u'] * F2un + self.Cmueqq_dict['C62d'] * F2dn + self.Cmueqq_dict['C62s'] * F2sn)\
                    + (1j / 2.) * m_minus * (self.Cmueqq_dict['C710u'] * F2un + self.Cmueqq_dict['C710d'] * F2dn + self.Cmueqq_dict['C710s'] * F2sn),

            'd15p': + self.Cmueqq_dict['C64u'] * FAup + self.Cmueqq_dict['C64d'] * FAdp + self.Cmueqq_dict['C64s'] * FAsp\
                    - 1j * m_minus * (self.Cmueqq_dict['C712u'] * FAup + self.Cmueqq_dict['C712d'] * FAdp + self.Cmueqq_dict['C712s'] * FAsp)\
                    # T-odd contribution (set to zero by default)
                    - 2 * (qSq / mp) * (self.Cmueqq_dict['C716u'] * FT3up + self.Cmueqq_dict['C716d'] * FT3dp + self.Cmueqq_dict['C716s'] * FT3sp),

            'd15n': + self.Cmueqq_dict['C64u'] * FAun + self.Cmueqq_dict['C64d'] * FAdn + self.Cmueqq_dict['C64s'] * FAsn\
                    - 1j * m_minus * (self.Cmueqq_dict['C712u'] * FAun + self.Cmueqq_dict['C712d'] * FAdn + self.Cmueqq_dict['C712s'] * FAsn)\
                    # T-odd contribution (set to zero by default)
                    - 2 * (qSq / mn) * (self.Cmueqq_dict['C716u'] * FT3un + self.Cmueqq_dict['C716d'] * FT3dn + self.Cmueqq_dict['C716s'] * FT3sn),

                    # T-odd contribution (set to zero by default)
            'd16p': - 1j * m_minus * (self.Cmueqq_dict['C712u'] * FA3up + self.Cmueqq_dict['C712d'] * FA3dp +self.Cmueqq_dict['C712s'] * FA3sp)\
                    + self.Cmueqq_dict['C64u'] * FA3up + self.Cmueqq_dict['C64d'] * FA3dp + self.Cmueqq_dict['C64s'] * FA3sp,

                    # T-odd contribution (set to zero by default)
            'd16n': - 1j * m_minus * (self.Cmueqq_dict['C712u'] * FA3un + self.Cmueqq_dict['C712d'] * FA3dn +self.Cmueqq_dict['C712s'] * FA3sn)\
                    + self.Cmueqq_dict['C64u'] * FA3un + self.Cmueqq_dict['C64d'] * FA3dn + self.Cmueqq_dict['C64s'] * FA3sn,

            'd17p': + (alpha / np.pi) * self.Cmueqq_dict['C52'] * (mL / qSq) * (qu * F1up + qd * F1dp + qs * F1sp)\
                    + mL * (self.Cmueqq_dict['C710u'] * F1up + self.Cmueqq_dict['C710d'] * F1dp + self.Cmueqq_dict['C710s'] * F1sp),

            'd17n': + (alpha / np.pi) * self.Cmueqq_dict['C52'] * (mL / qSq) * (qu * F1un + qd * F1dn + qs * F1sn)\
                    + mL * (self.Cmueqq_dict['C710u'] * F1un + self.Cmueqq_dict['C710d'] * F1dn + self.Cmueqq_dict['C710s'] * F1sn),

            'd18p': - (alpha / (2 * np.pi)) * self.Cmueqq_dict['C52'] * (mL / qSq) * (qu * F2up + qd * F2dp + qs * F2sp)\
                    - (mL / 2.) * (self.Cmueqq_dict['C710u'] * F2up + self.Cmueqq_dict['C710d'] * F2dp + self.Cmueqq_dict['C710s'] * F2sp),

            'd18n': - (alpha / (2 * np.pi)) * self.Cmueqq_dict['C52'] * (mL / qSq) * (qu * F2un + qd * F2dn + qs * F2sn)\
                    - (mL / 2.) * (self.Cmueqq_dict['C710u'] * F2un + self.Cmueqq_dict['C710d'] * F2dn + self.Cmueqq_dict['C710s'] * F2sn),

            'd19p': + mL * (self.Cmueqq_dict['C712u'] * FAup + self.Cmueqq_dict['C712d'] * FAdp + self.Cmueqq_dict['C712s'] * FAsp),

            'd19n': + mL * (self.Cmueqq_dict['C712u'] * FAun + self.Cmueqq_dict['C712d'] * FAdn + self.Cmueqq_dict['C712s'] * FAsn),

                    # T-odd contribution (set to zero by default)
            'd20p': + mL * (self.Cmueqq_dict['C712u'] * FA3up + self.Cmueqq_dict['C712d'] * FA3dp + self.Cmueqq_dict['C712s'] * FA3sp),

                    # T-odd contribution (set to zero by default)
            'd20n': + mL * (self.Cmueqq_dict['C712u'] * FA3un + self.Cmueqq_dict['C712d'] * FA3dn + self.Cmueqq_dict['C712s'] * FA3sn),

            'd21p': + self.Cmueqq_dict['C69u'] * FT0up + self.Cmueqq_dict['C69d'] * FT0dp + self.Cmueqq_dict['C69s'] * FT0sp,

            'd21n': + self.Cmueqq_dict['C69u'] * FT0un + self.Cmueqq_dict['C69d'] * FT0dn + self.Cmueqq_dict['C69s'] * FT0sn,

            'd22p': - (self.Cmueqq_dict['C69u'] * FT1up + self.Cmueqq_dict['C69d'] * FT1dp + self.Cmueqq_dict['C69s'] * FT1sp),

            'd22n': - (self.Cmueqq_dict['C69u'] * FT1un + self.Cmueqq_dict['C69d'] * FT1dn + self.Cmueqq_dict['C69s'] * FT1sn),

            'd23p': - (self.Cmueqq_dict['C69u'] * FT2up + self.Cmueqq_dict['C69d'] * FT2dp + self.Cmueqq_dict['C69s'] * FT2sp),

            'd23n': - (self.Cmueqq_dict['C69u'] * FT2un + self.Cmueqq_dict['C69d'] * FT2dn + self.Cmueqq_dict['C69s'] * FT2sn),

                    # T-odd contribution (set to zero by default)
            'd24p': - (self.Cmueqq_dict['C69u'] * FT3up + self.Cmueqq_dict['C69d'] * FT3dp + self.Cmueqq_dict['C69s'] * FT3sp),

                    # T-odd contribution (set to zero by default)
            'd24n': - (self.Cmueqq_dict['C69u'] * FT3un + self.Cmueqq_dict['C69d'] * FT3dn + self.Cmueqq_dict['C69s'] * FT3sn),

            'd25p': + self.Cmueqq_dict['C610u'] * FT0up + self.Cmueqq_dict['C610d'] * FT0dp + self.Cmueqq_dict['C610s'] * FT0sp,

            'd25n': + self.Cmueqq_dict['C610u'] * FT0un + self.Cmueqq_dict['C610d'] * FT0dn + self.Cmueqq_dict['C610s'] * FT0sn,

            'd26p': - (self.Cmueqq_dict['C610u'] * FT1up + self.Cmueqq_dict['C610d'] * FT1dp + self.Cmueqq_dict['C610s'] * FT1sp),

            'd26n': - (self.Cmueqq_dict['C610u'] * FT1un + self.Cmueqq_dict['C610d'] * FT1dn + self.Cmueqq_dict['C610s'] * FT1sn),

            'd27p': - (self.Cmueqq_dict['C610u'] * FT2up + self.Cmueqq_dict['C610d'] * FT2dp + self.Cmueqq_dict['C610s'] * FT2sp),

            'd27n': - (self.Cmueqq_dict['C610u'] * FT2un + self.Cmueqq_dict['C610d'] * FT2dn + self.Cmueqq_dict['C610s'] * FT2sn),

                    # T-odd contribution (set to zero by default)
            'd28p': - (self.Cmueqq_dict['C610u'] * FT3up + self.Cmueqq_dict['C610d'] * FT3dp + self.Cmueqq_dict['C610s'] * FT3sp),

                    # T-odd contribution (set to zero by default)
            'd28n': - (self.Cmueqq_dict['C610u'] * FT3un + self.Cmueqq_dict['C610d'] * FT3dn + self.Cmueqq_dict['C610s'] * FT3sn),

            'd29p': - mL * (self.Cmueqq_dict['C713u'] * (FT0up - (qSq / mp**2) * FT2up) + self.Cmueqq_dict['C713d'] * (FT0dp - (qSq / mp**2) * FT2dp) + self.Cmueqq_dict['C713s'] * (FT0sp - (qSq / mp**2) * FT2sp)),

            'd29n': - mL * (self.Cmueqq_dict['C713u'] * (FT0un - (qSq / mn**2) * FT2un) + self.Cmueqq_dict['C713d'] * (FT0dn - (qSq / mn**2) * FT2dn) + self.Cmueqq_dict['C713s'] * (FT0sn - (qSq / mn**2) * FT2sn)),

            'd30p': - mL * (self.Cmueqq_dict['C714u'] * (FT0up - (qSq / mp**2) * FT2up) + self.Cmueqq_dict['C714d'] * (FT0dp - (qSq / mp**2) * FT2dp) + self.Cmueqq_dict['C714s'] * (FT0sp - (qSq / mp**2) * FT2sp)),

            'd30n': - mL * (self.Cmueqq_dict['C714u'] * (FT0un - (qSq / mn**2) * FT2un) + self.Cmueqq_dict['C714d'] * (FT0dn - (qSq / mn**2) * FT2dn) + self.Cmueqq_dict['C714s'] * (FT0sn - (qSq / mn**2) * FT2sn)),

            'd31p': + (1 / 4.) * mL * (self.Cmueqq_dict['C716u'] * FT0up + self.Cmueqq_dict['C716d'] * FT0dp + self.Cmueqq_dict['C716s'] * FT0sp),

            'd31n': + (1 / 4.) * mL * (self.Cmueqq_dict['C716u'] * FT0un + self.Cmueqq_dict['C716d'] * FT0dn + self.Cmueqq_dict['C716s'] * FT0sn), 

            'd32p': + (1 / 4.) * mL * (self.Cmueqq_dict['C715u'] * FT0up + self.Cmueqq_dict['C715d'] * FT0dp + self.Cmueqq_dict['C715s'] * FT0sp),

            'd32n': + (1 / 4.) * mL * (self.Cmueqq_dict['C715u'] * FT0un + self.Cmueqq_dict['C715d'] * FT0dn + self.Cmueqq_dict['C715s'] * FT0sn),
            }

        # Convert to the isospin basis (written in two-component form {isoscalar, isovector})
        RET_isospin_dict = {
            'd1': np.array([RET_nucleon_dict['d1p'] + RET_nucleon_dict['d1n'], RET_nucleon_dict['d1p'] - RET_nucleon_dict['d1n']]) / 2.,
            'd2': np.array([RET_nucleon_dict['d2p'] + RET_nucleon_dict['d2n'], RET_nucleon_dict['d2p'] - RET_nucleon_dict['d2n']]) / 2.,
            'd3': np.array([RET_nucleon_dict['d3p'] + RET_nucleon_dict['d3n'], RET_nucleon_dict['d3p'] - RET_nucleon_dict['d3n']]) / 2.,
            'd4': np.array([RET_nucleon_dict['d4p'] + RET_nucleon_dict['d4n'], RET_nucleon_dict['d4p'] - RET_nucleon_dict['d4n']]) / 2.,
            'd5': np.array([RET_nucleon_dict['d5p'] + RET_nucleon_dict['d5n'], RET_nucleon_dict['d5p'] - RET_nucleon_dict['d5n']]) / 2.,
            'd6': np.array([RET_nucleon_dict['d6p'] + RET_nucleon_dict['d6n'], RET_nucleon_dict['d6p'] - RET_nucleon_dict['d6n']]) / 2.,
            'd7': np.array([RET_nucleon_dict['d7p'] + RET_nucleon_dict['d7n'], RET_nucleon_dict['d7p'] - RET_nucleon_dict['d7n']]) / 2.,
            'd8': np.array([RET_nucleon_dict['d8p'] + RET_nucleon_dict['d8n'], RET_nucleon_dict['d8p'] - RET_nucleon_dict['d8n']]) / 2.,
            'd9': np.array([RET_nucleon_dict['d9p'] + RET_nucleon_dict['d9n'], RET_nucleon_dict['d9p'] - RET_nucleon_dict['d9n']]) / 2.,
            'd10': np.array([RET_nucleon_dict['d10p'] + RET_nucleon_dict['d10n'], RET_nucleon_dict['d10p'] - RET_nucleon_dict['d10n']]) / 2.,
            'd11': np.array([RET_nucleon_dict['d11p'] + RET_nucleon_dict['d11n'], RET_nucleon_dict['d11p'] - RET_nucleon_dict['d11n']]) / 2.,
            'd12': np.array([RET_nucleon_dict['d12p'] + RET_nucleon_dict['d12n'], RET_nucleon_dict['d12p'] - RET_nucleon_dict['d12n']]) / 2.,
            'd13': np.array([RET_nucleon_dict['d13p'] + RET_nucleon_dict['d13n'], RET_nucleon_dict['d13p'] - RET_nucleon_dict['d13n']]) / 2.,
            'd14': np.array([RET_nucleon_dict['d14p'] + RET_nucleon_dict['d14n'], RET_nucleon_dict['d14p'] - RET_nucleon_dict['d14n']]) / 2.,
            'd15': np.array([RET_nucleon_dict['d15p'] + RET_nucleon_dict['d15n'], RET_nucleon_dict['d15p'] - RET_nucleon_dict['d15n']]) / 2.,
            'd16': np.array([RET_nucleon_dict['d16p'] + RET_nucleon_dict['d16n'], RET_nucleon_dict['d16p'] - RET_nucleon_dict['d16n']]) / 2.,
            'd17': np.array([RET_nucleon_dict['d17p'] + RET_nucleon_dict['d17n'], RET_nucleon_dict['d17p'] - RET_nucleon_dict['d17n']]) / 2.,
            'd18': np.array([RET_nucleon_dict['d18p'] + RET_nucleon_dict['d18n'], RET_nucleon_dict['d18p'] - RET_nucleon_dict['d18n']]) / 2.,
            'd19': np.array([RET_nucleon_dict['d19p'] + RET_nucleon_dict['d19n'], RET_nucleon_dict['d19p'] - RET_nucleon_dict['d19n']]) / 2.,
            'd20': np.array([RET_nucleon_dict['d20p'] + RET_nucleon_dict['d20n'], RET_nucleon_dict['d20p'] - RET_nucleon_dict['d20n']]) / 2.,
            'd21': np.array([RET_nucleon_dict['d21p'] + RET_nucleon_dict['d21n'], RET_nucleon_dict['d21p'] - RET_nucleon_dict['d21n']]) / 2.,
            'd22': np.array([RET_nucleon_dict['d22p'] + RET_nucleon_dict['d22n'], RET_nucleon_dict['d22p'] - RET_nucleon_dict['d22n']]) / 2.,
            'd23': np.array([RET_nucleon_dict['d23p'] + RET_nucleon_dict['d23n'], RET_nucleon_dict['d23p'] - RET_nucleon_dict['d23n']]) / 2.,
            'd24': np.array([RET_nucleon_dict['d24p'] + RET_nucleon_dict['d24n'], RET_nucleon_dict['d24p'] - RET_nucleon_dict['d24n']]) / 2.,
            'd25': np.array([RET_nucleon_dict['d25p'] + RET_nucleon_dict['d25n'], RET_nucleon_dict['d25p'] - RET_nucleon_dict['d25n']]) / 2.,
            'd26': np.array([RET_nucleon_dict['d26p'] + RET_nucleon_dict['d26n'], RET_nucleon_dict['d26p'] - RET_nucleon_dict['d26n']]) / 2.,
            'd27': np.array([RET_nucleon_dict['d27p'] + RET_nucleon_dict['d27n'], RET_nucleon_dict['d27p'] - RET_nucleon_dict['d27n']]) / 2.,
            'd28': np.array([RET_nucleon_dict['d28p'] + RET_nucleon_dict['d28n'], RET_nucleon_dict['d28p'] - RET_nucleon_dict['d28n']]) / 2.,
            'd29': np.array([RET_nucleon_dict['d29p'] + RET_nucleon_dict['d29n'], RET_nucleon_dict['d29p'] - RET_nucleon_dict['d29n']]) / 2.,
            'd30': np.array([RET_nucleon_dict['d30p'] + RET_nucleon_dict['d30n'], RET_nucleon_dict['d30p'] - RET_nucleon_dict['d30n']]) / 2.,
            'd31': np.array([RET_nucleon_dict['d31p'] + RET_nucleon_dict['d31n'], RET_nucleon_dict['d31p'] - RET_nucleon_dict['d31n']]) / 2.,
            'd32': np.array([RET_nucleon_dict['d32p'] + RET_nucleon_dict['d32n'], RET_nucleon_dict['d32p'] - RET_nucleon_dict['d32n']]) / 2.
            }

        # Return in list format compatible with Mu2e_NRET
        RET_isospin_array = [
            [int(1), RET_isospin_dict['d1'][0], RET_isospin_dict['d1'][1]],
            [int(2), RET_isospin_dict['d2'][0], RET_isospin_dict['d2'][1]],
            [int(3), RET_isospin_dict['d3'][0], RET_isospin_dict['d3'][1]],
            [int(4), RET_isospin_dict['d4'][0], RET_isospin_dict['d4'][1]],
            [int(5), RET_isospin_dict['d5'][0], RET_isospin_dict['d5'][1]],
            [int(6), RET_isospin_dict['d6'][0], RET_isospin_dict['d6'][1]],
            [int(7), RET_isospin_dict['d7'][0], RET_isospin_dict['d7'][1]],
            [int(8), RET_isospin_dict['d8'][0], RET_isospin_dict['d8'][1]],
            [int(9), RET_isospin_dict['d9'][0], RET_isospin_dict['d9'][1]],
            [int(10), RET_isospin_dict['d10'][0], RET_isospin_dict['d10'][1]],
            [int(11), RET_isospin_dict['d11'][0], RET_isospin_dict['d11'][1]],
            [int(12), RET_isospin_dict['d12'][0], RET_isospin_dict['d12'][1]],
            [int(13), RET_isospin_dict['d13'][0], RET_isospin_dict['d13'][1]],
            [int(14), RET_isospin_dict['d14'][0], RET_isospin_dict['d14'][1]],
            [int(15), RET_isospin_dict['d15'][0], RET_isospin_dict['d15'][1]],
            [int(16), RET_isospin_dict['d16'][0], RET_isospin_dict['d16'][1]],
            [int(17), RET_isospin_dict['d17'][0], RET_isospin_dict['d17'][1]],
            [int(18), RET_isospin_dict['d18'][0], RET_isospin_dict['d18'][1]],
            [int(19), RET_isospin_dict['d19'][0], RET_isospin_dict['d19'][1]],
            [int(20), RET_isospin_dict['d20'][0], RET_isospin_dict['d20'][1]],
            [int(21), RET_isospin_dict['d21'][0], RET_isospin_dict['d21'][1]],
            [int(22), RET_isospin_dict['d22'][0], RET_isospin_dict['d22'][1]],
            [int(23), RET_isospin_dict['d23'][0], RET_isospin_dict['d23'][1]],
            [int(24), RET_isospin_dict['d24'][0], RET_isospin_dict['d24'][1]],
            [int(25), RET_isospin_dict['d25'][0], RET_isospin_dict['d25'][1]],
            [int(26), RET_isospin_dict['d26'][0], RET_isospin_dict['d26'][1]],
            [int(27), RET_isospin_dict['d27'][0], RET_isospin_dict['d27'][1]],
            [int(28), RET_isospin_dict['d28'][0], RET_isospin_dict['d28'][1]],
            [int(29), RET_isospin_dict['d29'][0], RET_isospin_dict['d29'][1]],
            [int(30), RET_isospin_dict['d30'][0], RET_isospin_dict['d30'][1]],
            [int(31), RET_isospin_dict['d31'][0], RET_isospin_dict['d31'][1]],
            [int(32), RET_isospin_dict['d32'][0], RET_isospin_dict['d32'][1]],
            ]

        return RET_isospin_array

    def Cmueqq_to_RET_map(self):
        """ 
        Generate a mapping of the hadronization step between WET and RET consistent with NetworkX
        """
        # Initialize tuple array
        had_map = []
        # Initialize iterator
        iterator = 1
        for key in self.wc_keys:
            # Generate a zero dictionary
            for wc_name in self.wc_keys:
                self.Cmueqq_dict[wc_name] = 0.
            # Set the key to 1 (only interested in generated non-zero entries in d)
            self.Cmueqq_dict[key] = 1.0
            # Compute d coefficients with dummy q-vector
            d_isospin = self.hadronize_WET(np.array([1,1,1,1]))
            # Record non-zero values
            for i in range(len(d_isospin)):
                # Append non-zero values to map
                if d_isospin[i][1] != 0:
                    had_map.append((iterator,i))
            iterator += 1
        return had_map

if __name__ == "__main__":
    foobar_dict = {}
    foobar_input_dict = None
    had = Hadronization(foobar_dict, foobar_input_dict)
    had_map = had.hadronization_map()
    print(had_map)
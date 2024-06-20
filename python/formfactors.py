####################################################################################
# ------------------------------------ Vector ------------------------------------ #
####################################################################################

class F1:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor F1
        Return the nuclear form factor F1
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F1up_0'] == None:
                    return 2
                else: return self.ip['F1up_0']
            if self.quark == 'd':
                if self.ip['F1dp_0'] == None:
                    return 1
                else: return self.ip['F1dp_0']
            if self.quark == 's':
                if self.ip['F1sp_0'] == None:
                    return 0
                else: return self.ip['F1sp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F1un_0'] == None:
                    return 1
                else: return self.ip['F1un_0']
            if self.quark == 'd':
                if self.ip['F1dn_0'] == None:
                    return 2
                else: return self.ip['F1dn_0']
            if self.quark == 's':
                if self.ip['F1sn_0'] == None:
                    return 0
                else: return self.ip['F1sn_0']

    def first_deriv_zero_mom(self):
        """ 
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F1up_prime'] == None:
                    return (1 / 6) * (2 * self.ip['rEp2'] + self.ip['rEn2'] + self.ip['rEs2']) - (1 / (4 * mN**2)) * (2 * self.ip['mup'] + self.ip['mun'] + self.ip['mus'] - 2)
                else: return self.ip['F1up_prime']
            if self.quark == 'd':
                if self.ip['F1dp_prime'] == None:
                    return (1 / 6) * (self.ip['rEp2'] + 2 * self.ip['rEn2'] + self.ip['rEs2']) - (1 / (4 * mN**2)) * (self.ip['mup'] + 2 * self.ip['mun'] + self.ip['mus'] - 1)
                else: return self.ip['F1dp_prime']
            if self.quark == 's':
                if self.ip['F1sp_prime'] == None:
                    return (1 / 6) * self.ip['rEs2'] - (self.ip['mus'] / (4 * mN**2))
                else: return self.ip['F1sp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F1un_prime'] == None:
                    return (1 / 6) * (self.ip['rEp2'] + 2 * self.ip['rEn2'] + self.ip['rEs2']) - (1 / (4 * mN**2)) * (self.ip['mup'] + 2 * self.ip['mun'] + self.ip['mus'] - 1)
                else: return self.ip['F1un_prime']
            if self.quark == 'd':
                if self.ip['F1dn_prime'] == None:
                    return (1 / 6) * (2 * self.ip['rEp2'] + self.ip['rEn2'] + self.ip['rEs2']) - (1 / (4 * mN**2)) * (2 * self.ip['mup'] + self.ip['mun'] + self.ip['mus'] - 2)
                else: return self.ip['F1dn_prime']
            if self.quark == 's':
                if self.ip['F1sn_prime'] == None:
                    return (1 / 6) * self.ip['rEs2'] - (self.ip['mus'] / (4 * mN**2))
                else: return self.ip['F1sn_prime']

    def form_factor(self, q):
        """
        Return the full q-dependent form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F1up'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1up']
            if self.quark == 'd':
                if self.ip['F1dp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1dp']
            if self.quark == 's':
                if self.ip['F1sp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F1un'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1un']
            if self.quark == 'd':
                if self.ip['F1dn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1dn']
            if self.quark == 's':
                if self.ip['F1sn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F1sn']
        
class F2:
    def __init__(self, quark, nucleon, input_dict):
        """
        The nuclear form factor F2
        Return the nuclear form factor F2
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """
        Return the value of the form factor at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F2up_0'] == None:
                    return 2 * (self.ip['mup'] - 1) + self.ip['mun'] + self.ip['mus']
                else: return self.ip['F2up_0']
            if self.quark == 'd':
                if self.ip['F2dp_0'] == None:
                    return 2 * self.ip['mun'] + (self.ip['mup'] - 1) + self.ip['mus']
                else: return self.ip['F2dp_0']
            if self.quark == 's':
                if self.ip['F2sp_0'] == None:
                    return self.ip['mus']
                else: return self.ip['F2sp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F2un_0'] == None:
                    return 2 * self.ip['mun'] + (self.ip['mup'] - 1) + self.ip['mus']
                else: return self.ip['F2un_0']
            if self.quark == 'd':
                if self.ip['F2dn_0'] == None:
                    return 2 * (self.ip['mup'] - 1) + self.ip['mun'] + self.ip['mus']
                else: return self.ip['F2dn_0']
            if self.quark == 's':
                if self.ip['F2sn_0'] == None:
                    return self.ip['mus']
                else: return self.ip['F2sn_0']

    def first_deriv_zero_mom(self):
        """ 
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F2up_prime'] == None:
                    return (1 / 6) * (2 * (self.ip['rMp2'] - self.ip['rEp2']) + self.ip['rMn2'] - self.ip['rEn2'] + self.ip['rMs2'] - self.ip['rEs2']) + (1 / (4 * mN**2)) * (2 * self.ip['mup'] + self.ip['mun'] + self.ip['mus'] - 2)
                else: return self.ip['F2up_prime']
            if self.quark == 'd':
                if self.ip['F2dp_prime'] == None:
                    return (1 / 6) * (self.ip['rMp2'] - self.ip['rEp2'] + 2 * (self.ip['rMn2'] - self.ip['rEn2']) + self.ip['rMs2'] - self.ip['rEs2']) + (1 / (4 * mN**2)) * (self.ip['mup'] + 2 * self.ip['mun'] + self.ip['mus'] - 1)
                else: return self.ip['F2dp_prime']
            if self.quark == 's':
                if self.ip['F2sp_prime'] == None:
                    return (1 / 6) * (self.ip['rMs2'] - self.ip['rEs2']) + (self.ip['mus'] / (4 * mN**2))
                else: return self.ip['F2sp_prime']
    
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F2un_prime'] == None:
                    return (1 / 6) * (self.ip['rMp2'] - self.ip['rEp2'] + 2 * (self.ip['rMn2'] - self.ip['rEn2']) + self.ip['rMs2'] - self.ip['rEs2']) + (1 / (4 * mN**2)) * (self.ip['mup'] + 2 * self.ip['mun'] + self.ip['mus'] - 1)
                else: return self.ip['F2un_prime']
            if self.quark == 'd':
                if self.ip['F2dn_prime'] == None:
                    return (1 / 6) * (2 * (self.ip['rMp2'] - self.ip['rEp2']) + self.ip['rMn2'] - self.ip['rEn2'] + self.ip['rMs2'] - self.ip['rEs2']) + (1 / (4 * mN**2)) * (2 * self.ip['mup'] + self.ip['mun'] + self.ip['mus'] - 2)
                else: return self.ip['F2dn_prime']
            if self.quark == 's':
                if self.ip['F2sn_prime'] == None:
                    return (1 / 6) * (self.ip['rMs2'] - self.ip['rEs2']) + (self.ip['mus'] / (4 * mN**2))
                else: return self.ip['F2sn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F2up'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2up']
            if self.quark == 'd':
                if self.ip['F2dp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2dp']
            if self.quark == 's':
                if self.ip['F2sp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F2un'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2un']
            if self.quark == 'd':
                if self.ip['F2dn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2dn']
            if self.quark == 's':
                if self.ip['F2sn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['F2sn']

class F3:
    def __init__(self, quark, nucleon, input_dict):
        """
        The nuclear form factor F3
        Return the nuclear form factor F3
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)

        * This form factor is set to zero by default (T-odd) *
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def form_factor(self, q):
        """
        Return the full form factor
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['F3up'] == None:
                    return 0.
                else: return self.ip['F3up']
            if self.quark == 'd':
                if self.ip['F3dp'] == None:
                    return 0.
                else: return self.ip['F3dp']
            if self.quark == 's':
                if self.ip['F3sp'] == None:
                    return 0.
                else: return self.ip['F3sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['F3un'] == None:
                    return 0.
                else: return self.ip['F3un']
            if self.quark == 'd':
                if self.ip['F3dn'] == None:
                    return 0.
                else: return self.ip['F3dn']
            if self.quark == 's':
                if self.ip['F3sn'] == None:
                    return 0.
                else: return self.ip['F3sn']

####################################################################################
# ---------------------------------- Axial-vector -------------------------------- #
####################################################################################

class FA:
    def __init__(self, quark, nucleon, input_dict):
        """
        The nuclear form factor FA at zero momentum transfer
        Return the nuclear form factor FA, evaluated at zero momentum transfer.
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FAup_0'] == None:
                    return (1 / 2) * (self.ip['gA'] + self.ip['DeltaSigmaud'])
                else: return self.ip['FAup_0']
            if self.quark == 'd':
                if self.ip['FAdp_0'] == None:
                    return (1 / 2) * (- self.ip['gA'] + self.ip['DeltaSigmaud'])
                else: return self.ip['FAdp_0']
            if self.quark == 's':
                if self.ip['FAsp_0'] == None:
                    return self.ip['Deltas']
                else: return self.ip['FAsp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FAun_0'] == None:
                    return (1 / 2) * (- self.ip['gA'] + self.ip['DeltaSigmaud'])
                else: return self.ip['FAun_0']
            if self.quark == 'd':
                if self.ip['FAdn_0'] == None:
                    return (1 / 2) * (self.ip['gA'] + self.ip['DeltaSigmaud'])
                else: return self.ip['FAdn_0']
            if self.quark == 's':
                if self.ip['FAsn_0'] == None:
                    return self.ip['Deltas']
                else: return self.ip['FAsn_0']

    def first_deriv_zero_mom(self):
        """ 
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FAup_prime'] == None:
                    return (1 / 12) * (self.ip['gA'] * self.ip['rA2uMinusd'] + self.ip['DeltaSigmaud'] * self.ip['rA2uPlusd'])
                else: return self.ip['FAup_prime']
            if self.quark == 'd':
                if self.ip['FAdp_prime'] == None:
                    return (1 / 12) * (- self.ip['gA'] * self.ip['rA2uMinusd'] + self.ip['DeltaSigmaud'] * self.ip['rA2uPlusd'])
                else: return self.ip['FAdp_prime']
            if self.quark == 's':
                if self.ip['FAsp_prime'] == None:
                    return (1 / 6) * self.ip['Deltas'] * self.ip['rA2s']
                else: return self.ip['FAsp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FAun_prime'] == None:
                    return (1 / 12) * (- self.ip['gA'] * self.ip['rA2uMinusd'] + self.ip['DeltaSigmaud'] * self.ip['rA2uPlusd'])
                else: return self.ip['FAun_prime']
            if self.quark == 'd':
                if self.ip['FAdn_prime'] == None:
                    return (1 / 12) * (self.ip['gA'] * self.ip['rA2uMinusd'] + self.ip['DeltaSigmaud'] * self.ip['rA2uPlusd'])
                else: return self.ip['FAdn_prime']
            if self.quark == 's':
                if self.ip['FAsn_prime'] == None:
                    return (1 / 6) * self.ip['Deltas'] * self.ip['rA2s']
                else: return self.ip['FAsn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FAup'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAup']
            if self.quark == 'd':
                if self.ip['FAdp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAdp']
            if self.quark == 's':
                if self.ip['FAsp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAsp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FAun'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAun']
            if self.quark == 'd':
                if self.ip['FAdn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAdn']
            if self.quark == 's':
                if self.ip['FAsn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FAsn']

class FPprimed:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor FPprimed
        Return the nuclear form factor FPprimed
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FPprimedup_0']
            if self.quark == 'd':
                return self.ip['FPprimeddp_0']
            if self.quark == 's':
                return self.ip['FPprimedsp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FPprimedun_0']
            if self.quark == 'd':
                return self.ip['FPprimeddn_0']
            if self.quark == 's':
                return self.ip['FPprimedsn_0']

    def value_pion_pole(self):
        """ Return the coefficient of the pion pole
        The pion pole is given, in terms of the spatial momentum q, by 1 / (q^2 + mpi0^2)
        """
        # Define the residua of the pion pole contributions
        self.aPprimepiu = 2 * self.ip['gA']
        self.aPprimepid = - 2 * self.ip['gA']

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPprimedup_pion'] == None:
                    return self.aPprimepiu
                else: return self.ip['FPprimedup_pion']
            if self.quark == 'd':
                if self.ip['FPprimeddp_pion'] == None:
                    return self.aPprimepid
                else: return self.ip['FPprimeddp_pion']
            if self.quark == 's':
                if self.ip['FPprimedsp_pion'] == None:
                    return 0
                else: return self.ip['FPprimedsp_pion']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPprimedun_pion'] == None:
                    return self.aPprimepid
                else: return self.ip['FPprimedun_pion']
            if self.quark == 'd':
                if self.ip['FPprimeddn_pion'] == None:
                    return self.aPprimepiu
                else: return self.ip['FPprimeddn_pion']
            if self.quark == 's':
                if self.ip['FPprimedsn_pion'] == None:
                    return 0
                else: return self.ip['FPprimedsn_pion']

    def value_eta_pole(self):
        """ 
        Return the coefficient of the pion pole
        The eta pole is given, in terms of the spatial momentum q, by 1 / (q^2 + meta^2)
        """
        # Define the residua of the eta pole contributions
        self.aPprimeetau = (2 / 3) * (self.ip['DeltaSigmaud'] + 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])
        self.aPprimeetad = (2 / 3) * (self.ip['DeltaSigmaud'] + 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])
        self.aPprimeetas = - (1 / 3) * (self.ip['DeltaSigmaud'] + 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPprimedup_eta'] == None:
                    return self.aPprimeetau
                else: return self.ip['FPprimedup_eta']
            if self.quark == 'd':
                if self.ip['FPprimeddp_eta'] == None:
                    return self.aPprimeetad
                else: return self.ip['FPprimeddp_eta']
            if self.quark == 's':
                if self.ip['FPprimedsp_eta'] == None:
                    return self.aPprimeetas
                else: return self.ip['FPprimedsp_eta']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPprimedun_eta'] == None:
                    return self.aPprimeetad
                else: return self.ip['FPprimedun_eta']
            if self.quark == 'd':
                if self.ip['FPprimeddn_eta'] == None:
                    return self.aPprimeetau
                else: return self.ip['FPprimeddn_eta']
            if self.quark == 's':
                if self.ip['FPprimedsn_eta'] == None:
                    return self.aPprimeetas
                else: return self.ip['FPprimedsn_eta']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        pion_pole = (mN**2 / (self.ip['mpi0']**2 - qSq)) * self.value_pion_pole()
        eta_pole = (mN**2 / (self.ip['meta']**2 - qSq)) * self.value_eta_pole()
        b = self.value_zero_mom() - mN**2 * ((self.value_pion_pole() / self.ip['mpi0']**2) + (self.value_eta_pole() / self.ip['meta']**2))

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPprimedup'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimedup']
            if self.quark == 'd':
                if self.ip['FPprimeddp'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimeddp']
            if self.quark == 's':
                if self.ip['FPprimedsp'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimedsp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPprimedun'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimedun']
            if self.quark == 'd':
                if self.ip['FPprimeddn'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimeddn']
            if self.quark == 's':
                if self.ip['FPprimedsn'] == None:
                    return pion_pole + eta_pole + b
                else: return self.ip['FPprimedsn']

class FA3:
    def __init__(self, quark, nucleon, input_dict):
        """
        The nuclear form factor FA3 at zero momentum transfer
        Return the nuclear form factor FA3, evaluated at zero momentum transfer. 
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)

        * This form factor is set to zero by default (T-odd) *
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def form_factor(self, q):
        """
        Return the full form factor
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FA3up'] == None:
                    return 0.
                else: return self.ip['FA3up']
            if self.quark == 'd':
                if self.ip['FA3dp'] == None:
                    return 0.
                else: return self.ip['FA3dp']
            if self.quark == 's':
                if self.ip['FA3sp'] == None:
                    return 0.
                else: return self.ip['FA3sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FA3un'] == None:
                    return 0.
                else: return self.ip['FA3un']
            if self.quark == 'd':
                if self.ip['FA3dn'] == None:
                    return 0.
                else: return self.ip['FA3dn']
            if self.quark == 's':
                if self.ip['FA3sn'] == None:
                    return 0.
                else: return self.ip['FA3sn']

####################################################################################
# ------------------------------------ Scalar ------------------------------------ #
####################################################################################

class FS:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor FS
        Return the nuclear form factor FS
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        xi = (1 - self.ip['rud']) / (1 + self.ip['rud'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FSup_0'] == None:
                    sigmaup = (1 / 2) * self.ip['sigmapiNtilde'] * (1 - xi) + self.ip['c5hat'] * (1 - (1 / xi))
                    return sigmaup
                else: return self.ip['FSup_0']
            if self.quark == 'd':
                if self.ip['FSdp_0'] == None: 
                    sigmadp = (1 / 2) * self.ip['sigmapiNtilde'] * (1 + xi) + self.ip['c5hat'] * (1 + (1 / xi))
                    return sigmadp
                else: return self.ip['FSdp_0']
            if self.quark == 's': 
                if self.ip['FSsp_0'] == None:
                    return self.ip['sigmas']
                else: return self.ip['FSsp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FSun_0'] == None:
                    sigmaun = (1 / 2) * self.ip['sigmapiNtilde'] * (1 - xi) - self.ip['c5hat'] * (1 - (1 / xi))
                    return sigmaun
                else: return self.ip['FSun_0']
            if self.quark == 'd':
                if self.ip['FSdn_0'] == None:
                    sigmadn = (1 / 2) * self.ip['sigmapiNtilde'] * (1 + xi) - self.ip['c5hat'] * (1 + (1 / xi))
                    return sigmadn
                else: return self.ip['FSdn_0']
            if self.quark == 's':
                if self.ip['FSsn_0'] == None:
                    return self.ip['sigmas']
                else: return self.ip['FSsn_0']

    def first_deriv_zero_mom(self):
        """ 
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FSup_prime']
            if self.quark == 'd':
                return self.ip['FSdp_prime']
            if self.quark == 's':
                return self.ip['FSsp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FSun_prime']
            if self.quark == 'd':
                return self.ip['FSdn_prime']
            if self.quark == 's':
                return self.ip['FSsn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FSup'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSup']
            if self.quark == 'd':
                if self.ip['FSdp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSdp']
            if self.quark == 's':
                if self.ip['FSsp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSsp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FSun'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSun']
            if self.quark == 'd':
                if self.ip['FSdn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSdn']
            if self.quark == 's':
                if self.ip['FSsn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FSsn']

####################################################################################
# --------------------------------- Psuedoscalar --------------------------------- #
####################################################################################

class FP:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor FP
        Return the nuclear form factor FP
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)

        nucleon = 'p', 'n' -- the nucleon (proton or neutron)

        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_pion_pole(self):
        """ 
        Return the coefficient of the pion pole
        The pion pole is given, in terms of the spatial momentum q, by 1 / (q^2 + mpi0^2)
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        self.aPpiu = (self.ip['mpi0']**2 / mN) * (1 / (1 + (1 / self.ip['rud']))) * self.ip['gA']
        self.aPpid = (self.ip['mpi0']**2 / mN) * (1 / (1 + self.ip['rud'])) * self.ip['gA']

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPup_pion'] == None:
                    return self.aPpiu
                else: return self.ip['FPup_pion']
            if self.quark == 'd':
                if self.ip['FPdp_pion'] == None:
                    return self.aPpid
                else: return self.ip['FPdp_pion']
            if self.quark == 's':
                if self.ip['FPsp_pion'] == None:
                    return 0
                else: return self.ip['FPsp_pion']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPun_pion'] == None:
                    return self.aPpid
                else: return self.ip['FPun_pion']
            if self.quark == 'd':
                if self.ip['FPdn_pion'] == None:
                    return self.aPpiu
                else: return self.ip['FPdn_pion']
            if self.quark == 's':
                if self.ip['FPsn_pion'] == None:
                    return 0
                else: return self.ip['FPsn_pion']

    def value_eta_pole(self):
        """ 
        Return the coefficient of the pion pole
        The eta pole is given, in terms of the spatial momentum q, by 1 / (q^2 + meta^2)
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        self.aPetau = (self.ip['meta']**2 / mN) * (1 / (1 + (1 / self.ip['rud']))) * (1 / (1 + 2 * self.ip['rs'])) * (self.ip['DeltaSigmaud'] - 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])
        self.aPetad = (self.ip['meta']**2 / mN) * (1 / (1 + self.ip['rud'])) * (1 / (1 + 2 * self.ip['rs'])) * (self.ip['DeltaSigmaud'] - 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])
        self.aPetas = - (self.ip['meta']**2 / mN) * (1 / (2 + (1 / self.ip['rs']))) * (self.ip['DeltaSigmaud'] - 2 * self.ip['Deltas']) * (1 + self.ip['DeltaGT8'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPup_eta'] == None:
                    return self.aPetau
                else: return self.ip['FPup_eta']
            if self.quark == 'd':
                if self.ip['FPdp_eta'] == None:
                    return self.aPetad
                else: return self.ip['FPdp_eta']
            if self.quark == 's':
                if self.ip['FPsp_eta'] == None:
                    return self.aPetas
                else: return self.ip['FPsp_eta']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPun_eta'] == None:
                    return self.aPetad
                else: return self.ip['FPun_eta']
            if self.quark == 'd':
                if self.ip['FPdn_eta'] == None:
                    return self.aPetau
                else: return self.ip['FPdn_eta']
            if self.quark == 's':
                if self.ip['FPsn_eta'] == None:
                    return self.aPetas
                else: return self.ip['FPsn_eta']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        mtilde = ((1 / self.ip['mu_at_2GeV']) + (1 / self.ip['md_at_2GeV']) + (1 / self.ip['ms_at_2GeV']))**(-1)
        Deltau = (1 / 2) * (self.ip['gA'] + self.ip['DeltaSigmaud'])
        Deltad = (1 / 2) * (- self.ip['gA'] + self.ip['DeltaSigmaud'])
        pion_pole = (mN**2 / (self.ip['mpi0']**2 - qSq)) * self.value_pion_pole()
        eta_pole = (mN**2 / (self.ip['meta']**2 - qSq)) * self.value_eta_pole()

        bPud = (mN / 3) * (self.ip['DeltaSigmaud'] - 2 * self.ip['Deltas']) * ((6 * self.ip['rs'] - 3 * self.ip['DeltaGT8']) / (4 * self.ip['rs'] + 2))\
               + mN * (- (1 / 3) * self.ip['gA'] + self.ip['Deltas'])\
               - mN * mtilde * ((Deltau / self.ip['mu_at_2GeV']) + (Deltad / self.ip['md_at_2GeV']) + (self.ip['Deltas'] / self.ip['ms_at_2GeV']))

        bPs = bPud + (1 / 2) * mN * self.ip['DeltaGT8'] * (Deltau + Deltad - 2 * self.ip['Deltas'])

        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FPup'] == None:
                    return pion_pole + eta_pole + bPud
                else: return self.ip['FPup']
            if self.quark == 'd':
                if self.ip['FPdp'] == None:
                    return pion_pole + eta_pole + bPud
                else: return self.ip['FPdp']
            if self.quark == 's':
                if self.ip['FPsp'] == None:
                    return pion_pole + eta_pole + bPs
                else: return self.ip['FSsp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FPun'] == None:
                    return pion_pole + eta_pole + bPud
                else: return self.ip['FPun']
            if self.quark == 'd':
                if self.ip['FPdn'] == None:
                    return pion_pole + eta_pole + bPud
                else: return self.ip['FPdn']
            if self.quark == 's':
                if self.ip['FPsn'] == None:
                    return pion_pole + eta_pole + bPs
                else: return self.ip['FSsn']

####################################################################################
# ------------------------------------ Tensor ------------------------------------ #
####################################################################################

class FT0:
    def __init__(self, quark, nucleon, input_dict):
        """
        The nuclear form factor FT0
        Return the nuclear form factor FT0
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT0up_0']
            if self.quark == 'd':
                return self.ip['FT0dp_0']
            if self.quark == 's':
                return self.ip['FT0sp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT0un_0']
            if self.quark == 'd':
                return self.ip['FT0dn_0']
            if self.quark == 's':
                return self.ip['FT0sn_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT0up_prime']
            if self.quark == 'd':
                return self.ip['FT0dp_prime']
            if self.quark == 's':
                return self.ip['FT0sp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT0un_prime']
            if self.quark == 'd':
                return self.ip['FT0dn_prime']
            if self.quark == 's':
                return self.ip['FT0sn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FT0up'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0up']
            if self.quark == 'd':
                if self.ip['FT0dp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0dp']
            if self.quark == 's':
                if self.ip['FT0sp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FT0un'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0un']
            if self.quark == 'd':
                if self.ip['FT0dn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0dn']
            if self.quark == 's':
                if self.ip['FT0sn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT0sn']

class FT1(object):
    def __init__(self, quark, nucleon, input_dict):
        """ The nuclear form factor FT1

        Return the nuclear form factor FT1

        Arguments
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)

        nucleon = 'p', 'n' -- the nucleon (proton or neutron)

        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT1up_0']
            if self.quark == 'd':
                return self.ip['FT1dp_0']
            if self.quark == 's':
                return self.ip['FT1sp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT1un_0']
            if self.quark == 'd':
                return self.ip['FT1dn_0']
            if self.quark == 's':
                return self.ip['FT1sn_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT1up_prime']
            if self.quark == 'd':
                return self.ip['FT1dp_prime']
            if self.quark == 's':
                return self.ip['FT1sp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT1un_prime']
            if self.quark == 'd':
                return self.ip['FT1dn_prime']
            if self.quark == 's':
                return self.ip['FT1sn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FT1up'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1up']
            if self.quark == 'd':
                if self.ip['FT1dp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1dp']
            if self.quark == 's':
                if self.ip['FT1sp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FT1un'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1un']
            if self.quark == 'd':
                if self.ip['FT1dn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1dn']
            if self.quark == 's':
                if self.ip['FT1sn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT1sn']

class FT2:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor FT2
        Return the nuclear form factor FT2
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT2up_0']
            if self.quark == 'd':
                return self.ip['FT2dp_0']
            if self.quark == 's':
                return self.ip['FT2sp_0']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT2un_0']
            if self.quark == 'd':
                return self.ip['FT2dn_0']
            if self.quark == 's':
                return self.ip['FT2sn_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                return self.ip['FT2up_prime']
            if self.quark == 'd':
                return self.ip['FT2dp_prime']
            if self.quark == 's':
                return self.ip['FT2sp_prime']
        if self.nucleon == 'n':
            if self.quark == 'u':
                return self.ip['FT2un_prime']
            if self.quark == 'd':
                return self.ip['FT2dn_prime']
            if self.quark == 's':
                return self.ip['FT2sn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FT2up'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2up']
            if self.quark == 'd':
                if self.ip['FT2dp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2dp']
            if self.quark == 's':
                if self.ip['FT2sp'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FT2un'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2un']
            if self.quark == 'd':
                if self.ip['FT2dn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2dn']
            if self.quark == 's':
                if self.ip['FT2sn'] == None:
                    return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
                else: return self.ip['FT2sn']

class FT3:
    def __init__(self, quark, nucleon, input_dict):
        """ 
        The nuclear form factor FT3
        Return the nuclear form factor FT3
        Arguments:
        ---------
        quark = 'u', 'd', 's' -- the quark flavor (up, down, strange)
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        * This form factor is set to zero by default (T-odd) *
        """
        self.quark = quark
        self.nucleon = nucleon
        self.ip = input_dict

    def form_factor(self, q):
        """
        Return the full form factor
        """
        if self.nucleon == 'p':
            if self.quark == 'u':
                if self.ip['FT3up'] == None:
                    return 0.
                else: return self.ip['FT3up']
            if self.quark == 'd':
                if self.ip['FT3dp'] == None:
                    return 0.
                else: return self.ip['FT3dp']
            if self.quark == 's':
                if self.ip['FT3sp'] == None:
                    return 0.
                else: return self.ip['FT3sp']
        if self.nucleon == 'n':
            if self.quark == 'u':
                if self.ip['FT3un'] == None:
                    return 0.
                else: return self.ip['FT3un']
            if self.quark == 'd':
                if self.ip['FT3dn'] == None:
                    return 0.
                else: return self.ip['FT3dn']
            if self.quark == 's':
                if self.ip['FT3sn'] == None:
                    return 0.
                else: return self.ip['FT3sn']


####################################################################################
# ----------------------------------- Gluonic ------------------------------------ #
####################################################################################

class FG:
    def __init__(self, nucleon, input_dict):
        """ 
        The nuclear form factor FG
        Return the nuclear form factor FG
        Arguments:
        ---------
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        if self.nucleon == 'p':
            return self.ip['FGp_0']
        if self.nucleon == 'n':
            return self.ip['FGn_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            return self.ip['FGp_prime']
        if self.nucleon == 'n':
            return self.ip['FGn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.ip['FGp'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FGp']
        if self.nucleon == 'n':
            if self.ip['FGn'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FGn']

class FGtilde:
    def __init__(self, nucleon, input_dict):
        """ 
        The nuclear form factor FGtilde
        Return the nuclear form factor FGtilde
        Arguments:
        ---------
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        mtilde = ((1 / self.ip['mu_at_2GeV']) + (1 / self.ip['md_at_2GeV']) + (1 / self.ip['ms_at_2GeV']))**(-1)
        Deltau = (1 / 2) * (self.ip['gA'] + self.ip['DeltaSigmaud'])
        Deltad = (1 / 2) * (- self.ip['gA'] + self.ip['DeltaSigmaud'])

        if self.nucleon == 'p':
            if self.ip['FGtildep_0'] == None:
                return - mN * mtilde * (Deltau/self.ip['mu_at_2GeV']\
                                        + Deltad/self.ip['md_at_2GeV']\
                                        + self.ip['Deltas']/self.ip['ms_at_2GeV'])
            else: return self.ip['FGtildep_0']
        if self.nucleon == 'n':
            if self.ip['FGtilden_0'] == None:
                return - mN * mtilde * (Deltau/self.ip['mu_at_2GeV']\
                                        + Deltad/self.ip['md_at_2GeV']\
                                        + self.ip['Deltas']/self.ip['ms_at_2GeV'])
            else: return self.ip['FGtilden_0']

    def value_pion_pole(self):
        """ 
        Return the coefficient of the pion pole
        The pion pole is given, in terms of the spatial momentum q, by q^2 / (q^2 + mpi0^2)
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        mtilde = ((1 / self.ip['mu_at_2GeV']) + (1 / self.ip['md_at_2GeV']) + (1 / self.ip['ms_at_2GeV']))**(-1)

        if self.nucleon == 'p':
            if self.ip['FGtildep_pion'] == None:
                return - mN * mtilde * (1 / 2) * self.ip['gA'] * (1 / self.ip['mu_at_2GeV'] - 1 / self.ip['md_at_2GeV'])
            else: return self.ip['FGtildep_pion']
        if self.nucleon == 'n':
            if self.ip['FGtilden_pion'] == None:
                return - mN * mtilde * (1 / 2) * self.ip['gA'] * (1 / self.ip['mu_at_2GeV'] - 1 / self.ip['md_at_2GeV'])
            else: return self.ip['FGtilden_pion']

    def value_eta_pole(self):
        """ 
        Return the coefficient of the eta pole
        The eta pole is given, in terms of the spatial momentum q, by q^2 / (q^2 + meta^2)
        """
        mN = (1 / 2) * (self.ip['mproton'] + self.ip['mneutron'])
        mtilde = ((1 / self.ip['mu_at_2GeV']) + (1 / self.ip['md_at_2GeV']) + (1 / self.ip['ms_at_2GeV']))**(-1)
        Deltau = (1 / 2) * (self.ip['gA'] + self.ip['DeltaSigmaud'])
        Deltad = (1 / 2) * (- self.ip['gA'] + self.ip['DeltaSigmaud'])

        if self.nucleon == 'p':
            if self.ip['FGtildep_eta'] == None:
                return - mN * mtilde * (1 / 6) * (Deltau + Deltad - 2 * self.ip['Deltas'])\
                       * (1 / self.ip['mu_at_2GeV'] + 1 / self.ip['md_at_2GeV'] - 2 / self.ip['ms_at_2GeV'])
            else: return self.ip['FGtildep_eta']
        if self.nucleon == 'n':
            if self.ip['FGtilden_eta'] == None:
                return - mN * mtilde * (1 / 6) * (Deltau + Deltad - 2 * self.ip['Deltas'])\
                       * (1 / self.ip['mu_at_2GeV'] + 1 / self.ip['md_at_2GeV'] - 2 / self.ip['ms_at_2GeV'])
            else: return self.ip['FGtilden_eta']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        pion_pole = (qSq / (self.ip['mpi0']**2 - qSq)) * self.value_pion_pole()
        eta_pole = (qSq / (self.ip['meta']**2 - qSq)) * self.value_eta_pole()
        c = self.value_zero_mom()

        if self.nucleon == 'p':
            if self.ip['FGtildep'] == None:
                return pion_pole + eta_pole + c
            else: return self.ip['FGtildep']
        if self.nucleon == 'n':
            if self.ip['FGtilden'] == None:
                return pion_pole + eta_pole + c
            else: return self.ip['FGtilden']

####################################################################################
# ----------------------------------- Rayleigh ----------------------------------- #
####################################################################################

class Fgamma:
    def __init__(self, nucleon, input_dict):
        """ 
        The nuclear form factor Fgamma
        Return the nuclear form factor Fgamma
        Arguments:
        ---------
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        # T-odd - set to zero by default
        if self.nucleon == 'p':
            return self.ip['FFp_0']
        if self.nucleon == 'n':
            return self.ip['FFn_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.ip['FFp_prime'] == None:
                return 0.
            else: return self.ip['FFp_prime']
        if self.nucleon == 'n':
            if self.ip['FFn_prime'] == None:
                return 0.
            else: return self.ip['FFn_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.ip['FFp'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FFp']
        if self.nucleon == 'n':
            if self.ip['FFn'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FFn']

class Fgammatilde:
    def __init__(self, nucleon, input_dict):
        """ 
        The nuclear form factor Fgammatilde
        Return the nuclear form factor Fgammatilde
        Arguments:
        ---------
        nucleon = 'p', 'n' -- the nucleon (proton or neutron)
        input_dict (optional) -- a dictionary of hadronic input parameters
                                 (default is Num_input().input_parameters)
        """
        self.nucleon = nucleon
        self.ip = input_dict

    def value_zero_mom(self):
        """ 
        Return the value of the form factor at zero momentum transfer 
        """
        # T-odd - set to zero by default
        if self.nucleon == 'p':
            return self.ip['FFtildep_0']
        if self.nucleon == 'n':
            return self.ip['FFtilden_0']

    def first_deriv_zero_mom(self):
        """
        Return the value of the first derivative of the form factor
        w.r.t. q^2 at zero momentum transfer
        """
        if self.nucleon == 'p':
            if self.ip['FFtildep_prime'] == None:
                return 0.
            else: return self.ip['FFtildep_prime']
        if self.nucleon == 'n':
            if self.ip['FFtilden_prime'] == None:
                return 0.
            else: return self.ip['FFtilden_prime']

    def form_factor(self, q):
        """
        Return the full form factor
        """
        qSq = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3]
        if self.nucleon == 'p':
            if self.ip['FFtildep'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FFtildep']
        if self.nucleon == 'n':
            if self.ip['FFtilden'] == None:
                return self.value_zero_mom() + self.first_deriv_zero_mom() * qSq
            else: return self.ip['FFtilden']

if __name__ == '__main__':
    from parameters import NumericalInput
    input_dict = NumericalInput().input_parameters
    
    # Check numerical values of form factors for Al
    q = [0, 0, 0, 0.11081]

    # Vector F1 (proton)
    F1up = F1('u', 'p', input_dict)
    print('F1up @ q:', F1up.form_factor(q))
    F1dp = F1('d', 'p', input_dict)
    print('F1dp @ q:', F1dp.form_factor(q))
    F1sp = F1('s', 'p', input_dict)
    print('F1sp @ q:', F1sp.form_factor(q))
    print()

    # Vector F1 (neutron)
    F1un = F1('u', 'n', input_dict)
    print('F1un @ q:', F1un.form_factor(q))
    F1dn = F1('d', 'n', input_dict)
    print('F1dn @ q:', F1dn.form_factor(q))
    F1sn = F1('s', 'n', input_dict)
    print('F1sn @ q:', F1sn.form_factor(q))
    print()

    # Vector F2 (proton)
    F2up = F2('u', 'p', input_dict)
    print('F2up @ q:', F2up.form_factor(q))
    F2dp = F2('d', 'p', input_dict)
    print('F2dp @ q:', F2dp.form_factor(q))
    F2sp = F2('s', 'p', input_dict)
    print('F2sp @ q:', F2sp.form_factor(q))
    print()

    # Vector F2 (neutron)
    F2un = F2('u', 'n', input_dict)
    print('F2un @ q:', F2un.form_factor(q))
    F2dn = F2('d', 'n', input_dict)
    print('F2dn @ q:', F2dn.form_factor(q))
    F2sn = F2('s', 'n', input_dict)
    print('F2sn @ q:', F2sn.form_factor(q))
    print()

    # Vector F3 (proton) -- set to zero by default
    F3up = F3('u', 'p', input_dict)
    print('F3up @ q:', F3up.form_factor(q))
    F3dp = F3('d', 'p', input_dict)
    print('F3dp @ q:', F3dp.form_factor(q))
    F3sp = F3('s', 'p', input_dict)
    print('F3sp @ q:', F3sp.form_factor(q))
    print()

    # Vector F3 (neutron) -- set to zero by default
    F3un = F3('u', 'n', input_dict)
    print('F3un @ q:', F3un.form_factor(q))
    F3dn = F3('d', 'n', input_dict)
    print('F3dn @ q:', F3dn.form_factor(q))
    F3sn = F3('s', 'n', input_dict)
    print('F3sn @ q:', F3sn.form_factor(q))
    print()

    # Axial FA vector (proton)
    FAup = FA('u', 'p', input_dict)
    print('FAup @ q:', FAup.form_factor(q))
    FAdp = FA('d', 'p', input_dict)
    print('FAdp @ q:', FAdp.form_factor(q))
    FAsp = FA('s', 'p', input_dict)
    print('FAsp @ q:', FAsp.form_factor(q))
    print()

    # Axial FA vector (neutron)
    FAun = FA('u', 'n', input_dict)
    print('FAun @ q:', FAun.form_factor(q))
    FAdn = FA('d', 'n', input_dict)
    print('FAdn @ q:', FAdn.form_factor(q))
    FAsn = FA('s', 'n', input_dict)
    print('FAsn @ q:', FAsn.form_factor(q))
    print()

    # Axial vector FPprimed (proton)
    FPprimedup = FPprimed('u', 'p', input_dict)
    print('FPprimedup @ q:', FPprimedup.form_factor(q))
    FPprimeddp = FPprimed('d', 'p', input_dict)
    print('FPprimeddp @ q:', FPprimeddp.form_factor(q))
    FPprimedsp = FPprimed('s', 'p', input_dict)
    print('FPprimedsp @ q:', FPprimedsp.form_factor(q))
    print()

    # Axial vector FPprimed (neutron)
    FPprimedun = FPprimed('u', 'n', input_dict)
    print('FPprimedun @ q:', FPprimedun.form_factor(q))
    FPprimeddn = FPprimed('d', 'n', input_dict)
    print('FPprimeddn @ q:', FPprimeddn.form_factor(q))
    FPprimedsn = FPprimed('s', 'n', input_dict)
    print('FPprimedsn @ q:', FPprimedsn.form_factor(q))
    print()

    # Scalar FS (proton)
    FSup = FS('u', 'p', input_dict)
    print('FSup @ q:', FSup.form_factor(q))
    FSdp = FS('d', 'p', input_dict)
    print('FSdp @ q:', FSdp.form_factor(q))
    FSsp = FS('s', 'p', input_dict)
    print('FSsp @ q:', FSsp.form_factor(q))
    print()

    # Scalar FS (neutron) # DOES NOT AGREE WITH APPENDIX -- RECHECK
    FSun = FS('u', 'n', input_dict)
    print('FSun @ q:', FSun.form_factor(q))
    FSdn = FS('d', 'n', input_dict)
    print('FSdn @ q:', FSdn.form_factor(q))
    FSsn = FS('s', 'n', input_dict)
    print('FSsn @ q:', FSsn.form_factor(q))
    print()

    # Psuedoscalar FP (proton)
    FPup = FP('u', 'p', input_dict)
    print('FPup @ q:', FPup.form_factor(q))
    FPdp = FP('d', 'p', input_dict)
    print('FPdp @ q:', FPdp.form_factor(q))
    FPsp = FP('s', 'p', input_dict)
    print('FPsp @ q:', FPsp.form_factor(q))
    print()

    # Psuedoscalar FP (neutron)
    FPun = FP('u', 'n', input_dict)
    print('FPun @ q:', FPun.form_factor(q))
    FPdn = FP('d', 'n', input_dict)
    print('FPdn @ q:', FPdn.form_factor(q))
    FPsn = FP('s', 'n', input_dict)
    print('FPsn @ q:', FPsn.form_factor(q))
    print()

    # CP-even gluonic FG
    FGp = FG('p', input_dict)
    print('FGp @ q:', FGp.form_factor(q))
    FGn = FG('n', input_dict)
    print('FGn @ q:', FGn.form_factor(q))
    print()

    # CP-odd gluonic FGtilde
    FGtildep = FGtilde('p', input_dict)
    print('FGtildep @ q:', FGtildep.form_factor(q))
    FGtilden = FGtilde('n', input_dict)
    print('FGtilden @ q:', FGtilden.form_factor(q))
    print()

    # Tensor FT0 (proton)
    FT0up = FT0('u', 'p', input_dict)
    print('FT0up @ q:', FT0up.form_factor(q))
    FT0dp = FT0('d', 'p', input_dict)
    print('FT0dp @ q:', FT0dp.form_factor(q))
    FT0sp = FT0('s', 'p', input_dict)
    print('FT0sp @ q:', FT0sp.form_factor(q))
    print()

    # Tensor FT0 (neutron)
    FT0un = FT0('u', 'n', input_dict)
    print('FT0un @ q:', FT0un.form_factor(q))
    FT0dn = FT0('d', 'n', input_dict)
    print('FT0dn @ q:', FT0dn.form_factor(q))
    FT0sn = FT0('s', 'n', input_dict)
    print('FT0sn @ q:', FT0sn.form_factor(q))
    print()

    # Tensor FT1 (proton)
    FT1up = FT1('u', 'p', input_dict)
    print('FT1up @ q:', FT1up.form_factor(q))
    FT1dp = FT1('d', 'p', input_dict)
    print('FT1dp @ q:', FT1dp.form_factor(q))
    FT1sp = FT1('s', 'p', input_dict)
    print('FT1sp @ q:', FT1sp.form_factor(q))
    print()

    # Tensor FT1 (neutron)
    FT1un = FT1('u', 'n', input_dict)
    print('FT1un @ q:', FT1un.form_factor(q))
    FT1dn = FT1('d', 'n', input_dict)
    print('FT1dn @ q:', FT1dn.form_factor(q))
    FT1sn = FT1('s', 'n', input_dict)
    print('FT1sn @ q:', FT1sn.form_factor(q))
    print()

    # Tensor FT2 (proton)
    FT2up = FT2('u', 'p', input_dict)
    print('FT2up @ q:', FT2up.form_factor(q))
    FT2dp = FT2('d', 'p', input_dict)
    print('FT2dp @ q:', FT2dp.form_factor(q))
    FT2sp = FT2('s', 'p', input_dict)
    print('FT2sp @ q:', FT2sp.form_factor(q))
    print()

    # Tensor FT2 (neutron)
    FT2un = FT2('u', 'n', input_dict)
    print('FT2un @ q:', FT2un.form_factor(q))
    FT2dn = FT2('d', 'n', input_dict)
    print('FT2dn @ q:', FT2dn.form_factor(q))
    FT2sn = FT2('s', 'n', input_dict)
    print('FT2sn @ q:', FT2sn.form_factor(q))
    print()

    # Tensor FT3 (proton) -- set to zero by default
    FT3up = FT3('u', 'p', input_dict)
    print('FT3up @ q:', FT3up.form_factor(q))
    FT3dp = FT3('d', 'p', input_dict)
    print('FT3dp @ q:', FT3dp.form_factor(q))
    FT3sp = FT3('s', 'p', input_dict)
    print('FT3sp @ q:', FT3sp.form_factor(q))
    print()

    # Tensor FT3 (neutron) -- set to zero by default
    FT3un = FT3('u', 'n', input_dict)
    print('FT3un @ q:', FT3un.form_factor(q))
    FT3dn = FT3('d', 'n', input_dict)
    print('FT3dn @ q:', FT3dn.form_factor(q))
    FT3sn = FT3('s', 'n', input_dict)
    print('FT3sn @ q:', FT3sn.form_factor(q))
    print()

    # Gluonic-even (proton)
    FGp = FG('p', input_dict)
    print('FGp @ q:', FGp.form_factor(q))
    # Gluonic-even (neutron)
    FGn = FG('n', input_dict)
    print('FGn @ q:', FGn.form_factor(q))
    print()

    # Gluonic-odd (proton)
    FGtp = FGtilde('p', input_dict)
    print('FGtp @ q:', FGtp.form_factor(q))
    # Gluonic-odd (neutron)
    FGtn = FGtilde('n', input_dict)
    print('FGtn @ q:', FGtn.form_factor(q))
    print()

    # Rayleigh-even (proton)
    FFp = Fgamma('p', input_dict)
    print('FFp @ q:', FFp.form_factor(q))
    # Rayleigh-even (neutron)
    FFn = Fgamma('n', input_dict)
    print('FFn @ q:', FFn.form_factor(q))
    print()

    # Rayleigh-odd (proton)
    FFtp = Fgammatilde('p', input_dict)
    print('FFtp @ q:', FFtp.form_factor(q))
    # Rayleigh-odd (neutron)
    FFtn = Fgammatilde('n', input_dict)
    print('FFtn @ q:', FFtn.form_factor(q))
    print()
import sys
import numpy as np
import pathlib
from hadronization import Hadronization
from parameters import NumericalInput

# We have to figure out how to import Mu2eNRET.py
# from ../../Mu2E_NRET/v2/python

# we are in <top>/MuonConverter/python
# __file__ is the current file path
topdir=pathlib.Path(__file__).parent.parent.parent
mu2edir=topdir / "Mu2e_NRET" / "v2" / "python"
mu2edir=str(mu2edir.resolve())
# Add to search path - this is a bit of hack, but works
sys.path.insert(1, mu2edir) # insert after script path

# now import of Mu2eNRET will work
import Mu2eNRET

quit()

class Interface:
	def __init__(self, WET3_dict, q, input_dict = None):
		"""
		Provide an interface between external EFT codes (e.g. wilson) and the Mu2e_NRET software

		Input:
		WET3_dict (dictionary): Dictionary composed of mue sector numerical Wilson coefficients in the WET-3 JMS basis.
		"""
		self.WET3_dict = WET3_dict
		self.q = q
		# The dictionary of input parameters
		if input_dict == None:
			self.ip = NumericalInput().input_parameters
		else:
			self.ip = NumericalInput(my_input_dict = input_dict).input_parameters
	
	#########################################################################################################################################
	# ------------------------------------------------- Interface with external EFT codes ------------------------------------------------- #
	#########################################################################################################################################

	def JMS_to_MuonConverter(self):
		"""
		Given a Wilson coefficient dictionary in the WET-3 JMS basis (WET3), convert
		to the WCxf MuonConverter (MC) RET basis.
	
		Output:
		MC_WCxf_dict (dictionary): Dictionary containing translated coefficients in the WCxf MuonConverter nucleon RET basis
		"""
		# Initialize a MuonConverter dictionary
		MC_WCxf_dict = {'Tegamma_12': 0., 'ATegamma_12': 0., 'VVeu_1211': 0., 'VVed_1211': 0., 'VVed_1222': 0., 'AVVeu_1211': 0., 'AVVed_1211': 0., 'AVVed_1222': 0.,
						'VAVeu_1211': 0., 'VAVed_1211': 0., 'VAVed_1222': 0., 'AVAVeu_1211': 0., 'AVAVed_1211': 0., 'AVAVed_1222': 0., 'SSeu_1211': 0., 'SSed_1211': 0.,
						'SSed_1222': 0., 'ASeu_1211': 0., 'ASed_1211': 0., 'ASed_1222': 0., 'SAeu_1211': 0., 'SAed_1211': 0., 'SAed_1222': 0., 'AAeu_1211': 0.,
						'AAed_1211': 0., 'AAed_1222': 0., 'TTeu_1211': 0., 'TTed_1211': 0., 'TTed_1222': 0., 'ATTeu_1211': 0., 'ATTed_1211': 0., 'ATTed_1222': 0.}
	
		MC_WCxf_dict['Tegamma_12'] = 4. * np.pi**2 * self.WET3_dict['egamma_12']
		MC_WCxf_dict['ATegamma_12'] = - 4. * np.pi**2 * 1j * self.WET3_dict['egamma_12']
		
		MC_WCxf_dict['VVeu_1211'] = (self.WET3_dict['VeuLL_1211'] + self.WET3_dict['VeuRR_1211'] + self.WET3_dict['VeuLR_1211'] + self.WET3_dict['VueLR_1112']) / 4.
		MC_WCxf_dict['VVed_1211'] = (self.WET3_dict['VedLL_1211'] + self.WET3_dict['VedRR_1211'] + self.WET3_dict['VedLR_1211'] + self.WET3_dict['VdeLR_1112']) / 4.
		MC_WCxf_dict['VVed_1222'] = (self.WET3_dict['VedLL_1222'] + self.WET3_dict['VedRR_1222'] + self.WET3_dict['VedLR_1222'] + self.WET3_dict['VdeLR_2212']) / 4.
		
		MC_WCxf_dict['AVVeu_1211'] = (-self.WET3_dict['VeuLL_1211'] + self.WET3_dict['VeuRR_1211'] - self.WET3_dict['VeuLR_1211'] + self.WET3_dict['VueLR_1112']) / 4.
		MC_WCxf_dict['AVVed_1211'] = (-self.WET3_dict['VedLL_1211'] + self.WET3_dict['VedRR_1211'] - self.WET3_dict['VedLR_1211'] + self.WET3_dict['VdeLR_1112']) / 4.
		MC_WCxf_dict['AVVed_1222'] = (-self.WET3_dict['VedLL_1222'] + self.WET3_dict['VedRR_1222'] - self.WET3_dict['VedLR_1222'] + self.WET3_dict['VdeLR_2212']) / 4.
		
		MC_WCxf_dict['VAVeu_1211'] = (-self.WET3_dict['VeuLL_1211'] + self.WET3_dict['VeuRR_1211'] + self.WET3_dict['VeuLR_1211'] - self.WET3_dict['VueLR_1112']) / 4.
		MC_WCxf_dict['VAVed_1211'] = (-self.WET3_dict['VedLL_1211'] + self.WET3_dict['VedRR_1211'] + self.WET3_dict['VedLR_1211'] - self.WET3_dict['VdeLR_1112']) / 4.
		MC_WCxf_dict['VAVed_1222'] = (-self.WET3_dict['VedLL_1222'] + self.WET3_dict['VedRR_1222'] + self.WET3_dict['VedLR_1222'] - self.WET3_dict['VdeLR_2212']) / 4.
	
		MC_WCxf_dict['AVAVeu_1211'] = (self.WET3_dict['VeuLL_1211'] + self.WET3_dict['VeuRR_1211'] - self.WET3_dict['VeuLR_1211'] - self.WET3_dict['VueLR_1112']) / 4.
		MC_WCxf_dict['AVAVed_1211'] = (self.WET3_dict['VedLL_1211'] + self.WET3_dict['VedRR_1211'] - self.WET3_dict['VedLR_1211'] - self.WET3_dict['VdeLR_1112']) / 4.
		MC_WCxf_dict['AVAVed_1222'] = (self.WET3_dict['VedLL_1222'] + self.WET3_dict['VedRR_1222'] - self.WET3_dict['VedLR_1222'] - self.WET3_dict['VdeLR_2212']) / 4.
	
		MC_WCxf_dict['SSeu_1211'] = (self.WET3_dict['SeuRL_1211'] + self.WET3_dict['SeuRR_1211']) / 4.
		MC_WCxf_dict['SSed_1211'] = (self.WET3_dict['SedRL_1211'] + self.WET3_dict['SedRR_1211']) / 4.
		MC_WCxf_dict['SSed_1222'] = (self.WET3_dict['SedRL_1211'] + self.WET3_dict['SedRR_1222']) / 4.
	
		MC_WCxf_dict['ASeu_1211'] = -1j * (self.WET3_dict['SeuRL_1211'] + self.WET3_dict['SeuRR_1211']) / 4.
		MC_WCxf_dict['ASed_1211'] = -1j * (self.WET3_dict['SedRL_1211'] + self.WET3_dict['SedRR_1211']) / 4.
		MC_WCxf_dict['ASed_1222'] = -1j * (self.WET3_dict['SedRL_1222'] + self.WET3_dict['SedRR_1222']) / 4.
	
		MC_WCxf_dict['SAeu_1211'] = 1j * (self.WET3_dict['SeuRL_1211'] - self.WET3_dict['SeuRR_1211']) / 4.
		MC_WCxf_dict['SAed_1211'] = 1j * (self.WET3_dict['SedRL_1211'] - self.WET3_dict['SedRR_1211']) / 4.
		MC_WCxf_dict['SAed_1222'] = 1j * (self.WET3_dict['SedRL_1222'] - self.WET3_dict['SedRR_1222']) / 4.
	
		MC_WCxf_dict['AAeu_1211'] = (self.WET3_dict['SeuRL_1211'] - self.WET3_dict['SeuRR_1211']) / 4.
		MC_WCxf_dict['AAed_1211'] = (self.WET3_dict['SedRL_1211'] - self.WET3_dict['SedRR_1211']) / 4.
		MC_WCxf_dict['AAed_1222'] = (self.WET3_dict['SedRL_1222'] - self.WET3_dict['SedRR_1222']) / 4.
	
		MC_WCxf_dict['TTeu_1211'] = self.WET3_dict['TeuRR_1211'] / 4.
		MC_WCxf_dict['TTed_1211'] = self.WET3_dict['TedRR_1211'] / 4.
		MC_WCxf_dict['TTed_1222'] = self.WET3_dict['TedRR_1222'] / 4.
	
		MC_WCxf_dict['ATTeu_1211'] = -1j * self.WET3_dict['TeuRR_1211'] / 4.
		MC_WCxf_dict['ATTed_1211'] = -1j * self.WET3_dict['TedRR_1211'] / 4.
		MC_WCxf_dict['ATTed_1222'] = -1j * self.WET3_dict['TedRR_1222'] / 4.
	
		return MC_WCxf_dict
	
	def WCxf_to_Cmueqq(self):
		"""
		Convert from WCxf naming conventions to the internal C^{(d)}_{e,mu,q,q} naming convention
		"""
		MC_WCxf_dict = self.JMS_to_MuonConverter()
		# Create one final conversion dictionary for input
		MC_Cmueqq_dict = {'C51': MC_WCxf_dict['Tegamma_12'], 'C52': MC_WCxf_dict['ATegamma_12'], 
						  'C61u': MC_WCxf_dict['VVeu_1211'], 'C61d': MC_WCxf_dict['VVed_1211'], 'C61s': MC_WCxf_dict['VVed_1222'], 
						  'C62u': MC_WCxf_dict['AVVeu_1211'], 'C62d': MC_WCxf_dict['AVVed_1211'], 'C62s': MC_WCxf_dict['AVVed_1222'], 
						  'C63u': MC_WCxf_dict['VAVeu_1211'], 'C63d': MC_WCxf_dict['VAVed_1211'], 'C63s': MC_WCxf_dict['VAVed_1222'], 
						  'C64u': MC_WCxf_dict['AVAVeu_1211'], 'C64d': MC_WCxf_dict['AVAVed_1211'], 'C64s': MC_WCxf_dict['AVAVed_1222'],
						  'C65u': MC_WCxf_dict['SSeu_1211'], 'C65d': MC_WCxf_dict['SSed_1211'], 'C65s': MC_WCxf_dict['SSed_1222'],
						  'C67u': MC_WCxf_dict['ASeu_1211'], 'C67d': MC_WCxf_dict['ASed_1211'], 'C67s': MC_WCxf_dict['ASed_1222'],
						  'C66u': MC_WCxf_dict['SAeu_1211'], 'C66d': MC_WCxf_dict['SAed_1211'], 'C66s': MC_WCxf_dict['SAed_1222'],
						  'C68u': MC_WCxf_dict['AAeu_1211'], 'C68d': MC_WCxf_dict['AAed_1211'], 'C68s': MC_WCxf_dict['AAed_1222'],
						  'C69u': MC_WCxf_dict['TTeu_1211'], 'C69d': MC_WCxf_dict['TTed_1211'], 'C69s': MC_WCxf_dict['TTed_1222'],
						  'C610u': MC_WCxf_dict['ATTeu_1211'], 'C610d': MC_WCxf_dict['ATTed_1211'], 'C610s': MC_WCxf_dict['ATTed_1222']}
		return MC_Cmueqq_dict
	
	#########################################################################################################################################
	# ----------------------------------------------------- Interface with Mu2e_NRET ------------------------------------------------------ #
	#########################################################################################################################################
	def JMS_to_RET(self):
		"""
		Compose a series of translations and return RET "d"-coefficients for input into Mu2e_NRET
		"""
		# Convert to the internal Cmueqq basis
		MC_WET_dict = self.WCxf_to_Cmueqq()
		# Hadronize
		hadronize = Hadronization(MC_WET_dict, self.ip)
		MC_RET_isospin_array = hadronize.hadronize_WET(self.q)
		return MC_RET_isospin_array

	def HMMRZ_to_RET(self):
		"""
		Compose a series of translations and return RET "d"-coefficients for input into Mu2e_NRET
		"""
		# Hadronize
		hadronize = Hadronization(self.WET3_dict, self.ip)
		MC_RET_isospin_array = hadronize.hadronize_WET(self.q)
		return MC_RET_isospin_array

	def JMS_to_RET_map(self):
		"""
		Generate a map from the JMS WET-3 basis to the RET basis (d-coefficients)
		consistent with the NetworkX syntax
		"""
		# Initialize hadronization map
		had_map = []
		# Initialize WET keys
		WET_keys = ['egamma_12', 'VeuLL_1211', 'VedLL_1211', 'VedLL_1222', 'VeuRR_1211', 'VedRR_1211', 'VedRR_1222', 'VeuLR_1211',
					'VedLR_1211', 'VedLR_1222','VueLR_1112', 'VdeLR_1112', 'VdeLR_2212', 'SeuRL_1211', 'SedRL_1211', 'SedRL_1222', 
					'SeuRR_1211', 'SedRR_1211', 'SedRR_1222', 'TeuRR_1211', 'TedRR_1211', 'TedRR_1222']

		# Initialize iterator
		iterator = 1
		# Make a copy of WET dictionary
		WET_dict_clone = self.WET3_dict.copy()
		for keys in WET_keys:
			# Initialize zero dictionary
			for key in WET_keys:
				WET_dict_clone[key] = 0.
			# Set the key to 1 (only interested in generated non-zero entries in d)
			WET_dict_clone[keys] = 1.0
			# Compute d coefficients with dummy q-vector
			d_isospin = self.JMS_to_RET()
			# Record non-zero values
			for i in range(len(d_isospin)):
				# Append non-zero values to map
				if d_isospin[i][1] != 0:
					had_map.append((iterator,i))
			iterator += 1
		return had_map

	def compute_rate(self, element, isotope, interaction, oscb, isochar, input_basis, plots = "None", write_yaml = False, path = './Mu2e_NRET_data.yaml'):
		"""
		Compute the relative capture ratio for muon to electron conversion in the field of an "element" nucleus
	
		element (string): ------- target element -- can use UI element spelling Sulfur or it's symbol S (also Britsh sp)
		isotope (int): ---------- target element isotope -- an input of 0 averages over all isotopes
		interaction (string): --- type of interaction model.
							  	  Options:
							  	  'ck': ['ck', "Cohen and Kurath"],
    						  	  'bw': ['bw', "Brown-Wildenthal"],
    						  	  'usda': ['usda', "USDA"],
    						  	  'usdb': ['usdb', "USDB"],
    						  	  '4hw':  ['4hw', "4hw(16O)"],
    						  	  '2hw':  ['2hw', "2hw(18O)"],
    						  	  'kbp':  ['kbp', "KBP"],
    						  	  'gx1a': ['GX1A', "GX1A"],
    						  	  'kb3g': ['KB3G', "KB3G"],
    						  	  'gcn2850': ['GCN2850', "GCN2850"],
    						  	  'jj44b':   ['jj44b', "jj44b"],
    						  	  'jun45':   ['JUN45', "JUN45"],
		oscb:(int): ------------- Oscillator length scale (in fm). Omit or set to 0 to calc from average A (Abar) over isotopes
							  	  using:  b=Sqrt[41.467/(45*Abar^(-1/3) -25*Abar^(-2/3))]. Automatically fixed for some interactions.
		isochar (string or list): [1, 0]/isoscalar: isoscalar (1)
								  [0, 1]/isovector: isovector (tau_3)
								  [0.5, 0.5]/proton proton
								  [0.5,-0.5]/neutron neutron
		plots (list, string): --- vcrm: Vector Charge (or S.I.) Response (M)
							  	  alsr: Axial Longitudinal Spin Response (Sigma'') 
							  	  alsr: Axial Transverse Spin Response (Sigma')
							  	  ssd:  Standard Spin-Dependent (or S.D.) Response
							  	  vtmr: Vector Transverse Magnetic Response (Delta)
							  	  vlr:  Vector Longitudinal Response (Phi'')
							  	  vter: Vector Transverse Electric Response (PhiT')
							  	  all:  All of the above
							  	  none: None of the above
	
		path (string): ---------- path to where the data file should be saved. Defaults to './Mu2e_NRET_data.yaml'
		"""
		# Given WET3-JMS coefficients generate nuclone RET coefficients in the isospin basis
		if input_basis == 'JMS':
			ds = self.JMS_to_RET()
		elif input_basis == 'HMMRZ':
			ds = self.HMMRZ_to_RET()
		
		# Format the Wilson coefficients for input into the yaml file
		data = {'Element': element,
				'Isotope': isotope,
				'Interaction': interaction,
				'oscb': oscb,
				'isochar': isochar,
				'plots': plots,
				'mL': self.ip['mL'] * 1e3,
				'ds': ds}

		if write_yaml:
			yaml_data = 'Element: ' + element + '\n'\
						+ 'Isotope: ' + str(isotope) + '\n'\
						+ 'Interaction: ' + interaction + '\n'\
						+ 'oscb: ' + str(oscb) + '\n'\
						+ 'isochar: ' + isochar + '\n'\
						+ 'plots: ' + plots + '\n'\
						+ 'mL: ' + str(self.ip['mL'] * 1e3) + '\n'\
						+ 'ds: ' + str(ds) + '\n'
			import os.path
			output_file = str(os.path.expanduser(path))
			# Write the data file
			with open(output_file,'w') as f:
				f.write(yaml_data)

		# Compute the capture ratio
		ndata = Mu2eNRET.processdata(data)

		# Account for the normalization used in Mu2e_NRET code
		v = 246.200 # GeV
		# Return the branching ratio [dimensionless] and the decay rate [1/s]
		return v**4 * ndata['BranchingRatio'], v**4 * ndata['DecayRate'], ndata

if __name__ == "__main__":
	WET = {'egamma_12': -1.82228281e-14+1.75403582e-26j, 'VeuLL_1211': -9.56894111e-07+3.8245183e-21j, 'VedLL_1211': -9.89568678e-07+3.95510768e-21j, 
		   'VedLL_1222': -1.14666907e-10+4.53279135e-25j, 'VeuRR_1211': 2.66061e-20-1.10778543e-34j, 'VedRR_1211': -1.69273191e-20+6.9945152e-35j, 
		   'VedRR_1222': -1.69273184e-20+6.7305892e-35j, 'VeuLR_1211': 2.57046892e-10-1.02325324e-24j, 'VedLR_1211': -5.8115279e-11+2.30273539e-25j, 
		   'VedLR_1222': -5.81152727e-11+2.77900618e-25j, 'VueLR_1112': -5.28241239e-17+2.11732057e-31j, 'VdeLR_1112': -5.4624664e-17+2.18948049e-31j, 
		   'VdeLR_2212': -1.47378766e-20+6.01355236e-35j, 'SeuRL_1211': 0, 'SedRL_1211': 2.47508174e-16-9.89414374e-31j, 'SedRL_1222': -1.58467801e-17-1.53126562e-30j, 
		   'SeuRR_1211': -1.16478067e-16+2.17068834e-28j, 'SedRR_1211': 0, 'SedRR_1222': 0, 'TeuRR_1211': 3.16614544e-17-1.38571535e-29j, 'TedRR_1211': 0, 'TedRR_1222': 0}
	q = [1, 1, 1, 1]
	mu2e = Interface(WET, q)
	element     = 'Al'
	isotope     = 0
	interaction = 'bw'
	oscb        = 0
	isochar     = 'proton'
	CR = mu2e.compute_rate(element, isotope, interaction, oscb, isochar, input_basis = 'JMS')
	#print(CR)
	#had_map = mu2e.JMS_to_RET_map()
	#print(had_map)

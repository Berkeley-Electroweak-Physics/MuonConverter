# MuonConverter
The main purpose of ```MuonConverter``` is to provide an interface between external EFT software such as ```wilson``` ([https://wilson-eft.github.io/](https://wilson-eft.github.io/)), ```DsixTools``` ([https://dsixtools.github.io/](https://dsixtools.github.io/)), etc, and ```Mu2e_NRET```. This interface extends the functionality of the original ```Mu2e_NRET``` code and allows for full top-down (or bottom-up) phenomenological studies of muon-to-electron conversion in the field of a target nucleus. Explicitly, in conjunction with external EFT software, ```MuonConverter``` can be used to compute the influence of UV charged-lepton-flavor-violating operators on the predictions for branching and capture ratios reported by experimental collaborations.

Both versions of ```MuonConverter``` (Python and Mathematica) are comprised of four modular components:
- Numerical inputs --- all numerical inputs are stored within an associative array that can be modified by the user upon intialization of the ```MuonConverter``` class. The parameters and their default values can be found in ```parameters.py```(```.wl```).
- Form factors --- the form factor expressions required for the WET to NRET matching and whose numerical values are derived in App. D can be found in ```form_factors.py```(```.wl```). For maximum flexibility, the default form factor values may be manually overwritten within ```parameters.py```(```.wl```). 
- Matching --- to facilitate the WET to NRET matching, ```MuonConverter``` utilizes derived matching expressions for the relativistic $d_i$ coefficients (the $d_i$ coefficients are automatically translated to the nonrelativistic $c_i$, $b_i$ coefficients within ```Mu2e_NRET```). The matching expressions, as well as their translation to the isospin basis, can be found in ```hadronization.py```(```.wl```).
- Interfacing --- given an array of WET coefficients (in units of 1/GeV^2), the interface with ```Mu2e_NRET```, utilizing external and internal basis translations as well as the matching expressions implemented in ```hadronization.py```(```.wl```), can be found in ```MuonConverter.py```(```.wl```).

See the respecitve example notebooks within the desired language directory for usage examples - including interfacing with external codes.

# Dependencies
The Python version of the software was developed in Python 3.10 - a detailed list of software requirements/dependencies can be found in requirements.txt.

If using conda, the environment can be cloned via the command
```shell
conda create --name <env> --file requirements.txt
```
where ```<env>``` is the desired environment name.

The Mathematica version of the software was developed in Mathematica 14.0.
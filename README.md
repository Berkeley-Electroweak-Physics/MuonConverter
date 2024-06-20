# MuonConverter
The main purpose of \texttt{MuonConverter} is to provide an interface between external EFT software such as \texttt{wilson} \cite{Aebischer:2018bkb}, \texttt{DsixTools} \cite{Celis:2017hod,Fuentes-Martin:2020zaz}, etc, and the $\texttt{Mu2e\_NRET}$ software developed in \cite{Haxton:2022piv}. This interface extends the functionality of the original $\texttt{Mu2e\_NRET}$ code and allows for full top-down (or bottom-up) phenomenological studies of muon-to-electron conversion in the field of a target nucleus. Explicitly, in conjunction with external EFT software, $\texttt{MuonConverter}$ can be used to compute the influence of UV charged-lepton-flavor-violating operators on the predictions for branching and capture ratios reported by experimental collaborations.

Both versions of \texttt{MuonConverter} (Python and Mathematica) are comprised of four modular components:
    - Numerical inputs --- all numerical inputs are stored within an associative array that can be modified by the user upon intialization of the \texttt{MuonConverter} class. The parameters and their default values can be found in $\texttt{parameters.py}$($\texttt{.wl}$).
    - Form factors --- the form factor expressions required for the WET to NRET matching and whose numerical values are derived in App. D can be found in $\texttt{form\_factors.py}$($\texttt{.wl}$). For maximum flexibility, the default form factor values may be manually overwritten within $\texttt{parameters.py}$($\texttt{.wl}$). 
    - Matching --- to facilitate the WET to NRET matching, $\texttt{MuonConverter}$ utilizes derived matching expressions for the relativistic $d_i$ coefficients (the $d_i$ coefficients are automatically translated to the nonrelativistic $c_i$, $b_i$ coefficients within $\texttt{Mu2e\_NRET}$). The matching expressions, as well as their translation to the isospin basis, can be found in $\texttt{hadronization.py}$($\texttt{.wl}$).
    - Interfacing --- given an array of WET coefficients (in units of GeV$^{-2}$), the interface with $\texttt{Mu2e\_NRET}$, utilizing external and internal basis translations as well as the matching expressions implemented in $\texttt{hadronization.py}$($\texttt{.wl}$), can be found in $\texttt{MuonConverter.py}$($\texttt{.wl}$).

See the respecitve example notebooks within the desired language directory for usage examples - including interfacing with external codes.

# Dependencies
The Python version of the software was developed in Python 3.10 - a detailed list of required software can be found in dependencies.txt.

If using conda, the environment can be replicated via the command
```shell
conda create --name <env> --file dependencies.txt
```
where ```<env>``` is the desired environment name.

The Mathematica version of the software was developed in Mathematica 14.0.
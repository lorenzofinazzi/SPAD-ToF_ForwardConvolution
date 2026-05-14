# SPAD-ToF_ForwardConvolution
MATLAB code for forward-convolution dispersion analysis through TCSPC

InGaAs folder (results for 1550 nm range) and Si folder (results for 780 nm range).

Indications are valid for both folders.
DataPreparation.m code uploads in MATLAB environment most of the parameters. Manipulation is performed to make data usable. The script produce DataForAnalysis.mat file. This contains manipulated data for dispersion analysis.
Dispersion_DataAnalysis.m code needs DataForAnalysis.mat file. This code implements dispersive effect and produced simulated TCSPC curves. It introduces set up's Instrument Response Function effect on simulated curves. The code can perform sweep in dispersive medium length to analyze the minimum resolvable length.

Provided Txt files include raw data from experiments and simulations (source spectrum, experimental TCSPC curves, Lumerical simulated D curves, set up's IRF).

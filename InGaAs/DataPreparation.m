
%%
% Title: SPAD-based Time-of-Flight for dispersion measurement of integrated waveguides
% Author: Lorenzo Finazzi
% Year: 2026
% 
% Licensed under the MIT License.
% Please cite: SPAD-based Time-of-Flight for dispersion measurement of integrated waveguides
%%


clear all;
clc;

%%%% Upload laser spectrum (with PDE correction applied) %%%%

S = readmatrix('LaserSpectrum1550nm_PDEnorm.txt');
% Save wavelength (m → nm conversion)
lambda_nm = S(:,1)*1e9;
% Save intensity
intensity_spectrum = S(:,2);

% Wavelength → frequency
c = 299792458;                      % m/s
lambda_m = lambda_nm * 1e-9;
% Frequency ascending order
f = flip(c ./ lambda_m);            

%%% Manage spectrum to respect Nyquist for FFT %%%
% Number of interpolation point
N =2^21;
% Reduced number of interpolation point (help with computing)
N_reduced = 2^17;
% Spectrum zero-padding
f(3:end+2) = f;
f(1) = 0.92e14; %1 THz before central spectrum peak
f(2) = 1.7e14; %auxiliary point to help smooth interpolation
f(end+1) = 2.35e14; %auxiliary point to help smooth interpolation
f(end+1) = 2.92e14; %1 THz after central spectrum peak
% Intensity zero-padding
intensity_spectrum(3:end+2) = intensity_spectrum;
intensity_spectrum(1) = min(intensity_spectrum);
intensity_spectrum(2) = min(intensity_spectrum);
intensity_spectrum(end+1) = min(intensity_spectrum);
intensity_spectrum(end+1) = min(intensity_spectrum);

% Frequency interpolation
f_min = min(f);
f_max = max(f);
f_eq = linspace(f_min, f_max, N);   % Freq interp
f_eq_reduced = linspace(f_min, f_max, N_reduced); % Freq interp (reduced point)

% Spectrum intensity → electric field + interpolation
Ew_mag = sqrt(flip(intensity_spectrum));   % Electric field (magnitude)
Ew_mag_interp = interp1(f, Ew_mag, f_eq, 'pchip');
Ew_mag_interp = Ew_mag_interp-min(Ew_mag_interp); % Remove baseline
Ew_mag_interp = Ew_mag_interp/max(Ew_mag_interp); % Normalize

%%%% Upload experimental TCSPC curve after dispersive medium %%%%

T = readmatrix('Output_TCSPC_1550nmWG.txt');

time = T(:,1)*1000;   % nm to ps
intensity = T(:,2);
t_min = min(time);
t_max = max(time);
t_eq = linspace(t_min, t_max, length(f_eq)); % same number of point of f_eq
Et_meas = sqrt(intensity);          % measured electric field 
Et_meas = sgolayfilt(Et_meas, 3, 11);  % smoothing curves (pay attention not to change pulse FWHM!)
Et_meas = Et_meas-min(Et_meas); % remove baseline
Et_interp = interp1(time, Et_meas, t_eq, 'pchip');  % interpolation on fitted time
Et_interp = flip(Et_interp/max(Et_interp)); % normalization

t_eq = t_eq-t_eq(ceil(end/2)); %time axis centered on 0

%%%% Dispersion parameters

D = readmatrix('D_LumericalMODE_1550nmWG.txt'); %Lumerical simulation dispersion D curve of "1550 nm waveguide" (fundamental mode)
wl = D(:,1)*1e-6;   % save wavelength in m

dispersion = D(:,2);

f_D = flip(c ./ wl); % wl to freq
D_interp_f = interp1(f_D, flip(dispersion), f_eq, 'pchip'); % interpolation of D on f_eq


beta2_lambda = - (wl.^2) .* dispersion * 1e-6 ./ (2*pi*c);  % [s^2/m]1e-6 to convert to std unit   % beta2 formula from D
beta2_f = flip(beta2_lambda);


% extend frequency and beta2 to match f range from spectrum (see zero-padding above)
f_D(2:end+1) = f_D;
f_D(1) = min(f_eq); 
f_D(end+1) = f_D(end)+1; 

beta2_f(2:end+1) = beta2_f;
beta2_f(1) = beta2_f(2);
beta2_f(end+1) = beta2_f(end);

beta2_f_interp = interp1(f_D, beta2_f, f_eq, 'pchip');  % interpolate beta 2 on f_eq
beta2_f_reduced = interp1(f_D, beta2_f, f_eq_reduced, 'pchip'); % interpolate beta 2 on f_eq (reduced points)


% save data for next analysis
save('DataForAnalysis.mat', 'Ew_mag_interp' , 'f_eq', 'f_eq_reduced' ,'Et_interp', 't_eq', 'beta2_f_interp', 'beta2_f_reduced');


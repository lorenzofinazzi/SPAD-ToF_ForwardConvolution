
%%
% Title: SPAD-based Time-of-Flight for dispersion measurement of integrated waveguides
% Author: Lorenzo Finazzi
% Year: 2026
% 
% Licensed under the MIT License.
% Please cite: SPAD-based Time-of-Flight for dispersion measurement of integrated waveguides
%%


clear; clc;

%% Upload data from DataForAnalysis.mat
[file, path] = uigetfile('DataForAnalysis.mat'); 
if file == 0, return; end
load(fullfile(path, file));   
i = 0;

%% Data

% beta2_f_i = beta2_f_interp(:);    %beta 2 interpolated on full number of points -- high computing effort!
beta2_f_r = beta2_f_reduced(:); %beta 2 interpolated on reduced number of points

f = f_eq(:);
f_reduced = f_eq_reduced(:);
f_min = min(f);
f_max = max(f);

N = length(f);

% frequency to angular frequency 
omega_reduced = 2*pi*f_reduced;
omega_eq = 2*pi*f_eq;
omega_reduced = omega_reduced(:);

%% Waveguide parameters
%either single L or cycle on L can be performed (uncomment for and "end"
%for the cycle)

L = 10.06;

%L cycle
% for L = 0.2:0.2:10




%% Dispersive GVD phase calculation and IFFT

phi_reduced_2 = zeros(size(omega_reduced)); % initialize phase array
phi_reduced_2(1) = 0;
for idx = 2:length(omega_reduced)
    omega_int = omega_reduced(1:idx);
    beta2_int = beta2_f_r(1:idx);
    integrate_2 = beta2_int .* (omega_reduced(idx) - omega_int);
    phi_reduced_2(idx) = L * trapz(omega_int, integrate_2);
    idx;
end

phi_interp_2 = interp1(omega_reduced, phi_reduced_2, omega_eq, 'pchip'); % interpolate phase on omega

H_2 = exp(-1i * phi_interp_2(:));   % GVD phase

phi0 = zeros(size(Ew_mag_interp)); % input pulse's phase = 0 (if transform-limited)

Ew_in = Ew_mag_interp .* exp(1i * phi0);  % complex input spectrum
Et_in_t = fftshift(ifft(ifftshift(Ew_in)));   % input pulse (time domain)
Et_in_t = Et_in_t/max(abs(Et_in_t));

Ew_out_2 = Ew_in'.*H_2; % complex output spectrum with GVD effect
Et_out_t_2 = fftshift(ifft(ifftshift(Ew_out_2)));   % output pulse (time domain)
Et_out_t_2 = Et_out_t_2/max(abs(Et_out_t_2));

% time axis from IFFT
dt = 1/(f_max-f_min);   %time step
t_out = (-N/2:N/2-1)*dt*1e12;   %time axis


% abs value
Et_out_t_2 = abs(Et_out_t_2);
Et_interp = abs(Et_interp);

% electric field to intensity
I_simulated_2 = (Et_out_t_2).^2;
I_measured = (Et_interp).^2;

%% Apply setup IRF 

Ew_out_2 = Ew_out_2(:);  
t_out = t_out(:)/1e12; %t_out ps to s
t_eq = t_eq(:);  %experiment time
Et_interp = Et_interp(:);   %measured electric field
Et_out_t_2 = Et_out_t_2(:);     %simulated output without spad response only GVD
I_2 = abs(Et_out_t_2).^2;

% Upload IRF
T = readmatrix('IRF_SPAD_InGaAs.txt');
time_IRF = T(:,1)*1e-9; %s
time_IRF = max(time_IRF)-time_IRF; %start time axis from 0
IRF_data = T(:,2)-2250; %-2250 is intended to remove DCR, modify accordingly

t_IRF = (min(time_IRF):5e-15:max(time_IRF)); %time axis for convolution
IRF_interp = interp1(time_IRF, IRF_data, t_IRF);
IRF_interp = IRF_interp / trapz(t_IRF, IRF_interp); %normalization for convolution

dt = mean(diff(t_out));

disp('Start convolution between simulated signal and IRF')

M_full_2 = conv(abs(Et_out_t_2).^2, IRF_interp, 'full') * dt;
M_full_2 = M_full_2/max(M_full_2);
M_full_2 = M_full_2(1.5e6:1.5e6+length(t_out)-1); %!!the range must be changed accordly to the pulse position in the array...
%... the array dimention must match t_out dimention

%% Compare signal with IRF - useful to understand if the pulse broadening is detectable against setup's IRF

% normalize IRF
IRF_norm = IRF_data/max(IRF_data);
time_IRF = flip(time_IRF*1e12); %ps ascending order

t_out = t_out*1e12+abs(min(t_out*1e12)); % t_out in ps and first point at 0

hm_val = 0.5; %this indicates I want to find values at FWHM (i.e., 0.5 for normalized data)...
%... (if you want to find values at 10% hm_val= 0.1)

%FWHM signal+IRF
idx_above_hm = M_full_2 >= hm_val; %find array index for value > 0.5
% find first and last index above 0.5
first_idx = find(idx_above_hm, 1, 'first');
last_idx = find(idx_above_hm, 1, 'last');

%refine index for FWHM calculation
x1 = interp1(M_full_2([first_idx-1, first_idx]), t_out([first_idx-1, first_idx]), hm_val);
x2 = interp1(M_full_2([last_idx, last_idx+1]), t_out([last_idx, last_idx+1]), hm_val);
fwhm_value_signal = x2 - x1;

%FWHM IRF (same process as above)
idx_above_hm = IRF_norm >= hm_val;
first_idx = find(idx_above_hm, 1, 'first');
last_idx = find(idx_above_hm, 1, 'last');

x1 = interp1(IRF_norm([first_idx-1, first_idx]), time_IRF([first_idx-1, first_idx]), hm_val);
x2 = interp1(IRF_norm([last_idx, last_idx+1]), time_IRF([last_idx, last_idx+1]), hm_val);
fwhm_value_IRF = x2 - x1;

% difference between signal and IRF FWHM
delta_FWHM = fwhm_value_signal-fwhm_value_IRF;

threshold = 1.1; %empirical threshold to assess widening detection (e.g., 10% difference between...
%... measured signal and IRF's FWHM)

if (fwhm_value_signal > fwhm_value_IRF*threshold)
    disp('WG length can be reduced')
else
    disp('Cannot measure signal, broadening less than 10%')
end


% save FWHM as a function of WG's length
i = i+1;
output_FWHMdiff (i, 1) = L;
output_FWHMdiff (i, 2) = delta_FWHM %ps

% end %%uncomment if you want to cycle wg's L

%% Compare simulated and measured TCSPC curves

plot(t_out, M_full_2) %simulated
hold on
plot(t_eq+3900, I_measured) %measured %3900 is an empirical number to align the curves, which are based on...
%... two different time axis, change accordingly)

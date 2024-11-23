%% IDEAL - Kernel of metabolite separation method
% Input:
% - S = 5-dimensional array, containing image signals in the order:
% [columns, lines, slices, echoes, dynamic repetitions]
% - echo_times = echo times (TE)
% - delta_f = metabolite frequencies to reconstruct
% - it_max = max number of iterations
%
% Output:
% - metabolites = 5-dimensional array, containing metabolite intensities:
% [columns, lines, slices, metabolite, dynamic repetition]

%------------SETTINGS-------------

close all;
clear all;
clc;

path = "C:\Users\menze\Desktop\Matlab\IDEAL\ME bssfp data\2024_11_11\Reconstructed\32";

it_max = 10;

maxFmap = inf;
maxR2s = inf;
minR2s = -inf;
maxTheta = inf;
minT2s = -inf;

delta_f = [402.49 596.92]/1000; %kHz, fat water

load(path+"\par.mat");

TR = par.inf.rep_time; %ms 
echo_times = par.inf.te_time; %ms

%-------------END OF SETTINGS-------------

% swap echoes and repetitions
S=load(path+"\data.mat");
S=S.Data;
S = permute(S,[1 2 3 4 5 7 6]);
[Nx,Ny,Nz,Nr,N,Ne] = size(S);

% place echoes in 1st dimension
S = permute(reshape(S,[Nx*Ny*Nz*Nr Ne]), [2 1]);

% total number of pixels to calculate in all slices and repetetions
nPix = size(S,2);
warning off;

% number of metabolites to reconstruct
M = size(delta_f,2);

% construct matrix A for inversion
t = echo_times';
c = 1i*2*pi;
A = exp(c*t*delta_f);

% for bipolar correction
MOD = (1:Ne)';
MOD = (-1).^MOD;

% preallocations
% initial guess for off-resonance frequency
dFe = zeros(1,nPix) ;
ddFe = zeros(1,nPix);

% initial guess for 1/T2*
dR2e = ones(1,nPix)./minT2s;
ddR2e = zeros(1,nPix);

% initial guess for bipolar correction
dtheta = zeros(1,nPix);
ddtheta = zeros(1,nPix);

% IDEAL-LOOP
for iter = 1:it_max

    % correct signal for off-resonance and T2*, and bipolar correction
    S_c = S ./exp(c*t*dFe) ./exp(-abs(t)*dR2e) ./exp(MOD*dtheta);

    % estimated amplitudes
    AE = A\S_c;

    % estimated signal
    S_e = A*AE;

    % estimate new values for R2e, dfe and theta
    % based on signal in each voxel
    for k = 1:nPix
        % Taylor expansion matrix B of the signal equation
        B = [c*t.*S_e(:,k), -abs(t).*S_e(:,k), MOD.*S_e(:,k), A];
        % calculate the diff between corrected S_c and estimate S_e
        S_d(:,k) = (S_c(:,k)-S_e(:,k));
        % solve for taylor expansion terms
        dC = B\S_d(:,k);
        % use solution to correct dFe and R2e
        ddFe = real(dC(1));
        ddR2e = real(dC(2));
        ddtheta = dC(3);
        % adjust dFe and R2e for next iteration
        dFe(k) = max(min(dFe(k) + ddFe, maxFmap), -maxFmap);
        dR2e(k) = max(min(dR2e(k) + ddR2e, maxR2s),minR2s);
        dtheta(k) = max(min(real(dtheta(k)) + real(ddtheta), +maxTheta),-maxTheta) + 1i* max(min(imag(dtheta(k)) + imag(ddtheta),+maxTheta),-maxTheta);
    end

    % print current values
    fprintf('%4d%11.4f%11.4f%15.4f%15.4f%11.4f%11.4f\n', iter, mean(dFe), mean(ddFe), mean(dR2e), mean(ddR2e), mean(abs(dtheta)), mean(abs(ddtheta)));
end 

F_map = reshape(dFe, [Nx, Ny, Nz, Nr]);
R2s_map = reshape(dR2e, [Nx, Ny, Nz, Nr]);
Theta_map = reshape(dtheta, [Nx, Ny, Nz, Nr]);

%Final calculation steps:
% correction of signal with final correction values
S_c_final = S./exp(c*t*F_map(:)')./exp(-abs(t)*R2s_map(:)')./exp(MOD*Theta_map(:).');

% calculate final metabolite amplitudes ae
AE = A\S_c_final;

% calculate final estimate of signal
S_e_final = A*AE .*exp(c*t*F_map(:)').*exp(-abs(t)*R2s_map(:)').*exp(MOD*Theta_map(:).');

% calculate residual between S and Se for analysis
Residual = norm(S - S_e_final);

% create metabolite maps
metabolites = reshape(AE, [M,Nx,Ny,Nz,Nr]);
metabolites = permute(metabolites, [2,3,4,1,5]);

% turn on warning again
warning on;
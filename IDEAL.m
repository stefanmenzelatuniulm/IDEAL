close all;
clear all;
clc;

%------------SETTINGS-------------

dAlpha = 1/2; %alpha/2 spacing is dAlpha*TR
w = [402.49 596.92]; %kHz, fat water

nIt = 10^9;

path = "C:\Users\menze\Desktop\Matlab\IDEAL\ME bssfp data\2024_11_11\Reconstructed\32"; %Path to data
%TODO: Fabian mit Rekorechner um Hilfe bitten + Daten anschauen +
%vektorisieren mit kron/repmat/permute etc. sodass keine Schleife für
%Pixel, Slices, Averages nötig ist

T2sInitial = 50; %ms
phiInitial = 0;
omegaBInitial = 0;

%-------------END OF SETTINGS-------------

%Measured signal
load(path+"\data.mat");
S = squeeze(Data); %NS x NL x NE
clear Data;

load(path+"\par.mat");

TR = par.inf.rep_time; %ms 
TE = par.inf.te_time; %ms

NE = length(TE); %number of echoes
M = length(w); %number of metabolites
NS = length(S(:,1,1)); %number of samples per line in k-space
NL = length(S(1,:,1)); %number of lines in k-space

%Known chemical shifts
TEw = reshape(kron(w, TE), [NE M]);
C = exp(1i*TEw); %NE x M

%Write raw images at different echo times
SF = fftshift(fft2(S),1);
for k = 1:NE
    imwrite(rescale(abs(SF(:,:,k))),"TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
end

%This layout of S is better suited for the algorithm
S = permute(S, [3 1 2]); %NE x NS x NL

%IDEAL loop
for it = 1:nIt+1

    R = diag(exp(-abs(TE-dAlpha*TR)/T2sInitial)); 
    K = linspace(1,NE,NE);
    G = diag(exp(1i*(-1).^K*phiInitial));
    Q = diag(exp(1i*omegaBInitial*TE));

    %Encoding matrix
    E = R*G*Q*C; %NE x M
    
    %Estimate for metabolite signals
    rho = zeros(M,NS,NL);
    for k = 1:NL
        rho(:,:,k) = E\S(:,:,k); %A*x=B solved by A\B; dim of rho: M x NS x NL
    end
    
    %Signal that would have been measured with estimated parameters
    Se = permute(pagemtimes(repmat(E, [1 1 NL]),rho),[1 3 2]); %NE x NS x NL %TODO: Insert into algorithm in MA_TEX
    
    %Difference between measured signal and signal that would have been
    %measured with estimated parameters
    Sdiff = permute(S,[1 3 2])-Se; %NE x NS x NL
    
    %1st order Taylor expansion matrix of the signal equation
    %B = zeros(4,NE);
    %X = repmat(abs(TE-dAlpha*TR))
    %for k = 1:NL
    %B(1,:) = (1/T2sInitial^2)*X.*Se;
    B1 = permute((1/T2sInitial^2)*abs(permute(repmat(TE,[NL 1 NS]),[2 1 3])-dAlpha*TR).*Se,[1 3 2]);
    B2 = permute(permute(repmat(1i*((-1).^linspace(1,NE,NE)),[NL 1 NS]),[2 1 3]).*Se,[1 3 2]);
    B3 = permute(1i*permute(repmat(TE,[NL 1 NS]),[2 1 3]).*Se,[1 3 2]);
    B4 = permute(repmat(E,[1 1 NS NL]),[2 1 3 4]);
    %end
    %TODO: continue, not bad dimensions, estimates T2* MAP for all pixels,
    %not T* scalar
    dTaylor = B\Sdiff;
    
    T2sInitial = T2sInitial+dTaylor(1);
    phiInitial = phiInitial+dTaylor(2);
    omegaBInitial = omegaBInitial+dTaylor(3);
    rho = rho+dTaylor(4);

end

SFFT = fftshift(fft2(S),1);
rhoFFT = fftshift(fft2(rho),1);

%imshow(abs(SF(:,:,1,1,1,1,1)),[]);

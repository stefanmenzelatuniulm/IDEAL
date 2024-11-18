%------------SETTINGS-------------

close all;
clear all;
clc;

TR = 20; %ms
dAlpha = 1/2; %alpha/2 spacing is dAlpha*TR
TE = [1 2 3 4 5]; %ms
w = [0.316 -0.375 0.099]; %kHz

nIt = 10^9;

path = "C:\Users\menze\Desktop\Matlab\MR_Data\105"; %Path to data
%TODO: Fabian mit Rekorechner um Hilfe bitten + Daten anschauen +
%vektorisieren mit kron/repmat/permute etc. sodass keine Schleife für
%Pixel, Slices, Averages nötig ist

T2sInitial = 50; %ms
phiInitial = 0;
omegaBInitial = 0;

%-------------END OF SETTINGS-------------

%Measured signal
load(path+"\data.mat");
S = data;
clear data;

%Known chemical shifts
M = length(w);

C = zeros(NE, M);
for k = 1:NE
    for l = 1:M
        C(k,l) = exp(1i*w(l)*TE(k));
    end
end

%IDEAL loop
NE = length(TE);
for it = 1:nIt+1

    R = diag(exp(-abs(TE-dAlpha*TR)/T2sInitial)); 
    K = linspace(1,NE,NE);
    G = diag(exp(1i*(-1).^K*phiInitial));
    Q = diag(1i*omegaBInitial*TE);

    %Encoding matrix
    E = R*G*Q*C;
    
    %Estimate for metabolite signals
    rho = E\S; %A*x=B solved by A\B
    
    %Signal that would have been measured with estimated parameters
    Se = E*rho;
    
    %Difference between measured signal and signal that would have been
    %measured with estimated parameters
    Sdiff = S-Se;
    
    %1st order Taylor expansion matrix of the signal equation
    B = zeros(4,NE);
    %X = repmat(abs(TE-dAlpha*TR))
    for k = 1:NE
        %B(1,:) = (1/T2sInitial^2)*X.*Se;
        B(1,k) = (1/T2sInitial^2)*abs(TE(k)-dAlpha*TR).*Se(k,:);
        B(2,k) = 1i*(-1)^k*Se(k,:);
        B(3,k) = 1i*TE(k)*Se(k,:);
        B(4,k) = rho;
    end
    
    dTaylor = B\Sdiff;
    
    T2sInitial = T2sInitial+dTaylor(1);
    phiInitial = phiInitial+dTaylor(2);
    omegaBInitial = omegaBInitial+dTaylor(3);
    rho = rho+dTaylor(4);

end
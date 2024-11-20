close all;
clear all;
clc;

%------------SETTINGS-------------

dAlpha = 1/2; %alpha/2 spacing is dAlpha*TR
w = [402.49 596.92]/1000; %kHz, fat water

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
S_ = permute(squeeze(Data),[3 1 2]); %NE x NS x NL
clear Data;

load(path+"\par.mat");

TR = par.inf.rep_time; %ms 
TE = par.inf.te_time; %ms

NE = length(TE); %number of echoes
M = length(w); %number of metabolites
NS = length(S_(1,:,1)); %number of samples per line in k-space
NL = length(S_(1,1,:)); %number of lines in k-space

%Write raw images at different echo times
for k = 1:NE
    %imwrite(permute(rescale(abs(fftshift(fft2(S_(k,:,:))))),[3 2 1]),"Imag_TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
    imwrite(permute(rescale(abs(S_(k,:,:))),[3 2 1]),"kSpace_TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
    SF = permute(fftshift(fft2(S_(k,:,:))),[3 2 1]);
    imwrite(rescale(abs(SF)),"Imag_TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
end

T2sMap = zeros(NS, NL);
PhiMap = zeros(NS, NL);
OmegaBMap = zeros(NS, NL);
rhoSeparated = zeros(M, NS, NL);

for pixX = 1:NS

    for pixY = 1:NL
        
        %For every voxel

        S = S_(:,pixX,pixY); %NE

        %Known chemical shifts
        TEw = reshape(kron(w, TE), [NE M]);
        C = exp(1i*TEw); %NE x M
        
        %IDEAL loop
        for it = 1:nIt+1
        
            R = diag(exp(-abs(TE-dAlpha*TR)/T2sInitial)); %NE x NE
            K = linspace(1,NE,NE); %NE
            G = diag(exp(1i*(-1).^K*phiInitial)); %NE x NE
            Q = diag(exp(1i*omegaBInitial*TE)); %NE x NE
        
            %Encoding matrix
            E = R*G*Q*C; %NE x M
            
            %Estimate for metabolite signals
            %A*x=B solved by A\B or pagemldivide(A,B) for 3 D matrices
            rho = E\S; %M
            
            %Signal that would have been measured with estimated parameters
            Se = E*rho; %NE
            
            %Difference between measured signal and signal that would have been
            %measured with estimated parameters
            Sdiff = S-Se; %NE
            
            %1st order Taylor expansion matrix of the signal equation
            B1 = (1/T2sInitial^2)*abs(TE'-dAlpha*TR).*Se; %NE
            B2 = 1i*((-1).^linspace(1,NE,NE)').*Se; %NE
            B3 = 1i*TE'.*Se; %NE
            B4 = E; %NE x M
        
            B=[B1 B2 B3 B4]; %NE x (M + 3)

            dT = B\Sdiff; %M+3

            dT2s = dT(1); 
            dPhi = dT(2);
            dOmegaB = dT(3);
            dRho = dT(4:end);
            
            T2sInitial = T2sInitial+dT2s; %scalar
            phiInitial = phiInitial+dPhi; %scalar
            omegaBInitial = omegaBInitial+dOmegaB; %scalar
            rho = rho+dRho; %M
        
        end

        T2sMap(pixX, pixY) = T2sInitial;
        PhiMap(pixX, pixY) = phiInitial;
        OmegaBMap(pixX, pixY) = omegaBInitial;

        for met = 1:M
            rhoSeparated(met, pixX, pixY) = rho(met);
        end

    end

end

imwrite(rescale(abs(T2sMap)),"T2sMap.png");
imwrite(rescale(abs(PhiMap)),"PhiMap.png");
imwrite(rescale(abs(OmegaBMap)),"OmegaBMap.png");
imwrite(rescale(abs(OmegaBMap)),"OmegaBMap.png");
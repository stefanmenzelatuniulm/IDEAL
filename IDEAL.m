close all;
clear all;
clc;

%------------SETTINGS-------------

dAlpha = 1/2; %alpha/2 spacing is dAlpha*TR
w = [402.49 596.92]/1000; %kHz, fat water

nIt = 10;

path = "C:\Users\menze\Desktop\Matlab\IDEAL\ME bssfp data\2024_11_11\Reconstructed\32"; %Path to data

%Initial values
RT2s = 1/50; %1/ms
phi = 0;
omegaB = 0;

%-------------END OF SETTINGS-------------

%State of the Art Model

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
    img = permute(S_(k,:,:),[3 2 1]);
    imwrite(rescale(abs(img)),"kSpace_TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
    imwrite(rescale(abs(fftshift(fft2(img),2))),"Imag_TE_"+strrep(num2str(TE(k)),".","_")+"_ms.png");
end

RT2sMap = zeros(NS, NL);
PhiMap = zeros(NS, NL);
OmegaBMap = zeros(NS, NL);
rhoSeparated = zeros(M, NS, NL);

for pixY = 1:NL

    for pixX = 1:NS

        %For every voxel

        S = S_(:,pixX,pixY); %NE

        %Known chemical shifts
        TEw = reshape(kron(w, TE), [NE M]);
        C = exp(1i*TEw); %NE x M

        RT2sInitial = RT2s;
        phiInitial = phi;
        omegaBInitial = omegaB;
        
        %IDEAL loop
        for it = 1:nIt
            
            %if mod(it, 1000) == 0 
            disp("Pixel "+num2str((pixY-1)*NL+pixX)+"/"+num2str(NS*NL)+", iteration "+num2str(it)+"/"+num2str(nIt));
            %end

            R = diag(exp(-RT2sInitial*abs(TE-dAlpha*TR))); %NE x NE
            %R = diag(exp(-RT2sInitial*abs(TE-TR))+exp(-RT2sInitial*abs(TE))); %NE x NE
            K = linspace(1,NE,NE); %NE
            G = diag(exp(1i*(-1).^K*phiInitial)); %NE x NE
            Q = diag(exp(1i*omegaBInitial*TE)); %NE x NE
        
            %Encoding matrix
            E = R*G*Q*C; %NE x M
            
            %Estimate for metabolite signals
            %A*x=B solved by A\B or pagemldivide(A,B) for 3 D matrices
            lastwarn("");
            rho = E\S; %M
            msg = lastwarn;
            if msg ~= ""
                a = 1;
            end
            
            %Signal that would have been measured with estimated parameters
            Se = E*rho; %NE
            
            %Difference between measured signal and signal that would have been
            %measured with estimated parameters
            Sdiff = S-Se; %NE
            
            %1st order Taylor expansion matrix of the signal equation
            B1 = -abs(TE'-dAlpha*TR).*Se; %NE
            B2 = 1i*((-1).^linspace(1,NE,NE)').*Se; %NE
            B3 = 1i*TE'.*Se; %NE
            B4 = E; %NE x M
        
            B=[B1 B2 B3 B4]; %NE x (M + 3), [pSE/pRT2sInitial_TE1 pSE/pPhi_TE1 pSE/pOmegaB_TE1 pSE/prho_M1_TE1 ... pSE/prho_MM_TE1; ... (for all TE in the next rows)] 

            lastwarn("");
            dT = B\Sdiff; %M+3
            msg = lastwarn;
            if msg ~= ""
                a = 1;
            end

            %if it == nIt
            disp(num2str(norm(abs(dT))));
            %end

            dRT2s = real(dT(1)); 
            dPhi = real(dT(2));
            dOmegaB = real(dT(3));
            dRho = dT(4:end);
           
            continueFlag = false;

            RT2sInitial = RT2sInitial+dRT2s; %scalar            
            % if 1/real(RT2sInitial) > 1000
            %     RT2sInitial = 1/50;
            %     continueFlag = true;
            % end

            phiInitial = phiInitial+dPhi; %scalar
            % if real(phiInitial) > pi       
            %     phiInitial = 0;
            %     continueFlag = true;
            % end

            omegaBInitial = omegaBInitial+dOmegaB; %scalar
            % if real(omegaBInitial) > max(w)
            %     omegaBInitial = 0;
            %     continueFlag = true;
            % end

            % if continueFlag
            %     continue;
            % end

            rho = rho+dRho; %M
        
        end

        RT2sMap(pixX, pixY) = RT2sInitial;
        PhiMap(pixX, pixY) = phiInitial;
        OmegaBMap(pixX, pixY) = omegaBInitial;

        for met = 1:M
            rhoSeparated(met, pixX, pixY) = rho(met);
        end

    end

end

RT2sMap(isnan(RT2sMap))=0;
PhiMap(isnan(PhiMap))=0;
OmegaBMap(isnan(OmegaBMap))=0;
imwrite(rescale(RT2sMap),"RT2sMap_SoA.png");
imwrite(rescale(PhiMap),"PhiMap_SoA.png");
imwrite(rescale(OmegaBMap),"OmegaBMap_SoA.png");

for met = 1:M
    img = permute(rhoSeparated(met,:,:),[3 2 1]);
    img(isnan(img))=0;
    imwrite(rescale(abs(img)),"kSpace_met_"+strrep(num2str(met),".","_")+".png");
    imwrite(rescale(abs(fftshift(fft2(img),2))),"Imag_met_"+strrep(num2str(met),".","_")+".png");
end
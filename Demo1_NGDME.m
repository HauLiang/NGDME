% -------------- The superiority of exploiting sparsity -------------------
%
% This is a simple example to test the Nonlinear Group Delay Mode Estimation (NGDME) algorithm 
%
% Author: Hao Liang and Xiaotong Tu
%
% Last modified by: 21/08/23
%

clc; clear; close all

fs = 100;   % sample frequency
T = 15;     % time duration
Nt = 1500;  % the number of samples in time-domain
Nf = floor(Nt/2)+1; % the number of samples in frequency-domain
f = (0 : Nf-1)/T;   % frequency variables
t = (0 : Nt-1)/fs;  % time variables

% Group delay of the signal modes
gd1 = 0.1*f + 1;
gd2 = 0.2*f + 5;

% Amplitude of the signal modes
a1 = sawtooth(2*pi*f,0.8)+1;
a2 = (1 + 0.2*sin(2*pi*2/50*f)); 

% 'NGDm1' denotes a mode of nonlinear group delay signal, 'iNGDFs1' denotes 
% its bilateral spectrums, and 'ifftSig1' denotes its time-domain representation
NGDm1 = a1.*exp(-1j*2*pi*(0.1/2*f.^2 + 1*f + 0.1));  iNGDFs1 = [NGDm1,conj(fliplr(NGDm1(2:ceil(Nt/2))))];  ifftSig1 = ifft(iNGDFs1);
NGDm2 = a2.*exp(-1j*2*pi*(0.2/2*f.^2 + 5*f + 0.8));  iNGDFs2 = [NGDm2,conj(fliplr(NGDm2(2:ceil(Nt/2))))];  ifftSig2 = ifft(iNGDFs2);

% time-domain signal
Sig = real(ifftSig1 + ifftSig2);

% frequency-domain signal
NGDmodes = NGDm1 + NGDm2; 


%% NGDME
gamma = 1e7; lambda = 5e-3; tol = 5e-5;
iniGD = [3*ones(1,length(f));10*ones(1,length(f))];

tic;
[eGDest, Desest] = NGDME(NGDmodes,T,iniGD,lambda,gamma,tol);
toc;

NGDMEeD1 = Desest(1,:,end); NGDMEeD2 = Desest(2,:,end); 

% Obtain the time-domain signal by inverse FFT
ieADFs1 = [NGDMEeD1,conj(fliplr(NGDMEeD1(2:ceil(Nt/2))))]; iffteASig1 = ifft(ieADFs1); 
ieADFs2 = [NGDMEeD2,conj(fliplr(NGDMEeD2(2:ceil(Nt/2))))]; iffteASig2 = ifft(ieADFs2);

%% Reconstructed amplitudes
% Amplitude 1
figure; x1 = 8.5; y1 = -0.2; x2 = 9.5; y2 = 2.2;
plot(f,abs(NGDm1),'b','linewidth',1);
hold on;
plot(f,abs(NGDMEeD1),'r','linewidth',1);
set(gcf,'Position',[20 100 900 350]);	
set(gcf,'Color','w'); 
xlabel('Frequency/Hz','FontName','Times New Roman');
ylabel('Amplitude','FontName','Times New Roman');
ylim([-1.5 3.5]);
set(gca,'FontSize',24)
set(gca,'linewidth',2);
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);
legend('Original','Estimated','FontSize',17,'FontName','Times New Roman');
h1 = axes('position',[0.6 0.3 0.3 0.3]);
axis(h1);
plot(f,abs(NGDm1),'b','linewidth',1.5);
hold on;
plot(f,abs(NGDMEeD1),'r','linewidth',1.5);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12,'linewidth',1)

% Amplitude 2
figure;
set(gcf,'Position',[20 100 900 350]);	    
set(gcf,'Color','w'); 
plot(f,abs(NGDm2),'b','linewidth',1);   % time-domain signal modes
hold on
plot(f,abs(NGDMEeD2),'r','linewidth',1);
xlabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
legend('Original','Estimated','FontSize',17,'FontName','Times New Roman');
ylim([0.7 1.3])


%% Estimated group delay by NGDME
figure
plot([gd1;gd2],f,'b','linewidth',3) % true group delay
hold on
plot(eGDest(:,:,end),f,'r--','linewidth',3) % estimated group delay
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
ylabel('Time/s','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	

%% EQF results
Amp1_EQF_NGDME = 20*log10(norm(abs(NGDm1) - abs(NGDMEeD1),2)/norm(abs(NGDm1),2))
Amp2_EQF_NGDME = 20*log10(norm(abs(NGDm2) - abs(NGDMEeD2),2)/norm(abs(NGDm2),2))

GD1_EQF_NGDME = 20*log10(norm(eGDest(1,:,end) - gd1,2)/norm(gd1,2))
GD2_EQF_NGDME = 20*log10(norm(eGDest(2,:,end) - gd2,2)/norm(gd2,2))

Mode1_EQF_NGDME = 20*log10(norm(real(ifftSig1) - real(iffteASig1),2)/norm(real(ifftSig1),2))
Mode2_EQF_NGDME = 20*log10(norm(real(ifftSig2) - real(iffteASig2),2)/norm(real(ifftSig2),2))



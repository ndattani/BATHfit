%% THIS CODE IS PROTECTED BY COPYRIGHT LAW.
%% The owners of this code are Nike S. Dattani and David Mark Wilkins
%% The copyright is currently owned by Nike S. Dattani.

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT. LEGAL CONSEQUENCES WILL APPLY 
function leastSquaresFitToBathResponseFunction(filenameContainingSpectralDensityParameters,T,numberOfInitialTrialExponentialTerms,maxNumberOfFittedExponentialTerms)
% '2011Olbrich_JPCL_2_1771_tableS1_BChl1.mat' -> input file, given over here since matlab's ls command isn't giving the full filename ! 
%% 1 Conversion factors
joules2electronVolts=6.24150974*10^18;
fs=1e15;
%% 2 Fundamental constants
kb=1.3806504*10^(-23); % Joules / Kelvin
kb=kb*joules2electronVolts; % eV / Kelvin

hbar=1.054571628*10^(-34); % Joules*seconds
hbar=4.135667516*10^(-15); % eV*seconds
hbar=4.135667516*10^(-15)*fs; % eV*femtoseconds
%% 3 Physical parameters
beta=1/(kb*T); % in Kelvins/eV
spectralDensityParameters=flipud(struct2array(load(filenameContainingSpectralDensityParameters))); % filename is different from variable name because matlab doesn't let us begin a variable name with a number .  The flipud is because the high frequency modes (the ones we don't want) seem to be given FIRST !
eta=spectralDensityParameters(:,1)*1e-5; % in eV^2
omegaTilde=2*pi./spectralDensityParameters(:,2); % in fs
gamma=1./spectralDensityParameters(:,3); % in fs
%% 4 Calculate the spectral density
w=0:0.001:2.5;w=w.'/hbar;J=zeros(length(w),length(eta)); % in fs-1

for ii=1:length(w);
J(ii,:)=(2/(pi*hbar))*tanh(beta*hbar*w(ii)/2).*cumsum((eta.*gamma./(2*(gamma.^2+(w(ii)-omegaTilde).^2))+(eta.*gamma./(2*(gamma.^2+(w(ii)+omegaTilde).^2))))); % Equation 4 of 2011 Olbrich et al. JPCL ,2, 1771-1776. 
end

plot(hbar*w,J(:,length(eta)),'b','LineWidth',2);
%% 5 Labels
title('BChl 1, convert $\omega$ from fs$^{-1}$ to cm for correct units ?','FontSize',32,'Interpreter','latex')
xlabel('$\hbar\omega$ [eV]','FontSize',32,'Interpreter','latex');ylabel('$J(\omega)$ [eV]','FontSize',32,'Interpreter','latex')
%% 6 Plot spectral density with varying number of Lorentzian terms
jetVariable=jet;jetVariable=jetVariable(find(mod(1:length(jetVariable),2)),:);set(gca,'ColorOrder',jetVariable(find(mod(1:length(jetVariable),2)),:));hold('all'); % varycolor.m on the FEX is the ideal way to plot many lines, each a different color. Instead we've used jet(find(mod(1:length(jet),2)),:), which reduces the number of colors in jet to 32 by removing every other row.
plot(hbar*w,J,'LineWidth',2);
legendHandle=legend(arrayfun(@(i) num2str(i), 1:length(eta), 'Uniform', 0));
set(get(legendHandle,'title'),'string','Number of terms');
set(legendHandle,'Interpreter','latex','FontSize',16,'LineWidth',3)
%% 7 Bath response function
tMesh=0.1;finalTime=1000;t=0:tMesh:finalTime; alpha=t;% in fs
numberOfLorentzianTermsInSpectralDensity=length(eta);
p=expand(eta(1:numberOfLorentzianTermsInSpectralDensity)/2,[2,1]);Omega=p; % initialize omega's length to be that of p
Omega(1:2:end)=1i*omegaTilde(1:numberOfLorentzianTermsInSpectralDensity)-gamma(1:numberOfLorentzianTermsInSpectralDensity);
Omega(2:2:end)=conj(1i*omegaTilde(1:numberOfLorentzianTermsInSpectralDensity)-gamma(1:numberOfLorentzianTermsInSpectralDensity));

for ii=1:length(t)
alpha(ii)=complex(sum(p.*exp(Omega*t(ii))),0);
end; alpha=alpha.'; %alpha needs to be a row for the least-squares fitting program
figure(2);plot(t,alpha);
%% 8 non-linear Least Squares Fit
clearvars -except alpha tMesh finalTime t J w eta omegaTilde gamma p Omega  numberOfInitialTrialExponentialTerms maxNumberOfFittedExponentialTerms % clearing these variables at the beginning of the cell allows us to run the cell many times with different fitting characteristics. It may be appropriate to keep more of the useful variables that don't affect the fit.
weights=ones(size(alpha));weights(find(t==400):find(t==401))=2.5;weights([find(t==10):find(t==40) find(t==30):find(t==400)])=1.5; % If all weights are 1, we have a non-weighted least-squares fit. 
%options=optimset('Display','iter','TolFun',1e-20,'Algorithm','levenberg-marquardt','TolX',1e-18,'DiffMinChange',1e-7,'FinDiffType','central','MaxFunEvals',4000);
options=optimset('Display','iter','TolFun',1e-20,'Algorithm',[],'TolX',1e-18,'DiffMinChange',1e-7,'FinDiffType','central','MaxFunEvals',50000,'MaxIter',1000);

for ii=numberOfInitialTrialExponentialTerms:maxNumberOfFittedExponentialTerms
switch ii
    case numberOfInitialTrialExponentialTerms;trialParameters=[p(1:numberOfInitialTrialExponentialTerms) ; Omega(1:numberOfInitialTrialExponentialTerms)].'; % complex initial conditions
    otherwise trialParameters=[fittedParameters(1:length(fittedParameters)/2) randn*(p(1)) fittedParameters(length(fittedParameters)/2+1:length(trialParameters)) randn*(Omega(1))];  
end
% trialParameters=rand*fittedParameters; % sometimes the parameters need a bit of a nudge to escape a local minimum ?

x0=[real(trialParameters) imag(trialParameters)]; % Convert initial conditions to a real part, then Concatenate the imaginary part to the real part. I think this is needed because lsqnonlin treats the real and imaginary parts as separate fitting parameters, and all fitting parameters need to be in a row.
%% 8.1 Main fitting computation !
%x=lsqnonlin(@(x)fit2(x,t,alpha),x0,[],[],options); %t is a row and alpha's a column
x=lsqnonlin(@(x)fitComplex(x,t,alpha,weights),x0,[],[],options); %t is a row and alpha's a column
fittedParameters = complex(x(1:length(trialParameters)),x(length(trialParameters)+1:2*length(trialParameters))); % Convert answer to complex. x and x0 are rows
%% 8.2 Plot the fitted result (before plotting the original points that were being fitted to, so that originals aren't covered)
alphaFitted=exp(repmat(fittedParameters(length(fittedParameters)/2+1:length(fittedParameters)),length(t),1).*repmat(t.',1,length(fittedParameters)/2))*fittedParameters(1:length(fittedParameters)/2).';
figure(ii);hold('on');
plot(t,real(alphaFitted),'r','LineWidth',4);plot(t,imag(alphaFitted),'r','LineWidth',4);
axis([0,finalTime,min([real(alpha) ; imag(alpha)]),max([real(alpha) ; imag(alpha)])])
%% 8.3 Plot the results on a finer mesh (if points are sparse). Why not always do this ? Because the points in between the fitted points weren't constrained by the fit and might not match well.
% splineMeshSize=tMesh/1000;splineMesh=0:splineMeshSize:finalTime;% 
% alphaFittedOnSplineMesh=exp(repmat(fittedParameters(length(fittedParameters)/2+1:length(fittedParameters)),length(splineMesh),1).*repmat(splineMesh.',1,length(fittedParameters)/2))*fittedParameters(1:length(fittedParameters)/2).';
% plot(splineMesh,real(alphaFittedOnSplineMesh),'m');plot(splineMesh,imag(alphaFittedOnSplineMesh),'m');
%% 8.4 Plot the spline of the results (if the above two don't look good)
% alphaFittedSplined=complex(spline(t,real(alphaFitted),splineMesh),spline(t,imag(alphaFitted),splineMesh));
% plot(splineMesh,real(alphaFittedSplined),'g','LineWidth',4);plot(splineMesh,spline(alphaFittedSplined),'g','LineWidth',4);
%% 8.5 Plot the original points that were fitted to. Do this at the end so that the points are seen.
plot(t,real(alpha),'.');plot(t,imag(alpha),'.') ;hold('off')
end % loop over number of expoenntial terms fitted
end % function
%% THIS CODE IS PROTECTED BY COPYRIGHT LAW.
%% The owners of this code are Nike S. Dattani and David M. Wilkins
%% The copyright is currently owned by Nike S. Dattani.

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT. LEGAL CONSEQUENCES WILL APPLY 
function fittedParametersInSIunits=leastSquaresFitToBathResponseFunction(t,T,w,J,numberOfInitialTrialExponentialTerms,maxNumberOfFittedExponentialTerms,figureTitle)
%% 1 Fundamental constants
kb=1.3806504*10^(-23); % Joules / Kelvin
hbar=1.054571628*10^(-34); % Joules*seconds
c=299792458; % meters / second
%% 2 Conversion factors
joules2electronVolts=6.24150974*10^18;
joules2wavenumbers=1/(hbar*2*pi*c*100); % wavenumbers2joules is h*c , were c is in cm/s
pico=1e12;femto=1e15;
%% 3 Physical parameters
beta=1/(kb*T); % dimensionless
%% 4 Obtain the spectral density and Bath Response Function
[J alpha trialParameters]=boseEinsteinBathResponseFunctionByNumericalIntegration(hbar,t,beta,w,J); 
scaledAlpha=alpha.'/max(real(alpha));scaledTime=t/max(t); %  alpha needs to be a column for the least-squares fitting program % the fitting is easier when the ordinate and abscissa variables are the same size (I think this is so that diffMinChange and diffMaxChange can be the same for all fitting parameters)
%% 5 non-linear Least Squares Fit
clearvars -except beta hbar alpha scaledAlpha tMesh finalTime t scaledTime J w figureTitle eta omegaTilde gamma p Omega  trialParameters numberOfInitialTrialExponentialTerms maxNumberOfFittedExponentialTerms % clearing these variables at the beginning of the cell allows us to run the cell many times with different fitting characteristics. It may be appropriate to keep more of the useful variables that don't affect the fit.
weights=ones(size(scaledAlpha));%weights(find(t==400):find(t==600))=2.5;weights([find(t==10):find(t==40) find(t==30):find(t==400)])=1.5; % If all weights are 1, we have a non-weighted least-squares fit. 
%options=optimset('Display','iter','TolFun',1e-20,'Algorithm','levenberg-marquardt','TolX',1e-18,'DiffMinChange',1e-7,'FinDiffType','central','MaxFunEvals',4000);
options=optimset('Display','iter','TolFun',1e-20,'Algorithm',[],'TolX',1e-18,'DiffMinChange',1e-7,'FinDiffType','central','MaxFunEvals',50000,'MaxIter',1000);

for ii=numberOfInitialTrialExponentialTerms:maxNumberOfFittedExponentialTerms
if ii>numberOfInitialTrialExponentialTerms;trialParameters=[fittedParametersScaled(1:length(fittedParametersScaled)/2) randn*(trialParameters(1)) fittedParametersScaled(length(fittedParametersScaled)/2+1:length(trialParameters)) randn*(trialParameters(end))];end  
% trialParameters=rand*fittedParameters; % sometimes the parameters need a bit of a nudge to escape a local minimum ?

x0=[real(trialParameters) imag(trialParameters)]; % Convert initial conditions to a real part, then Concatenate the imaginary part to the real part. I think this is needed because lsqnonlin treats the real and imaginary parts as separate fitting parameters, and all fitting parameters need to be in a row.
%% 5.1 Main fitting computation !
%x=lsqnonlin(@(x)fit2(x,t,alpha),x0,[],[],options); %t is a row and alpha's a column
x=lsqnonlin(@(x)fitComplex(x,scaledTime,scaledAlpha,weights),x0,[],[],options); %t is a row and alpha's a column
fittedParametersScaled = complex(x(1:length(trialParameters)),x(length(trialParameters)+1:2*length(trialParameters))); % Convert answer to complex. x and x0 are rows
fittedParametersScaled(1:length(fittedParametersScaled)/2)=fittedParametersScaled(1:length(fittedParametersScaled)/2)*max(real(alpha)); %scale the p's so that the ordinate values of alpha are back in SI units
fittedParametersInSIunits=fittedParametersScaled;
fittedParametersInSIunits(length(fittedParametersScaled)/2+1:end)=fittedParametersScaled(length(fittedParametersScaled)/2+1:end)/max(t); % scale the Omega's so that the abscissa values of alpha are back in SI units
alphaFitted=exp(repmat(fittedParametersInSIunits(length(fittedParametersInSIunits)/2+1:length(fittedParametersInSIunits)),length(t),1).*repmat(t.',1,length(fittedParametersInSIunits)/2))*fittedParametersInSIunits(1:length(fittedParametersInSIunits)/2).';
%% 5.2 Plot the fitted result for alpha (before plotting the original points that were being fitted to, so that originals aren't covered)
figure(ii);subplot(2,1,1);hold('on');title(figureTitle,'FontSize',32,'Interpreter','latex'); % would be nice if we could make a title for the whole thing though, not just the subplot
plotHandle(1)=plot(t,real(alphaFitted),'r','LineWidth',4);plot(t,imag(alphaFitted),'r','LineWidth',4);
%% 5.3 Plot the results on a finer mesh (if points are sparse). Why not always do this ? Because the points in between the fitted points weren't constrained by the fit and might not match well.
% % splineMeshSize=tMesh/1000;splineMesh=0:splineMeshSize:finalTime;% 
% % alphaFittedOnSplineMesh=exp(repmat(fittedParametersInSIunits(length(fittedParametersInSIunits)/2+1:length(fittedParametersInSIunits)),length(splineMesh),1).*repmat(splineMesh.',1,length(fittedParametersInSIunits)/2))*fittedParametersInSIunits(1:length(fittedParametersInSIunits)/2).';
% % plot(splineMesh,real(alphaFittedOnSplineMesh),'m');plot(splineMesh,imag(alphaFittedOnSplineMesh),'m');
%% 5.4 Plot the spline of the results (if the above two don't look good)
% % alphaFittedSplined=complex(spline(t,real(alphaFitted),splineMesh),spline(t,imag(alphaFitted),splineMesh));
% % plot(splineMesh,real(alphaFittedSplined),'g','LineWidth',4);plot(splineMesh,spline(alphaFittedSplined),'g','LineWidth',4);
%% 5.5 Plot the original points that were fitted to. Do this at the end so that the points are seen.
plotHandle(2)=plot(t,real(alpha),'.');plot(t,imag(alpha),'.');hold('off')
%% 5.6 Axes and labels
axis('tight');box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',3,'FontSize',16);
%yLabelHandle=get(gca,'YLabel');set(yLabelHandle,'Position',get(yLabelHandle,'Position') - [0.000125 0 0]);
%xLabelHandle=get(gca,'XLabel');set(xLabelHandle,'Position',get(xLabelHandle,'Position') - [0 0.000125 0]);
ylabel('$\alpha(t)$ [Joules $\cdot$ seconds]','FontSize',24,'Interpreter','latex')
xlabel('Time [seconds]','FontSize',32,'Interpreter','latex')
legendHandle=legend(plotHandle,{'$\alpha(t)$ fitted' '$\alpha(t)$ original'});set(legendHandle,'Interpreter','latex','FontSize',16,'LineWidth',3,'Position',[0.747672758188061 0.753355153875044 0.134706814580032 0.139318885448916])
%% 5.7 Plot J(w) corresponding to fitted alpha(t) and original J(w)
subplot(2,1,2);hold('on');
Jfitted=(1/pi)*(1-exp(-beta*hbar*w.')).*real(sum(bsxfun(@rdivide,1i*fittedParametersInSIunits(1:length(fittedParametersInSIunits)/2).',bsxfun(@minus,w.',1i*fittedParametersInSIunits(length(fittedParametersInSIunits)/2+1:length(fittedParametersInSIunits)).'))));
plotHandle2(1)=plot(w,Jfitted,'r','LineWidth',4); 
plotHandle2(2)=plot(w,J,'b','LineWidth',4); 
axis('tight');box('on');
set(gca,'XMinorTick','on','YMinorTick','on','LineWidth',3,'FontSize',16);
ylabel('$J(\omega)$ [Joules]','FontSize',32,'Interpreter','latex')
xlabel('$\omega$ [seconds$^{-1}$]','FontSize',32,'Interpreter','latex')
legendHandle2=legend(plotHandle2,{'$J(\omega)$ fitted' '$J(\omega)$ original'});set(legendHandle2,'Interpreter','latex','FontSize',16,'LineWidth',3,'Location','NorthEast')
end % loop over number of expoenntial terms fitted

function realAndImaginaryPartsOfErrorInOneColumnVector = fitComplex(fittedParametersInOneColumnVector,X,Y,weights)
n=length(fittedParametersInOneColumnVector);
fittedParameters=complex(fittedParametersInOneColumnVector(1:n/2),fittedParametersInOneColumnVector(n/2+1:n)); % Remember took the complex trial parameters and made a column vector composed half of the real parts of the trial parameters, and half of the imag parts of the trial parameters. I think this is needed because lsqnonlin treats the real and imaginary parts as separate fitting parameters, and all fitting parameters need to be in a row.
%% Evaluate the error (difference between the fitted function and the desired function to which we're fitting)
error= exp(repmat(fittedParameters((n/4)+1:n/2),length(X),1).*repmat(X.',1,(n/4)))*(fittedParameters(1:n/4).')- Y; %[exp(omega1)AtTimeT1 exp(omega2)AtTimeT1 ... ; exp(omega1)AtTimeT2 exp(omega2)AtTimeT2]*[p1 ; p2 ; ..]
%errsor = (1/pi)*(1-exp(-beta*X.')).*real((1./(repmat(-1i*fittedParameters((n/4)+1:n/2)*3.33563759345165*10^-11,length(X),1)+repmat(X.',1,(n/4))))*(1i*fittedParameters(1:n/4).'))- Y;
%%
error=error.*weights;
realAndImaginaryPartsOfErrorInOneColumnVector=[real(error);imag(error)]; % Convert the answer back into a column vector such that the first half is the real part, and the second half is the imaginary part. 

%%%%%%%%%%%%%%%%%%%  Five Different types of spectral densities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J alpha trialParameters]=boseEinsteinBathResponseFunctionByNumericalIntegration(hbar,t,beta,w,J)
%alpha=trapz(w,repmat(J.*(coth(beta*hbar*w/2)),length(t),1).*cos(repmat(w/(2*pi),length(t),1).*repmat(t',1,length(w)))-1i*repmat(J,length(t),1).*sin(repmat(w/(2*pi),length(t),1).*repmat(t',1,length(w))),2); %very memory demanding and slow. It's also wrong because the sin is never multiplied by J !
alpha=zeros(size(t));for ii=1:length(t);alpha(ii)=trapz(w,J.*(coth(beta*hbar*w/2).*cos(w*t(ii)/(2*pi))-1i*sin(w*t(ii)/(2*pi))));end
trialParameters=[complex(0.6563,1.004) complex(1.179,6.772) complex(6.257,-7.812) complex(-0.518,49.15) complex(- 4.666,43.84) complex(-81.81,0.4144)];

function LDD(hbar,t,beta,w,wTilde,lambdaVector,gamma,numberOfExponentialTermsInAlpha,typeOfSeries) % we need to call it lambdaVector because lambda is the reorganization energy
function Kleinekathoefer(hbar,t,beta,w,wTilde,lambdaVector,gamma,numberOfExponentialTermsInAlpha,typeOfSeries)
function MeierTannor(hbar,t,beta,w,wTilde,lambdaVector,gamma,numberOfExponentialTermsInAlpha,typeOfSeries)
pi*w/2;
function [J alpha trialParameters]=Leggett(hbar,t,beta,w,s,gimel,wc,q);zt=complex(1,wc*t)/(beta*hbar*wc); % Hebrew gimel looks very much like lambda, but won't confuse you with th reorganization energy
J=gimel*w.^s.*exp(-(w/wc).^q);
alpha=gimel*complex(real((-1/(beta*hbar))^(s+1)*(psin(s,zt)+psin(s,1+zt))),imag(gamma(s+1)./(beta*hbar*zt).^(s+1)));
trialParameters=[complex(alpha(1),0) ; complex(-1/(wc*q),1/(wc*q))];

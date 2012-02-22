%% THIS CODE IS PROTECTED BY COPYRIGHT LAW.
%% The owners of this code are Nike S. Dattani and David Mark Wilkins
%% The copyright is currently owned by Nike S. Dattani.

%% NO PART OF THIS CODE IS TO BE MODIFIED OUTSIDE OF GIT. LEGAL
%% CONSEQUENCES WILL APPLY 
function realAndImaginaryPartsOfErrorInOneColumnVector = fitComplex(fittedParametersInOneColumnVector,X,Y,weights)
global beta; 
n=length(fittedParametersInOneColumnVector);
fittedParameters=complex(fittedParameters(1:n/2),fittedParametersInOneColumnVector(n/2+1:n)); % Remember took the complex trial parameters and made a column vector composed half of the real parts of the trial parameters, and half of the imag parts of the trial parameters. I think this is needed because lsqnonlin treats the real and imaginary parts as separate fitting parameters, and all fitting parameters need to be in a row.
%% Evaluate the error (difference between the fitted function and the desired function to which we're fitting)
error= exp(repmat(fittedParameters((n/4)+1:n/2),length(X),1).*repmat(X.',1,(n/4)))*(fittedParameters(1:n/4).')- Y; %[exp(omega1)AtTimeT1 exp(omega2)AtTimeT1 ... ; exp(omega1)AtTimeT2 exp(omega2)AtTimeT2]*[p1 ; p2 ; ..]
%error = (1/pi)*(1-exp(-beta*X.')).*real((1./(repmat(-1i*fittedParameters((n/4)+1:n/2)*3.33563759345165*10^-11,length(X),1)+repmat(X.',1,(n/4))))*(1i*fittedParameters(1:n/4).'))- Y;
%%
error=error.*weights;
realAndImaginaryPartsOfErrorInOneColumnVector=[real(error);imag(error)]; % Convert the answer back into a column vector such that the first half is the real part, and the second half is the imaginary part. 
end
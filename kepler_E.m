%% Solution of Keplers Equation (E - e*sin(E) = M) by Newtons Method
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function E = kepler_E(e, M)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%...Set an error tolerance:
error = 1.e-10; 

%...Select a starting value for E:
if M < pi
    E = M + e/2; 
else
    E = M - e/2;
end

%...Iterate on Equation 3.17 until E is determined to within
%...the error tolerance:
ratio = 1;

while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end
%fprintf("\n Converged! Error Tolerance: %g  \n",error)
%fprintf('––––––––––––––––––––––––––––––––––––––––––––––––––––– \n')
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Created by Cameron Edwards, Florida Institute of Technology (April - May 2026)
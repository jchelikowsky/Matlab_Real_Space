function [vxc, exc] = exc_nspn(Domain, rho, fid)
%% [rhoxc, exc] = function (Domain, rho)
%%      LDA variables:
%%      various recurring parameters in the exchange-correlation
%%      formulas:
%%      quantities common to all local density expressions:
%%      rs is the local value of the Wigner-Seitz radius
%%      rho temporarily stores the charge at a given grid point
%%      Ceperley-Alder parameters
%%      fid = output file id
global OPTIMIZATIONLEVEL enableMexFilesTest



g  =-0.2846;  b1 = 1.0529;
b2 = 0.3334; c1 = 0.0622;
c2 = 0.096;  c3 = 0.004;
c4 = 0.0232; c5 = 0.0192;
%%
zero = 0.0;     one = 1.0;    two = 2.0;     four = 4.0;  nine = 9.0;
third = 1/3;
%%      actual number of grid points
ndim = length(rho);
%%   ---------------------------------------------------------------
a0 = (four/(nine*pi))^third;
twovpia0 = two/(pi*a0);
p75vpi = 0.75/pi;
%%
vxc = rho;
vxcCopy=vxc;
%%      find max and min values of the charge density and report
dmax = max(vxc(1:ndim));
dmin = min(vxc(1:ndim));
%%      warn (but do not kill) if negative values found
%%
fprintf(fid,' max and min values of charge density [e/bohr^3]');
fprintf(fid,'   %10.5e   %10.5e  \n', dmax, dmin);
if (dmin < zero)
    fprintf(1,'warning in excorr.f: \n')
    fprintf(1,'NEGATIVE CHARGE DENSITY FOUND .\n')
end

%%      Ceperly-Alder exchange correlation

if (OPTIMIZATIONLEVEL~=0)
    [vxc exc]=Ceperly_Alder(vxc,ndim,twovpia0);
else
    %  initialize the total exchange-correlation energy to zero
    exc = zero;

    %     [vxc2 exc2]=Ceperly_Alder(vxc,ndim,twovpia0);

    for i=1:ndim
        rho = vxc(i);
        vxc(i) = zero;
        if (rho > zero)
            rs = (p75vpi/rho)^third;
            vxc(i) = -twovpia0/rs;
            exc = exc + 0.75*rho*vxc(i);
            if (rs >= one)
                sqrs = sqrt(rs);
                ec = g/(one + b1*sqrs + b2*rs);
                vxc(i) = vxc(i) + ...
                    ec*ec*(one+3.5*b1*sqrs*third+four*b2*rs*third)/g;
            else
                alpha = log(rs);
                ec = c1*alpha - c2 + (c3*alpha - c4)*rs;
                vxc(i) = vxc(i) + ...
                    ec - (c1 + (c3*alpha - c5)*rs)*third;
            end
            exc = exc + rho*ec;
        end
    end
    %
    if (enableMexFilesTest==1)
        [vxc2 exc2]=Ceperly_Alder(vxcCopy,ndim,twovpia0);
        if (any(abs(vxc-vxc2)>0.000001) || abs(exc-exc2)>0.000001)
            exception=struct('message','Mex file descreptency for Ceperly_Alder.c','identifier',[],'stack',{mfilename('fullpath') 'filler'});
            logError(exception);
        end
    end
end


%%
%%      scale the total energy integral by h^3 (factor necessary due to the;
%%      expression of an integral as a summation over grid points);
%%
exc  = exc * (Domain.h)^3;

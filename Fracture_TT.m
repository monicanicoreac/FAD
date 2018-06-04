function [Lr, Kr_a, Krcr, KI_a, rhoa, sigmaResM, sigmaResB, Ksb] = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT, fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Ma, Mma,  Mba, Kmat)
%% Plastic collape parameter Lr 
% Local collapse analysis:
% Cut-off value of Lr (to prevent plastic colapse)
if fyweldT == 0  % no weld
    Lrmax = (fyT+fuT)/(2*fyT);
else
    Lrmax = (fyweldT+fuweldT)/(2*fyweldT);
end
% P.5 in BS7910:2013 A1 2015 Reference stress and/or limit load solutions for deep notched plates
constraint = 2;
sigmaref = referenceStress(flaw, constraint, a, c, T, W, sigmaSdb, sigmaSdm);
if fyweldT == 0
    Lr = sigmaref/fyT;
else
    Lr = sigmaref/fyweldT;
end

%% Fracture collapse parameter Kr

% 7.1.8.2 in BS7910:2013 A1 2015 Residual stresses in as-welded structures:
% Membrane residual stress:
% structure in as-welded consition, with a flaw lying in a plane parallel
% to the welding direction (the stresses to be considered are perpendicular
% to the weld): lesser of the RT yield strengths of the weld or parent
% material
% structure in as-welded condition, with a flaw lying in a plane
% transversal to the welding direction (the stresses to be considered are 
% parallel to the weld): RT yield strength of the material in which the
% flaw is located

% 7.1.8.3 in BS7910:2013 A1 2015: Post-welded heat-treated (PWHT) structures
%  a flaw lying in a plane parallel to the welding direction (the stresses
%  to be considered are perpendicular to the weld): 20% of the lesser of
%  the yield strength of the weld or parent material
%  a flaw lying in a plane transverse to the welding direction (the stresses
%  to be considered are parallel to the weld): 30% of RT yield strength of 
% the material in which the flaw is located
if (weld == 0)
    sigmaResM = 0;
    sigmaResB = 0;
    Ksb       = 0;
elseif (weld == 11) || (weld == 21)
    % full penetration flaw, flaw perpendicular to weld direction,as-welded
    sigmaResM = min(fyT, fyT*(1.4 - 2*sigmaref/(fyT + fuT)));
    sigmaResB = 0;
    Ksb       = 0;
elseif (weld == 12) || (weld == 22)
    % full penetration flaw, flaw parallel to weld direction,as-welded
    fy = min(fyT, fyweldT);
    fu = min(fuT, fuweldT);
    sigmaResM = min(fy, fy*(1.4 - 2*sigmaref/(fy + fu)));
    sigmaResB = 0;
    Ksb       = 0;
elseif (weld == 13) || (weld == 23)
    % full penetration flaw, flaw perpendicular to weld direction, PWHT
    sigmaResM = 0.30 * fyT;
    sigmaResB = 0;
    Ksb       = 0;
elseif (weld == 14) || (weld == 24)
    % full penetration flaw, flaw parallel to weld direction, PWHT
    sigmaResM = 0.20 * min(fyT, fyweldT);
    sigmaResB = 0;
    Ksb       = 0;
end

% Annex R: Plasticity correction factor for secondary stresses
KIp_a = (Yma * sigmaSdm + Yba * sigmaSdb)*sqrt(pi*a);

KIs_a = (Mma * sigmaResM + Mba * sigmaResB)*sqrt(pi*a) + Ksb;

chia = KIs_a * Lr / KIp_a;

if KIs_a <= 0
    rhoa = 0;
elseif (chia > 4) && (KIs_a >0)
    disp('Error in Fracture: plasticity correction factor rhoa not implemnted for chia>4')
elseif (chia <= 4) && (Lr <= 0.8) && (KIs_a >0)
    rhoa = 0.1*chia^0.714 - 0.007*chia^2 + 3e-5*chia^5;
elseif (chia <= 4) && (Lr > 0.8) && (Lr < 1.05) && (KIs_a >0)
    rhoa = 4*(0.1*chia^0.714 - 0.007*chia^2 + 3e-5*chia^5)*(1.05 - Lr);
elseif (chia <= 4) && (Lr > 1.05) && (KIs_a >0)
    rhoa = 0;
end

KI_a = KIp_a + KIs_a;

Kr_a = KI_a/Kmat + rhoa;

% Failure assessment diagram
if fyT <= 1000
    deltaEpsilon = 0.037 * (1 - 0.001 * fyT);
else
    disp('Error in Fracture: Yield strength has to lower dan 1000MPa');
return
end
lambda = max(1, 1 + E * deltaEpsilon/fyT);
Nfrac  = 0.3 * (1 - fyT/fuT);
if (Lr >= Lrmax)
    Krcr = 0;
elseif (Lr < 1)
    Krcr = 1/sqrt(1 + 0.5*Lr^2);
elseif (Lr == 1)
    Krcr = 1/sqrt(lamda + 0.5/lambda);
elseif (Lr > 1) && (Lr < Lrmax)
    Krcr = 1/sqrt(lambda + 0.5/lambda) * Lr^((Nfrac-1)/(2*Nfrac));
end
end

function [Lr, Kr_a, Kr_c, Krcr, KI_a, KI_c,rhoa, rhoc] = Fracture_tubular(a, c, fyT, fuT, crackLocation, weld, GeometryTub, LoadsTub, JointType, Kmat, Kmati, Yma, Yba, Ymc, Ybc, sigmaSdm, sigmaSdb,Mma, Mba, Mmc, Mbc)
%% Plastic collape parameter Lr 
% Local collapse analysis:
% Cut-off value of Lr (to prevent plastic colapse)
Lrmax = (fyT+fuT)/(2*fyT);

% Lrmax = (fy+fu)/(2*fy);
% P.5 in BS7910:2013 A1 2015 Reference stress and/or limit load solutions for deep notched plates

% Global analysis:
% Flow strenght: 
fflow = (fyT + fuT) / 2;
[ Pcchord Ptchord Mcichord Mcochord Pcbrace Mcibrace Mcobrace ] = LrLimitLoad( a, c, GeometryTub, LoadsTub, fyT, JointType );
if LoadsTub(4,1) >= 0 % tension force in brace
    Lrchord = (fflow / fyT) * (abs(LoadsTub(4,1)/Ptchord) + (LoadsTub(5,1) / Mcichord)^2 + abs(LoadsTub(6,1)/Mcochord));
else  % compression force in brace
    Lrchord = (fflow / fyT) * (abs(LoadsTub(4,1)/Pcchord) + (LoadsTub(5,1) / Mcichord)^2 + abs(LoadsTub(6,1)/Mcochord));
end
Lrbrace = (fflow / fyT) * (abs(LoadsTub(4,1)/Pcbrace) + (LoadsTub(5,1) / Mcibrace)^2 + abs(LoadsTub(6,1)/Mcobrace));
if (crackLocation == 'chord')
    Lr = Lrchord;
elseif (crackLocation == 'brace')
    Lr = Lrbrace;
end

%% Fracture collapse LoadsTub(4,1)rameter Kr

% 7.1.8.2 in BS7910:2013 A1 2015 Residual stresses in as-welded structures:
% Membrane residual stress:
% structure in as-welded consition, with a flaw lying in a plane LoadsTub(4,1)rallel
% to the welding direction (the stresses to be considered are perpendicular
% to the weld): lesser of the RT yield strengths of the weld or LoadsTub(4,1)rent
% material
% structure in as-welded condition, with a flaw lying in a plane
% transversal to the welding direction (the stresses to be considered are 
% LoadsTub(4,1)rallel to the weld): RT yield strength of the material in which the
% flaw is located

% 7.1.8.3 in BS7910:2013 A1 2015: Post-welded heat-treated (PWHT) structures
%  a flaw lying in a plane LoadsTub(4,1)rallel to the welding direction (the stresses
%  to be considered are perpendicular to the weld): 20% of the lesser of
%  the yield strength of the weld or LoadsTub(4,1)rent material
%  a flaw lying in a plane transverse to the welding direction (the stresses
%  to be considered are LoadsTub(4,1)rallel to the weld): 30% of RT yield strength of 
% the material in which the flaw is located
if (weld == 0)
    sigmaResM = 0;
    sigmaResB = 0;
elseif (weld == 11)
    % tubular T-butt weld (full penetration flaw), flaw perpendicular 
    % to weld direction,as-welded (Annex Q)
    sigmaResM = 0.564 * fyT;
    sigmaResB = (0.307 - 0.614*a/GeometryTub(1,1)) * fyT;
    Ksb = 0.3 * fyT * sqrt(GeometryTub(1,1));
elseif (weld == 12)
    % tubular T-butt weld (full penetration flaw), flaw parallel to weld 
    % direction,as-welded (Annex Q)
    sigmaResM = 0.496 * fyT;
    sigmaResB = (0.216 - 0.432*a/GeometryTub(1,1))*fyT;
    Ksb = 0.21 * fyT * sqrt(GeometryTub(1,1));
elseif (weld == 13)
    % full penetration flaw, flaw perpendicular to weld direction, PWHT
    sigmaResM = 0.30 * fyT;
    sigmaResB = 0;
elseif (weld == 14) || (weld == 24)
    % full penetration flaw, flaw LoadsTub(4,1)rallel to weld direction, PWHT
    sigmaResM = 0.20 * min(fyT, fyweldT);
    sigmaResB = 0;
end

% Annex R: Plasticity correction factor for secondary stresses
KIp_a = (Yma * sigmaSdm + Yba * sigmaSdb)*sqrt(pi*a);
KIp_c = (Ymc * sigmaSdm + Ybc * sigmaSdb)*sqrt(pi*a);

KIs_a = (Mma * sigmaResM + Mba * sigmaResB) * sqrt(pi*a) + Ksb;
KIs_c = (Mmc * sigmaResM + Mbc * sigmaResB) * sqrt(pi*a) + Ksb;

chia = KIs_a * Lr / KIp_a;
chic = KIs_c * Lr / KIp_c;

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

if KIs_c <= 0
    rhoc = 0;
elseif (chic > 4) && (KIs_c >0)
    disp('Error in Fracture: plasticity correction factor rhoa not implemnted for chia>4')
elseif (chic <= 4) && (Lr <= 0.8) && (KIs_c >0)
    rhoc = 0.1*chic^0.714 - 0.007*chic^2 + 3e-5*chic^5;
elseif (chic <= 4) && (Lr > 0.8) && (Lr < 1.05) && (KIs_c >0)
    rhoc = 4*(0.1*chic^0.714 - 0.007*chic^2 + 3e-5*chic^5)*(1.05 - Lr);
elseif (chic <= 4) && (Lr > 1.05) && (KIs_c >0)
    rhoc = 0;
end
KI_a = KIp_a + KIs_a;
KI_c = KIp_c + KIs_c;

Kr_a = KI_a/Kmat + rhoa;
% For c (surface point) no correction for crack front length is needed:
Kr_c = KI_c/Kmati + rhoc;

% Failure assessment diagram
if fyT <= 1000
    deltaEpsilon = 0.037 * (1 - 0.001 * fyT);
else
    disp('Error in Fracture: Yield strength haqs to lower dan 1000MLoadsTub(4,1)');
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

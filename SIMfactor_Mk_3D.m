function [Mkma, Mkba, Mkmc, Mkbc] = SIMfactor_Mk_3D(a,c, T, L, phi, rho)
% M.11 Welded joints
% M.11.1.2 Solution based on 3D finite element analysis
% Value of crack depth to be used in calculations
ac = a/c;
ca = c/a;
aT = a/T;
LT = L/T;
% if (a <= 0.005*T)
%     disp('Error calculating weld geometry factor in Subroutine SIMfactor_Mk_3D: a/T < 0.005')
%     return
% end 
% if (a >= 0.9*T)
% 	disp('Error calculating weld geometry factor in Subroutine SIMfactor_Mk_3D: a/T > 0.9')
%     return
% end 
% % if (a < 0.1*c)
% % 	disp('Error calculating weld geometry factor in Subroutine SIMfactor_Mk_3D: a/c < 0.1')
% %     return
% % end
% if (a > c)
% 	disp('Error calculating weld geometry factor in Subroutine SIMfactor_Mk_3D: a/c > 1')
%     return
% end 
% if L < 0.5*T
%     disp('Error calculating weld geometry factor in Subroutine SIMfactor_Mk_3D: L/T < 0.5')
%     return
% end
if (L > 2.75*T) 
    LT = 2.75;
end
% Membrane - deepest point
g1 = -1.0343 * (ac^2) - 0.15657*(ac) + 1.3409;
g2 =  1.3218 * (ac^(-0.61153));
g3 = -0.87238 * (ac) + 1.2788;
g4 = -0.46190 *(ac^3) + 0.67090 *(ac^2) - 0.37571 *(ac) + 4.6511;
g5 = -0.015647 *(LT^3) + 0.090889 *(LT^2 ) - 0.17180 *(LT) - 0.24587;
g6 = -0.20136 *(LT^2) + 0.93311 *(LT) - 0.41496;
g7 = 0.20188 *(LT^2) - 0.97857 *(LT) + 0.068225;
g8 = -0.027338 *(LT^2) + 0.12551 *(LT) - 11.218;
f1 = 0.43358*(aT)^( g1 + (g2*aT)^g3 ) + 0.93163*exp( (aT)^(-0.050966) ) + g4;
f2 = -0.21521*((1-aT)^176.4199) + 2.8141*(aT^(-0.10740*aT));
f3 = 0.33994*(aT^g5) + 1.9493*(aT^0.23003) + (g6*(aT^2) + g7*aT + g8);
Mkma = max(f1 + f2 + f3, 1.0);

% Bending  - deepest point
if (aT < 0.5)
	Mkba = 1.;
elseif (aT > 0.9)
	Mkba = 1.;
else
	g1 = -0.014992*(ac^2) - 0.021401*(ac) - 0.23851;
	g2 = 0.61775*(ac^(-1.0278));
	g3 = 0.00013242*(ac) - 1.4744;
	g4 = -0.28783*(ac^3.) + 0.58706*(ac^2) - 0.37198*(ac) - 0.89887;
	g5 = -17.195*(aT^2.) + 12.468*(aT) - 0.51662;
	g6 = -0.059798*(LT^3.) + 0.38091*(LT^2) - 0.8022020*(LT) + 0.31906;
	g7 = -0.35848*(LT^2.) + 1.3975*(LT) - 1.7535;
	g8 = 0.31288*(LT^2.) - 1.3599*(LT) + 1.6611;
	g9 = -0.001470*(LT^2.) - 0.0025074*(LT) - 0.0089846;
	f1 = 0.065916*(aT)^( g1 + (g2*aT)^g3 ) + 0.52086*exp( (aT)^(-0.10364) ) + g4;
	f2 = -0.021950*((1-aT)^2.8086) + 0.021403*(aT^g5) ;
	f3 = 0.23344*(aT^g6) - 0.14827*(aT^(-0.20077)) + (g7*(aT^2) + g8*aT + g9);
	Mkba = max(f1 + f2 + f3, 1.0);
end
% Membrane - surface point
g1 = 0.0078157*((ca^2.)) - 0.070664*(ca) + 1.8508;
g2 = -0.000054546*(LT^2.) + 0.0013651*(LT) - 0.00047844;
g3 = 0.00049192*(LT^2) - 0.0013595*LT + 0.011400;
g4 = 0.0071654*(LT^2.) - 0.033399*(LT) - 0.25064;
g5 = -0.018640*(ca^2.) + 0.24311*(ca) - 1.7644;
g6 = -0.0016713*(LT^2.) + 0.0090620*(LT) - 0.016479;
g7 = -0.0031615*(LT^2.) - 0.010944*(LT) + 0.13967;
g8 = -0.045206*(LT^3.) + 0.32380*(LT^2.) - 0.68935*(LT) + 1.4954;
g9 = -0.25473*(ac^2.) + 0.40928*(ac) + 0.0021892;
g10 = 37.423*(ac^2.) - 15.741*(ac) + 64.903;
g11 = -0.10553*(LT^3.) + 0.59894*(LT^2.) - 1.0942*(LT) - 1.2650;
g12 = 0.043891*(LT^3.) - 0.24898*(LT^2.) + 0.44732*(LT) + 0.60136;
g13 = -0.011411*(ac^2.) + 0.004369*(ac) + 0.51732;
f1 = g1*((aT)^( g2*(ca^2.) + g3*(ca) +g4 )) + g5*((1-aT)^(g6*(ca^2.) + g7*(ca) +g8));
f2 = (-0.28639*(ac^2.) + 0.35411*(ac) + 1.6430)*(aT^g9) + 0.27349*((1-aT)^g10);
f3 = g11*(aT^0.75429) + g12*exp(aT^g13);
Mkmc = max(f1*f2*f3, 1.0);

% Bending - surface point
g1 = 0.0023232*((ca^2.)) - 0.00037156*(ca) + 4.5985;
g2 = -0.000044010*(LT^2.) + 0.00014425*(LT) - 0.00086706;
g3 = 0.00039951*(LT^2.) - 0.0013715*(LT) + 0.014251;
g4 = 0.0046169*(LT^2.) - 0.017917*(LT) - 0.16335;
g5 = -0.018524*((ca^2.)) + 0.27810*(ca) - 5.4253;
g6 = -0.00037981*(LT^2.) + 0.0025078*(LT) + 0.00014693;
g7 = -0.0038508*(LT^2.) + 0.0023212*(LT) - 0.026862;
g8 = -0.011911*(LT^3.) + 0.082625*(LT^2.) - 0.16086*(LT) + 1.2302;
g9 = 0.27798*(aT^3.) - 1.2144*(aT^2.) - 2.4680*(aT) + 0.099981;
g10 = -0.25922*(ac^2.) + 0.39566*(ac) + 0.011759;
g11 = 6.5964*(ac^2.) + 55.787*(ac) + 37.053;
g12 = -0.14895*(LT^3.) + 0.815*(LT^2.) - 1.4795*(LT) - 0.89808;
g13 = 0.055459*(LT^3.) - 0.30180*(LT^2.) + 0.54154*(LT) + 0.53433;
g14 = -0.01343*(ac^2.) + 0.0066702*(ac) + 0.75939;
f1 = g1*((aT)^( g2*((ca^2.)) + g3*(ca) +g4 )) + g5*((1-aT)^(g6*((ca^2.)) + g7*(ca) +g8)) + g9;
f2 = (-0.35006*(ac^2.) + 0.40768*(ac) + 1.7053)*(aT^g10) + 0.24988*((1-aT)^g11);
f3 = g12*(aT^0.94761) + g13*exp(aT^g14);
Mkbc = max(f1*f2*f3, 1.0);

% Correction factor for weld angle
if (phi >= 25) && (phi <= 65)
    phi_cor = 13.096e-3 + 28.199e-3*phi - 139.45e-6*phi^2;
else
    disp('Error in SIMfactor_Mk: phi has to be between 25 and 65');
    return
end
if (a/T <= 0.1) && (a/T >= 0.001)
    fphia = (10*a/T)^(-0.5*log10(phi_cor));
    fphic = (10*0.15/T)^(-0.5*log10(phi_cor));
elseif a/T > 0.1
    fphia = 1;
elseif a/T < 0.001
    disp('Error in in SIMfactor_Mk_3D: a/T must be larger than 0.001');
end
if (0.15/T <= 0.1) && (0.15/T >= 0.001)
    fphic = (10*0.15/T)^(-0.5*log10(phi_cor));
elseif 0.15/T > 0.1
    fphic = 1;
elseif 0.15/T < 0.001
    disp('Error in in SIMfactor_Mk: c/T must be larger than 0.001');
end

% Correction factor for weld toe radius
if rho == 0
    frhoma = 1;
    frhoba = 1;
    frhomc = 1;
    frhobc = 1;
elseif rho > 0.125*T
    disp('Error in in SIMfactor_Mk: rho should be smaller than T/8');
else
    arhom = 0.71032 - 0.024015/(rho/T + 0.028061);
    brhom = 105.29 - 1993.8*rho^2/(T^2);
    arhob = 0.70754 - 0.020160/(rho/T + 0.024502);
    brhob = 75.323 - 1541.7*rho^2/(T^2);
    frhoma = min(1, 1 - arhom * exp(-brhom*a/T));
    frhoba = min(1, 1 - arhob * exp(-brhob*a/T));
    frhomc = min(1, 1 - arhom * exp(-brhom*0.15/T));
    frhobc = min(1, 1 - arhob * exp(-brhob*0.15/T));
end
Mkma = Mkma * fphia * frhoma;
Mkba = Mkba * fphia * frhoba;
Mkmc = Mkmc * fphic * frhomc;
Mkbc = Mkbc * fphic * frhobc;
end
    
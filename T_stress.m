function [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W )
%% T-stress correction:
% Surface crack, tension, deepest point:
X0am = @(a,c) -0.581176 + 0.429125 * (a/c) - 0.670213 * (a/c)^2 + 0.290259 * (a/c)^3;
X1am = @(a,c) 0.795958 - 6.53119 * (a/c) + 12.1073 * (a/c)^2 - 6.47294 * (a/c)^3;
X2am = @(a,c) -4.70569 + 28.4706 * (a/c) - 51.1588 * (a/c)^2 + 27.3436 * (a/c)^3;
X3am = @(a,c) 6.16031 - 40.1692 * (a/c) + 68.7053 * (a/c)^2 - 35.4684 * (a/c)^3;

if (a/c >= 0.2) && (a/T <= 0.8)
    betaSref_am = sigmaSdm * (X0am(a,c) + X1am(a,c) * (a/T)^2 + X2am(a,c) * (a/T)^4 + X3am(a,c) * (a/T)^6);
elseif (a/c >= 0.2) && (a/T > 0.8)
    betaSref_am = sigmaSdm * (X0am(a,c) + X1am(a,c) * (0.8)^2 + X2am(a,c) * (0.8)^4 + X3am(a,c) * (0.8)^6);
elseif (a/c < 0.2) && (a/T <= 0.8)
    eq1 = sigmaSdm * (X0am(a,c) + X1am(a,c) * (a/T)^2 + X2am(a,c) * (a/T)^4 + X3am(a,c) * (a/T)^6);
    eq2 = sigmaSdm * (X0am(a,a/0.2) + X1am(a,a/0.2) * (a/T)^2 + X2am(a,a/0.2) * (a/T)^4 + X3am(a,a/0.2) * (a/T)^6);
    eq3 = sigmaSdm * (X0am(0.2*c,c) + X1am(0.2*c,c) * (a/T)^2 + X2am(0.2*c,c) * (a/T)^4 + X3am(0.2*c,c) * (a/T)^6);
    betaSref_am = min([eq1  eq2  eq3]);
elseif (a/c < 0.2) && (a/T > 0.8)
    eq1 = sigmaSdm * (X0am(a,c) + X1am(a,c) * (0.8)^2 + X2am(a,c) * (0.8)^4 + X3am(a,c) * (0.8)^6);
    eq2 = sigmaSdm * (X0am(a,a/0.2) + X1am(a,a/0.2) * (0.8)^2 + X2am(a,a/0.2) * (0.8)^4 + X3am(a,a/0.2) * (0.8)^6);
    eq3 = sigmaSdm * (X0am(0.2*c,c) + X1am(0.2*c,c) * (0.8)^2 + X2am(0.2*c,c) * (0.8)^4 + X3am(0.2*c,c) * (0.8)^6);
    betaSref_am = min([eq1  eq2  eq3]);
end
% Surface crack, tension, surface point:
X0cm = @(a,c) 0.062026 - 1.49614  * (a/c) + 0.94377 * (a/c)^2 - 0.15367 * (a/c)^3;
X1cm = @(a,c) 0.59687 - 4.2027 * (a/c) + 6.94022 * (a/c)^2 - 3.29765 * (a/c)^3;
X2cm = @(a,c) -0.020640 + 13.1854 * (a/c) - 30.022 * (a/c)^2  + 16.3561 * (a/c)^3;
X3cm = @(a,c)  -7.82399 + 19.532 * (a/c) - 14.4551 * (a/c)^2  + 3.36514 * (a/c)^3;

if (a/c >= 0.2) && (a/T <= 0.8)
    betaSref_cm = sigmaSdm * (X0cm(a,c) + X1cm(a,c) * (a/T)^2 + X2cm(a,c) * (a/T)^4 + X3cm(a,c) * (a/T)^6);
elseif (a/c >= 0.2) && (a/T > 0.8)
    betaSref_cm = sigmaSdm * (X0cm(a,c) + X1cm(a,c) * (0.8)^2 + X2cm(a,c) * (0.8)^4 + X3cm(a,c) * (0.8)^6);
elseif (a/c < 0.2) && (a/T < 0.8)
    eq1 = sigmaSdm * (X0cm(a,c) + X1cm(a,c) * (a/T)^2 + X2cm(a,c) * (a/T)^4 + X3cm(a,c) * (a/T)^6);
    eq2 = sigmaSdm * (X0cm(a,a/0.2) + X1cm(a,a/0.2) * (a/T)^2 + X2cm(a,a/0.2) * (a/T)^4 + X3cm(a,a/0.2) * (a/T)^6);
    eq3 = sigmaSdm * (X0cm(0.2*c,c) + X1cm(0.2*c,c) * (a/T)^2 + X2cm(0.2*c,c) * (a/T)^4 + X3cm(0.2*c,c) * (a/T)^6);
    betaSref_cm = min([eq1  eq2  eq3]);
elseif (a/c < 0.2) && (a/T > 0.8)
    eq1 = sigmaSdm * (X0cm(a,c) + X1cm(a,c) * (0.8)^2 + X2cm(a,c) * (0.8)^4 + X3cm(a,c) * (0.8)^6);
    eq2 = sigmaSdm * (X0cm(a,a/0.2) + X1cm(a,a/0.2) * (0.8)^2 + X2cm(a,a/0.2) * (0.8)^4 + X3cm(a,a/0.2) * (0.8)^6);
    eq3 = sigmaSdm * (X0cm(0.2*c,c) + X1cm(0.2*c,c) * (0.8)^2 + X2cm(0.2*c,c) * (0.8)^4 + X3cm(0.2*c,c) * (0.8)^6);
    betaSref_cm = min([eq1  eq2  eq3]);
end
% Surface crack, bending, deepest point:
X0ab = @(a,c) -0.37758 + 0.346823 * (a/c) - 0.539154 * (a/c)^2 + 0.221408 * (a/c)^3;
X1ab = @(a,c) 4.61176 - 5.35305 * (a/c) + 9.4601 * (a/c)^2 - 5.03237 * (a/c)^3;
X2ab = @(a,c) -9.85178 + 20.786 * (a/c) - 36.8869 * (a/c)^2  + 19.806 * (a/c)^3;
X3ab = @(a,c) 10.8584 - 30.3672 * (a/c) + 49.8859 * (a/c)^2  - 25.7138 * (a/c)^3;

if (a/c >= 0.2) && (a/T <= 0.8)
    betaSref_ab = sigmaSdb * (X0ab(a,c) + X1ab(a,c) * (a/T)^2 + X2ab(a,c) * (a/T)^4 + X3ab(a,c) * (a/T)^6);
elseif (a/c >= 0.2) && (a/T > 0.8)
    betaSref_ab = sigmaSdb * (X0ab(a,c) + X1ab(a,c) * (0.8)^2 + X2ab(a,c) * (0.8)^4 + X3ab(a,c) * (0.8)^6);
elseif (a/c < 0.2) && (a/T < 0.8)
    eq1 = sigmaSdb * (X0ab(a,c) + X1ab(a,c) * (a/T)^2 + X2ab(a,c) * (a/T)^4 + X3ab(a,c) * (a/T)^6);
    eq2 = sigmaSdb * (X0ab(a,a/0.2) + X1ab(a,a/0.2) * (a/T)^2 + X2ab(a,a/0.2) * (a/T)^4 + X3ab(a,a/0.2) * (a/T)^6);
    eq3 = sigmaSdb * (X0ab(0.2*c,c) + X1ab(0.2*c,c) * (a/T)^2 + X2ab(0.2*c,c) * (a/T)^4 + X3ab(0.2*c,c) * (a/T)^6);
    betaSref_ab = min([eq1  eq2  eq3]);
elseif (a/c < 0.2) && (a/T > 0.8)
    eq1 = sigmaSdb * (X0ab(a,c) + X1ab(a,c) * (0.8)^2 + X2ab(a,c) * (0.8)^4 + X3ab(a,c) * (0.8)^6);
    eq2 = sigmaSdb * (X0ab(a,a/0.2) + X1ab(a,a/0.2) * (0.8)^2 + X2ab(a,a/0.2) * (0.8)^4 + X3ab(a,a/0.2) * (0.8)^6);
    eq3 = sigmaSdb * (X0ab(0.2*c,c) + X1ab(0.2*c,c) * (0.8)^2 + X2ab(0.2*c,c) * (0.8)^4 + X3ab(0.2*c,c) * (0.8)^6);
    betaSref_ab = min([eq1  eq2  eq3]);
end
% Surface crack, bending, surface point:
X0cb = @(a,c) 0.0433442 - 1.4636 * (a/c) + 0.97 * (a/c)^2 - 0.182882 * (a/c)^3;
X1cb = @(a,c) 0.787135 - 5.06806 * (a/c) + 8.63486 * (a/c)^2 - 4.13922 * (a/c)^3;
X2cb = @(a,c) -1.66573 + 16.8231 * (a/c) - 31.5732 * (a/c)^2  + 15.831 * (a/c)^3;
X3cb = @(a,c) -2.80697 + 1.64157 * (a/c) + 6.29347 * (a/c)^2  - 4.54519 * (a/c)^3;

if (a/c >= 0.2) && (a/T <= 0.8)
    betaSref_cb = sigmaSdb * (X0cb(a,c) + X1cb(a,c) * (a/T)^2 + X2cb(a,c) * (a/T)^4 + X3cb(a,c) * (a/T)^6);
elseif (a/c >= 0.2) && (a/T > 0.8)
    betaSref_cb = sigmaSdb * (X0cb(a,c) + X1cb(a,c) * (0.8)^2 + X2cb(a,c) * (0.8)^4 + X3cb(a,c) * (0.8)^6);
elseif (a/c < 0.2) && (a/T < 0.8)
    eq1 = sigmaSdb * (X0cb(a,c) + X1cb(a,c) * (a/T)^2 + X2cb(a,c) * (a/T)^4 + X3cb(a,c) * (a/T)^6);
    eq2 = sigmaSdb * (X0cb(a,a/0.2) + X1cb(a,a/0.2) * (a/T)^2 + X2cb(a,a/0.2) * (a/T)^4 + X3cb(a,a/0.2) * (a/T)^6);
    eq3 = sigmaSdb * (X0cb(0.2*c,c) + X1cb(0.2*c,c) * (a/T)^2 + X2cb(0.2*c,c) * (a/T)^4 + X3cb(0.2*c,c) * (a/T)^6);
    betaSref_cb = min([eq1  eq2  eq3]);
elseif (a/c < 0.2) && (a/T > 0.8)
    eq1 = sigmaSdb * (X0cb(a,c) + X1cb(a,c) * (0.8)^2 + X2cb(a,c) * (0.8)^4 + X3cb(a,c) * (0.8)^6);
    eq2 = sigmaSdb * (X0cb(a,a/0.2) + X1cb(a,a/0.2) * (0.8)^2 + X2cb(a,a/0.2) * (0.8)^4 + X3cb(a,a/0.2) * (0.8)^6);
    eq3 = sigmaSdb * (X0cb(0.2*c,c) + X1cb(0.2*c,c) * (0.8)^2 + X2cb(0.2*c,c) * (0.8)^4 + X3cb(0.2*c,c) * (0.8)^6);
    betaSref_cb = min([eq1  eq2  eq3]);
end
%% Centre-cracked plate in uniaxial tension
sigmaref_TT = sqrt(3)*sigmaSdm/(2*(1-2*a/W));
betaSref_TT = sigmaref_TT * (-1.1547 + 1.1476 * (2*a/W) -2.4091 * (2*a/W)^2 + 4.059 * (2*a/W)^3 -1.9907 * (2*a/W)^4);
% betaSref_TT = sigmaref_TT * (-1.1547 + 1.1511 * (2*a/W) + -0.7826 * (2*a/W)^2 + 0.4751 * (2*a/W)^3 + -0.1761 * (2*a/W)^4);

%% Cracks at the edge of circular hole
sigmaref_EH = sqrt(3)*sigmaSdm/(2*(1-2*a/W));
betaSref_EH = sigmaref_EH * (-0.65);
%% Double edge-crack plate: uniaxial tension
sigmaref_DECT = sigmaSdm * sqrt(3)/(2*(1+log((0.5*W-a)/(0.5*W-2*a)))*(1-2*2*a/W));
betaSref_DECT = sigmaref_DECT * (-0.5889 - 0.0097*(4*a/W) + 1.1103 * (4*a/W)^2 - 1.3852 * (4*a/W)^3 + 0.6573 * (4*a/W)^4);

% Bending or tension
if (sigmaSdb <= 0.5*sigmaSdm) && (flaw == 2 || flaw == 5)
    Tstr_a = min(0, betaSref_am);
elseif (sigmaSdb > 0.5*sigmaSdm) && (flaw == 2 || flaw == 5)
    Tstr_a = min(0, betaSref_ab);
elseif (flaw == 1) || (flaw == 3) % through thickness or trough width 
    Tstr_a = min(0, betaSref_TT);
elseif (flaw == 9)
    Tstr_a = min(0, betaSref_EH);
elseif (flaw == 10)
    Tstr_a = min(0, betaSref_DECT);
else
    Tstr_a = 0;
end
if (sigmaSdb <= 0.5*sigmaSdm) && (flaw == 2 || flaw == 5)
    Tstr_c = min(0, betaSref_cm);
elseif (sigmaSdb > 0.5*sigmaSdm) && (flaw == 2 || flaw == 5)
    Tstr_c = min(0, betaSref_cb);
else
    Tstr_c = 0;
end
end


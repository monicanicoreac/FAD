function [Mkm, Mkb] = SIMfactor_Mk(weld, z, T, L, tw, phi, rho)
% M.11 Welded joints
% M.11.1.2 Solution based on 2D finite element analysis
% Value of crack depth to be used in calculations
zc = z;
if zc < 0.15
    zc = 0.15;
elseif zc > T
    zc = T;
end

% Check weld geometry
if (weld < 1) && (weld > 2)
    disp('Error in weld geometry definition: weld is 1 or 2');
    return
end

% Axial loading
if L/T <= 2
    if zc/T <= 0.05*(L/T)^0.55
        vm = 0.51 * (L/T)^0.27;
        wm = -0.31;
    else
        vm = 0.83;
        wm = -0.15 * (L/T)^0.46;
    end
elseif L/T > 2
    if zc/T <= 0.073
        vm = 0.615;
        wm = -0.31;
    else
        vm = 0.83;
        wm = -0.20;
    end
end

% Bending loading
if L/T <= 1
    if zc/T <= 0.03*(L/T)^0.55
        vb = 0.45 * (L/T)^0.21;
        wb = -0.31;
    else
        vb = 0.68;
        wb = -0.19 * (L/T)^0.21;
    end
elseif L/T > 1
    if zc/T <= 0.03
        vb = 0.45;
        wb = -0.31;
    else
        vb = 0.68;
        wb = -0.19;
    end
end

% Weld geometry 2: flaws at the toes of load-carrying or partial
% penetration welds (Fig. M.26) the values of v and w are those for L/T >2
% for axial and L/T >1 for bending
if weld == 2
    if zc/T <= 0.073
        vm = 0.615;
        wm = -0.31;
    else
        vm = 0.83;
        wm = -0.20;
    end
    if zc/T <= 0.03
        vb = 0.45;
        wb = -0.31;
    else
        vb = 0.68;
        wb = -0.19;
    end
    vm = vm * sqrt(T/tw);
    vb = vb * sqrt(T/tw);
end

fMkm = max(1, vm*(zc/T).^wm);
fMkb = max(1, vb*(zc/T).^wb);

% Correction factor for weld angle
if (phi >= 25) && (phi <= 65)
    phi_cor = 13.096e-3 + 28.199e-3*phi - 139.45e-6*phi^2;
else
    disp('Error in SIMfactor_Mk: phi has to be between 25 and 65');
    return
end
if (zc/T <= 0.1) && (zc/T >= 0.001)
    fphi = (10*zc/T)^(-0.5*log10(phi_cor));
elseif zc/T > 0.1
    fphi = 1;
elseif zc/T < 0.001
    disp('Error in in SIMfactor_Mk: a/T must be larger than 0.001');
end

% Correction factor for weld toe radius
if rho == 0
    frhom = 1;
    frhob = 1;
elseif rho > 0.125*T
    disp('Error in in SIMfactor_Mk: rho should be smaller than T/8');
else
    arhom = 0.71032 - 0.024015/(rho/T + 0.028061);
    brhom = 105.29 - 1993.8*rho^2/(T^2);
    arhob = 0.70754 - 0.020160/(rho/T + 0.024502);
    brhob = 75.323 - 1541.7*rho^2/(T^2);
    frhom = min(1, 1 - arhom * exp(-brhom*zc/T));
    frhob = min(1, 1 - arhob * exp(-brhob*zc/T));
end

Mkm = max(1, fMkm * fphi * frhom);
Mkb = max(1, fMkb * fphi * frhob);
end
    
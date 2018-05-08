function sigmaref = referenceStress(flaw, constraint, a, c, T, W, sigmaSdb, sigmaSdm)
if flaw == 2    % Surface flaws in plates: P.6.1 in BS7910:2013 A1 2015
    if W >= 2*(c+T)
        alfadacc = (a/T)/(1+T/c);
    else
        alfadacc = 2*(a/T)*(c/W);
    end
    if constraint == 1 % Normal bending restraint (uniform remote stress + bending)
        sigmaref = (sigmaSdb + sqrt(sigmaSdb^2 + 9 * sigmaSdm^2*(1-alfadacc)^2))./(3*(1-alfadacc)^2);
    elseif constraint == 2 % Neglijible bending restraint:
        sigmaref = (sigmaSdb + 3 * sigmaSdm * alfadacc + sqrt((sigmaSdb + 3 * sigmaSdm * alfadacc)^2 + (9 * sigmaSdm^2 * (1-alfadacc)^2))) ./ (3 * (1-alfadacc)^2);
    end
%%
elseif flaw == 1 % Through-thickness flaws in plates: P.5.1 in BS7910:2013 A1 2015
    sigmaref = (sigmaSdb + sqrt(sigmaSdb^2 + 9 * sigmaSdm^2))/(3*(1-2*a/W));
%%
elseif flaw == 3 % Through-width flaws in plates: P.6.2 in BS7910:2013 A1 2015
    alfadacc = a/T;
    if constraint == 1 % Normal bending restraint (uniform reote stress + bending)
        sigmaref = (sigmaSdb + sqrt(sigmaSdb^2 + 9 * sigmaSdm^2*(1-alfadacc)^2))./(3*(1-alfadacc)^2);
    elseif constraint == 2 % Neglijible bending restraint:
        sigmaref = (sigmaSdb + 3 * sigmaSdm * alfadacc + sqrt((sigmaSdb + 3 * sigmaSdm * alfadacc)^2 + (9 * sigmaSdm^2 * (1-alfadacc)^2))) ./ (3 * (1-alfadacc)^2);
    end
elseif flaw == 4 %semi-circular surface flaw in round bar
    disp('error in referenceStress: flaw 4 not implemented');
    return
elseif flaw == 5 % Corner flaws in plates: P.7.1 in BS7910:2013 A1 2015
    if W >= c + 2*T
        alfadacc = (a/T)/(1+2*T/c);
    else
        alfadacc = 2*(a/T)*(c/W);
    end
    sigmaref = (sigmaSdb + sqrt(sigmaSdb^2 + 9 * sigmaSdm^2*(1-alfadacc)^2))/(3*(1-alfadacc)^2);
elseif flaw == 6 % Semi-circular surface flaw in bolts, solution 2
    disp('error in referenceStress: flaw 6 not implemented');
    return
elseif flaw == 7 % Edge flaws in plates
    disp('error in referenceStress: flaw 7 not implemented');
    return
elseif flaw == 8 % Through-thickness deck plate crack (orthotropic bridge deck)
    disp('error in referenceStress: flaw 8 not implemented');
    return
elseif flaw == 9 % symmetric crack emanating from hole (taken as through thickness)
    sigmaref = (sigmaSdb + sqrt(sigmaSdb^2 + 9 * sigmaSdm^2))/(3*(1-2*a/W));
end
    
    
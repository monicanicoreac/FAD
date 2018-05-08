function err = SIMconditions(a, c, T, W, radius, flaw)
% Stress intensity magnification factor Mm
% -1 = Through-thickness flaws in plates
% -2 = Surface flaws in plates (= flaw 1 in mathcad)
% -3 = Through-width flaws in plates
% -4 = semi-circular surface flaw in round bar
% -5 = Corner flaws in plates
% 6 = Semi-circular surface flaw in bolts, solution 2
% -7 = Edge flaws in plates
% -8 = Through-thickness deck plate crack (orthotropic bridge deck)
err = 0;
if flaw == 2
    if a < 0
        err = 1;
    elseif c < 0
        err = 1;
    elseif a/T > 1
        err = 1;
    elseif a/c > 2
        err = 1;
%     elseif c/W > 0.4
%         err = 1;
    end
elseif flaw == 3
    if a/T > 0.6 
        err = 1;
    end
elseif flaw == 4
    if a/radius > 1.2
        err = 1;
    end
elseif flaw == 5
    if a/c < 0.2
        err = 1;
    elseif a/c > 2
        err = 1;
    elseif a/T > 1
        err = 1;
    elseif c/W > 0.5
        err =1;
    end
end
end
    
function [Ym, Yb] = SICfactorY (index, a, c, T, W, L, tw, phi, radius, radiusBar, weld, flaw, theta, hole)
% Stress intensity correction factor Annex M 
err = SIMconditions(a, c, T, W, radius, flaw);

if err == 0
    if weld == 0 % no weld
        Mkma = 1;
        Mkba = 1;
        Mkmc = 1;
        Mkbc = 1;
    elseif (index == 9) || (index == 10) || (index == 27) || (index == 28) || (index == 30) || (index == 32) || (index == 33)   % weld but no weld influence??
        Mkma = 1;
        Mkba = 1;
        Mkmc = 1;
        Mkbc = 1;
    elseif (flaw == 10)
        Mkma = 1;
        Mkba = 1;
        Mkmc = 1;
        Mkbc = 1;
%     elseif (a <= 0.005*T) || (a >= 0.9*T) || (a < 0.1*c) || (a > c) || (L < 0.5*T)
%         % Use 2D solution
%         [Mkma, Mkba] = SIMfactor_Mk(weld, a, T, L, tw, phi, radius);
%         [Mkmc, Mkbc] = SIMfactor_Mk(weld, 0.15, T, L, tw, phi, radius);
    else
        % Use 3D solution
       [Mkma, Mkba, Mkmc, Mkbc] = SIMfactor_Mk_3D(a, c, T, L, phi, radius);
    end
    [M, Mm, Mb, fw] = SIMfactor(flaw, a, c, W, T, theta, hole);
    if theta == pi/2
        Ym = M * fw * Mm * Mkma;
        Yb = M * fw * Mb * Mkba;
    elseif theta == 0
        Ym = M * fw * Mm * Mkmc;
        Yb = M * fw * Mb * Mkbc;
    end
end
end
    
function [M, Mm, Mb, fw] = SIMfactor(flaw, a, c, W, T, theta, hole)
%
% Stress intensity factor acc. to Annex M in BS 7910:2013
% 1 = Through-thickness flaws in plates
% 2 = Surface flaws in plates (= flaw 1 in mathcad)
% 3 = Through-width flaws in plates
% 4 = semi-circular surface flaw in round bar
% 5 = Corner flaws in plates
% 6 = Semi-circular surface flaw in bolts, solution 2
% 7 = Edge flaws in plates
% 8 = Through-thickness deck plate crack (orthotropic bridge deck)
    %% M.3.1 Through-thickness flaws in plates
if flaw == 1
    M = 1;
    Mm = 1;
    Mb = 1;
    fw = (sec((pi*a)/W))^0.5;
    %% M.4.2 Extended surface flaws in plates
elseif flaw == 3
    M = 1;
    fw = 1;
    if a/T > 0.6
        disp ('Outside validity range in SIMfactor: a/T > 0.6');
        return
    end
    Mm = 1.12 - 0.23*(a/T) + 10.16*(a/T)^2 -21.7*(a/T)^3 + 30.4*(a/T)^4;
    Mb = 1.12 - 1.39*(a/T) + 7.32*(a/T)^2 -13.1*(a/T)^3 +14*(a/T)^4;
    %% M.10.2 Semi-circular surface flaws in round bars
elseif flaw == 4
    M = 1;
    fw = 1;
    pir = pi*a/(4*radius);
    if a/2*radius > 0.6
        disp ('Outside validity range in SIMfactor: a/2r > 0.6');
        return
    end
    g = ((1.84/pi) * (tan(pir)/(pir))^0.5) / cos(pir);
    Mm = g * (0.752 + 2.02 * (a/(2*radius)) + 0.37*(1 - sin(pir))^3);
    Mb = g * (0.923 + 0.199 * (1-sin(pir))^4);
    %% M.3.2 Edge flaws in plates
elseif flaw == 7
    M = 1;
    Mm = 1.12-0.23*(a/W) + 10.16*(a/W)^2-21.7*(a/W)^3 + 30.4*(a/W)^4;
    Mb = Mm;
    fw = 1;
    %% M.4.1 Surface flaws in plates
elseif flaw == 2
     M = 1;
    fw = (sec((pi*c/W)*(a/T)^0.5))^0.5;
    % Validity conditions
    if a < 0
        disp('Incorrect input in SIMfactor: a<0');
        return
    end
    if c < 0
        disp('Incorrect input in SIMfactor: c<0');
        return
    end
    if a/c <= 0.2
        if a/T > 1.25*(a/c + 0.6);
            disp('Outside validity range in SIMfactor: a/2c<0.1');
            return
        end
    elseif (a/c > 0.2) && (a/c <= 2)
        if a/T >= 1
            disp('Outside validity range in SIMfactor: a/2c>0.1');
            return
        end
    end
    if a/c > 2
        as = a;
        cs = a/2;
    else 
        as = a;
        cs = c;
    end
    % Membrane loading Mm:
    if as/cs <= 1
        M1 = 1.13 - 0.09*(as/cs);
        M2 = (0.89/(0.2 + as/cs))-0.54;
        M3 = 0.5 - 1/(0.65 + as/cs) + 14*(1-as/cs)^24;
        g  = 1 + (0.1 + 0.35*(as/T)^2) * (1 - sin(theta))^2;
        ftheta = ((as/cs)^2 * (cos(theta))^2 + (sin(theta))^2)^0.25;
    elseif as/cs > 1
        M1 = (cs/as)^0.5*(1 + 0.04*(cs/as));
        M2 = 0.2*(cs/as)^4;
        M3 = -0.11*(cs/as)^4;
        g = 1 + (0.1 + 0.35*(cs/as)*(a/T)^2) * (1 - sin(theta))^2;
        ftheta = ((cs/as)^2 * (sin(theta))^2 + (cos(theta))^2)^0.25;
    end
    PHI = EllipticalIntegral(as, cs);
    Mm = (M1 + M2*(as/T)^2 + M3*(as/T)^4) * g * ftheta / PHI;
    
    % Bending loading Mb:
    if as/cs <= 1
        q  = 0.2 + (as/cs) + 0.6*(as/T);
        H1 = 1 - 0.34*(as/T) - 0.11*(as/cs)*(as/T);
        G1 = -1.22 - 0.12*(as/cs);
        G2 = 0.55 - 1.05*(as/cs)^0.75 + 0.47*(as/cs)^1.5;
    elseif as/cs > 1
        q  = 0.2 + (cs/as) + 0.6*(as/T);
        H1 = 1 - (0.04 + 0.41*(cs/as))*(as/T) + (0.55 - 1.93*(cs/as)^0.75 + 1.38*(cs/as)^1.5)*(as/T)^2;
        G1 = -2.11 + 0.77*(cs/as);
        G2 = 0.55 - 0.72*(cs/as)^0.75 + 0.14*(cs/as)^1.5;
    end
    H2 = 1 + G1*(as/T) + G2*(as/T)^2;
    H  = H1 + (H2 - H1)*(sin(theta))^q;
    Mb = H * Mm;
    
    %% Correction for deep cracks
    xrma = 5.738; xsma = -9.87; xrba = 5.431; xsba = -8.974;
    yrma = 0.3283; ysma = -0.4710; yrba = 0.1982; ysba = -0.2629;
    zrma = -0.02371; zsma = 0.04587; zrba = -0.007021; zsba = 0.014;
    
    xrmc = -1.587; xsmc = 0.249; xrbc = -1.029; xsbc = -0.196;
    yrmc = -0.1001; ysmc = -0.06314; yrbc = -0.07077; ysbc = -0.03291;
    zrmc = 0.006702; zsmc = -0.007079; zrbc = 0.003178; zsbc = -0.002074;
    
    Rma = xrma + yrma * (c/a -15) + zrma*(c/a - 15)^2;
    Sma = xsma + ysma * (c/a -15) + zsma*(c/a - 15)^2;
    Rba = xrba + yrba * (c/a -15) + zrba*(c/a - 15)^2;
    Sba = xsba + ysba * (c/a -15) + zsba*(c/a - 15)^2;
    Rmc = xrmc + yrmc * (c/a -15) + zrmc*(c/a - 15)^2;
    Smc = xsmc + ysmc * (c/a -15) + zsmc*(c/a - 15)^2;
    Rbc = xrbc + yrbc * (c/a -15) + zrbc*(c/a - 15)^2;
    Sbc = xsbc + ysbc * (c/a -15) + zsbc*(c/a - 15)^2;
    if a/T < 0.25
        Mm_cor = 0;
        Mb_cor = 0;
    else
        if theta == pi/2
            Mm_cor = Rma * (a/T - 0.25)^2 + Sma * (a/T - 0.25)^3;
            Mb_cor = Rba * (a/T - 0.25)^2 + Sba * (a/T - 0.25)^3;
        elseif theta == 0
            Mm_cor = Rmc * (a/T - 0.25)^2 + Smc * (a/T - 0.25)^3;
            Mb_cor = Rbc * (a/T - 0.25)^2 + Sbc * (a/T - 0.25)^3;
        end
    end
    Mm = Mm + Mm_cor;
    Mb = Mb + Mb_cor;    
        %% M.5.1 Corner flaws in plates
elseif flaw == 5 
    M = 1;
    lambda = (c/W)*sqrt(a/T);
    if c/W <= 0.5
        fw = 1 - 0.2*lambda + 0.94*lambda^2 - 19.4*lambda^3 + 27.1*lambda^4;
    else
        disp('Outside validity range in SIMfactor (flaw 5): c/W > 0.5');
        return
    end
    % Validity conditions
    if a < 0
        disp('Incorrect input in SIMfactor: a<0');
        return
    end
    if c < 0
        disp('Incorrect input in SIMfactor: c<0');
        return
    end
    if a/c <= 0.2
        as = c * 0.2;
        cs = c;
    elseif a/c > 2
        as = a;
        cs = a/2;
    else 
        as = a;
        cs = c;
    end
   if a/T > 1
       disp('Outside validity range in SIMfactor flaw 5: a/T>1');
       return;
   end
   if theta < 0
      disp('Incorrect input in SIMfactor:  theta<0');
      return
   elseif theta > pi/2
      disp('Incorrect input in SIMfactor:  theta>pi/2');
      return
   end
   % Membrane loading
   if as/cs <= 1
       M1 = 1.08 - 0.03*(a/c);
       M2 = (1.06/(0.3 + a/c)) - 0.44;
       M3 = -0.5 + 0.25*(a/c) + 14.8*(1 - a/c)^15;
       g1 = 1 + (0.08 + 0.4*(a/T)^2) * (1 - sin(theta))^3;
       g2 = 1 + (0.08 + 0.15*(a/T)^2) * (1 - cos(theta))^3;
       ftheta = ((a/c)^2 * (cos(theta))^2 + (sin(theta))^2)^0.25;
   elseif as/cs > 1
       M1 = (1.08 - 0.03*(a/c)) * (c/a)^0.5;
       M2 = 0.375 * (c/a)^2;
       M3 = -0.25 * (c/a)^2;
       g1 = 1 + (0.08 + 0.4*(c/T)^2) * (1 - sin(theta))^3;
       g2 = 1 + (0.08 + 0.15*(c/T)^2) * (1 - cos(theta))^3;
       ftheta = ((c/a)^2 * (sin(theta))^2 + (cos(theta))^2)^0.25;
   end
   PHI = EllipticalIntegral(as, cs);
   Mm = (M1 + M2*(as/T)^2 + M3*(as/T)^4)*g1*g2*ftheta/PHI;
   % Bending loading
   if as/cs <= 1
       q = 0.2 + (a/c) + 0.6*(a/T);
       H1 = 1 - 0.34*(a/T) - 0.11*(a/c)*(a/T);
       G1 = -1.22 - 0.12*(a/c);
       G2 = 0.64 - 1.05*(a/c)^0.75 + 0.47*(a/c)^1.5;
   elseif as/cs > 1
       q = 0.2 + (c/a) + 0.6*(a/T);
       H1 = 1 - (0.04 + 0.41*(c/a)) * (a/T) + (0.55 - 1.93*(c/a)^0.75 + 1.38*(c/a)^1.5) * (a/T)^2;
       G1 = -2.11 + 0.77*(c/a);
       G2 = 0.64 - 0.72*(c/a)^0.75 + 0.14*(c/a)^1.5;
   end
   H2 = 1 + G1*(a/T) + G2*(a/T)^2;
   H = H1 + (H2 - H1)*(sin(theta))^q;
   Mb = H*Mm;
   %% Through-thickness deck plate crack (orthotropic bridge deck)
elseif flaw == 8
    M = 1;
    Mb = 0.375;
    Mm = 1;
   %% Symmetrical crack emanating from hole
elseif flaw == 9
    lambda = 2*(hole/2 + a)/W;
    F1 = (1 - 0.025 * lambda^2 + 0.06*lambda^4)*1/sqrt(cos(pi*lambda/2));
    aR = a*2/hole;
    F2 = 1 + 1/(2*aR^2 + 1.93 * aR + 0.539) + 1/(2*aR+2);
    RW = hole / W;
    F3 = 1 - 0.08 * RW + RW^2;
    M = F1 * F2 * F3;
    Mm = 1;
    Mb = 1;
    fw = 1;           
end
%%
if Mm < 0
    Mm = 0;
end
if Mb < 0
	Mb = 0;
end
    %% Modified stress intensity concentrations factors for deep cracks
end
    
            
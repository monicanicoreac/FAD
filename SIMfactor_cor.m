function [Mm_cor, Mb_cor] = SIMfactor_cor(a, c, T)
% Modified stress intensity concentration factors for deep cracks
% X / Y / Z = [r_ma s_ma, r_ba, s_ba;...
%              r_mc, s_mc, r_bc, r_bc];
% Mm_cor = [Mma_cor; Mmc_cor];
% Mb_cor = [Mba_cor; Mbc_cor];
X = [5.738, -9.87, 5.431, -8.974;...
    -1.587, 0.249, -1.029, -0.196];
Y = [0.3282, -0.471, 0.1982, -0.2629;...
    -0.1001, -0.06314, -0.07077, -0.03291];
Z = [-0.02371, 0.04587, -0.007021, 0.014;...
    0.006702, -0.007079, 0.003178, -0.002074];

RS(1,:) = X(1,:) + Y(1,:)*(c/a - 15) + Z(1,:)*(c/a - 15)^2;
RS(2,:) = X(2,:) + Y(2,:)*(c/a - 15) + Z(2,:)*(c/a - 15)^2;

if a/T < 0.25
    Mm_cor(1:2,1) = 0;
    Mb_cor(1:2,1) = 0;
else
    Mm_cor = RS(:,1)*(a/T - 0.25)^2 + RS(:,2)*(a/T - 0.25)^3;
    Mb_cor = RS(:,3)*(a/T - 0.25)^2 + RS(:,4)*(a/T - 0.25)^3;
end
end

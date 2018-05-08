function PHI = EllipticalIntegral(a,c)
% Elliptical integral of the second kind
if a/c  <= 1
    PHI = (1 + 1.464*(a/c)^1.65)^0.5;
elseif (a/c > 1) && (a/c <= 2)
    PHI = (1 + 1.464*(c/a)^1.65)^0.5;
end
end
function [Kmat, Kmati] = FractureToughnessOSK(fyT, fuT, E, niu, ctod, Bctod, Lcrack, direction, Tstr)
% Material toughness measured by stress intensity factor
% Constraint factor: 7.1.4.6 in BS7910:2013 A1 2015
m_ct = 1.517 * (fyT/fuT)^(-0.3188);
% Fracture toughness from CTOD: 7.1.4.6 in BS7910:2013 A1 2015
% if condition == 1 % plane stress
%     Kmat0 = sqrt(m_ct * fyT * ctod * E);
% elseif condition == 2 % plane strain
    Kmat0 = sqrt(m_ct * fyT * ctod * E/(1 - niu^2));
% else
%     disp('Error in FractureToughnessOSK: condition has to be 1 (plane stress) or 2 (plane strain');
% end

KmatB = 632 + (Kmat0 - 632)*(Bctod/25)^0.25;
% Correction for three-dimensionality of stress
% Single edge-cracked plate: 3PB
% for LS: a/W=0.3
% for LT: a/W=0.5
if direction == 1 % LS direction
    betaCTOD = -0.7887 - 0.1795 * 0.3 +32.9014 * 0.3^2 - 153.45 * 0.3^3 + ...
    316.11 * 0.3^4 - 308.47 * 0.3^5 + 115.15 * 0.3^6;
%     ratio = 1 / (sqrt(3) * 1.22* (1-0.3)^2);
elseif direction == 2 % LT direction
    betaCTOD = -0.7887 - 0.1795 * 0.5 +32.9014 * 0.5^2 - 153.45 * 0.5^3 + ...
    316.11 * 0.5^4 - 308.47 * 0.5^5 + 115.15 * 0.5^6;
%     ratio = 1 / (sqrt(3) * 1.22* (1-0.5)^2);
else
  disp('Error in FractureToughnessOSK: direction has to be 1 (LS) or 2 (LT)');
end
% Based on earlier CTOD tests
C1 = 1.273;
C2 = 6.09;
SrefFy = C1 * exp(-(C2/(KmatB/sqrt(fyT))));
SrefCTOD = fyT * SrefFy;
% SrefCTOD = ratio*fyT;
T_stressCTOD = min(0,betaCTOD*SrefCTOD);
KmatCTOD = 948 + (KmatB - 948)*exp(0.019*(T_stressCTOD/10));

% Adjustment for specimen size and crack front length:
% L.9.2 and L.9.5 in BS7910:2013 A1 2015
% Find the maximum crack front length corresponding to the plane strain
% situation:
fun     = @(x) 632+(KmatCTOD - 632)*(25/x)^0.25 - fyT*(x/2.5)^0.5;

% constrained expansion of the fzero search interval
int     = [1e-2, 200];
ii      = 1;
while sign(fun(int(1))) == sign(fun(int(2))) && ii < 1e3
    int     = [int(1), int(2)*1.2];
    ii      = ii + 1;
end
Lmax    = fzero(fun, int);


% % polyroot
% aa  = 632;
% bb  = (Kmat0 - 632)*25^0.25;
% cc  = -fyT*(1/2.5)^0.5;
% r   = roots([bb, aa, 0, cc]);
% rr  = (r).^(-1/0.25);
% Lmax = rr(imag(rr) == 0 & real(rr) > 0);
% Lmax = min(Lmax);

Lcrack = min(Lcrack, Lmax);

Kmati = 632 + (316 + (KmatCTOD - 948)*max(1, exp(0.019*(-Tstr/10)))) * (25 / 25)^0.25;
Kmat = 632 + (316 + (KmatCTOD - 948)*max(1, exp(0.019*(-Tstr/10)))) * (25 / Lcrack)^0.25;

%{
% Annex J in BS 7910:2013: Use of Charpy V-notch impact tests to estimate
% fracture toughness
% J.2.1 If Charpy results are available at the temperature at which Kmat is
% required, use equation J.1
Kmat = (((12*sqrt(Charpy)-20)*(25/T)^0.25)+20)*sqrt(1000);
% J.2.4 To avoid overestimating fracture toughness at the service temperature in
% materials with potentially low upper shelf Charpy energy, Kmat should not
% exceed the value from equation J.5
Kmat = min(Kmat, (0.54*Charpy+55)*sqrt(1000));
%}
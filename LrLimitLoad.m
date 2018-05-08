function [ Pcchord, Ptchord, Mcichord, Mcochord, Pcbrace, Mcibrace, Mcobrace ] = LrLimitLoad( a, c, GeometryTub, LoadsTub, fy, JointType )
% P limit load solution (Annex P.14 - B.S.7910:2013 A1 2015)
% Calculate characteristic strength of welded tubular joint
lambdaPLN  = 0.030;   % for brace axial load
lambdaPLMi = 0.045;   % for brace in-plane moment load
lambdaPLMo = 0.021;   % for brace out-of-plane moment load
gammaPL    = GeometryTub(2,1)/(2*GeometryTub(1,1));
Ka         = 0.5*(1+1/sin(GeometryTub(5,1)*pi/180));
% Chord utilisation factor
UPL    = sqrt((0.23*LoadsTub(1,1)*GeometryTub(2,1))^2 + LoadsTub(2,1)^2 + LoadsTub(3,1)^2)/(0.72*GeometryTub(1,1)*fy*GeometryTub(2,1)^2); % this is a function of the applied loads in the chord
cond1  = sqrt(LoadsTub(2,1)^2 + LoadsTub(3,1)^2)/(23.*GeometryTub(2,1));

% Qf is a factor that allows for the presence of axial and moment loads in the chord
if (LoadsTub(1,1) >= cond1) 
	QfextrN  = 1;     % for extreme conditions
	QfextrMi = 1;
	QfextrMo = 1;
	QfopN  = 1;	%for operating conditions
	QfopMi = 1;
	QfopMo = 1;
else
	QfextrN  = (1-1.638*lambdaPLN*gammaPL*UPL^2);
	QfextrMi = (1-1.638*lambdaPLMi*gammaPL*UPL^2);
	QfextrMo = (1-1.638*lambdaPLMo*gammaPL*UPL^2);
	QfopN   = (1-2.89*lambdaPLN*gammaPL*UPL^2);
	QfopMi  = (1-2.89*lambdaPLMi*gammaPL*UPL^2);
	QfopMo  = (1-2.89*lambdaPLMo*gammaPL*UPL^2);
end
betaTub = GeometryTub(4,1)/GeometryTub(2,1);
if (betaTub <= 0.6)
	Qbeta = 1;
else
	Qbeta = 0.3/(betaTub*(1-0.833*betaTub));
end
% Coefficient Qu is a strength factor that varies with loading type and joint configuration
%
if (JointType == 1)
    QuNc = (2 + 20*betaTub)*sqrt(Qbeta);
elseif (JointType == 2)
    Qg = 1.7 - 0.9*(gap/GeometryTub(4,1))^0.5;
    Qg = max(Qg,1);
    QuNc = (2. + 20.*betaTub)*sqrt(Qbeta)*Qg;
elseif (JointType == 3)
    QuNc = (2.5+14.*betaTub)*Qbeta;
else
    disp('Incorrect input in LrLimitLoad')
return
end
%
if (JointType == 1)
    QuNt = (8. + 22.*betaTub);
elseif (JointType == 2) 
    QuNt = (8. + 22.*betaTub)*Qg;
elseif(JointType == 3) 
    QuNt = (7. + 17.*betaTub)*Qbeta;
else
    disp('Incorrect input in LrLimitLoad')
return
end
%
	QuMi = 5*betaTub*(gammaPL^0.5)*sin(GeometryTub(5,1)*pi/180);
%
if (JointType == 1)
    QuMo = (1.6 + 7.*betaTub)*Qbeta;
elseif (JointType == 2) 
    QuMo = (1.6 + 7.*betaTub)*Qbeta;
elseif (JointType == 3) 
    QuMo = (1.6 + 7.*betaTub*sqrt(Qbeta));
else
    disp('Incorrect input in LrLimitLoad')
    return
end

% Characteristic strength
%
Pkc  = QuNc*QfextrN*(fy*Ka*GeometryTub(1,1)^2)/(sin(GeometryTub(5,1)*pi/180));
Pkt  = QuNt*QfextrN*(fy*Ka*GeometryTub(1,1)^2)/(sin(GeometryTub(5,1)*pi/180));
Mki = QuMi*QfextrMi*(fy*GeometryTub(4,1)*GeometryTub(1,1)^2)/(sin(GeometryTub(5,1)*pi/180));
Mko = QuMo*QfextrMo*(fy*GeometryTub(4,1)*GeometryTub(1,1)^2)/(sin(GeometryTub(5,1)*pi/180));
% Correction factor for axial loading for cracked geometry
anorm = a / GeometryTub(1,1);
if (anorm < 1) 
	Ac = 0.5*pi*a*c;        % crack area for surface breaking flaw
else
	Ac = 2*a*GeometryTub(1,1);		% crack are for trough thickness flaw
end

% Lenght of weld 
hh = (pi/2-0)/100;
f1 = sqrt((sin(pi)^2)/(sin(GeometryTub(5,1)*pi/180)^2) + (cos(pi)^2)/((1-sin(pi)*betaTub^2)^2));
f2 = sqrt((sin(0)^2)/(sin(GeometryTub(5,1)*pi/180)^2) + (cos(0)^2)/((1-sin(0)*betaTub^2)^2));
total = 0.5 * (f1 + f2);
for i = 1: 99
	total = total + sqrt((sin(0 + i*hh)^2)/(sin(GeometryTub(5,1)*pi/180)^2) + (cos(0 + i*hh)^2)/((1-sin(0 + i*hh)*betaTub^2)^2));
end

Lweld = 2*betaTub * GeometryTub(2,1) * hh * total;

if (anorm < 1)
	mq = 0;	% for tubular joint containing part thickness flaws
else
	mq = 1;
end
% Reduction factor to allow for the loss of load-bearing cross-sectional area due to the presence of flaw
% For tubular joints containing a flaw in chord
FAR  = (1-(Ac/(GeometryTub(1,1)*Lweld)))*(1/Qbeta)^mq;    
%
PI=4*atan(1);
fct = Lweld / (PI*GeometryTub(4,1));
PHIC = 2 * c * fct / GeometryTub(4,1) ;                    % angle around the crack length 2c
FBR = cos(PHIC/2)*(1-sin(PHIC/2));              % for in plane and out of plane bending

Pcchord  = Pkc*FAR;
Ptchord  = Pkt*FAR;
Mcichord = Mki * FBR;
Mcochord = Mko * FBR;

% For T and K tubular joints containing a flaw in the brace
GSBbrace = 0.25 * pi * (GeometryTub(4,1)^2 -(GeometryTub(4,1) - 2*GeometryTub(3,1))^2);
Wplbrace = (GeometryTub(4,1)^3 - (GeometryTub(4,1) - 2 * GeometryTub(3,1))^3)/6;
Anet = GSBbrace - a * 2 * c;
if (JointType == 1)
	Pcbrace  = fy * Anet;
	Mcibrace = FBR * Wplbrace * fy;
	Mcobrace = FBR * Wplbrace * fy;
else
	disp( 'Incorrect input in subroutine LrLimitLoad: procedure for flaw in brace is implemented only for T and Y joints');
    return
end

end


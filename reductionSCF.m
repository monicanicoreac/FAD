function [SCF, PhiSCF] = reductionSCF(SCFp, c, Lscf)
% Cosinus reduction is used over an influence length Lscf
% SCFp is the stress concetration factor at hotspot 
% Lscf is the influence length

if AND((c>=-Lscf),(c<=Lscf))
    SCF = ((SCFp-1)/2)*(cos(PI*c)/Lscf+1)+1;
elseif c>Lscf
    SCF = 1;
elseif c<Lscf
    SCF = 1;
end
PhiSCF = SCF/SCFp;
end

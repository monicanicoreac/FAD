% function N = FractureAnalysis%, T, ai, ci, af, flaw, weld, tw, phi, L, rho, Charpy, sigmaSdm, CGL)
% global T  ai  ci  af  flaw  weld  tw  phi  L  rho  sigmaSdm CGL ecc
% Units: mm, N, C
clear all; clc;
global index ctodLT_tub ctodLS_tub a c T W flaw E niu

% ---------------------------------------------------------------------
% CONTROL & OPTIONS
% ---------------------------------------------------------------------
row_num         = 1;
shift_row       = 1;
specimen_type   = 'tubular';
% specimen_treat  = 'per_specimen';
% % specimen_treat  = 'average_ctod';
specimen_treat  = 'magic_average_ctod';
fname           = [specimen_type, '_', specimen_treat, '_to_Arpi.xlsx'];
if exist(fname, 'file')
    delete(fname)
end

[batch_name, CTOD, Lr_trunc_flag, Lr_sup] = get_ctod_data(specimen_type);

% get rid of the first specimen type, A
% idx                     = 1;
% batch_name(idx)         = [];
% CTOD(idx,:)             = [];
% Lr_trunc_flag(idx,:)    = [];
% Lr_sup(idx,:)           = [];

switch lower(specimen_treat)
    case 'per_specimen'
        % the calculation is completed for each specimen individually
        % do nothing
    case 'average_ctod'
        % proper averaging the ctod values
        CTOD = mean(CTOD,2,'omitnan');
    case 'magic_average_ctod'
        % who knows how the "averaged" ctod values are calculated
        CTOD = [0.17; 0.2100;    0.1530;    1.1100;   0.4700];
        CTOD = CTOD(:); % just in case
    otherwise
       error(['Unknown specimen treatment (specimen_treat): ', specimen_treat]) 
end

load('tubular.mat')

% ctodLT_tub(idx,:) = [];
% ctodLS_tub(idx,:) = [];
% forces_tub(idx,:) = [];
% ratio(idx,:)      = [];
% SCF_tub(idx,:)    = [];
% sigmaSd(idx,:)    = [];
% tubular(idx,:)    = [];

n_batch = length(batch_name);

% loop over the batches
for ii = 1:n_batch
    
    index = ii;
    
    % ctodLT = ctodLT_tub(index,1);
    % ctodLS = ctodLS_tub(index,1);
    
    % ---------------------------------------------------------------------
    % GEOMETRY INPUT PARAMETERS
    % ---------------------------------------------------------------------
    JointType = 1;
    % 1 = T and Y joints
    % 2 = K joints
    % 3 = X and DT joints
    % Chord thickness
    T = tubular(index,1);
    Tchord = tubular(index,1);
    % Chord diameter
    Dchord = tubular(index,2);
    % Brace thickness
    tbrace = tubular(index,3);
    % Diameter brace
    dbrace = tubular(index,4);
    % Angle between chord and brace
    thetachbr = tubular(index,14);
    % Gap between chords (for K and X joints)
    gap = tubular(index,15);
    % Plate width:
    W = tubular(index,6);
    % Initial defect size:
    a_i = tubular(index,7);
    c_i = tubular(index,8);
    % Final crack depth:
    a_f = T;
    % Flaw type:
    flaw = 2;
    % 1 = Through-thickness flaws in plates
    % 2 = Semi elliptical surface flaws in plates (= flaw 1 in mathcad)
    % 3 = Through-width flaws in plates
    % 4 = semi-circular surface flaw in round bar
    % 5 = Corner flaws in plates
    % 6 = Semi-circular surface flaw in bolts, solution 2
    % 7 = Edge flaws in plates
    % 8 = Through-thickness deck plate crack (orthotropic bridge deck)
    % 9 = symmetric crack emanating from hole
    % Type of weld
    weld            = tubular(index,9);
    % Weld angle
    Wangle          = tubular(index,10);
    % Attachment lenght
    L               = 1.5*T;
    % Weld toe radius
    radius          = 0;
    % Weld throat
    tw              = 1.2*tbrace;
    crackLocation   = 'chord';
    
    % ---------------------------------------------------------------------
    % MATERIAL INPUT PARAMETERS
    % ---------------------------------------------------------------------
    
    % Yield stress
    fy              = tubular(index,11);
    % Ultimate tensile strength
    fu              = tubular(index,12);
    % Modulus of elasticity
    E               = 210000;
    % Poisson ratio:
    niu             = 0.3;
    % Width CTOD specimen
    Bctod           = T-1;
    % Crack front lenght:
    Lcrack          = 2 * c_i;
    % Temperature assessment:
    Temp_a          = 0;
    % CTOD temperature:
    Temp_ctod       = 0;
    % Assessment temperature: 7.1.3.4 in BS7910:2013 A1 2015
    if Temp_a == Temp_ctod
        fyT = fy;
        fuT = fu;
        % elseif (Temp_a < 20) % room temperature 20C
        %     fyT = fy + 10E5 / (1.8 * Temp_a + 491) - 189;
        %     fuT = fu * (0.7857 + 0.2423 * exp(-Temp_a/170.646));
        % else
        %     fyT = fy;
        %     fuT = fu;
    end
    
    % ---------------------------------------------------------------------
    % PERMANENT STRESS ULS
    % ---------------------------------------------------------------------
    % Permanent stresses:
    sigmaSdm        = sigmaSd(index,1) * ratio(index,1)/(1 + ratio(index,1));
    sigmaSdb        = sigmaSd(index,1) * 1/(1 + ratio(index,1));
    % Forces in tubular joint: Pa axial tension/compression, Mai in plane
    % bending moment, Mao out-of-plane bending moment
    Pa              = forces_tub(index,1);
    Mai             = forces_tub(index,2);
    Mao             = forces_tub(index,3);
    Pachord         = forces_tub(index,4);
    Maichord        = forces_tub(index,5);
    Maochord        = forces_tub(index,6);
    %
    GeometryTub     = [Tchord; Dchord; tbrace; dbrace; thetachbr; gap];  % geometry of tubular joint
    LoadsTub        = [Pachord*1E3; Maichord*1E6; Maochord*1E6; Pa*1E3; Mai*1E6; Mao*1E6];  % forces in chord and brace
    % ---------------------------------------------------------------------
    % FRACTURE ANALYSIS
    % ---------------------------------------------------------------------
    a               = a_i;
    c               = c_i;
    
    n_specimen      = sum(~isnan(CTOD(ii,:)));
    specimen        = 1;
    Lr_tf           = Lr_trunc_flag(index);
    Lr_s            = Lr_sup(index);
    % write the header: column names
    if index == 1
        xlRange = [ 'A', num2str(1)];
%         xlswrite(fname, {'batch', 'specimen', 'CTOD', 'KI', 'Kmat', 'rho', 'Kr', 'Lr', 'Lr_trunc_flag', 'Lr_sup'}, 'Sheet1', xlRange);
    end
    
    % loop over the specimens in a given batch
    for jj = 1:n_specimen
        ctod            = CTOD(index,jj);
        ctodLT          = ctodLT_tub(index,1);
        ctodLS          = ctodLS_tub(index,1);
% %         ctodLT          = ctodLT_tub(index,1);
% %         ctodLS          = ctodLS_tub(index,1);
    
        % T-stress
        [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
        %
        [Kmat_a, Kmati_a] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLS, Bctod, Lcrack, 1, Tstr_a); % 2 is for plane strain
        [Kmat_c, Kmati_c] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_c); % 1 is for plane stress
        [Yma, Yba] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, 0);
        [Ymc, Ybc] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, 0, 0);
        [~, Mma, Mba, ~] = SIMfactor(flaw, a, c, W, T, pi/2, 0);
        [~, Mmc, Mbc, ~] = SIMfactor(flaw, a, c, W, T, 0, 0);
        % ---------------------------------------------------------------------
        % PLASTIC COLLAPSE PARAMETER Lr
        % ---------------------------------------------------------------------
        
        % Local collapse analysis:
        % Cut-off value of Lr (to prevent plastic colapse)
        Lrmax = (fyT+fuT)/(2*fyT);
        
        % P.5 in BS7910:2013 A1 2015 Reference stress and/or limit load solutions for deep notched plates
        
        % Global analysis:
        % Flow strenght:
        fflow = (fyT + fuT) / 2;
        [ Pcchord, Ptchord, Mcichord, Mcochord, Pcbrace, Mcibrace, Mcobrace ] = LrLimitLoad( a, c, GeometryTub, LoadsTub, fyT, JointType );
        if LoadsTub(4,1) >= 0 % tension force in brace
            Lrchord = (fflow / fyT) * (abs(LoadsTub(4,1)/Ptchord) + (LoadsTub(5,1) / Mcichord)^2 + abs(LoadsTub(6,1)/Mcochord));
        else  % compression force in brace
            Lrchord = (fflow / fyT) * (abs(LoadsTub(4,1)/Pcchord) + (LoadsTub(5,1) / Mcichord)^2 + abs(LoadsTub(6,1)/Mcochord));
        end
        Lrbrace = (fflow / fyT) * (abs(LoadsTub(4,1)/Pcbrace) + (LoadsTub(5,1) / Mcibrace)^2 + abs(LoadsTub(6,1)/Mcobrace));
        if (crackLocation == 'chord')
            Lr = Lrchord;
        elseif (crackLocation == 'brace')
            Lr = Lrbrace;
        end
        
        % ---------------------------------------------------------------------
        % FRACTURE COLLAPSE LoadsTub(4,1)rameter Kr
        % ---------------------------------------------------------------------
        
        % Annex Q in BS7910:2013 A1 2015 Residual stresses in as-welded structures:
        % and
        % 7.1.8.3 in BS7910:2013 A1 2015: Post-welded heat-treated (PWHT) structures
        %  a flaw lying in a plane LoadsTub(4,1)rallel to the welding direction (the stresses
        %  to be considered are perpendicular to the weld): 20% of the lesser of
        %  the yield strength of the weld or LoadsTub(4,1)rent material
        %  a flaw lying in a plane transverse to the welding direction (the stresses
        %  to be considered are LoadsTub(4,1)rallel to the weld): 30% of RT yield strength of
        % the material in which the flaw is located
        if (weld == 0)
            sigmaResM = 0;
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 11)
            % tubular T-butt weld (full penetration flaw), flaw perpendicular
            % to weld direction,as-welded (Annex Q)
            sigmaResM = 0.564 * fyT;
            sigmaResB = (0.307 - 0.614*a/GeometryTub(1,1)) * fyT;
            Ksb       = 0.3 * fyT * sqrt(GeometryTub(1,1));
        elseif (weld == 12)
            % tubular T-butt weld (full penetration flaw), flaw parallel to weld
            % direction,as-welded (Annex Q)
            sigmaResM = 0.496 * fyT;
            sigmaResB = (0.216 - 0.432*a/GeometryTub(1,1))*fyT;
            Ksb       = 0.21 * fyT * sqrt(GeometryTub(1,1));
        elseif (weld == 13)
            % full penetration flaw, flaw perpendicular to weld direction, PWHT
            sigmaResM = 0.30 * fyT;
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 14) || (weld == 24)
            % full penetration flaw, flaw LoadsTub(4,1)rallel to weld direction, PWHT
            sigmaResM = 0.20 * fyT;
            sigmaResB = 0;
            Ksb = 0;
        end
        
        % Annex R: Plasticity correction factor for secondary stresses
        KIp_a = (Yma * sigmaSdm + Yba * sigmaSdb)*sqrt(pi*a);
        KIp_c = (Ymc * sigmaSdm + Ybc * sigmaSdb)*sqrt(pi*a);
        
        KIs_a = (Mma * sigmaResM + Mba * sigmaResB) * sqrt(pi*a) + Ksb;
        KIs_c = (Mmc * sigmaResM + Mbc * sigmaResB) * sqrt(pi*a) + Ksb;
        
        chia = KIs_a * Lr / KIp_a;
        chic = KIs_c * Lr / KIp_c;
        
        if KIs_a <= 0
            rhoa = 0;
        elseif (chia > 4) && (KIs_a >0)
            disp('Error in Fracture: plasticity correction factor rhoa not implemnted for chia>4')
        elseif (chia <= 4) && (Lr <= 0.8) && (KIs_a >0)
            rhoa = 0.1*chia^0.714 - 0.007*chia^2 + 3e-5*chia^5;
        elseif (chia <= 4) && (Lr > 0.8) && (Lr < 1.05) && (KIs_a >0)
            rhoa = 4*(0.1*chia^0.714 - 0.007*chia^2 + 3e-5*chia^5)*(1.05 - Lr);
        elseif (chia <= 4) && (Lr > 1.05) && (KIs_a >0)
            rhoa = 0;
        end
        
        if KIs_c <= 0
            rhoc = 0;
        elseif (chic > 4) && (KIs_c >0)
            disp('Error in Fracture: plasticity correction factor rhoa not implemnted for chia>4')
        elseif (chic <= 4) && (Lr <= 0.8) && (KIs_c >0)
            rhoc = 0.1*chic^0.714 - 0.007*chic^2 + 3e-5*chic^5;
        elseif (chic <= 4) && (Lr > 0.8) && (Lr < 1.05) && (KIs_c >0)
            rhoc = 4*(0.1*chic^0.714 - 0.007*chic^2 + 3e-5*chic^5)*(1.05 - Lr);
        elseif (chic <= 4) && (Lr > 1.05) && (KIs_c >0)
            rhoc = 0;
        end
        %
        KI_a = KIp_a + KIs_a;
        KI_c = KIp_c + KIs_c;
        
        Kr_a = KI_a/Kmat_a + rhoa;
        Kr_c = KI_c/Kmati_c + rhoc;
        
        % Failure assessment diagram
        if fyT <= 1000
            deltaEpsilon = 0.037 * (1 - 0.001 * fyT);
        else
            disp('Error in Fracture: Yield strength haqs to lower dan 1000MLoadsTub(4,1)');
            return
        end
        lambda = max(1, 1 + E * deltaEpsilon/fyT);
        Nfrac  = 0.3 * (1 - fyT/fuT);
        if (Lr >= Lrmax)
            Krcr = 0;
        elseif (Lr < 1)
            Krcr = 1/sqrt(1 + 0.5*Lr^2);
        elseif (Lr == 1)
            Krcr = 1/sqrt(lamda + 0.5/lambda);
        elseif (Lr > 1) && (Lr < Lrmax)
            Krcr = 1/sqrt(lambda + 0.5/lambda) * Lr^((Nfrac-1)/(2*Nfrac));
        end
%         KI_a, KI_c, Kmat_a, Kmati_c, rhoa, rhoc, Kr_a, Kr_c, Lr
        
        rho                             = rhoa;
        Kr                              = Kr_a;
        KI                              = KI_a;
        Kmat                            = Kmat_a;
%         if Kr_a > Kr_c
%             rho                             = rhoa;
%             Kr                              = Kr_a;
%             KI                              = KI_a;
%             Kmat                            = Kmat_a;
%         else
%             rho                             = rhoc;
%             Kr                              = Kr_c;
%             KI                              = KI_c;
%             Kmat                            = Kmat_c;
%         end
        batch                           = batch_name{ii};
        xlRange                         = [ 'A', num2str(row_num + shift_row)];
%         xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
        row_num     = row_num + 1;
        specimen    = specimen + 1;
    end
OUT(index,1) = sigmaSdm+sigmaSdb;
OUT(index,2) = sigmaResM;
OUT(index,3) = sigmaResB;
OUT(index,4) = Ksb;
OUT(index,5) = KI;
OUT(index,6) = rho;
OUT(index,7) = Kmat;
OUT(index,8) = Kr;
OUT(index,9) = Lr;    
end
% function N = FractureAnalysis%, T, ai, ci, af, flaw, weld, tw, phi, L, rho, Charpy, sigmaSdm, CGL)
% global T  ai  ci  af  flaw  weld  tw  phi  L  rho  sigmaSdm CGL ecc
% Units: mm, N, C
clear all; clc;
global index ctod ctodLT ctodLS a c T W flaw E niu

% ---------------------------------------------------------------------
% CONTROL & OPTIONS
% ---------------------------------------------------------------------
% row shift in the excel sheet for saving the data (make space for
% metadata, e.g. header)
row_num         = 1;
shift_row       = 1;
specimen_type   = 'wide_plate_add';
% specimen_treat  = 'per_specimen';
specimen_treat  = 'average_ctod';
% specimen_treat  = 'magic_average_ctod';

fname           = [specimen_type, '_', specimen_treat, '_to_Arpi.xlsx'];
if exist(fname, 'file')
    delete(fname)
end

[batch_name, CTOD, Lr_trunc_flag, Lr_sup] = get_ctod_data(specimen_type);

switch lower(specimen_treat)
    case 'per_specimen'
        % the calculation is completed for each specimen individually
        % do nothing
    case 'average_ctod'
        % proper averaging the ctod values
        CTOD = mean(CTOD,2,'omitnan');
    case 'magic_average_ctod'
        % who knows how the "averaged" ctod values are calculated
        % nope
%         load('ctod.mat')
%         CTOD = ctod;
    otherwise
       error(['Unknown specimen treatment (specimen_treat): ', specimen_treat]) 
end

n_batch = length(batch_name);

% loop over the batches
for ii = 1:n_batch
    load('NIL_add.mat')
    index       = ii;
    FADtype = 2; % curve suitable for materials that exhibit a yield discontinuity
    % ---------------------------------------------------------------------
    % GEOMETRY INPUT PARAMETERS
    % ---------------------------------------------------------------------
    % Plate thickness
    T           = NIL_add(index,2);
    % Plate width:
    W           = NIL_add(index,3);
    % Initial defect size:
    a_i         = NIL_add(index,4);
    c_i         = NIL_add(index,5);
    % Final crack depth:
    a_f         = T;
    % Flaw type:
    flaw        = NIL_add(index,10);
    % 1 = Through-thickness flaws in plates
    % 2 = Semi elliptical surface flaws in plates (= flaw 1 in mathcad)
    % 3 = Through-width flaws in plates
    % 4 = semi-circular surface flaw in round bar
    % 5 = Corner flaws in plates
    % 6 = Semi-circular surface flaw in bolts, solution 2
    % 7 = Edge flaws in plates
    % 8 = Through-thickness deck plate crack (orthotropic bridge deck)
    % 9 = symmetric crack emanating from hole
    hole        = 0;
    % Type of weld
    weld        = NIL_add(index,11);
    % Weld angle
    Wangle      = NIL_add(index,16);
    % Attachment lenght
    L           = NIL_add(index,15);
    % Weld toe radius
    radius      = 0;
    % Weld throat
    tw          = 0;
    crackLocation = 'chord';

    % ---------------------------------------------------------------------
    % MATERIAL INPUT PARAMETERS
    % ---------------------------------------------------------------------
    % Yield stress
    fy          = NIL_add(index,6);
    fyweld      = NIL_add(index,17);
    % Ultimate tensile strength
    fu          = NIL_add(index,7);
    fuweld      = NIL_add(index,18);
    % Modulus of elasticity
    E           = NIL_add(index,12);
    % Poisson ratio:
    niu         = 0.3;
    % Width CTOD specimen
    Bctod       = NIL_add(index,13);
    % Crack front lenght:
    Lcrack      = NIL_add(index,14);
    % Temperature assessment:
    Temp_a      = NIL_add(index,1);
    % CTOD temperature:
    Temp_ctod   = NIL_add(index,1);
    % Assessment temperature: 7.1.3.4 in BS7910:2013 A1 2015
    % The strength parameters are determined at RT, the assessment
    % temperature is different than RT
    fyT = fy + 1E5 / (1.8 * Temp_a + 491) - 189;
    fuT = fu * (0.7857 + 0.2423 * exp(-Temp_a/170.646));
    fyweldT = fyweld + 1E5 / (1.8 * Temp_a + 491) - 189;
    fuweldT = fuweld * (0.7857 + 0.2423 * exp(-Temp_a/170.646));
    % Increase in strain (eq. 8 in BS7910)
    DeltaEpsilon = 0.0375 * (1 - 0.001 * fyT);
    % ---------------------------------------------------------------------
    % PERMANENT STRESS ULS
    % --------------------------------------------------------------------- 
    sigmaSdm = NIL_add(index,8)*1E3 / (W*T);
    sigmaSdb = 0;
    
    % ---------------------------------------------------------------------
    % FRACTURE ANALYSIS
    % ---------------------------------------------------------------------
    a = a_i;
    c = c_i;
    % write the header: column names
    if index == 1
        xlRange = [ 'A', num2str(1)];
%         xlswrite('NIL_add.xlsx', {'Kr_a', 'Kr_c', 'Lr', 'KI_a', 'KI_c', 'Kmat_a', 'Kmati_c', 'rhoa', 'rhoc'}, 'Sheet1', xlRange);
        xlswrite(fname, {'batch', 'specimen', 'E', 'fy', 'fu', '??', 'FADtype', 'Kr', 'Lr'}, 'Sheet1', xlRange);
        xlswrite(fname, {'batch', 'specimen', 'CTOD', 'KI', 'Kmat', 'rho', 'Kr', 'Lr', 'Lr_trunc_flag', 'Lr_sup'}, 'Sheet2', xlRange);
    end
    
    n_specimen      = sum(~isnan(CTOD(ii,:)));
    specimen        = 1;
    Lr_tf           = Lr_trunc_flag(index);
    Lr_s            = Lr_sup(index);
    % loop over the specimens in a given batch
    for jj = 1:n_specimen
        disp('')
        ctod    = CTOD(ii,jj);
%         ctodLT  = ctod % WARNING!!
%         ctodLS  = ctod;
        ctodLT  = ctod; % WARNING!!
        ctodLS  = ctod;
        % T-stress
        [ Tstr_a, Tstr_c ]              = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
        [Kmat, Kmati]                   = FractureToughnessOSK(fyweldT, fuweldT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_a);
        [Yma, Yba]                      = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
        [Ma, Mma, Mba, ~]               = SIMfactor(flaw, a, c, W, T, pi/2, 0);
        %% Plastic collape parameter Lr 
        % Local collapse analysis:
        % Cut-off value of Lr (to prevent plastic colapse)
        if fyweldT == 0  % no weld
            Lrmax = (fyT+fuT)/(2*fyT);
        else
            Lrmax = (fyweldT+fuweldT)/(2*fyweldT);
        end
        % P.5 in BS7910:2013 A1 2015 Reference stress and/or limit load solutions for deep notched plates
        constraint = 2;
        if flaw == 10
            % plane strain
            if (2*a/W <= 0.884) && (2*a/W >= 0)
                beta1 = 1 + log((W - a)/(W - 2*a));
            elseif (2*a/W < 1) && (2*a/W > 0.884)
                beta1 = 1 + PI()/2;
            end
            FeNpl_strain = beta1 * 4 * T * (0.5*W - a) * fyweldT / sqrt(3);
            % plane stress
            if (2*a/W <= 0.286) && (2*a/W >= 0)
                beta2 = 1 + 0.54*2*a/W;
            elseif (2*a/W < 1) && (2*a/W > 0.286)
                beta2 = 2/sqrt(3);
            end
            FeNpl_stress = beta2 * 2 * T * (0.5*W - a) * fyweldT;
            % Strength mis-match
            M = fyweldT / fyT; % mis-match ratio
            psi = (0.5*W - a) / (0.5*L);
            if M < 1
                FeNpl_stressP = FeNpl_stress * M;
                if (psi >= 0) && (psi <= 1)
                    FeNpl_strainP = FeNpl_strain * M;
                elseif (psi>1)
                    ratio5 = 1 - (1 - M)/psi;
                    if ((2*a/W) >0) && ((2*a/W)<=0.35)
                        A = 0.125 - (beta1 - 2.3422)/(psi - 1);
                        B = 0;
                    elseif ((2*a/W)>0.35)
                        A = 0.125 - (2*(beta1-2.3422))/(psi-1);
                        B = (beta1 - 2.3422)/((psi-1)^2);
                    end
                    psi4 = 32.6 - 70.4*(2*a/W) + 39.8*(2*a/W)^2;
                    if (psi>=1) && (psi<=psi4)
                        ratio6 = M * (beta1 + A*(psi-1) + B*(psi-1)^2)/beta1;
                    elseif (psi>psi4)
                        ratio6 = M * (0.125*psi + 2.2172)/beta;
                    end
                    FeNpl_strainP = min(ratio5,ratio6)*FeNpl_strain;
                end
            else
                FeNpl_stressP = FeNpl_stress;
                FeNpl_strainP = FeNpl_strain;
            end
            FeN = (FeNpl_strainP + FeNpl_stressP)/2;
            sigmaref = sigmaSdm * W * T * fyweldT/FeN;
        else
            sigmaref = referenceStress(flaw, constraint, a, c, T, W, sigmaSdb, sigmaSdm);
        end
        if flaw == 10
            Lr = sigmaSdm * W * T/FeN;
        else
            if fyweldT == 0
                Lr = sigmaref/fyT;
            else
                Lr = sigmaref/fyweldT;
            end
        end
        %% Fracture collapse parameter Kr

        % 7.1.8.2 in BS7910:2013 A1 2015 Residual stresses in as-welded structures:
        % Membrane residual stress:
        % structure in as-welded consition, with a flaw lying in a plane parallel
        % to the welding direction (the stresses to be considered are perpendicular
        % to the weld): lesser of the RT yield strengths of the weld or parent
        % material
        % structure in as-welded condition, with a flaw lying in a plane
        % transversal to the welding direction (the stresses to be considered are 
        % parallel to the weld): RT yield strength of the material in which the
        % flaw is located

        % 7.1.8.3 in BS7910:2013 A1 2015: Post-welded heat-treated (PWHT) structures
        %  a flaw lying in a plane parallel to the welding direction (the stresses
        %  to be considered are perpendicular to the weld): 20% of the lesser of
        %  the yield strength of the weld or parent material
        %  a flaw lying in a plane transverse to the welding direction (the stresses
        %  to be considered are parallel to the weld): 30% of RT yield strength of 
        % the material in which the flaw is located
        if (weld == 0)
            sigmaResM = 0;
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 11) || (weld == 21)
            % full penetration flaw, flaw perpendicular to weld direction,as-welded
            sigmaResM = min(fyT, fyT*(1.4 - 2*sigmaref/(fyT + fuT)));
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 12) || (weld == 22)
            % full penetration flaw, flaw parallel to weld direction,as-welded
            fy = min(fyT, fyweldT);
            fu = min(fuT, fuweldT);
            sigmaResM = min(fy, fy*(1.4 - 2*sigmaref/(fy + fu)));
%             sigmaResM = 0;
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 13) || (weld == 23)
            % full penetration flaw, flaw perpendicular to weld direction, PWHT
            sigmaResM = 0.30 * fyT;
            sigmaResB = 0;
            Ksb       = 0;
        elseif (weld == 14) || (weld == 24)
            % full penetration flaw, flaw parallel to weld direction, PWHT
            sigmaResM = 0.20 * min(fyT, fyweldT);
            sigmaResB = 0;
            Ksb       = 0;
        end

        % Annex R: Plasticity correction factor for secondary stresses
        KIp_a = (Yma * sigmaSdm + Yba * sigmaSdb)*sqrt(pi*a);

        KIs_a = (Mma * sigmaResM + Mba * sigmaResB)*sqrt(pi*a) + Ksb;

        chia = KIs_a * Lr / KIp_a;

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

        KI_a = KIp_a + KIs_a;

        Kr_a = KI_a/Kmat + rhoa;

        % Failure assessment diagram
        if fyT <= 1000
            deltaEpsilon = 0.037 * (1 - 0.001 * fyT);
        else
            disp('Error in Fracture: Yield strength has to lower dan 1000MPa');
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
%         [Lr, Kr_a, Krcr, KI_a, rhoa, sigmaResM, sigmaResB, Ksb]    = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT, fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Ma, Mma, Mba, Kmat);
         %     Kr_a, Lr, KI_a, Kmat, rhoa
%             xlRange                         = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL_add.xlsx', [Kr_a, 0, Lr, KI_a, 0, Kmat, 0, rhoa, 0], 'Sheet1', xlRange);
        rho                             = rhoa;
        Kr                              = Kr_a;
        KI                              = KI_a;
        batch                           = batch_name(ii);
        xlRange                         = [ 'A', num2str(row_num + shift_row)];
        xlswrite(fname, {batch, specimen, E, fyT, fuT, DeltaEpsilon, FADtype, Kr, Lr}, 'Sheet1', xlRange);
        xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet2', xlRange);
       
        row_num     = row_num + 1;
        specimen    = specimen + 1;
    end
end


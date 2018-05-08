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
specimen_type   = 'wide_plate';
specimen_treat  = 'per_specimen';
% specimen_treat  = 'average_ctod';
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
    load('NIL.mat')
%     load('ctod.mat')
%     load('ctodLT.mat')
%     load('ctodLS.mat')
    load('SCF.mat')
    
    index       = ii;
%     ctodLT      = ctodLT(index,1);
%     ctodLS      = ctodLS(index,1);
%     ctod        = ctod(index,1);

    % ---------------------------------------------------------------------
    % GEOMETRY INPUT PARAMETERS
    % ---------------------------------------------------------------------
    % Plate thickness
    T           = NIL(index,2);
    % Plate width:
    W           = NIL(index,3);
    % Initial defect size:
    a_i         = NIL(index,4);
    c_i         = NIL(index,5);
    % Final crack depth:
    a_f         = T;
    % Flaw type:
    flaw        = NIL(index,10);
    % 1 = Through-thickness flaws in plates
    % 2 = Semi elliptical surface flaws in plates (= flaw 1 in mathcad)
    % 3 = Through-width flaws in plates
    % 4 = semi-circular surface flaw in round bar
    % 5 = Corner flaws in plates
    % 6 = Semi-circular surface flaw in bolts, solution 2
    % 7 = Edge flaws in plates
    % 8 = Through-thickness deck plate crack (orthotropic bridge deck)
    % 9 = symmetric crack emanating from hole
    hole        = 200;
    % Type of weld
    weld        = NIL(index,11);
    % Weld angle
    Wangle      = NIL(index,16);
    % Attachment lenght
    L           = NIL(index,15);
    % Weld toe radius
    radius      = 0;
    % Weld throat
    tw          = 0;
    crackLocation = 'chord';

    % ---------------------------------------------------------------------
    % MATERIAL INPUT PARAMETERS
    % ---------------------------------------------------------------------
    % Yield stress
    fy          = NIL(index,6);
    fyweld      = NIL(index,17);
    % Ultimate tensile strength
    fu          = NIL(index,7);
    fuweld      = NIL(index,18);
    % Modulus of elasticity
    E           = NIL(index,12);
    % Poisson ratio:
    niu         = 0.3;
    % Width CTOD specimen
    Bctod       = NIL(index,13);
    % Crack front lenght:
    Lcrack      = NIL(index,14);
    % Temperature assessment:
    Temp_a      = NIL(index,1);
    % CTOD temperature:
    Temp_ctod   = NIL(index,1);
    % Assessment temperature: 7.1.3.4 in BS7910:2013 A1 2015
    if Temp_a == Temp_ctod
        fyT         = fy;
        fuT         = fu;
        fyweldT     = fyweld;
        fuweldT     = fuweld;
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
    % only for fatigue specimens:
    GS = W * (762^2 - 717^2) / (2*762);
    if (index <= 18)
        sigmaSdm = NIL(index,8)*1E3 / (W*T);
    elseif (index > 18) && (index <= 22)
        sigmaSdm = NIL(index,8)*1E3 / (GS);
    elseif (index > 22) && (index <=25)
        sigmaSdm = NIL(index,8)*1E3 / (W*T);
    else
        sigmaSdm = SCF(index-25,1) * NIL(index,8) * 1E3 / (W*T);
    end
    sigmaSdb = 0;
    
    % ---------------------------------------------------------------------
    % FRACTURE ANALYSIS
    % ---------------------------------------------------------------------
    a = a_i;
    c = c_i;
    % write the header: column names
    if index == 1
        xlRange = [ 'A', num2str(1)];
%         xlswrite('NIL.xlsx', {'Kr_a', 'Kr_c', 'Lr', 'KI_a', 'KI_c', 'Kmat_a', 'Kmati_c', 'rhoa', 'rhoc'}, 'Sheet1', xlRange);
        xlswrite(fname, {'batch', 'specimen', 'CTOD', 'KI', 'Kmat', 'rho', 'Kr', 'Lr', 'Lr_trunc_flag', 'Lr_sup'}, 'Sheet1', xlRange);
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
        if (index <= 6)
            % T-stress
            [ Tstr_a, Tstr_c ]              = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
            [Kmat, Kmati]                   = FractureToughnessOSK(fyT, fuT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_a);
            [Yma, Yba]                      = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
            [~, Mma, Mba, ~]                = SIMfactor(flaw, a, c, W, T, pi/2, 0);
           [Lr, Kr_a, Krcr, KI_a, rhoa, ~, ~, ~]    = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat);
            %     Kr_a, Lr, KI_a, Kmat, rhoa
%             xlRange                         = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL.xlsx', [Kr_a, 0, Lr, KI_a, 0, Kmat, 0, rhoa, 0], 'Sheet1', xlRange);
            rho                             = rhoa;
            Kr                              = Kr_a;
            KI                              = KI_a;
            batch                           = batch_name{ii};
            xlRange                         = [ 'A', num2str(row_num + shift_row)];
            xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
% disp('_______')
        elseif (index==7) || (index == 8)
            [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
            [Kmat, Kmati] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_a);
            [Yma, Yba] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
            [Ma, Mma, Mba, ~] = SIMfactor(flaw, a, c, W, T, pi/2, hole);
            [Lr, ~, ~, ~, ~] = Fracture_TT(flaw, weld, a+hole/2, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat);
            [~, Kr_a, Krcr, KI_a, rhoa,  ~, ~, ~] = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat);
            %     Kr_a, Lr, KI_a, Kmat, rhoa
%             xlRange = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL.xlsx', [Kr_a, 0, Lr, KI_a, 0, Kmat, 0, rhoa, 0], 'Sheet1', xlRange);
            
            rho                             = rhoa;
            Kr                              = Kr_a;
            KI                              = KI_a;
            batch                           = batch_name{ii};
            xlRange                         = [ 'A', num2str(row_num + shift_row)];
            xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
% disp('_______')            
        elseif (index == 9) || (index == 10)
            [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
            [Kmat, Kmati] = FractureToughnessOSK(fyweldT, fuweldT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_a);
            [Yma, Yba] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
            [Ma, Mma, Mba, ~] = SIMfactor(flaw, a, c, W, T, pi/2, hole);
            [Lr, ~, ~, ~, ~] = Fracture_TT(flaw, weld, a+hole/2, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat);
            [~, Kr_a, Krcr, KI_a, rhoa,  ~, ~, ~] = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat);
            %     Kr_a, Lr, KI_a, Kmat, rhoa
%             xlRange = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL.xlsx', [Kr_a, 0, Lr, KI_a, 0, Kmat, 0, rhoa, 0], 'Sheet1', xlRange);

            rho                             = rhoa;
            Kr                              = Kr_a;
            KI                              = KI_a;
            batch                           = batch_name{ii};
            xlRange                         = [ 'A', num2str(row_num + shift_row)];
            xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
% disp('_______')
        elseif (index >= 15) && (index <= 18)
            [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
            [Kmat, Kmati] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLS, Bctod, Lcrack, 1, Tstr_a);
            [Yma, Yba] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
            [~, Mma, Mba, ~] = SIMfactor(flaw, a, c, W, T, pi/2, 0);
            [Lr, Kr_a, Krcr, KI_a, rhoa,  ~, ~, ~] = Fracture_TT(flaw, weld, a, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Yba, Mma, Mba, Kmat)    ;
            %     Kr_a, Lr, KI_a, Kmat, rhoa
%             xlRange = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL.xlsx', [Kr_a, 0, Lr, KI_a, 0, Kmat, 0, rhoa, 0], 'Sheet1', xlRange);

            rho                             = rhoa;
            Kr                              = Kr_a;
            KI                              = KI_a;
            batch                           = batch_name{ii};
            xlRange                         = [ 'A', num2str(row_num + shift_row)];
            xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
% disp('_______')
        else
            [ Tstr_a, Tstr_c ] = T_stress( a, c, T, sigmaSdm, sigmaSdb, flaw, W );
            if fyweld == 0
                [Kmat_a, Kmati_a] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLS, Bctod, Lcrack, 1, Tstr_a);
                [Kmat_c, Kmati_c] = FractureToughnessOSK(fyT, fuT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_c);
            else
                [Kmat_a, Kmati_a] = FractureToughnessOSK(fyweldT, fuweldT, E, niu, ctodLS, Bctod, Lcrack, 1, Tstr_a);
                [Kmat_c, Kmati_c] = FractureToughnessOSK(fyweldT, fuweldT, E, niu, ctodLT, Bctod, Lcrack, 2, Tstr_c);
            end
            [Yma, Yba] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, pi/2, hole);
            [Ymc, Ybc] = SICfactorY (index, a, c, T, W, L, tw, Wangle, radius, 0, weld, flaw, 0, hole);
            [~, Mma, Mba, ~] = SIMfactor(flaw, a, c, W, T, pi/2, 0);
            [~, Mmc, Mbc, ~] = SIMfactor(flaw, a, c, W, T, 0, 0);
            [Lr, Kr_a, Kr_c, Krcr, KI_a, KI_c,rhoa, rhoc,  ~, ~, ~] = Fracture(index, flaw, weld, a, c, T, W, fyT, fuT,fyweldT, fuweldT, E, sigmaSdm, sigmaSdb, Yma, Ymc, Yba, Ybc, Mma, Mmc, Mba, Mbc, Kmat_a, Kmati_c);
            %     Kr_a, Kr_c, Lr, KI_a, KI_c, Kmat_a, Kmati_c, rhoa, rhoc
%             xlRange = [ 'A', num2str(index + shift_row)];
%             xlswrite('NIL.xlsx', [Kr_a, Kr_c, Lr, KI_a, KI_c, Kmat_a, Kmati_c, rhoa, rhoc], 'Sheet1', xlRange);

            rho     = rhoa;
            Kr      = Kr_a;
            KI      = KI_a;
            Kmat    = Kmat_a;
%             if Kr_a > Kr_c
%                 rho     = rhoa;
%                 Kr      = Kr_a;
%                 KI      = KI_a;
%                 Kmat    = Kmat_a;
%             else
%                 rho     = rhoc;
%                 Kr      = Kr_c;
%                 KI      = KI_c;
%                 Kmat    = Kmat_c;
%             end
            batch                           = batch_name{ii};
            xlRange                         = [ 'A', num2str(row_num + shift_row)];
            xlswrite(fname, {batch, specimen, ctod, KI, Kmat, rho, Kr, Lr, Lr_tf, Lr_s}, 'Sheet1', xlRange);
        end
% disp('_______')        
        row_num     = row_num + 1;
        specimen    = specimen + 1;
    end
end


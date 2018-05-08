% Get experimental CTOD  data
%
%SYNOPSYS
% [batch_name, ctod, Lr_trunc_flag] = GET_CTOD_DATA(type)
%
% Timestamp: 2017 november 3
%
function [batch_name, CTOD, Lr_trunc_flag, Lr_sup] = get_ctod_data(type)

switch lower(type)
    case 'wide_plate'
        fname           = 'experimental_ctod_wide_plate.csv';
    case 'tubular'
        fname           = 'experimental_ctod_tubular.csv';
    otherwise
        error(['Unknown type: ', type])
end

% grab the data
A               = importdata(fname, ',', 1);

header          = A.textdata(1,:);
batch_name      = A.textdata(2:end,1);
Lr_trunc_flag   = A.data(:,1);
Lr_sup          = A.data(:,2);
CTOD            = A.data(:,3:end);

end



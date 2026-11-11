function [z,T,S] = load_ctd_nc(f,zname,Tname,Sname,prof)
% LOAD_CTD_NC  Read one CTD profile from a NetCDF file and clean it.
%   [z,T,S] = load_ctd_nc(file, depthVar, tempVar, salVar, profileIndex)
%   Returns column vectors, sorted by depth, NaNs and duplicate depths removed.
%   No unit conversion: check that T is conservative and S absolute first.
if nargin < 5, prof = 1; end

z = readvar(f,zname,prof);
T = readvar(f,Tname,prof);
S = readvar(f,Sname,prof);

ok = ~isnan(z) & ~isnan(T) & ~isnan(S);
z = z(ok); T = T(ok); S = S(ok);
[z,ui] = unique(z);          % sorts and removes duplicate depths
T = T(ui); S = S(ui);
end

function v = readvar(f,name,prof)
v = double(ncread(f,name));
if ~isvector(v)              % 2-D: depth x profile (check orientation with ncdisp)
    v = v(:,prof);
end
v = v(:);
v(abs(v) > 1e30 | v == -999) = NaN;   % catch common fill values
end
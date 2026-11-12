function [Q, tt, Qts] = outlet_discharge_mean(f_nc, outlet_id, cast_file, ndays)
% OUTLET_DISCHARGE_MEAN  Mean Mankoff outlet discharge over the ndays before a CTD cast.
%   f_nc       path to the yearly Mankoff MAR discharge .nc file
%   outlet_id  outlet ID(s) to extract (several IDs are summed)
%   cast_file  CTD filename containing yyyymmddHHMMSS (OMG naming)
%   ndays      number of whole days before the cast day to average
%   Q = mean (m^3/s); tt, Qts = full daily series for checking.

info = ncinfo(f_nc);
vnames = {info.Variables.Name};

% --- (1) extract the daily series for the outlet ---
idvar = first_present(vnames,{'station','outlet','cat','id'});
qvar  = first_present(vnames,{'discharge','runoff','Q'});
ids   = double(ncread(f_nc,idvar));
q     = double(ncread(f_nc,qvar));

idx = find(ismember(ids,outlet_id));
if isempty(idx)
    error('Outlet ID %s not found in "%s" (%d IDs, e.g. %s).', ...
        mat2str(outlet_id), idvar, numel(ids), mat2str(ids(1:min(5,end))'))
end
fprintf('Outlet %s found: %d match(es) in "%s"\n', mat2str(outlet_id), numel(idx), idvar)

if size(q,1) == numel(ids)            % outlet x time
    Qts = sum(q(idx,:),1)';
elseif size(q,2) == numel(ids)        % time x outlet
    Qts = sum(q(:,idx),2);
else
    error('Cannot match discharge dimensions %s to %d outlets.', mat2str(size(q)), numel(ids))
end

% --- time axis from the file's own units attribute ---
tv    = double(ncread(f_nc,'time'));
units = ncreadatt(f_nc,'time','units');                  % e.g. 'days since 1950-01-01'
tok   = regexp(units,'^(\w+)\s+since\s+(.+)$','tokens','once');
epoch = datetime(regexp(tok{2},'\d{4}-\d{1,2}-\d{1,2}','match','once'),'InputFormat','yyyy-MM-dd');
switch lower(tok{1}(1:3))
    case 'day', tt = epoch + days(tv);
    case 'hou', tt = epoch + hours(tv);
    case 'sec', tt = epoch + seconds(tv);
    otherwise,  error('Unknown time unit: %s', units)
end
tt = tt(:);
fprintf('Time axis: %s to %s (%d steps)\n', string(tt(1)), string(tt(end)), numel(tt))

% --- (2) cast date from the filename ---
tok = regexp(cast_file,'_(\d{14})_','tokens','once');
cast_time = datetime(tok{1},'InputFormat','yyyyMMddHHmmss');
fprintf('Cast time: %s\n', string(cast_time))

% --- (3) mean over the ndays whole days before the cast day ---
d_end = dateshift(cast_time,'start','day');
win   = tt >= d_end - days(ndays) & tt < d_end;
if nnz(win) < ndays
    warning('Only %d of %d days in the window.', nnz(win), ndays)
end
Q = mean(Qts(win),'omitnan');
fprintf('%d-day window %s to %s: mean %.2f, min %.2f, max %.2f m^3/s\n', ndays, ...
    string(min(tt(win))), string(max(tt(win))), Q, min(Qts(win)), max(Qts(win)))

% quick look: series with the window shaded
figure; plot(tt,Qts,'k'); hold on
xregion(d_end - days(ndays), d_end)
xline(cast_time,'r--','CTD cast'); ylabel('Discharge (m^3/s)'); title('Outlet discharge')
end

function v = first_present(vnames,candidates)
for c = candidates
    if any(strcmp(vnames,c{1})), v = c{1}; return, end
end
error('None of %s found. Variables in file: %s', strjoin(candidates,', '), strjoin(vnames,', '))
end


%------------- Call the Script -----------------------------
%f_nc = 'C:\Users\efc4\OneDrive - University of St Andrews\Desktop\PhD\umq_fjords\data\raw\mankoff2020_v8.0\discharge\MAR_2016.nc';
%cast_file = 'OMG_Ocean_CTD_L2_20160922074809_gsw.nc';     % near-glacier cast
%[Q, tt, Qts] = outlet_discharge_mean(f_nc, 80710, cast_file, 14)
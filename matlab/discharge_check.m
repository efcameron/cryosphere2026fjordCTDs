%% Mean IKS discharge during 14 days before CTD cast

S = load(qfile);

% Identify IKS glaciers
iks_idx = find(S.fjord_id == "iks");

% Time and discharge
tt = S.tAll;
Q  = S.QAll;

% CTD cast time
cast_time = datetime(2016,9,15,16,44,46);

% Define 14-day window before cast day
d_end   = dateshift(cast_time,'start','day');
d_start = d_end - days(14);

win = tt >= d_start & tt < d_end;

% Extract individual glacier discharge
Q_store = Q(win,iks_idx(1));
Q_lille = Q(win,iks_idx(2));

% Combined discharge
Q_total = Q_store + Q_lille;

% Calculate means
mean_store = mean(Q_store,'omitnan');
mean_lille = mean(Q_lille,'omitnan');
mean_iks   = mean(Q_total,'omitnan');

% Display results
fprintf('\nIKS discharge before CTD cast\n');
fprintf('Cast: %s\n',string(cast_time));
fprintf('Window: %s to %s\n',string(d_start),string(d_end));
fprintf('Number of days: %d\n',nnz(win));

fprintf('Store mean discharge: %.2f m3/s\n',mean_store);
fprintf('Lille mean discharge: %.2f m3/s\n',mean_lille);
fprintf('Combined IKS discharge: %.2f m3/s\n',mean_iks);

%% MAR: IKS discharge throughout 2016

win2016 = tt >= datetime(2016,1,1) & ...
    tt < datetime(2017,1,1);

t2016 = tt(win2016);

Q_kuj = Q(win2016,iks_idx(1));
Q_ava = Q(win2016,iks_idx(2));

Q_iks = Q_kuj + Q_ava;

figure('Color','w')

plot(t2016,Q_iks,'k','LineWidth',1.3)
hold on

xline(datetime(2016,9,15),'r--','CTD cast');

ylabel('Discharge (m^3 s^{-1})')
xlabel('Date')
title('IKS – MAR-derived daily discharge, 2016')
grid on
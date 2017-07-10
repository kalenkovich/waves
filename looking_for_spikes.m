% In bst: jordy_spont03 -> import to db -> export to Matlab as X

% Bandpass-filter the data
Fs = fix(1/(X.Time(2)-X.Time(1)));
[b,a] = fir1(128, [1 50]/(Fs/2), 'bandpass');
channel_idx = 32:size(X.F, 1); % MEG channels start with 32
DataF = filtfilt(b, a, X.F(channel_idx, :)')';

% N = length( X.Events(2).samples);
% i=1;
% Av = zeros(size(DataF(:,X.Events(2).samples(i)-Fs/2:X.Events(2).samples(i)+Fs/2)));
% for i=1:N
%    Av = Av+DataF(:,X.Events(2).samples(i)-Fs/2:X.Events(2).samples(i)+Fs/2);
% end;
% close all
% figure
% plot(Av')

% Run ICA on filtered data
% W - matrix of ICA weights
% S - sphering (whitening matrix)
[W, S] = runica(DataF);
% Actual unmixing matrix Q is W*S
Q = W*S;
% The mixing matrix is the inverse of Q
% The topographies of the ICA components are columns of iQ
iQ = inv(Q);

% Let's sort the components in the descending order of kurtoses of their
% corresponding time series Z
Z = Q*DataF;
kZ = kurtosis(Z');
[~, kkey] = sort(kZ, 'descend');
iQs = iQ(:, kkey);

% Let's create a dataset that will contain only the sorted topographies as
% data
ICA_topo = X;
channel_cnt = length(channel_idx);
IC_cnt = channel_cnt;
ICA_topo.F = X.F(:, 1:channel_cnt);
ICA_topo.F(channel_idx, 1:IC_cnt) = iQs;
ICA_topo.Time = 1:channel_cnt;
ICA_topo.Comment = [ICA_topo.Comment ' ICA topographies sorted'];
% Import to Brainstorm

% In bst Right click -> MEG -> 2D Disc
% Right click -> Snapshot -> Time contact sheet
% Look for components with focal topographies that are close to each other
IC_picked = [6 10];

% Now unmix, then mix back again with only the picked components
spike_ICs = X;
Zs = Z(kkey, :);
spike_ICs.F(32:end, :) = iQs(:, IC_picked) * Zs(IC_picked, :);
spike_ICs.Comment = [spike_ICs.Comment ' only IC ' num2str(IC_picked) ' retained'];
% Import to Brainstorm

% Export to Matlab

DataF6_12 = iQs(:,[6,12])*Zs([6,12],:);
Zs = Q(kkey,:)*DataF1;
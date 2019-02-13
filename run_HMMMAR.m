 

if ~exist('data','var') || ~exist('T','var')
    error('You need to load the data (data and T - see Documentation)')
end

data_modality = 'fMRI' ; % one of: 'fMRI', 'M/EEG', 'M/EEG power' or 'LFP' 
no_states = 25; % the number of states depends a lot on the question at hand
Hz = 1; % the frequency of the data
stochastic_inference = 1; % set to 1 if a normal run is too computationally expensive (memory or time)
N = length(T); % number of subjects

% getting the number of channels
if iscellstr(data) 
    dfilenames = data;
    if ~isempty(strfind(dfilenames{1},'.mat')), load(dfilenames{1},'X');
    else X = dlmread(dfilenames{1});
    end
elseif iscell(data)
    X = data{1};
end
ndim = size(X,2); 

% Setting the options

options = struct();
options.K = no_states;
options.standardise = 1;
options.verbose = 1;
options.Fs = Hz;

if iscell(T), sumT = 0; for j = 1:N, sumT = sumT + sum(T{j}); end
else, sumT = sum(T); 
end

if strcmp(data_modality,'fMRI') % Gaussian observation model
    options.order = 1;
    options.zeromean = 0;
    options.covtype = 'full';    
end

% stochastic options
if stochastic_inference
    options.BIGNbatch = max(round(N/30),5);
    options.BIGtol = 1e-7;
    options.BIGcyc = 500;
    options.BIGundertol_tostop = 5;
    options.BIGforgetrate = 0.7;
    options.BIGbase_weights = 0.9;
end

% HMM computation
[hmm, Gamma, Xi, vpath] = hmmmar(data,T,options);

% HMM TEST FOR Autism

figure; subplot(3,1,1)
plot(Gamma(1:1000,:)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')
subplot(3,1,2)
plot(vpath(1:1000)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')


% figure
% 
% subplot(2,4,1), imagesc(getFuncConn(hmm,1)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,2), imagesc(getFuncConn(hmm,2)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,3), imagesc(getFuncConn(hmm,3)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,4), imagesc(getFuncConn(hmm,4)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
figure
% 
% subplot(2,4,1), imagesc(getFuncConn(hmm,1)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,2), imagesc(getFuncConn(hmm,2)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,3), imagesc(getFuncConn(hmm,3)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
% subplot(2,4,4), imagesc(getFuncConn(hmm,4)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))










%

% Some useful information about the dynamics
maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM 
            % is capturing dynamics or grand between-subject 
            % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats

% W = getMARmodel(hmm,k)
% m = getMean(hmm,k)
TP = getTransProbs (hmm)
%figure; subplot(3,1,1)
%plot(FO(1:1000,:)), set(gca,'Title',text('String','True state path'))
%set(gca,'ylim',[-0.2 1.2]); ylabel('state #')
% subplot(3,1,2)
% plot(vpath(1:1000)), set(gca,'Title',text('String','True state path'))
% set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')







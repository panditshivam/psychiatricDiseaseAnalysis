addpath(genpath('.')) % assuming we are in the HMM-MAR directory
K = 4; % number of states
ndim = 30;%3; % number of channels
T = 146 * ones(1,1);
%T = 146;
N = 1; % number of trials
Fs = 200;
X = data;

%T = 10000 * ones(N,1); % number of data points
hmmtrue = struct();
hmmtrue.K = K;
hmmtrue.state = struct();
hmmtrue.train.covtype = 'full';
hmmtrue.train.zeromean = 0;
hmmtrue.train.order = 0;
hmmtrue.train.orderoffset = 0;
hmmtrue.train.timelag = 2;
hmmtrue.train.exptimelag = 0;
hmmtrue.train.S = ones(ndim);
hmmtrue.train.Sind = ones(1,ndim);
hmmtrue.train.multipleConf = 0;

for k = 1:K
    hmmtrue.state(k).W.Mu_W = rand(1,ndim);
    r = randn(ndim);
    hmmtrue.state(k).Omega.Gam_rate = 0.1 * 1000 * r' * r + eye(ndim);
    hmmtrue.state(k).Omega.Gam_shape = 1000;
end

hmmtrue.P = rand(K) + 100 * eye(K);  
for j=1:K
    hmmtrue.P(j,:) = hmmtrue.P(j,:) ./ sum(hmmtrue.P(j,:));
end
hmmtrue.Pi = ones(1,K); %rand(1,K);
hmmtrue.Pi = hmmtrue.Pi./sum(hmmtrue.Pi);
%[X,T,Gammatrue] = simhmmmar(T,hmmtrue,[]);
%Flips = binornd(1,0.2,N,3);
%unflipped_X = X; 
%$for n = 1:N
  %  if mean(Flips(n,:))>0.5 % we force to be more unflipped than flipped
   %     Flips(n,:) = 1 - Flips(n,:);
    %end
    %for d = 1:ndim
     %   ind = (1:T(n)) + sum(T(1:n-1));
      %  if Flips(n,d) 
       %     X(ind,d) = -X(ind,d);
        %end
   % end
%end 





options = struct();
options.K = K; 
options.Fs = Fs; 
options.covtype = 'full';
options.order = 0;
options.DirichletDiag = 2; 
options.zeromean = 0;
options.verbose = 1;

[hmm, Gamma, ~, vpath] = hmmmar(data,T,options);

figure; subplot(3,1,1)
plot(Gammatrue(1:1000,:)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')


subplot(3,1,2)
plot(Gamma(1:1000,:)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[-0.2 1.2]); ylabel('state #')
subplot(3,1,3)
plot(vpath(1:1000)), set(gca,'Title',text('String','True state path'))
set(gca,'ylim',[0 hmm.K+1]); ylabel('state #')


figure
subplot(2,4,1), imagesc(getFuncConn(hmmtrue,1)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,2), imagesc(getFuncConn(hmmtrue,2)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,3), imagesc(getFuncConn(hmmtrue,3)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))
subplot(2,4,4), imagesc(getFuncConn(hmmtrue,4)), colormap('gray'), set(gca,'Title',text('String','Simulated covariance'))


subplot(2,4,5), imagesc(getFuncConn(hmm,1)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,6), imagesc(getFuncConn(hmm,2)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,7), imagesc(getFuncConn(hmm,3)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))
subplot(2,4,8), imagesc(getFuncConn(hmm,4)), colormap('gray'), set(gca,'Title',text('String','Inferred covariance'))


options = struct();
options.K = K; 
options.Fs = Fs; 
options.covtype = 'full';
options.order = 2;
options.DirichletDiag = 2; 
options.zeromean = 1;
options.verbose = 1;

[hmm, Gamma,Xi] = hmmmar(data,T,options);


[mcv,cv] = cvhmmmar(data,T,options);
fe = hmmfe(X,T,hmm,Gamma,Xi);
sum(fe)

options.Fs = Fs; 
options.completelags = 1;
options.MLestimation = 1; 
options.order = 20; % increase the order
options.p = 0.01;
spectral_info = hmmspectramar(X,T,[],Gamma,options);


figure
for k=1:2
    subplot(1,2,k)
    plot(spectral_info.state(k).f,spectral_info.state(k).psd(:,1,1),'k', 'LineWidth',3)
end


options_mt = struct('Fs',Fs);
options_mt.fpass = [1 48];
options_mt.tapers = [4 7];
options_mt.p = 0;
options_mt.win = 500;
options_mt.Gamma = Gamma;
spectral_info = hmmspectramt(X,T,Gamma,options_mt);


figure
for k=1:2
    subplot(1,2,k)
    plot(spectral_info.state(k).f,spectral_info.state(k).psd(:,1,1),'k')
end


FO = getFractionalOccupancy (Gamma,T); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T); % state life times
Intervals = getStateIntervalTimes (Gamma,T); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T); % rate of switching between stats




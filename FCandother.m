maxFO = getMaxFractionalOccupancy(Gamma,T,options); % useful to diagnose if the HMM 
            % is capturing dynamics or grand between-subject 
            % differences (see Wiki)
FO = getFractionalOccupancy (Gamma,T,options); % state fractional occupancies per session
LifeTimes = getStateLifeTimes (Gamma,T,options); % state life times
Intervals = getStateIntervalTimes (Gamma,T,options); % interval times between state visits
SwitchingRate =  getSwitchingRate(Gamma,T,options); % rate of switching between stats



for k = 1:10
    W = getMARmodel(hmm,k);
    %m = getMean(hmm,k);
end

TP = getTransProbs (hmm);

% for k = 1:10
%     p = cholcov(C);
%     if p ~= 0
%         %figure(k);
%         plot(2,4), imagesc(getFuncConn(hmm,k)), set(gca,'Title',text('String','Inferred covariance'))
%     end
%     
% end

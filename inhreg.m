function [spikes, stimes, isi] = inhreg(t, dt, f)
% function [spikes, stimes isi] = inhreg(t, dt, f)
% Inhomogenous regular distributed ISIs
% t - time vector
% dt - time step
% f - instantaneous rate vector (per timebase)
% Basic assumptions are:
% (1) constant rate (frequency) over a time step
% (2) only a single arrival possible in a time step
% (so time step should be small relative to the rate of change in
% frequency and arrival rate)
% BPG 14-1-08

spikes=zeros(1,length(t));    
spikes(1)=0;
tp=0;   % index of previous spike
for i=2:length(t)
    cisi=1/f(i);    % current ISI (seconds)
    if (i-tp)*dt >= cisi
        spikes(i)=1;
        tp=i;       % index of new spike
    end;
end;
stimes=t(spikes==1);
isi=stimes(2:length(stimes))-stimes(1:length(stimes)-1);

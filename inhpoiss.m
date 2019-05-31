function [spikes, stimes, isi] = inhpoiss(t, dt, f)
% function [spikes, stimes isi] = inhpoiss(t, dt, f)
% Inhomogenous Poisson distributed ISIs
% t - time vector
% dt - time step
% f - instantaneous rate vector (per timebase)
% Basic assumptions are:
% (1) constant rate (frequency) over a time step
% (2) only a single arrival possible in a time step
% (so time step should be small relative to the rate of change in
% frequency and arrival rate)
% BPG 11-1-08

pr=1-exp(-f*dt);    % spike occurrence probabilities
rn=rand(1,length(t));    % uniform random numbers
spikes=rn<=pr;
stimes=t(spikes==1);
isi=stimes(2:length(stimes))-stimes(1:length(stimes)-1);

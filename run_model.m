% Run stochastic CoH synapse model
% Ref: Yang et al, Neural Computation, in press
% Z. Yang, M. Hennig and B. Graham, University of Stirling, 2008

% Stimulus parameters - set to what you want
fre=[10 20 50 100];     % frequencies (Hz)
stimtime=1;             % stimulation time (s)
fstimtype=1;            % Type: (1) regular ISIs, (2) Poisson ISIs

% Time step (no need to change this)
dt = 0.0001; % time step for spike train generation (secs)

% For plotting
syms = ['.', '*', '+', 'o', 's', 'd']';
lines = ['k-','k--','k-.','k:']';
colors = ['k','r','b','m','y','c']';
lwidth = 1;

% Load and plot experimental data
e10 = load('expdata/Ca2mM_10Hz_norm.dat');
plot(e10(:,1), e10(:,2)/100, 'k-');
hold on;
e20 = load('expdata/Ca2mM_20Hz_norm.dat');
plot(e20(:,1), e20(:,2)/100, 'k-');
e50 = load('expdata/Ca2mM_50Hz_norm.dat');
plot(e50(:,1), e50(:,2)/100, 'k-');
e100 = load('expdata/Ca2mM_100Hz_norm.dat');
plot(e100(:,1), e100(:,2)/100, 'k-');

% Do simulations
for i=1:length(fre)
  
  fvec = fre(i)*ones(1, stimtime/dt);
  tvec = dt:dt:stimtime;
  % Stimulus type can be regular or Poisson distributed ISIs
  if fstimtype == 1   % regular ISIs
    [spikes, stimes, isi] = inhreg(tvec, dt, fvec); 
  elseif fstimtype == 2 % Poisson ISIs
    [spikes, stimes, isi] = inhpoiss(tvec, dt, fvec);
  end;
  num = length(isi); 

  % Canonical synapse model
  [psr, npsr] = coh_stoch_mod(isi);

  xtime = stimes(1:num-1);
  resps = npsr(1:num-1);
  
  p=plot(xtime, resps, syms(1,:));
%  p=plot(xtime, resps, syms(mod(i-1,length(syms))+1,:));
  set(p,'Color',colors(mod(i-1,length(colors))+1,:),'LineWidth',lwidth);
  hold on;
end;

function [resps, nresps] = coh_stoch_mod(isi)
% Stochastic model of vesicle recycling and release at the calyx of Held
% isi - sequence of interspike intervals
% resps - postsynaptic responses
% nresps - normalised responses
% Ref: Yang et al, Neural Computation, in press
% Z. Yang, M. Hennig and B. Graham, University of Stirling, 2008

rand('state', sum(100*clock));  % seed random numbers

% Full model parameter values
gke   = 0.058;      % re
gkd   = 0.4;        % rp (1/s)
gfac  = 0.091;      % nf
gfrel = 0.0252;     % tf (s)
grd   = 4;          % nd
grr   = 0.043;      % td (s)
gcai    = 0.003;    % ni
gcairel = 8;        % ti (s)
ggli    = 0.21;     % nb
gglrel  = 0.6;      % tb (s)
gcares   = 0;       % residual calcium (not used)   
gcarel   = 10;      % decay time (s)

ca = 10;            % Ca++ transient amplitude (uM)
k=0.00001628;       % factor to convert (ca*[0:1])^4 into release prob.
exponent = 4;       % power law for release probability

vesiNum = 5;        % number of vesicles per active zone
indexx = 550;       % number of AZs
totalNum = vesiNum*indexx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = 1:length(isi);
pulse = length(time);

% initial values
ns = zeros(pulse, indexx);
for index = 1 : indexx
    ns(1,index) = 5;
end
nrels = zeros(pulse, 1);
noldd = zeros(pulse, 1);
nold = zeros(pulse, 1);


rdess(1) = 0;
rdess2(1) = 0;
cares(1) = 0;
caress(1) = 0;
pp = 1;
ppb = 1;
ppin = 0;
ppmglu  = 0;
prel(1)  = 1-exp(-k*(ca)^exponent);
pprel(1) = 1-exp(-k*(ca)^exponent);

ca_amp = zeros(pulse, 1);
for index = 1 : indexx
    ca_amp(1) = pp;
end

ppbase(1) = ca_amp(1);
ppfac(1) = ca_amp(1);
  
% Enhanced Recycling
kee = 1-exp(-gke);

for t=time,

  dt = isi(t);
  labNum = 0;
  
  % Recycling
  kd  = 1-exp(-dt*gkd);

  % Facilitation
  fac = gfac;
  if gfrel == 0
      frel = 1;
  else
      frel  = 1-exp(-dt/gfrel);
  end
      
  % Ca++ dynamics
  % inactivation
  cairel  = 1-exp(-dt/gcairel);
  % mGluR
  glrel = 1-exp(-dt/gglrel);
  % residual Ca++
  carel = 1-exp(-dt/gcarel);
  
  % Desensitisation
  rd = grd;
  if grr == 0
      rr = 1;
  else
      rr = 1-exp(-dt/grr);
  end
  
  
  for index = 1:indexx

    nr = ns(t, index);
 
    lab = 0;
    noldd(t) = noldd(t) + ns(t, index);
    if nr >= 1
        for jj = 1:nr
        
           if ((rand < pprel(t)) && (t <=(pulse-1)))   
                % fuse a vesicle
  
               ns(t, index) = ns(t, index) - 1; 
             
               labNum = labNum + 1;
               ns(t+1, index) = ns(t, index);
        
               lab = 1;
            end
        end 
    end

    if lab == 0
        ns(t+1, index) = ns(t, index);  
       
    else
        lab = 0;
    end
    
    %if t <= (pulse-1)
    %   nrel(t) = nrel(t) + labNum;
    %end
    
    
    if (vesiNum - nr) >= 1
          %passive recycling
        for kk = 1:(vesiNum - nr)
            if ((rand < gkd*dt ) && ( t <= (pulse-1))) 
              
                ns(t, index) = ns(t, index) + 1;
        
                ns(t+1, index) = ns(t, index);
           
            end
        end
    end
    
    
    nen = ns(t, index);
    
    if (vesiNum - nen) >= 1
          %enhanced recycling
        for kkk = 1:(vesiNum - nen)
            if ((rand < gke ) && ( t <= (pulse-1))) 
               
                ns(t, index) = ns(t, index) + 1;
                
                ns(t+1, index) = ns(t, index);
               
            end
        end
    end            
    
    nd = ns(t, index);
    if nd >= 1
       %passive depletion
        for kkkk = 1:nd
            
            if ((rand < (1000*dt*0.01*1000*(ca*cares(t)))) && (t <= (pulse-1)))
                ns(t+1, index) = ns(t, index) - 1;
                ns(t+1, index) = ns(t, index);
            end
        end
    end
    
  end
    
  nold(t) = noldd(t) / totalNum;
  nrels(t) = labNum / totalNum;  %nold(t) * pprel(t);
    
  % DESENSITISATION
  desstmp = 1-rdess(t);
  dess(t) = desstmp;

  rdess(t+1) = rdess(t) + nrels(t) * rd * desstmp;
  rdess(t+1) = rdess(t+1) -  rr * rdess(t+1);
   
  % CALCULATE RESPONSE
  resps(t) = nrels(t) * desstmp;

  % FACILITATION
  pp = pp + fac * ca_amp(t);

  % temporary variables
  % mGluR activation
  pglut = ggli * nrels(t) * ppb;
  % Ca++ channel inactivation  
  pcait = gcai *  ppb;
    
  carest =  ca_amp(t) * gcares * (0.05-cares(t));

  % BASE PROBABILITY
  ppb = ppb - pcait - pglut;
    
  % mGluR ACTIVATION
  ppmglu = ppmglu + pglut;

  % CA++ CHANNEL INACT
  ppin = ppin + pcait;
  
  % FINISH BASE PROBABILITY
  ppb = ppb + cairel * ppin + glrel * ppmglu;

  % RECOVERY FROM mGluR ACTIVATION
  ppmglu = ppmglu - glrel * ppmglu;

  % RECOVERY FROM CA++ CHANNEL INACT
  ppin = ppin - cairel * ppin;
    
  % RECOVERY FROM FACILITATION
  pp = pp - frel * (pp - ppb);

  % RESIDUAL CA++
  cares(t+1) = cares(t) + carest;
  cares(t+1) = cares(t+1) - carel * cares(t+1);
    
  % NEW RELEASE PROB
  if t < pulse
     ca_amp(t+1) = pp+cares(t+1);
     ppbase(t+1) = ppb;
     ppfac(t+1) = ca_amp(t+1);
     pprel(t+1) = 1-exp(-k*(ca * ca_amp(t+1))^exponent);
  end

end

nresps = resps/resps(1);    % normalised responses

function [Isyn,s] = synapse(s,V_pre,V_post)
% To caculate synaptic current injection at a particular time
% s,V_pre,syn.E_rev - (1,nsyn) V_post - (1,1)
global syn;

NT = syn.Tmax./( 1+exp(-syn.k_pre*(V_pre'-syn.E_pre)) );
s = s + syn.alpha*NT.*(1-s) - syn.beta*s;
s(s>1)=1;
s(s<0)=0;
Isyn = sum(syn.Gsyn*s.*(V_post-syn.E_rev));

end


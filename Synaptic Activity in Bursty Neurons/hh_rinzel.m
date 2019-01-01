    %HODGKIN_HUXLEY_RINZEL Rinzel's simplified Hodgkin-Huxley neuron model.
%   [V, W] = HODGKIN_HUXLEY_RINZEL(T,I_EXT) simulates the membrane
%   voltage V and recovery variable W of Rinzel's simplified
%   Hodgkin-Huxley neuron model in response to an external current
%   I_EXT over times T. The lengths of T and I_EXT must be the same.

%   Authors: Lev Givon and Yevgeniy Slutskiy
%   Copyright 2009-2010 Lev Givon and Yevgeniy Slutskiy

function [V, R] = hh_rinzel(t,I_ext,V_pre)
% V_pre is a 2D array of presynaptic potentials of all neurons projecting
% onto Rinzel neuron (for generality) - (nsyn, length(t))

% Assume that the time is given in seconds, and convert it to
% number of milliseconds:

dt = 1e3*diff(t(1:2));

V = zeros(size(t));
R = zeros(size(t));
V(1) = 36.7794;
R(1) = 0.7330;
s=zeros(1,size(V_pre,1));            % Initial synapse dynamics
Isyn=0;

global syn;
if(~isempty(V_pre) && length(syn.alpha)>1)
    alpha=syn.alpha;
end
    
for i = 2:length(t)
    if(~isempty(V_pre))
        if(length(syn.alpha)>1)
            syn.alpha=alpha(i);
        end
        [Isyn,s] = synapse(s,V_pre(:,i-1),V(i-1));
    end
    dx = rinzel_ode([V(i-1),R(i-1)], I_ext(i) - Isyn);
    V(i) = V(i-1) + dt*dx(1);
    R(i) = R(i-1) + dt*dx(2);
end

return

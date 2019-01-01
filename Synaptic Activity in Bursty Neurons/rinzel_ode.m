%% RINZEL_ODE
function dx = rinzel_ode(x, I_ext)
    v = x(1);

    a = exp(-(v+55)/10)-1;
    alpha_n = ([a==0 a~=0])*[0.01; -0.01*(v+55)/a];
    a = exp(-(v+40)/10)-1;
    alpha_m = ([a==0 a~=0])*[1; -0.1*(v+40)/a];
    alpha_h = 0.07*exp(-(v+65)/20) ;                     

    beta_n = 0.125*exp(-(v+65)/80); 
    beta_m = 4*exp(-(v+65)/18);
    beta_h = 1/(exp(-(v+35)/10)+1);                  

    n_infty = alpha_n/(alpha_n + beta_n);
    m_infty = alpha_m/(alpha_m + beta_m);    
    h_infty = alpha_h/(alpha_h + beta_h);

    % compute dR
    S = 1.2714;  % s = (1-h_inf(0)) / n_inf(0) 
    R_infty = S/(1+S^2)*(n_infty + S*(1-h_infty));
    tau_R = 1 + 5*exp(-(v+55)^2/55^2);   
    dR = -3*(x(2)-R_infty) / tau_R;

    % Update the ionic currents and membrane voltage:
    dV = I_ext -120*m_infty^3*(1-x(2))*(v-50) - 36*(x(2)/S).^4*(v+77) - 0.3*(v+54.387);
    dx = [dV,dR];
return
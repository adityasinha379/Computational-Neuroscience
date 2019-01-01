function [x_e,x_r,x_i] = STPlast(x_eo,x_ro,x_io,t_re,t_ei,t_ir,spike)
    x_r=x_ro-spike*x_ro/t_re+x_io/t_ir;
    x_e=x_eo+spike*x_ro/t_re-x_eo/t_ei;
    x_i=x_io+x_eo/t_ei-x_io/t_ir;
end
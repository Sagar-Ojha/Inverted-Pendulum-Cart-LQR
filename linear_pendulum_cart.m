%% Drives any IC to linearized equilibrium point, i.e., [0 0 0 0]
function xdot = linear_pendulum_cart(t,x,Ka_pair,Aa_0,Ba_0,Aa_pi,Ba_pi,Br,up,x_eq_up,x_eq_down)
    if (t < 10)
        ref_q = 0;
    elseif (t < 20)
        ref_q = 5;
    elseif (t < 30)
        ref_q = -2;
    else
        ref_q = -0;
    end

    if (up == 1)
        % x = x - x_eq_up;
        f_ext = -Ka_pair(2) * x;
        xdot = Aa_pi * x + Ba_pi .* f_ext + Br .* ref_q;
    else
        % x = x - x_eq_down;
        f_ext = -Ka_pair(2) * x;
        xdot = Aa_0 * x + Ba_0 .* f_ext + Br .* ref_q;
    end
end
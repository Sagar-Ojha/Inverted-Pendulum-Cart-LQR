%% State derivatives
function xdot = pendulum_cart_tracking(t, x, Ka_pair, x_eq_down, x_eq_up, up, m1, m2, l, k, kappa, g)
    if (t < 10)
        ref_q = 0;
    elseif (t < 30)
        ref_q = 1;
    elseif (t < 50)
        ref_q = -1;
    else
        ref_q = -0;
    end
    % State variables
    q = x(1); qdot = x(2); theta = x(3); thetadot = x(4); error_q = ref_q - q;

    x_eq = x_eq_down;
    if (3*pi/4 < x(3)) && (x(3) < 5*pi/4) && (up == 1)
        x_eq = x_eq_up;
    end
    del_x = x - x_eq;

    if (up == 1) && (3*pi/4 < theta) && (theta < 5*pi/4)
        f_ext = -Ka_pair(2,1:4) * del_x(1:4) - Ka_pair(2,5) * del_x(5);
    elseif (up == 1)
        f_ext = -Ka_pair(1,1:4) * del_x(1:4) - Ka_pair(1,5) * del_x(5);
    elseif (theta == pi)
        f_ext = -Ka_pair(1,1:4) * del_x(1:4) - Ka_pair(1,5) * del_x(5);
    else
        f_ext = -Ka_pair(2,1:4) * del_x(1:4) - Ka_pair(2,5) * del_x(5);
    end

    xdot_1 = qdot;
    xdot_3 = thetadot;
    qddot = (1/(m1 + m2 - m2 * (cos(theta))^2)) *...
            (f_ext + m2 * g * sin(theta) * cos(theta) +...
            (kappa/l) * cos(theta) * theta + k*q +...
             m2 * l * sin(theta) * thetadot^2);
    xdot_2 = qddot;
    xdot_4 = -(cos(theta)/l) * qddot - (g/l) * sin(theta) -...
             (kappa/(m2*(l^2))) * theta;
    xdot = [xdot_1; xdot_2; xdot_3; xdot_4; error_q];
end
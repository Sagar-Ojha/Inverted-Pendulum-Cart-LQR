function animate_pendulum_cart(anim_pause, x, t, m1, m2, l)
    % System Parameters
    pendulum_length = l;
    pendulum_radius = 35 * m2/1.5;
    cart_length = 0.1*m1; cart_height = 0.05*m1;

    cart_x_pos = -x(:, 1);
    cart_y_pos = zeros(length(cart_x_pos));
    pendulum_angle = x(:, 3);
    pendulum_x_pos = -pendulum_length*sin(pendulum_angle) + cart_x_pos;
    pendulum_y_pos = -pendulum_length*cos(pendulum_angle) + cart_y_pos;

    animVideo = VideoWriter("pendulumCart.avi");
    animVideo.Quality = 100;
    animVideo.FrameRate = 120;
    open(animVideo);

    figure(1);
    set(gcf, 'WindowState', 'maximized');
    % First Frame
    pendulum = plot([cart_x_pos(1) pendulum_x_pos(1)], [cart_y_pos(1) pendulum_y_pos(1)], 'Color','b', 'LineWidth', 1.5);
    xlim([min(cart_x_pos)-1 max(cart_x_pos)+1]);
    ylim([-0.6 0.6]);
    hold on;
    daspect([1 1 1]);
    pendulum_bob = plot(pendulum_x_pos(1), pendulum_y_pos(1), 'b.', 'MarkerSize', pendulum_radius);
    % cart = plot(cart_x_pos(1), cart_y_pos(1), 'Marker', 'square', 'MarkerSize', 30, 'MarkerFaceColor', [0.4940 0.1840 0.5560]);

    yline(-cart_height/2, "LineWidth", 1);
    cart_ends = [cart_x_pos(1) - cart_length/2, cart_y_pos(1) - cart_height/2, cart_length, cart_height];

    cart = rectangle('Position', cart_ends, 'Curvature', 0.2);
    cart.FaceColor = [0.4940 0.1840 0.5560];

    frame = getframe(gcf);
    writeVideo(animVideo, frame);

    % All other frames
    for i = 2:3:length(t)
        set(pendulum, 'XData', [cart_x_pos(i) pendulum_x_pos(i)]);
        set(pendulum, 'YData', [cart_y_pos(i) pendulum_y_pos(i)]);
        set(pendulum_bob, 'XData', pendulum_x_pos(i));
        set(pendulum_bob, 'YData', pendulum_y_pos(i));
        % set(cart, 'XData', cart_x_pos(i));
        % set(cart, 'YData', cart_y_pos(i));

        delete(cart);
        cart_ends = [cart_x_pos(i) - cart_length/2, cart_y_pos(i) - cart_height/2, cart_length, cart_height];

        cart = rectangle('Position', cart_ends, 'Curvature', 0.2);
        cart.FaceColor = [0.4940 0.1840 0.5560];
        pause(anim_pause);
        drawnow;
        frame = getframe(gcf);
        writeVideo(animVideo, frame);
    end
    close(animVideo);

end
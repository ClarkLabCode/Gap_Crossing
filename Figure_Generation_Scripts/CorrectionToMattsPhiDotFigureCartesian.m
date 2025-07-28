% Equation for phiDot when following curved path
phiDot = @(v_x, D, R, phi) v_x*(-D + R*sin(phi))/(D*R);
% Equation for phiDot when following straight path (R -> Inf)
phiDotStraight = @(v_x, D, phi) v_x*sin(phi)/D;

% Set parameters in arbitrary units
v_x = 0.5;
R = 1;

% Sweep x and y to make a grid
x = [-1.5*R:0.125:1.5*R];
y = [-1.5*R:0.125:1.5*R];

% Now make the grid
[x_grid, y_grid] = meshgrid(x,y);

% Initialize phiDot in grid form
phiDotOnCartesianGridCurved = zeros(size(x_grid));
phiDotOnCartesianGridStraight = zeros(size(x_grid));

% Go through all the grid points and compute phiDot for both straight and
% curved paths
for i_x = 1:length(x)
    for i_y = 1:length(y)
        % Skip the origin since it blows up the calculation
        if (x(i_x) == 0) && (y(i_y) == 0)
            continue
        end
        phiDotOnCartesianGridCurved(i_y, i_x) = phiDot(v_x, sqrt(x(i_x)^2 + y(i_y)^2), R, atan2(y(i_y),x(i_x)));
        phiDotOnCartesianGridStraight(i_y, i_x) = phiDotStraight(v_x, sqrt(x(i_x)^2 + y(i_y)^2), atan2(y(i_y),x(i_x)));
    end
end

% Reshape matrices to be vectors so it's easier to work with
x_grid = reshape(x_grid, [], 1);
y_grid = reshape(y_grid, [], 1);
phiDotOnCartesianGridCurved = reshape(phiDotOnCartesianGridCurved, [], 1);
phiDotOnCartesianGridStraight = reshape(phiDotOnCartesianGridStraight, [], 1);

% Remove all things that are supposed to be 0 but weren't because of pi
% rounding errors in the trig functions
phiDotOnCartesianGridCurved(abs(phiDotOnCartesianGridCurved) < 10^(-10)) = 0;
phiDotOnCartesianGridStraight(abs(phiDotOnCartesianGridStraight) < 10^(-10)) = 0;

% Separate out progressive and regressive motion for curved path
x_grid_prog_Curved = x_grid(phiDotOnCartesianGridCurved > 0);
y_grid_prog_Curved = y_grid(phiDotOnCartesianGridCurved > 0);
phiDot_grid_prog_Curved = phiDotOnCartesianGridCurved(phiDotOnCartesianGridCurved > 0);
x_grid_reg_Curved = x_grid(phiDotOnCartesianGridCurved < 0);
y_grid_reg_Curved = y_grid(phiDotOnCartesianGridCurved < 0);
phiDot_grid_reg_Curved = phiDotOnCartesianGridCurved(phiDotOnCartesianGridCurved < 0);

% Separate out progressive and regressive motion for straight path
x_grid_prog_Straight = x_grid(phiDotOnCartesianGridStraight > 0);
y_grid_prog_Straight = y_grid(phiDotOnCartesianGridStraight > 0);
phiDot_grid_prog_Straight = phiDotOnCartesianGridStraight(phiDotOnCartesianGridStraight > 0);
x_grid_reg_Straight = x_grid(phiDotOnCartesianGridStraight < 0);
y_grid_reg_Straight = y_grid(phiDotOnCartesianGridStraight < 0);
phiDot_grid_reg_Straight = phiDotOnCartesianGridStraight(phiDotOnCartesianGridStraight < 0);

% Make a figure with points on the grid and color coded by direction of
% motion for a curved trajectory
figure
hold on
% First plot the dashed lines that show the radius of curvature
plot([R + R*sinpi(3/2+1/6),R], [R*cospi(3/2+1/6),0], 'k--','LineWidth',0.5); 
plot([0,R], [0,0], 'k--','LineWidth',0.5);
% Plot the direction of motion by just plotting the sign (prog = red, reg = blue)
scatter(y_grid, x_grid, 10, sign(phiDotOnCartesianGridCurved), 'filled');
colorbar
colormap(redblue());
% Plot a small segment of the trajectory arc
plot(R + R*sinpi(3/2+[0:0.0001:1/6]), R*cospi(3/2+[0:0.0001:1/6]), 'k-','LineWidth',2);
% Replace the line above with this one below if you want to plot the entire trajectory
% plot(R + R*sinpi(0:0.0001:2), R*cospi(0:0.0001:2), 'k--');
% Plot the boundary of where there is no motion and do it in purple
plot(R/2 + R/2*sinpi(0:0.0001:2), R/2*cospi(0:0.0001:2), '--', 'Color', [1, 0, 1]/2,'LineWidth',2);
% Plot the center of the radius of curvature
plot(R,0,'ko','MarkerSize',6,'MarkerFaceColor','k')
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
title('Prog vs Reg Motion in Space for Curved Trajectory')

% Make a figure with points on the grid and color coded by direction of
% motion for a straight trajectory
figure
hold on
% Plot the direction of motion by just plotting the sign (prog = red, reg = blue)
scatter(y_grid(abs(phiDotOnCartesianGridStraight) > 10^(-10)), x_grid(abs(phiDotOnCartesianGridStraight) > 10^(-10)), 10, sign(phiDotOnCartesianGridStraight(abs(phiDotOnCartesianGridStraight) > 10^(-10))), 'filled');
colorbar
colormap(redblue());
% Plot a small segment of the trajectory
plot([0,0],[0,R/2], 'k-','LineWidth',2);
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
title('Prog vs Reg Motion in Space for Straight Trajectory')

% Make a figure with arrows on the grid to show direction of phi dot fpr a
% curved path
figure
hold on
% First plot the dashed lines that show the radius of curvature
plot([R + R*sinpi(3/2+1/6),R], [R*cospi(3/2+1/6),0], 'k--','LineWidth',0.5); 
plot([0,R], [0,0], 'k--','LineWidth',0.5);
% Plot the arrows showing the direction of phi dot
q_reg = quiver(y_grid_reg_Curved,...
       x_grid_reg_Curved,...
       phiDot_grid_reg_Curved.*sin(pi/2+atan2(y_grid_reg_Curved,x_grid_reg_Curved)),...
       phiDot_grid_reg_Curved.*cos(pi/2+atan2(y_grid_reg_Curved,x_grid_reg_Curved)),...
       1.5,'b','Alignment','center');
q_prog = quiver(y_grid_prog_Curved,...
       x_grid_prog_Curved,...
       phiDot_grid_prog_Curved.*sin(pi/2+atan2(y_grid_prog_Curved,x_grid_prog_Curved)),...
       phiDot_grid_prog_Curved.*cos(pi/2+atan2(y_grid_prog_Curved,x_grid_prog_Curved)),...
       1.5,'r','Alignment','center');
% Plot a small segment of the trajectory arc
plot(R + R*sinpi(3/2+[0:0.0001:1/6]), R*cospi(3/2+[0:0.0001:1/6]), 'k-','LineWidth',2);
% Replace the line above with this one below if you want to plot the entire trajectory
% plot(R + R*sinpi(0:0.0001:2), R*cospi(0:0.0001:2), 'k--');
% Plot the boundary of where there is no motion and do it in purple
plot(R/2 + R/2*sinpi(0:0.0001:2), R/2*cospi(0:0.0001:2), '--', 'Color', [1, 0, 1]/2,'LineWidth',2);
% Plot the center of the radius of curvature
plot(R,0,'ko','MarkerSize',6,'MarkerFaceColor','k')
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
hold off
title('Phi Dot in Space for Curved Trajectory')

% Make a figure with arrows on the grid to show direction of phi dot fpr a
% straight path
figure
hold on
% Plot the arrows showing the direction of phi dot
q_reg = quiver(y_grid_reg_Straight,...
       x_grid_reg_Straight,...
       phiDot_grid_reg_Straight.*sin(pi/2+atan2(y_grid_reg_Straight,x_grid_reg_Straight)),...
       phiDot_grid_reg_Straight.*cos(pi/2+atan2(y_grid_reg_Straight,x_grid_reg_Straight)),...
       1.5,'b','Alignment','center');
q_prog = quiver(y_grid_prog_Straight,...
       x_grid_prog_Straight,...
       phiDot_grid_prog_Straight.*sin(pi/2+atan2(y_grid_prog_Straight,x_grid_prog_Straight)),...
       phiDot_grid_prog_Straight.*cos(pi/2+atan2(y_grid_prog_Straight,x_grid_prog_Straight)),...
       1.5,'r','Alignment','center');
% Plot a small segment of the trajectory
plot([0,0],[0,R/2], 'k-','LineWidth',2);
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
hold off
title('Phi Dot in Space for Straight Trajectory')

% Make a figure with arrows on the grid to show direction of motion fpr a
% curved path
figure
hold on
% First plot the dashed lines that show the radius of curvature
plot([R + R*sinpi(3/2+1/6),R], [R*cospi(3/2+1/6),0], 'k--','LineWidth',0.5); 
plot([0,R], [0,0], 'k--','LineWidth',0.5); 
% Plot the arrows showing the direction of motion
q_reg = quiver(y_grid_reg_Curved,...
       x_grid_reg_Curved,...
       -v_x/R*x_grid_reg_Curved,...
       -v_x + v_x/R*y_grid_reg_Curved,...
       'b','Alignment','center');
q_prog = quiver(y_grid_prog_Curved,...
       x_grid_prog_Curved,...
       -v_x/R*x_grid_prog_Curved,...
       -v_x + v_x/R*y_grid_prog_Curved,...
       'r','Alignment','center');
% Plot a small segment of the trajectory arc
plot(R + R*sinpi(3/2+[0:0.0001:1/6]), R*cospi(3/2+[0:0.0001:1/6]), 'k-','LineWidth',2);
% Replace the line above with this one below if you want to plot the entire trajectory
% plot(R + R*sinpi(0:0.0001:2), R*cospi(0:0.0001:2), 'k--');
% Plot the boundary of where there is no motion and do it in purple
plot(R/2 + R/2*sinpi(0:0.0001:2), R/2*cospi(0:0.0001:2), '--', 'Color', [1, 0, 1]/2,'LineWidth',2);
% Plot the center of the radius of curvature
plot(R,0,'ko','MarkerSize',6,'MarkerFaceColor','k')
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
hold off
title('Velocity in Space for Curved Trajectory')

% Make a figure with arrows on the grid to show direction of motion fpr a
% straight path
figure
hold on
% Plot the arrows showing the direction of motion
q_reg = quiver(y_grid_reg_Straight,...
       x_grid_reg_Straight,...
       zeros(size(x_grid_reg_Straight)),...
       -v_x*ones(size(x_grid_reg_Straight)),...
       'b','Alignment','center');
q_prog = quiver(y_grid_prog_Straight,...
       x_grid_prog_Straight,...
       zeros(size(x_grid_reg_Straight)),...
       -v_x*ones(size(x_grid_reg_Straight)),...
       'r','Alignment','center');
% Plot a small segment of the trajectory
plot([0,0],[0,R/2], 'k-','LineWidth',2);
hold off
% Remove the axes colors to make it prettier
set(gca,'XColor', 'none','YColor','none')
set(gca, 'color', 'none');
% Fix the aspect ratio
pbaspect([1,1,1])
hold off
title('Velocity in Space for Straight Trajectory')
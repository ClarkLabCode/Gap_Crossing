tic

% Make a function that computes phiDot given all parameters
phiDot = @(v_x, phi, omega, D) -omega + v_x*sin(phi)/D;

% Option that allows you to choose to only compute R^2 for the ratios of
% phiDots to save time (since we know this is the main one in parallax)
computeOnlyRatioTerms = 1;

% Option that allows you to choose whether to normalize occurrenceFreqDens
% to have a sum of 1 when plotting
% Note that if this is set to 0, then the fracInfoCaptured is weighed by
% the size of the bins in log-log space (i.e., just taking the average in
% log-log space of the R^2 value in the quadrant)
% If this option is toggled to 1, the frequency of occurrences will be used
% as the weight metric for computing the info captured in each quadrant
normOccurrenceFreqDensTo1 = 1;

% Number of times to simulate
% numSamples = 100000000;
numSamples = 10000000;

% Choose the number of histogram bins before and after ratios of 1
orderMagBelow1 = 2;
orderMagAbove1 = 1;

% Size of the receptive field of neuron of interest (in radians)
sizeOfRF = pi/6;

% Minimum distance from center line to keep phi (because simulations too
% close to the center line on either side lead to nasty behavior
minAngleFromCenter = pi/36;

% Option to control whether to fix phi1 = phi2 = phi_RF
% If this option is toggled to 1, then phi1 = phi2 = phi_RF
% Otherwise, phi1 & phi2 each randomly pulled from phi_RF +/- 1/2*sizeOfRF
phi1and2atCenterOfRF = 1;

% Make a vector of all the different parameters we want to plot
min_v_x = [5];
v_x_meanVec = [20];
phi_constVec = -pi*[1/2];
minD1Vec = [0];
min_D2_from_D1 = [0];
max_D1Vec = [25];
max_D2_from_D1Vec = [25];
omega_sigmaVec = [5*pi/6];

% Lots of meaningless warnings pop-up in this script, so temporarily
% disable them (they are turned back on at end of script)
warning('off','all')

% Loop through all the different parameter vectors
for i_v_x = 1:length(v_x_meanVec)
for i_phi = 1:length(phi_constVec)
for i_D1 = 1:length(max_D1Vec)
for i_D2 = 1:length(max_D2_from_D1Vec)
for i_omega = 1:length(omega_sigmaVec)
    % Choose the parameters that are set as constants
    v_x_mean = v_x_meanVec(i_v_x);
    phi_const = phi_constVec(i_phi);
    max_D1 = max_D1Vec(i_D1);
    max_D2_from_D1 = max_D2_from_D1Vec(i_D2);
    omega_sigma = omega_sigmaVec(i_omega);

    % Randomly select distances for object 1 and 2
    D1 = max_D1*rand(numSamples,1) + minD1Vec;
    D2 = (1+min_D2_from_D1)*D1 + max_D2_from_D1*rand(numSamples,1);

    % Randomly select angular velocities (omega) for each sample
    % omega = omega_sigma*(2*(rand(numSamples,1)-0.5)); % Pull omega uniformly
    omega = sign(phi_const)*abs(normrnd(0,omega_sigma,numSamples,1)); % Pull omega from a Gaussian
    
    % Initialize phiDot vectors for object 1 and 2
    phiDot1 = zeros(size(D1));
    phiDot2 = zeros(size(D2));

    % Initialize phi vectors for object 1 and 2
    phi1 = zeros(size(D1));
    phi2 = zeros(size(D2));

    % Randomly select v_x and phi_RF vectors within reasonable ranges
    % v_x = v_x_mean*ones(numSamples,1); % Pull v_x as just a constant
    % v_x = v_x_mean/2 + v_x_mean*rand(numSamples,1); % Pull v_x uinformly
    v_x = normrnd((v_x_mean-min_v_x),v_x_mean/4,numSamples,1); % Pull v_x from a Gaussian
    v_x(v_x<0) = 0; % Only keep this uncommented if pulling from Gaussian
    v_x = v_x + min_v_x; % Only keep this uncommented if pulling from Gaussian
    phi_RF = phi_const*ones(numSamples,1); % Pull phi_RF as just a constant

    Rc = abs(v_x./omega);
    
    % Compute phiDot for each object given the random parameters
    for i_sample = 1:numSamples
        % Let's first randomly select phi1 and phi2 for each random
        % sample so that we can make sure the randomly chosen angle is
        % outside of the very noisy region along the center line.
        % First choose phi2 at random within the receptive field
        % centered at the chosen phi_const param and keep randomly 
        % selecting until we get a phi2 value sufficiently away from 
        % the center line
        % Note the (1-phi1and2atCenterOfRF) term which controls whether 
        % phi2 gets selected at the center of the receptive field or not
        while abs(phi2(i_sample)) < minAngleFromCenter
            phi2(i_sample) = sign(phi_const)*abs(phi_RF(i_sample) + (1-phi1and2atCenterOfRF)*sizeOfRF*(rand(1)-0.5));
        end
        % Now choose phi1 values and keep doing so sufficiently far
        % from the center line
        while abs(phi1(i_sample)) < minAngleFromCenter
            % If phi1and2atCenterOfRF was selected, then choose phi1 = phi2
            % since we already set phi2 = phi_RF
            if (phi1and2atCenterOfRF == 1)
                phi1(i_sample) = phi2(i_sample);
            % Otherwise, randomly choose phi1 within the receptive field
            else
                phi1(i_sample) = sign(phi_const)*abs(phi_RF(i_sample) + sizeOfRF*(rand(1)-0.5));
            end
        end
        % Now compute phiDot for object 2
        phiDot2(i_sample) = ...
            phiDot(v_x(i_sample), phi2(i_sample), omega(i_sample), D2(i_sample));
        % Now compute phiDot for object 1
        phiDot1(i_sample) = ...
            phiDot(v_x(i_sample), phi1(i_sample), omega(i_sample), D1(i_sample));
    end

    % Compute D1/Rc and D2/Rc to use for later when plotting since these
    % form the boundaries of the different cases of prog / reg motion of
    % obj 1 / obj 2
    D1overRc = D1./Rc;
    D2overRc = D2./Rc;

    % Compute the max values of D1/Rc and D2/Rc because these will be used
    % as bounds for the histogram being made later
    maxD1overRc = max(D1overRc);
    maxD2overRc = max(D2overRc);
    
    % Now form the edges of the histogram bins
    xBinLims = 10.^[-orderMagBelow1:0.1:orderMagAbove1]; % Log-log bins
    yBinLims = 10.^[-orderMagBelow1:0.1:orderMagAbove1]; % Log-log bins
    % xBinLims = 10^(-orderMagBelow1):10^(-orderMagBelow1):10^(orderMagAbove1); % Linear bins
    % yBinLims = 10^(-orderMagBelow1):10^(-orderMagBelow1):10^(orderMagAbove1); % Linear bins

    % Initialize matrices that will hold the R^2 values which measure the
    % variance of D1/D2 captured by {phiDot1, phiDot2, phiDot2/phiDot1}
    if computeOnlyRatioTerms == 0
        phiDot1ToRatioRSquared = zeros(length(yBinLims),length(xBinLims));
        phiDot2ToRatioRSquared = zeros(length(yBinLims),length(xBinLims));
    end
    ratioToRatioRSquared = zeros(length(yBinLims),length(xBinLims));

    % Initialize matrices that will hold the R^2 values which measure the
    % variance of D1 captured by {phiDot1, phiDot2, phiDot2/phiDot1}
    if computeOnlyRatioTerms == 0
        phiDot1ToD1RSquared = zeros(length(yBinLims),length(xBinLims));
        phiDot2ToD1RSquared = zeros(length(yBinLims),length(xBinLims));
    end
    ratioToD1RSquared = zeros(length(yBinLims),length(xBinLims));

    % Initialize a matrix that holds the frequency with which each bin of
    % the histogram occurs and a matrix that holds this frequency density
    % The frequency matrix is useful if using linear bins, whereas the
    % frequency density matrix is useful when using log-log bins
    occurrenceFreq = zeros(length(yBinLims),length(xBinLims));
    occurrenceFreqDens = zeros(length(yBinLims),length(xBinLims));

    % Initialize a matrix that holds the size of the bins in log space
    sizeOfBinsInLogSpace = zeros(size(meshgrid(xBinLims,yBinLims)));

    % Now go through and compute the R^2 values for all the points within
    % each bin of the 2D histogram
    for i_D1overRc = 1:(length(xBinLims)-1)
        for i_D2overRc = 1:(length(yBinLims)-1)
            % Create a logical array for which values are in the hist bin
            inHistBin = ...
                (D1overRc > xBinLims(i_D1overRc)) & (D1overRc < xBinLims(i_D1overRc+1)) & ... % D1/Rc within range
                (D2overRc > yBinLims(i_D2overRc)) & (D2overRc < yBinLims(i_D2overRc+1));      % D2/Rc within range
            % Now grab all the things we need for values in the hist bin
            % and hold them in temporary variables
            tempD1 = D1(inHistBin);
            tempD2 = D2(inHistBin);
            tempPhiDot1 = phiDot1(inHistBin);
            tempPhiDot2 = phiDot2(inHistBin);

            % If there are fewer than ten points in the bin, skip it
            % (otherwise you'll always be able to draw a "very good" line)
            if sum(inHistBin) < 10
                continue
            end

            % Compute R^2 of the best fit lines
            % Fit a line for D1/D2 vs {phiDot1, phiDot2, phiDot1/phiDot2}
            if computeOnlyRatioTerms == 0
                [~,temp_mdl_phiDot1_to_ratio] = fit(tempD1./tempD2,tempPhiDot1,'poly1');
                [~,temp_mdl_phiDot2_to_ratio] = fit(tempD1./tempD2,tempPhiDot2,'poly1');
            end
            [~,temp_mdl_ratio_to_ratio] = fit(tempD1./tempD2,tempPhiDot2./tempPhiDot1,'poly1');

            % Temporarily store R^2 for D1/D2 captured by {phiDot1, phiDot2, phiDot1/phiDot2}
            if computeOnlyRatioTerms == 0
                tempPhiDotToRatio1RSquared = temp_mdl_phiDot1_to_ratio.rsquare;
                tempPhiDot2ToRatioRSquared = temp_mdl_phiDot2_to_ratio.rsquare;
            end
            tempRatioToRatioRSquared = temp_mdl_ratio_to_ratio.rsquare;

            % Compute R^2 of the best fit lines for the inverses
            % Now fit a line for D1/D2 vs {1/phiDot1, 1/phiDot2, phiDot2/phiDot1}
            % if computeOnlyRatioTerms == 0    
            %     [~,temp_mdl_phiDot1_to_ratio] = fit(tempD1./tempD2,1./(tempPhiDot1),'poly1');
            %     [~,temp_mdl_phiDot2_to_ratio] = fit(tempD1./tempD2,1./(tempPhiDot2),'poly1');
            % end
            % [~,temp_mdl_ratio_to_ratio] = fit(tempD1./tempD2,1./(tempPhiDot2./tempPhiDot1),'poly1');

            % Now take the better of the two R^2 (linear vs inverse)
            % if computeOnlyRatioTerms == 0
            %     tempPhiDotToRatio1RSquared = max(temp_mdl_phiDot1_to_ratio.rsquare,tempPhiDotToRatio1RSquared);
            %     tempPhiDot2ToRatioRSquared = max(temp_mdl_phiDot2_to_ratio.rsquare,tempPhiDot2ToRatioRSquared);
            % end
            % tempRatioToRatioRSquared = max(temp_mdl_ratio_to_ratio.rsquare,tempRatioToRatioRSquared);

            % Fill in the matrices appropriately now
            if computeOnlyRatioTerms == 0
                phiDot1ToRatioRSquared(i_D2overRc,i_D1overRc) = tempPhiDotToRatio1RSquared;
                phiDot2ToRatioRSquared(i_D2overRc,i_D1overRc) = tempPhiDot2ToRatioRSquared;
            end
            ratioToRatioRSquared(i_D2overRc,i_D1overRc) = tempRatioToRatioRSquared;





            if (xBinLims(i_D1overRc+1) == 0.1) && (yBinLims(i_D2overRc) == 0.1)
                figure
                scatter(tempD1./tempD2,tempPhiDot2./tempPhiDot1,10,[0.7,0.7,0.7],'filled')
            elseif (xBinLims(i_D1overRc+1) == 1) && (yBinLims(i_D2overRc) == 1)
                figure
                scatter(tempD1./tempD2,tempPhiDot2./tempPhiDot1,10,[0.7,0.7,0.7],'filled')
            end






            % Now fit a line for D1 vs {phiDot1, phiDot2, phiDot1/phiDot2}
            if computeOnlyRatioTerms == 0
                [~,temp_mdl_phiDot1_to_D1] = fit(tempD1,tempPhiDot1,'poly1');
                [~,temp_mdl_phiDot2_to_D1] = fit(tempD1,tempPhiDot2,'poly1');
            end
            [~,temp_mdl_ratio_to_D1] = fit(tempD1,tempPhiDot2./tempPhiDot1,'poly1');

            % Temporarily store R^2 for D1 captured by {phiDot1, phiDot2, phiDot1/phiDot2}
            if computeOnlyRatioTerms == 0
                tempPhiDotToD11RSquared = temp_mdl_phiDot1_to_D1.rsquare;
                tempPhiDot2ToD1RSquared = temp_mdl_phiDot2_to_D1.rsquare;
            end
            tempRatioToD1RSquared = temp_mdl_ratio_to_D1.rsquare;

            % Now fit a line for D1 vs {1/phiDot1, 1/phiDot2, phiDot2/phiDot1}
            if computeOnlyRatioTerms == 0
                [~,temp_mdl_phiDot1_to_D1] = fit(tempD1,1./(tempPhiDot1),'poly1');
                [~,temp_mdl_phiDot2_to_D1] = fit(tempD1,1./(tempPhiDot2),'poly1');
            end
            [~,temp_mdl_ratio_to_D1] = fit(tempD1,1./(tempPhiDot2./tempPhiDot1),'poly1');

            % Now take the better of the two R^2 (linear vs inverse)
            if computeOnlyRatioTerms == 0
                tempPhiDotToD11RSquared = max(temp_mdl_phiDot1_to_D1.rsquare,tempPhiDotToD11RSquared);
                tempPhiDot2ToD1RSquared = max(temp_mdl_phiDot2_to_D1.rsquare,tempPhiDot2ToD1RSquared);
            end
            tempRatioToD1RSquared = max(temp_mdl_ratio_to_D1.rsquare,tempRatioToD1RSquared);

            % Fill in the matrices appropriately now
            if computeOnlyRatioTerms == 0
                phiDot1ToD1RSquared(i_D2overRc,i_D1overRc) = tempPhiDotToD11RSquared;
                phiDot2ToD1RSquared(i_D2overRc,i_D1overRc) = tempPhiDot2ToD1RSquared;
            end
            ratioToD1RSquared(i_D2overRc,i_D1overRc) = tempRatioToD1RSquared;

            % Compute the density of samples that are in this histogram bin
            % (Note that this is normalized by the area of the bins because
            % the binning being used for the histogram is not square and
            % the area of each bin varies largely across the phase space)
            % When along the diagonal, the bins are trapezoids
            if i_D1overRc == i_D2overRc
                occurrenceFreqDens(i_D2overRc,i_D1overRc) = (sum(inHistBin)/length(inHistBin))/...
                             ((1/2)*(xBinLims(i_D1overRc+1)-xBinLims(i_D1overRc))*... % Length in x of trapezoid (1/2 from trapezoid formula)
                                    ((yBinLims(i_D2overRc+1)-xBinLims(i_D1overRc+1))+... % Length in y of right side of trapezoid
                                     (yBinLims(i_D2overRc+1)-xBinLims(i_D1overRc)))); % Length in y of left side of trapezoid
                % Also compute the size of the bins in log-log space
                sizeOfBinsInLogSpace(i_D2overRc,i_D1overRc) = ...
                             ((1/2)*(log10(xBinLims(i_D1overRc+1))-log10(xBinLims(i_D1overRc)))*... % Length in x of trapezoid (1/2 from trapezoid formula)
                                   ((log10(yBinLims(i_D2overRc+1))-log10(xBinLims(i_D2overRc+1)))+... % Length in y of right side of trapezoid
                                    (log10(yBinLims(i_D1overRc+1))-log10(xBinLims(i_D2overRc))))); % Length in y of left side of trapezoid
            % When not along the diagonal, the bins are rectangles
            else
                occurrenceFreqDens(i_D2overRc,i_D1overRc) = (sum(inHistBin)/length(inHistBin))/...
                             ((yBinLims(i_D2overRc+1)-yBinLims(i_D2overRc))*... % Length in y of rectangle
                              (xBinLims(i_D1overRc+1)-xBinLims(i_D1overRc))); % Length in x of rectangle
                % Also compute the size of the bins in log-log space
                sizeOfBinsInLogSpace(i_D2overRc,i_D1overRc) = ...
                             ((log10(yBinLims(i_D2overRc+1))-log10(yBinLims(i_D2overRc)))*... % Length in y of rectangle (in log)
                              (log10(xBinLims(i_D1overRc+1))-log10(xBinLims(i_D1overRc)))); % Length in x of rectangle (in log)
            end
            % Now compute the occurrence frequency (not the density) to be
            % used for the information integration metric (only if bins are
            % linear rather than log-log)
            occurrenceFreq(i_D2overRc,i_D1overRc) = (sum(inHistBin)/length(inHistBin));
        end
    end

    % Zero out the size in log space of the bins outside the requested ranges
    sizeOfBinsInLogSpace((yBinLims >= 10^orderMagAbove1),:) = 0;
    sizeOfBinsInLogSpace(:,(xBinLims >= 10^orderMagAbove1)) = 0;

    % If option was selected to normalize the sum of occurrenceFreqDens to
    % 1 for plotting, then do that
    if normOccurrenceFreqDensTo1 == 1
        % Normalize occurrenceFreq such that its sum is 1
        occurrenceFreqDensForPlot = occurrenceFreqDens/sum(occurrenceFreqDens,'all');
    % Otherwise, just set the occurrenceFreqDensForPlot equal to the
    % non-normalized one
    else
        occurrenceFreqDensForPlot = occurrenceFreqDens;
        % occurrenceFreqDensForPlot = occurrenceFreq; % Useful for linear bins
    end

    % First we plot the histogram for frequency of occurrences
    % Open a new full screen figure window
    figure('WindowState','maximized')
    hold on
    % Go through the bins of the 2D histogram and start plotting
    for i_x = 1:(length(xBinLims)-1)
        for i_y = 1:(length(yBinLims)-1)
            % Skip all the iterations below the line y = x
            if i_x > i_y
                continue
            end
            % Along the diagonal, the patches should be trapezoids
            % Note that y values along the line are chosen via the x
            % values to ensure that you get them along y = x
            if i_x == i_y
                patch([xBinLims(i_x),xBinLims(i_x),xBinLims(i_x+1),xBinLims(i_x+1)],...
                      [xBinLims(i_x),yBinLims(i_y+1),yBinLims(i_y+1),xBinLims(i_x+1)],...
                      [occurrenceFreqDensForPlot(i_y,i_x);occurrenceFreqDensForPlot(i_y+1,i_x);occurrenceFreqDensForPlot(i_y+1,i_x+1);occurrenceFreqDensForPlot(i_y+1,i_x+1)],... % 
                      'EdgeColor','none')
            % Away from the diagonals, the patches should be rectangles
            else
                patch([xBinLims(i_x),xBinLims(i_x),xBinLims(i_x+1),xBinLims(i_x+1)],...
                      [yBinLims(i_y),yBinLims(i_y+1),yBinLims(i_y+1),yBinLims(i_y)],...
                      [occurrenceFreqDensForPlot(i_y,i_x);occurrenceFreqDensForPlot(i_y+1,i_x);occurrenceFreqDensForPlot(i_y+1,i_x+1);occurrenceFreqDensForPlot(i_y,i_x+1)],... % 
                      'EdgeColor','none')
            end
        end
    end
    
    % Temporarily set any bins at the edges of the plot outside of our 
    % region to 0 to avoid picking up any color values we're not plotting
    tempOccurrenceFreqDensForPlot = occurrenceFreqDensForPlot;
    tempOccurrenceFreqDensForPlot((yBinLims >= 10^orderMagAbove1),:) = 0;
    tempOccurrenceFreqDensForPlot(:,(xBinLims >= 10^orderMagAbove1)) = 0;
    
    % Add a color bar that ranges from min to max of occurrenceFreq
    colorbar
    clim([min(tempOccurrenceFreqDensForPlot(tempOccurrenceFreqDensForPlot>0),[],'all'),...
          max(tempOccurrenceFreqDensForPlot(tempOccurrenceFreqDensForPlot>0),[],'all')])
    set(gca, 'ColorScale', 'log')
    colormap jet
    % Add dashed lines at D1/Rc = 1, D2/Rc = 1, and D1/Rc = D2/Rc
    dashLineVec = [10^(-orderMagBelow1):10^(orderMagAbove1)];
    plot(dashLineVec,ones(size(dashLineVec)),'k--','LineWidth',2);
    plot(1:maxD1overRc,ones(size(1:maxD1overRc)),'k--','LineWidth',2)
    plot(ones(size(dashLineVec)),dashLineVec,'k--','LineWidth',2);
    plot(ones(size(1:maxD2overRc)),1:maxD2overRc,'k--','LineWidth',2)
    plot(10^(-orderMagBelow1):maxD2overRc,10^(-orderMagBelow1):maxD2overRc,'k--','LineWidth',2);
    hold off
    % Take care of axes, limits, labels, titles, etc
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xticks(10.^[-orderMagBelow1:orderMagAbove1])
    xticklabels(10.^[-orderMagBelow1:orderMagAbove1])
    yticks(10.^[-orderMagBelow1:orderMagAbove1])
    yticklabels(10.^[-orderMagBelow1:orderMagAbove1])
    pbaspect([1,1,1])
    xlim([10^(-orderMagBelow1),10^(orderMagAbove1)])
    ylim([10^(-orderMagBelow1),10^(orderMagAbove1)])
    xlabel('D_1/R_C')
    ylabel('D_2/R_C')
    % Print the title of what the plot is showing
    title('Frequency of Occurrences','Interpreter','latex')
    % Print the subtitle with what the parameters used were
    if phi1and2atCenterOfRF == 1
        subtitle(['$\bar{v}_{x} = $', num2str(v_x_mean), ' mm/s, ',...
                  '$\phi_{RF, cent}$ = -(1/',num2str(-1/(phi_const/pi)),')*$\pi$, ',...
                  '$\omega_{\sigma}$ = ', num2str(omega_sigma*180/pi), ' $^{\circ}$/s, ',...
                  '$D_{1,max}$ = ',num2str(max_D1), ' mm, ',...
                  '$D_{2,max}$ = ',num2str(max_D2_from_D1 + max_D1), ' mm, '...
                  '$\phi_1 = \phi_2 = \phi_{RF, cent}$'],'Interpreter','latex')
    else
        subtitle(['$\bar{v}_{x} = $', num2str(v_x_mean), ' mm/s, ',...
                  '$\phi_{RF, cent}$ = -(1/',num2str(-1/(phi_const/pi)),')*$\pi$, ',...
                  '$\omega_{\sigma}$ = ', num2str(omega_sigma*180/pi), ' $^{\circ}$/s, ',...
                  '$D_{1,max}$ = ',num2str(max_D1), ' mm, ',...
                  '$D_{2,max}$ = ',num2str(max_D2_from_D1 + max_D1), ' mm, '...
                  '$\phi_1 \neq \phi_2 \neq \phi_{RF, cent}$'],'Interpreter','latex')
    end

    % Add text to label the different regions of prog/reg bar/BG in the plot
    text(10^(-orderMagBelow1/2),10^(-1/8),'Prog Bar / Prog BG','HorizontalAlignment','center','FontSize',18,'Color','k');
    text(10^(-orderMagBelow1/2),10^(1/4),'Prog Bar / Reg BG','HorizontalAlignment','center','FontSize',18,'Color','k');
    text(10^(2/5),10^(orderMagAbove1)*10^(-1/8),'Reg Bar / Reg BG','HorizontalAlignment','center','FontSize',18,'Color','w');

    % Next we plot the histogram for phiDot2/phiDot1 vs D1/D2 and print the
    % title of what the plot is showing
    makePlotForPhiDotSims(xBinLims,yBinLims,ratioToRatioRSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
    title('Variance of $D_1/D_2$ accounted for by $\dot{\phi_2}/\dot{\phi_1}$','Interpreter','latex')

    % If the option is toggled to not only compute the ratio terms, do the rest
    if computeOnlyRatioTerms == 0
        % Next we plot the histogram for phiDot1 vs D1/D2 and print the
        % title of what the plot is showing
        makePlotForPhiDotSims(xBinLims,yBinLims,phiDot1ToRatioRSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
        title('Variance of $D_1/D_2$ accounted for by $\dot{\phi_1}$','Interpreter','latex')
    
        % Next we plot the histogram for phiDot2 vs D1/D2 and print the
        % title of what the plot is showing
        makePlotForPhiDotSims(xBinLims,yBinLims,phiDot2ToRatioRSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
        title('Variance of $D_1/D_2$ accounted for by $\dot{\phi_2}$','Interpreter','latex')

        % Next we plot the histogram for phiDot2/phiDot1 vs D1 and print the
        % title of what the plot is showing
        makePlotForPhiDotSims(xBinLims,yBinLims,ratioToD1RSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
        title('Variance of $D_1$ accounted for by $\dot{\phi_2}/\dot{\phi_1}$','Interpreter','latex')

        % Next we plot the histogram for phiDot1 vs D1 and print the
        % title of what the plot is showing
        makePlotForPhiDotSims(xBinLims,yBinLims,phiDot1ToD1RSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
        title('Variance of $D_1$ accounted for by $\dot{\phi_1}$','Interpreter','latex')
    
        % Next we plot the histogram for phiDot2 vs D1 and print the
        % title of what the plot is showing
        makePlotForPhiDotSims(xBinLims,yBinLims,phiDot2ToD1RSquared,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot);
        title('Variance of $D_1$ accounted for by $\dot{\phi_2}$','Interpreter','latex')
    
    end
end
end
end
end
end

% Turn all the warnings back on
warning('on','all')

toc

% Define a function that does all the plotting for each of the plots we
% want to make (so that we don't have to copy this a ton of times)
function makePlotForPhiDotSims(xBinLims,yBinLims,RSquaredUsed,orderMagBelow1,orderMagAbove1,maxD1overRc,maxD2overRc,sizeOfBinsInLogSpace,phi1and2atCenterOfRF,v_x_mean,phi_const,omega_sigma,max_D1,max_D2_from_D1,normOccurrenceFreqDensTo1,occurrenceFreqDensForPlot)

% Open a new full screen figure window
    figure('WindowState','maximized')
    hold on
    % Go through the bins of the 2D histogram and start plotting
    for i_x = 1:(length(xBinLims)-1)
        for i_y = 1:(length(yBinLims)-1)
            % Skip all the iterations below the line y = x
            if i_x > i_y
                continue
            end
            % Along the diagonal, the patches should be trapezoids
            % Note that y values along the line are chosen via the x
            % values to ensure that you get them along y = x
            if i_x == i_y
                patch([xBinLims(i_x),xBinLims(i_x),xBinLims(i_x+1),xBinLims(i_x+1)],...
                      [xBinLims(i_x),yBinLims(i_y+1),yBinLims(i_y+1),xBinLims(i_x+1)],...
                      [RSquaredUsed(i_y,i_x)*occurrenceFreqDensForPlot(i_y,i_x)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y+1,i_x)*occurrenceFreqDensForPlot(i_y+1,i_x)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y+1,i_x+1)*occurrenceFreqDensForPlot(i_y+1,i_x+1)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y+1,i_x+1)*occurrenceFreqDensForPlot(i_y+1,i_x+1)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'))],... % 
                      'EdgeColor','none')
            % Away from the diagonals, the patches should be rectangles
            else
                patch([xBinLims(i_x),xBinLims(i_x),xBinLims(i_x+1),xBinLims(i_x+1)],...
                      [yBinLims(i_y),yBinLims(i_y+1),yBinLims(i_y+1),yBinLims(i_y)],...
                      [RSquaredUsed(i_y,i_x)*occurrenceFreqDensForPlot(i_y,i_x)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y+1,i_x)*occurrenceFreqDensForPlot(i_y+1,i_x)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y+1,i_x+1)*occurrenceFreqDensForPlot(i_y+1,i_x+1)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'));...
                       RSquaredUsed(i_y,i_x+1)*occurrenceFreqDensForPlot(i_y,i_x+1)/...
                       (((xBinLims(i_x)<1)&(yBinLims(i_y)<1))*sum(occurrenceFreqDensForPlot(yBinLims<1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)<1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims<1),'all')+...
                       ((xBinLims(i_x)>=1)&(yBinLims(i_y)>=1))*sum(occurrenceFreqDensForPlot(yBinLims>=1,xBinLims>=1),'all'))],... % 
                      'EdgeColor','none')
            end
        end
    end
    % Add a color bar that ranges from 0 to 1
    colorbar
    clim([0 1])
    % set(gca, 'ColorScale', 'log')
    colormap jet
    % Add dashed lines at D1/Rc = 1, D2/Rc = 1, and D1/Rc = D2/Rc
    dashLineVec = [10^(-orderMagBelow1):10^(orderMagAbove1)];
    plot(dashLineVec,ones(size(dashLineVec)),'w--','LineWidth',2);
    plot(1:maxD1overRc,ones(size(1:maxD1overRc)),'k--','LineWidth',2)
    plot(ones(size(dashLineVec)),dashLineVec,'k--','LineWidth',2);
    plot(ones(size(1:maxD2overRc)),1:maxD2overRc,'w--','LineWidth',2)
    plot(10^(-orderMagBelow1):maxD2overRc,10^(-orderMagBelow1):maxD2overRc,'k--','LineWidth',2);
    hold off
    % Take care of axes, limits, labels, titles, etc
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xticks(10.^[-orderMagBelow1:orderMagAbove1])
    xticklabels(10.^[-orderMagBelow1:orderMagAbove1])
    yticks(10.^[-orderMagBelow1:orderMagAbove1])
    yticklabels(10.^[-orderMagBelow1:orderMagAbove1])
    pbaspect([1,1,1])
    xlim([10^(-orderMagBelow1),10^(orderMagAbove1)])
    ylim([10^(-orderMagBelow1),10^(orderMagAbove1)])
    xlabel('D_1/R_C')
    ylabel('D_2/R_C')
    % Print the subtitle with what the parameters used were
    if phi1and2atCenterOfRF == 1
        subtitle(['$\bar{v}_{x} = $', num2str(v_x_mean), ' mm/s, ',...
                  '$\phi_{RF, cent}$ = -(1/',num2str(-1/(phi_const/pi)),')*$\pi$, ',...
                  '$\omega_{\sigma}$ = ', num2str(omega_sigma*180/pi), ' $^{\circ}$/s, ',...
                  '$D_{1,max}$ = ',num2str(max_D1), ' mm, ',...
                  '$D_{2,max}$ = ',num2str(max_D2_from_D1 + max_D1), ' mm, '...
                  '$\phi_1 = \phi_2 = \phi_{RF, cent}$'],'Interpreter','latex')
    else
        subtitle(['$\bar{v}_{x} = $', num2str(v_x_mean), ' mm/s, ',...
                  '$\phi_{RF, cent}$ = -(1/',num2str(-1/(phi_const/pi)),')*$\pi$, ',...
                  '$\omega_{\sigma}$ = ', num2str(omega_sigma*180/pi), ' $^{\circ}$/s, ',...
                  '$D_{1,max}$ = ',num2str(max_D1), ' mm, ',...
                  '$D_{2,max}$ = ',num2str(max_D2_from_D1 + max_D1), ' mm, '...
                  '$\phi_1 \neq \phi_2 \neq \phi_{RF, cent}$'],'Interpreter','latex')
    end

    % Now compute the max info captured in each quadrant
    % If renormalizing OccurrenceFreqDens to have sum 1, then use that for
    % the info metric (this corresponds to weighting by frequency)
    if normOccurrenceFreqDensTo1 == 1
        fracInfoCapturedProgProg = sum(occurrenceFreqDensForPlot(yBinLims < 1, xBinLims < 1).*RSquaredUsed(yBinLims < 1, xBinLims < 1),'all');
        fracInfoCapturedProgReg = sum(occurrenceFreqDensForPlot(yBinLims >= 1, xBinLims < 1).*RSquaredUsed(yBinLims >= 1, xBinLims < 1),'all');
        fracInfoCapturedRegReg = sum(occurrenceFreqDensForPlot(yBinLims >= 1, xBinLims >= 1).*RSquaredUsed(yBinLims >= 1, xBinLims >= 1),'all');
    % If not renormalizing OccurrenceFreqDens to have sum 1, then use the
    % size of the bins in log space as the relative weights
    else
        tempSizeMat = sizeOfBinsInLogSpace(yBinLims < 1, xBinLims < 1);
        fracInfoCapturedProgProg = sum((tempSizeMat/sum(tempSizeMat,'all')).*RSquaredUsed(yBinLims < 1, xBinLims < 1),'all');
        tempSizeMat = sizeOfBinsInLogSpace(yBinLims >= 1, xBinLims < 1);
        fracInfoCapturedProgReg = sum((tempSizeMat/sum(tempSizeMat,'all')).*RSquaredUsed(yBinLims >= 1, xBinLims < 1),'all');
        tempSizeMat = sizeOfBinsInLogSpace(yBinLims >= 1, xBinLims >= 1);
        fracInfoCapturedRegReg = sum((tempSizeMat/sum(tempSizeMat,'all')).*RSquaredUsed(yBinLims >= 1, xBinLims >= 1),'all');
    end

    % Add text to label the different regions of prog/reg bar/BG in the plot
    % Prog / Prog text
    text(10^(-orderMagBelow1/2),10^(-1/8),'Prog Bar / Prog BG','HorizontalAlignment','center','FontSize',18,'Color','w');
    if fracInfoCapturedProgProg < 0.01
        text(10^(-orderMagBelow1/2),10^(-1/4),[num2str(fracInfoCapturedProgProg,'%.2e')],'HorizontalAlignment','center','FontSize',18,'Color','w');
    else
        text(10^(-orderMagBelow1/2),10^(-1/4),[num2str(fracInfoCapturedProgProg,3)],'HorizontalAlignment','center','FontSize',18,'Color','w');
    end
    % Prog / Reg text
    text(10^(-orderMagBelow1/2),10^(1/4),'Prog Bar / Reg BG','HorizontalAlignment','center','FontSize',18,'Color','w');
    if fracInfoCapturedProgReg < 0.01
        text(10^(-orderMagBelow1/2),10^(1/8),[num2str(fracInfoCapturedProgReg,'%.2e')],'HorizontalAlignment','center','FontSize',18,'Color','w');
    else
        text(10^(-orderMagBelow1/2),10^(1/8),[num2str(fracInfoCapturedProgReg,3)],'HorizontalAlignment','center','FontSize',18,'Color','w');
    end
    % Reg / Reg text
    text(10^(2/5),10^(orderMagAbove1)*10^(-1/8),'Reg Bar / Reg BG','HorizontalAlignment','center','FontSize',18,'Color','w');
    if fracInfoCapturedRegReg < 0.01
        text(10^(2/5),10^(orderMagAbove1)*10^(-1/4),[num2str(fracInfoCapturedRegReg,'%.2e')],'HorizontalAlignment','center','FontSize',18,'Color','w');
    else
        text(10^(2/5),10^(orderMagAbove1)*10^(-1/4),[num2str(fracInfoCapturedRegReg,3)],'HorizontalAlignment','center','FontSize',18,'Color','w');
    end
    % If weighing by frequency, note it on plot
    if normOccurrenceFreqDensTo1 == 1
        text(10^(2/5),10^(-1/8),'Condition','HorizontalAlignment','center','FontSize',18,'Color','k');
        text(10^(2/5),10^(-1/4),'Weighted by Freq','HorizontalAlignment','center','FontSize',18,'Color','k');
    % If weighing by log-log average, note it on plot
    else
        text(10^(2/5),10^(-1/8),'Condition','HorizontalAlignment','center','FontSize',18,'Color','k');
        text(10^(2/5),10^(-1/4),'Log-Log Averaged','HorizontalAlignment','center','FontSize',18,'Color','k');
    end

end
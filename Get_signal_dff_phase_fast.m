function Get_signal_dff_phase_fast(layer, fstim, phase, N, data_half, F)
%%
%This function plot the signal stack of a selected region in the brain
%layer choosen.
%It also plot the stimulus and the tracking.

%%
% =========== Parameters
figure;
nlayer = length(F.sets);            % number of layers in stack
nframes = size(F.set.frames, 2);    % number of frames in one layer
f = 1/F.dt * 1000;                  % image acquisition frequency (at which the camera runs)
Time = [(1/f)*layer:(1/f)*nlayer:(nframes*(1/f)*nlayer)];             % time of experiment sampled at image aquisition frequency
Time_stim = [0:(1/f)*nlayer:(nframes*(1/f)*nlayer)-((1/f)*nlayer)];             % time of experiment sampled at image aquisition frequency
a = [1:3:N;(1:3:N)+1];
Lin_first = size(a,2)+1;
Lin = Lin_first+N;
Col = 3;
%

pstim = size(Time, 2) * F.dt*0.001 * length(F.sets) * fstim;
if mod(nframes/pstim,1)==0
    pstim = pstim;
else
    pstim = pstim/2;
end

% =========== Load phasestack image
try
    im_PhaseMap = imread([ F.Files 'Phase_map_normalized/PhaseMap_DFF_pix_fstim', num2str(fstim),'/Images/Images0_' num2str(layer,'%02d') '.tif']);
    if N > 3
        ax(1) = subplot(Lin,Col,reshape(a,1,[]));
    else
        ax(1) = subplot(Lin,Col,[1:2,4:5]);
    end
    imshow(rescalegd2(im_PhaseMap));
    axis off
    hold(ax(1), 'on');
catch
    warning('Phase stack doesn t exit (gray stack default)');
    im_GreyStack = imread([ F.Files 'grey_stack/Image_' num2str(layer,'%02d') '.tif']);
    if N > 3
        ax(1) = subplot(Lin,Col,reshape(a,1,[]));
    else
        ax(1) = subplot(Lin,Col,[1:2,4:5]);
    end
    imshow(rescalegd2(im_GreyStack));
    axis off
    hold(ax(1), 'on');
end

% =========== Load greystack image
im_GreyStack = imread([ F.Files 'grey_stack/Image_' num2str(layer,'%02d') '.tif']);
if N > 3
    ax(2) = subplot(Lin,Col,(1:3:N)+2);
else
    ax(2) = subplot(Lin,Col,[3,6]);
end
h_im = imshow(rescalegd2(im_GreyStack));
axis off
hold(ax(2), 'on');

% =========== Load Tracking's trace
try
    [Tracking_sig, Tracking_time] = tracking_signal(F);
    ax(3) = subplot(Lin,Col,Lin_first*Col-2:Lin_first*Col-1);
    hold(ax(3), 'on');
    plot(Tracking_time, -Tracking_sig, 'r');
    [Tracking_int] = tracking_interpolate(Time_stim, F);
    Tracking_int(Tracking_int==0) = nan;
    Tracking_sig_mean = nanmean(reshape(Tracking_int,[size(Tracking_int,2)/pstim, pstim]),2);
    ax(6) = subplot(Lin,Col,Lin_first*Col);
    hold(ax(6), 'on');
    plot(Time_stim(1:size(Tracking_sig_mean,1)), -(Tracking_sig_mean-min(Tracking_sig_mean)) , 'r');
catch
    warning('No Tracking trace');
end

% =========== Load Motor's trace
try
    [Motor_sig, Motor_time] = motor_signal(F);
    ax(3) = subplot(Lin,Col,Lin_first*Col-2:Lin_first*Col-1);
    [Motor_int] = motor_interpolate(Time_stim, F);
    if data_half==1
        plot(Time_stim(1:size(Time_stim,2)/2), Motor_int(1:size(Time_stim,2)/2), 'k');
    elseif data_half==2
        plot(Time_stim((size(Time_stim,2)/2)+1:end), Motor_int((size(Time_stim,2)/2)+1:end), 'k');
    else
        plot(Time_stim, Motor_int, 'k');
    end
    ylabel(ax(3),'Amplitude [°]')
    Motor_sig_mean = reshape(Motor_int,[size(Motor_int,2)/pstim, pstim]);
    Motor_sig_mean = Motor_sig_mean(:,2);
    ax(6) = subplot(Lin,Col,Lin_first*Col);
    plot(Time_stim(1:size(Motor_sig_mean,1)), Motor_sig_mean, 'k');
catch
    warning('No Motor trace');
end

% =========== Load data
S = matfile([F.Files 'signal_stacks' filesep num2str(layer) filesep 'DFF_bg.mat']); % DFF per pixel
index = S.index;

% % =========== Mean Data selected
% ROI = imrect(ax(1));                 % get user specified ROI
% mask_ROI = createMask(ROI,h_im);         % make a black and white mask of the ROI
% delete(ROI);                             % delete the ROI overlay from the image
% ind_ROI = find(mask_ROI);                % get indices of pixels inside the ROI
% [~,ind_ROI] = intersect(index,ind_ROI);  % keep only indices of ROI-pixels that are inside the brain contour
%
% tmp_dff = zeros(size(ind_ROI(:),1), nframes);
% for i = 1:size(ind_ROI(:),1)
%     tmp_dff(i,:) = S.DFF_pix(ind_ROI(i,:),:);
% end
% mean_Data = nanmean(tmp_dff);
% clear tmp_dff
%
% ax(4) = subplot(3,3,4:5);
% hold(ax(4), 'on');
% plot(Time, mean_Data, '.') % Plot DFF of ROI (first point is shifted to zero)

%N = 10;             % number of colors used in plot
colours = jet(N);   % set color map used in plot
count = 1;

while count <= N
    % =========== Select ROI and get pixel indices
    ROI = imfreehand(ax(1));                 % get user specified ROI
    mask_ROI = createMask(ROI,h_im);         % make a black and white mask of the ROI
    delete(ROI);                             % delete the ROI overlay from the image
    ind_ROI = find(mask_ROI);                % get indices of pixels inside the ROI
    [~,ind_ROI] = intersect(index,ind_ROI);  % keep only indices of ROI-pixels that are inside the brain contour
    [x,y] = find(mask_ROI);                  % get x,y coordinates of all pixels in the ROI
    plot(ax(1),y,x,'color',colours(count,:)) % plot all pixels of ROI as overlay in ax(1); color is defined by color map and the color counter
    plot(ax(2),y,x,'color',colours(count,:)) % plot all pixels of ROI as overlay in ax(1); color is defined by color map and the color counter
    
    tmp_dff = zeros(size(ind_ROI(:),1), nframes);
    for i = 1:size(ind_ROI(:),1)
        tmp_dff(i,:) = S.DFF_pix(ind_ROI(i,:),:);
    end
    Data_ROI = nanmean(tmp_dff);
    clear tmp_dff
    
    if data_half==1
        Data_ROI = Data_ROI(1:(size(Data_ROI,2)/2));
        Data_ROI_mean = mean(reshape(Data_ROI,[size(Data_ROI,2)/(pstim/2), (pstim/2)]),2);
    elseif data_half==2
        Data_ROI = Data_ROI((size(Data_ROI,2)/2)+1:end);
        Data_ROI_mean = mean(reshape(Data_ROI,[size(Data_ROI,2)/(pstim/2), (pstim/2)]),2);
    else
        Data_ROI_mean = mean(reshape(Data_ROI,[size(Data_ROI,2)/pstim, pstim]),2);
    end
    ax(4) = subplot(Lin,Col,(Col*Lin_first+((count-1)*3)+1):(Col*Lin_first+((count-1)*3)+2));
    plot(Time(1:size(Data_ROI,2)), Data_ROI,'color',colours(count,:),'LineWidth',1) % Plot DFF of ROI (first point is shifted to zero)
    ax(5) = subplot(Lin,Col,(Col*Lin_first+((count-1)*3)+3));
    ylabel(ax(4),'DFF')
    yyaxis left
    plot(Time(1:size(Data_ROI_mean,1)),Data_ROI_mean, 'color',colours(count,:),'LineWidth',1)
    set(ax(5),'ycolor',colours(count,:))
    ylabel(ax(5),'DFF')
    set(ax(5),'ycolor','k')
    yyaxis right
    plot(Time_stim(1:size(Motor_sig_mean,1)), Motor_sig_mean, 'k');
    set(ax(5),'ycolor','k')

    if count == N
        xlabel(ax(4),'Time [s]')
        xlabel(ax(5),'Time [s]')
        %grid(ax(3:6), 'on')
    else
        set(ax(3),'XTick',[])
        set(ax(6),'XTick',[])
        set(ax(4),'XTick',[])
        set(ax(5),'XTick',[])
    end
    
    % Save
    save_ROI{count, 1} = Data_ROI;
    save_ROI{count, 2} = mask_ROI;
    save_ROI{count, 3} = Data_ROI_mean;
    count = count+1;
end

%% =========== Save
mkdir([F.Files 'analysis' filesep phase])
save([F.Files 'analysis' filesep phase filesep 'ROI_layer' num2str(layer) '.mat'],'save_ROI', 'Time');



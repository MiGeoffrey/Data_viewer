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
%Time = [0:(1/f)*nlayer:(nframes*(1/f)*nlayer)-((1/f)*nlayer)];             % time of experiment sampled at image aquisition frequency

pstim = size(Time, 2) * F.dt*0.001 * length(F.sets) * fstim;
if mod(nframes/pstim,1)==0
    pstim = pstim;
else
    pstim = pstim/2;
end

% =========== Load phasestack image
im_PhaseMap = imread([ F.Files 'Phase_map_normalized/PhaseMap_DFF_pix_fstim', num2str(fstim),'/Images/Images0_' num2str(layer,'%02d') '.tif']);
ax(1) = subplot(3,3,1);
h_im = imshow(rescalegd2(im_PhaseMap));
hold(ax(1), 'on');

% =========== Load greystack image
im_GreyStack = imread([ F.Files 'grey_stack/Image_' num2str(layer,'%02d') '.tif']);
ax(2) = subplot(3,3,2);
h_im = imshow(rescalegd2(im_GreyStack));
hold(ax(2), 'on');

% =========== Load Tracking's trace
try
    [Tracking_sig, Tracking_time] = tracking_signal(F);
    ax(3) = subplot(3,3,7:8);
    hold(ax(3), 'on');
    plot(Tracking_time, -Tracking_sig, 'r');
    [Tracking_int] = tracking_interpolate(Time, F);
    Tracking_int(Tracking_int==0) = nan;
    Tracking_sig_mean = nanmean(reshape(Tracking_int,[size(Tracking_int,2)/pstim, pstim]),2);
    ax(6) = subplot(3,3,9);
    hold(ax(6), 'on');
    plot(Time(1:size(Tracking_sig_mean,1)), -(Tracking_sig_mean-min(Tracking_sig_mean)) , 'r');
catch
    warning('No Tracking trace');
end

% =========== Load Motor's trace
try
    [Motor_sig, Motor_time] = motor_signal(F);
    ax(3) = subplot(3,3,7:8);
    plot(Motor_time, Motor_sig, 'b');
    axis([0 600 -20 20]);
    [Motor_int] = motor_interpolate(Time, F);
    Motor_sig_mean = reshape(Motor_int,[size(Motor_int,2)/pstim, pstim]);
    Motor_sig_mean = Motor_sig_mean(:,2);
    ax(6) = subplot(3,3,9);
    plot(Time(1:size(Motor_sig_mean,1)), Motor_sig_mean, 'b');
    ax(5) = subplot(3,3,6);
    hold(ax(5), 'on');
    yyaxis right
    plot(Time(1:size(Motor_sig_mean,1)), Motor_sig_mean, '-r');
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
    ax(4) = subplot(3,3,4:5);
    plot(Time(1:size(Data_ROI,2)), Data_ROI,'color',colours(count,:),'LineWidth',1) % Plot DFF of ROI (first point is shifted to zero)
    ax(5) = subplot(3,3,6);
    hold(ax(5), 'on');
    yyaxis left
    plot(Time(1:size(Data_ROI_mean,1)),Data_ROI_mean, 'color',colours(count,:),'LineWidth',1)
    
    % Save
    save_ROI{count, 1} = Data_ROI;
    save_ROI{count, 2} = mask_ROI;
    save_ROI{count, 3} = Data_ROI_mean;
    count = count+1;
end

%% =========== Save
mkdir([F.Files 'analysis' filesep phase])
save([F.Files 'analysis' filesep phase filesep 'ROI_layer' num2str(layer) '.mat'],'save_ROI', 'Time');



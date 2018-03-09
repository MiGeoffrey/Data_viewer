function [Tracking_int] = tracking_interpolate(Time, F)

%% =========== Load tracking signal
file = dir([F.Data,'*Tracking.dat']);
fileID = fopen([F.Data file.name]);
Data = fread(fileID,[ 2, Inf ],'single=>single','ieee-le');
fclose(fileID);

Tracking = Data(2,2:2:end);
TimeTracking = Data(1,2:2:end);

TimeOffset = TimeTracking(1);
TimeTracking = TimeTracking - TimeOffset;

%% =========== Out put
Tracking_int = interp1(TimeTracking(1:end),Tracking(1:end),Time); % interpolate tracking signal at time points where images were acquired

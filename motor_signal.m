function [Motor_sig, Motor_time] = motor_signal(F)

%% =========== Load motor signal
filename = [F.Data,'Stimulus.txt']; % data path and name of where the motor movement is saved
delimiter = '\t';                    
formatSpec = '%f%f%f%[^\n\r]';      % formate of the text file
fileID = fopen(filename,'r');       % open text file
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false); % scan text file
fclose(fileID);                     % close text file

TimeMotor = dataArray{:, 2};   
TimeOffset = TimeMotor(1);
TimeMotor = TimeMotor - TimeOffset;
Motor = dataArray{:, 3};

%% =========== Out put
Motor_sig = Motor;
Motor_time = TimeMotor;

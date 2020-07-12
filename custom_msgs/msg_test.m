%examplePackages = fullfile(fileparts('../../ibeacon_ros/'));
userFolder = 'msg/msg';
%copyfile(examplePackages, userFolder);

rosgenmsg(userFolder)
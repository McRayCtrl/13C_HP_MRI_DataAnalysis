%% HP Pyruvate/Lactate Movie Maker

%Go to folder that has 90 HP pyruvate/lactate map png files and then run
%this code

clear;close all
% Specify the folder containing PNG images
image_folder = pwd;

for i = 1:90
    image_name = num2str(i);
    file_pattern = fullfile(image_folder, [image_name,'.png']);
    png_files(i) = dir(file_pattern);
end

% Specify the output video file name
output_video = 'output_video.mp4';

% Create a VideoWriter object
writerObj = VideoWriter(output_video, 'MPEG-4');

% Set frame rate (adjust as needed)
writerObj.FrameRate = 18; % frames per second

% Open the VideoWriter object
open(writerObj);

% Loop through each PNG file and write it to the video
for i = 1:numel(png_files)
    % Read the PNG file
    filename = fullfile(image_folder, png_files(i).name);
    frame = imread(filename);
    
    % Write the frame to the video
    writeVideo(writerObj, frame);
end

% Close the VideoWriter object
close(writerObj);

% Display message
disp('Video created successfully.');

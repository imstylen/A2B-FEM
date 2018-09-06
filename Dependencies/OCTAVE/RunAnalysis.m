clear; clc; close all;


disp("____________________________________")
disp("A2B Finite Element Modeling Software")
disp("____________________________________")
disp("")
disp("Gathering input files")
disp("")


inputFileList = dir('../../0. INPUT/*.node');
nModels = length(inputFileList);

fprintf('Performing (%i) analyses \n',nModels)
disp("____________________________________")
disp("")
disp("")


for i = 1:nModels
  
    fprintf("Analysis %i of %i \n",i,nModels);
    modelName = inputFileList(i).name;
    modelName = modelName(1:length(modelName)-5);
  
    Analyze(modelName)
    disp('Complete!');
    disp('');
  
  
endfor

fprintf("Completed (%i) anlayses \n",nModels);

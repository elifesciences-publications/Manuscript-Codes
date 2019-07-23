%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Find membrane voltage over time %%%%%%
%%%%%%%%%%% Created: 09-20-2016 %%%%%%%%%%%%%
%%%%%%%%%%%% Edited: 10-18-2016 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAR PREV WORK
close all
clearvars


%% INIT VARS
%sample rate
sample_rate=10000;

%% LOAD FILES
Vm_dir=uigetdir;

%get three minutes of baseline
cd([Vm_dir '/baseline'])
contents = dir('*.txt');
filenames.bl = {contents.name}';
files.bl = fullfile(cd,filenames.bl);
for ii = 1:length(files.bl)
    fid=fopen(files.bl{ii});
    tempbaseData{ii}=textscan(fid,'%f'); %#ok<*SAGROW>
end
baseData=[tempbaseData{1}; tempbaseData{2}; tempbaseData{3}];
baseData=[baseData{1}; baseData{2}; baseData{3}];

% get ten minutes of drug
cd([Vm_dir '/cno']);
contents = dir('*.txt');
filenames.drug = {contents.name}';
files.drug = fullfile(cd,filenames.drug);
for ii = 1:length(files.drug)
    fid=fopen(files.drug{ii});
    tempdrugData{ii} = textscan(fid,'%f');
end
drugData=[tempdrugData{1}; tempdrugData{2}; tempdrugData{3}; tempdrugData{4}; tempdrugData{5}; ...
    tempdrugData{6}; tempdrugData{7}; tempdrugData{8}; tempdrugData{9}; tempdrugData{10}];
drugData=[drugData{1}; drugData{2}; drugData{3}; drugData{4}; drugData{5}; ...
    drugData{6}; drugData{7}; drugData{8}; drugData{9}; drugData{10}];

%% DATA ANALYSIS

%get mean baseline Vm
mrawbaseVm = mean(baseData);
srawbaseVm = baseData-mrawbaseVm; %gives 3min baseline subtracted baseline Vm

%get baseline subtracted drug Vm
srawdVm = drugData-mrawbaseVm; %gives baseline subtracted drug Vm
mean_rawdVm=mean(srawdVm(end-3*sample_rate*60+1:end)); %mean change Vm cause by drug defined as last 3min of drug

%get mean values every twenty seconds
t=(1/2):(1/2):13;
for ii = 1:6 %6 sets of 30s within 3 minute baseline
    baseVm(ii)=mean(srawbaseVm((ii-1)*(length(srawbaseVm)/6)+1:ii*(length(srawbaseVm)/6)));
end
for ii = 1:20 %20 sets of 30s within 10 minute drug period
    dVm(ii)=mean(srawdVm((ii-1)*(length(srawdVm)/20)+1:ii*(length(srawdVm)/20)));
end

%% DISPLAY DATA
traceFig=figure(1);
traceFig.Position = [725 450 725 350];
subplot(2,1,1)
hold on
scatter(t,[baseVm dVm],'filled')
line([t(10) t(end)],[5 5],'linewidth',3,'color','k')
line([0 t(end)],[0 0],'linewidth',1.5,'color','k','linestyle','--')
xlim([0 t(end)])
ylim([1.2*min(srawdVm) 6])
subplot(2,1,2)
hold on
plot([srawbaseVm; srawdVm])
line([1+length(srawbaseVm) length([srawbaseVm; srawdVm])],[5 5],'linewidth',3,'color','k')
line([1 length([srawbaseVm; srawdVm])],[0 0],'linewidth',1.5,'color','k','linestyle','--')
xlim([0 length([srawbaseVm; srawdVm])])
ylim([1.2*min(srawdVm) 6])

%% SAVE DATA
button = questdlg('Save data');
if strcmp(button,'Yes')
    [save_file,path_save]=uiputfile;
    full_save_file=fullfile(path_save,save_file);
    save(full_save_file);
end

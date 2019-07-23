%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Ra analysis over mult files %%%%%
%%%%%%%%% Created: 08-15-2016 %%%%%%%%%%
%%%%%%%%%% Edited: 08-15-2016 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars

%init vars
sample_rate=10000; %sample rate
injstart = 5625; %voltage injection start
injend = 6624; %voltage injection end
%get number of cells
file_no=inputdlg('How many files?');
file_no=str2double(file_no{1});
all_dataRa=[];

%load data files
for aa = 1:file_no
    dirRa=uigetdir;
    cd(dirRa);
    contents = dir('*.abf');
    
    %load
    filenames = {contents.name}';
    fileRa = fullfile(cd,filenames);
    dataRa = abfload(fileRa{1},'sweeps','a');
    
    if isempty(all_dataRa)
        all_dataRa = dataRa;
    else
        all_dataRa = cat(3,all_dataRa,dataRa);
    end
end

%check Ra
[Ra, meanRa, keep] = getRa(all_dataRa,injstart,injend,sample_rate) %#ok<NOPTS>

%plot Ra to view
Rafig=figure;
plot(Ra,'-o')
line([0 length(Ra)+1],[meanRa * 1.2, meanRa * 1.2],'LineStyle','--','Color','k')
line([0 length(Ra)+1],[meanRa * .8, meanRa * .8],'LineStyle','--','Color','k')
ylim([min(Ra)-5 max(Ra)+5])
xlim([0 length(Ra)+1])

%save data
button = questdlg('Save data');
if strcmp(button,'Yes')
    [save_file,path_save]=uiputfile;
    full_save_file=fullfile(path_save,save_file);
    save(full_save_file);
end
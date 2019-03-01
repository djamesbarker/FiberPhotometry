%% Code by David J. Barker for Morales Laboratory
% Beta Release 7-6-2018; Contact djamesbarker@gmail.com
% Code should be cited as "Barker, D. J. et al., (2017). Lateral preoptic 
% control of the lateral habenula through convergent glutamate and GABA 
% transmission. Cell reports, 21(7), 1757-1769."

clear all; close all; %Clear all variables and close all figures
GroupFold=uigetdir; %Opens prompt to select the folder with the TDT Tanks for a whole group of subjects
cd(GroupFold); %Change directory to group folder
Subjects=dir(GroupFold); % Create a list of files in the folder
Master=struct; % Establish summary variable


for s=3:length(Subjects) %Main loop (file 3 is first real file, first two directories are '.' and '..').
cd(GroupFold) 
cd(Subjects(s).name); % Move to first subject folder
SubjectID=strtok(Subjects(s).name,'-'); %Use filename to extract subject ID number for data storage
BlockDir=dir('*Photometry*'); %Find TDT block with photometry data
BlockDir=BlockDir.name; %Get folder name
cd(BlockDir); %Change directory to photometry folder
[Tank,Block,~]=fileparts(cd); %List full directory for Tank and Block
data=TDT2mat2(Tank,Block); %Use TDT2Mat to extract data.

%% SYNAPSE- Extract Relevant Data from Data file
%Create Variables for each Photometry Channel and timestamps
Ch490=data.streams.x90R.data; %GCaMP
Ch405=data.streams.x05R.data; %Isosbestic Control
Ts = ((1:numel(data.streams.x90R.data(1,:))) / data.streams.x90R.fs)'; % Get Ts for samples based on Fs
StartTime=5000; %Set the starting sample(recommend eliminating a few seconds for photoreceiver/LED rise time).
EndTime=length(Ch490)-1000; %Set the ending sample (again, eliminate some).
Fs=data.streams.x90R.fs; %Variable for Fs
Ts=Ts(StartTime:EndTime); % eliminate timestamps before starting sample and after ending.
Ch490=Ch490(StartTime:EndTime);
Ch405=Ch405(StartTime:EndTime);

%% Behavioral Events
% Events are in data.epocs variable, usually in PAB_ or a similar variable
% data.epocs.PAB_.data contains the bit the data came in on
% data.epocs.PAB_.onset is the time the bit turned ON and offset is the
% time it turned OFF.
% The rows of data correspond directly to the onset rows.

%In the example:
% TONE=16
% SHOCK+TONE=2048
% SHOCK+TONE OFF=0

Beh=[];Tone=[];Shock=[];End=[];
Beh=data.epocs.PAB_; %Make a variable with all behavioral data.

%find the rows of data with each specific bit and then extract the
%timestamp for those specific bits. Copy and change these lines for any
%bits that you need data from.
Tone=Beh.onset(Beh.data(:,1)==16,:); 
Shock=Beh.onset(Beh.data(:,1)==2048,:);
End=Beh.onset(Beh.data(:,1)==0,:);

%% Function to get DeltaF/F using regression fit
%Function to get DF/F for whole session
Delta490=DeltaF(Ch490,Ch405); 
%Variant of whole session that adds ticks for events. Two events possible.
[Delta490] = DeltaFevent(Ch490,Ch405,Ts,[Tone Tone+5],'Tone-Shock',[Tone+5 Tone+10]);
%Function to downsample the whole session trace for saving
[DSDeltaF, DTs, DFs] =DownSample(Delta490,Fs);

%Save CSV and Matlab file with downsampled whole session trace.
save('Delta490','DSDeltaF');
csvwrite('Delta490.csv',DSDeltaF);
close all


%% Photometry RASTER-PETH (Tone and Shock);
Pre=5; %Time to sample for Raster-PETH before event
Post=10; %Time to sample for raster-PETH after event.
BL_Width=5; %Set the duration of the baseline
BL_Start=5; %Time before raster to start taking baseline

%Code to create Raster-PETH. Delta F/F is calculated for each trial.
%ToneTrace is average trace across all trials.
%TrialTrace is individual trial data.
[ToneTrace, TrialTrace]=TrialPETH( Fs, Ts, Pre, Post, Tone,Ch490,Ch405, 'TONE PETH',100,0,5, BL_Width,BL_Start );

%% Summary Data
%Set up and/or clear variables.
AUC=[];MaxValue=[];
Master.(['S_' num2str(SubjectID)]).Tone=ToneTrace;

%Example AUC. Take AUC Between specific timestamps from PETH data.
AUC(1,1)=trapz(ToneTrace(ToneTrace(:,1)<=-2.5,2));
AUC(1,2)=trapz(ToneTrace(ToneTrace(:,1)>=0 & ToneTrace(:,1)<=2.5,2));
AUC(1,3)=trapz(ToneTrace(ToneTrace(:,1)>=5 & ToneTrace(:,1)<=7.5,2));

% Example Peak Data- Take peak between specific timestamps from PETH data.
MaxValue(1,1)=max(ToneTrace(ToneTrace(:,1)<=-2.5,2));
MaxValue(1,2)=max(ToneTrace(ToneTrace(:,1)>=0 & ToneTrace(:,1)<=2.5,2));
MaxValue(1,3)=max(ToneTrace(ToneTrace(:,1)>=5 & ToneTrace(:,1)<=7.5,2));

%Wilcoxon stats example. Compare Tone or Shock to Baseline to determine if
% 'responsive'. NOT USED FOR FINAL STATS
WilcoxTone=signrank(ToneTrace(ToneTrace(:,1)<=-2.5,2),ToneTrace(ToneTrace(:,1)>=0 & ToneTrace(:,1)<=2.5,2));
WilcoxShock=signrank(ToneTrace(ToneTrace(:,1)<=-2.5,2),ToneTrace(ToneTrace(:,1)>=5 & ToneTrace(:,1)<=7.5,2));
    
%Example Master Data format used to save individual subject data into matlab structure.
Master.(['S_' num2str(SubjectID)]).AUC=AUC;
Master.(['S_' num2str(SubjectID)]).Peak=MaxValue;
Master.(['S_' num2str(SubjectID)]).TrialTrace=TrialTrace;
Master.(['S_' num2str(SubjectID)]).WilcoxTone=WilcoxTone;
Master.(['S_' num2str(SubjectID)]).WilcoxShock=WilcoxShock;

% Matrix for all group data; Written into a matlab structure after all
% subjects data is collected.
GroupZTrace(:,1)=ToneTrace(:,1);
GroupZTrace(:,s-1)=ToneTrace(:,2);

GroupDFTrace(:,1)=ToneTrace(:,1);
GroupDFTrace(:,s-1)=ToneTrace(:,3);

GroupAUC(s-2,:)=AUC;
GroupPeakMax(s-2,:)=MaxValue;
end 

cd(GroupFold)
Master.GroupZTrace=GroupZTrace;
Master.GroupDFTrace=GroupDFTrace;
Master.GroupAUC=GroupAUC;
Master.GroupPeakMax=GroupPeakMax;


% Very quick example for plotting summary data.
save('Master','Master')
figure
%UN-COMMENT SHADED ERROR BAR WHEN USING MORE THAN ONE FILE! 
%shadedErrorBar(GroupZTrace(:,1),mean(GroupZTrace(:,2:end),2),(std(GroupZTrace(:,2:end),0,2)/sqrt(size(GroupZTrace(:,2:end),2))));
%ONLY MEAN PLOTTED FOR SINGLE SAMPLE-DATA FOLDER.
plot(GroupZTrace(:,1),mean(GroupZTrace(:,2:end),2));
xlabel('Time');ylabel('Z-Score')
title('Average Data')

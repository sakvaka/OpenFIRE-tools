% plot_stack.m
% Created by S. Vakeva (2017) for the OpenFIRE project.
% Ported from the legacy F77 code su_plot_fire.f by Mr J. Keskinen (Uni Helsinki).
%
% Plots a seismic nmo or dmo stack. The input file must be provided in .su format.
% The input file must also be migrated and depth-converted.
%
% Tested on Matlab r2014b.

clear ; clc; close

% requests an .su file
filename=uigetfile('*.su');

fp=fopen(filename,'r');

% requests if the section is to be plotted from left to right or vice versa
bflip=questdlg('Flip horizontal axis?','Flip','Yes','No','No');

% single precision is sufficient
traces=single([]);
envtraces=[];
envtraces_stacked=[];
envtraces_averaged=[];

% half-width of the window used for moving average calculation 
% within traces (in samples) after Hilbert transform
iwin=10;
% 1/3 width of the window used for lateral moving average calculation (in samples)
istack=5;
% plot only every nth sample
ispar=10;
% coefficient to multiply with before envelope calculation
rother=0.65;
% coefficient to multiply with before logarithm calculation
rcontur=0.16;

depthsamprate=12.5; % IN METERS
cdpsamprate=25; % IN METERS

% don't touch below this line
depthsamprate=depthsamprate/1000;
cdpsamprate=cdpsamprate/1000;

curtrace=0;
rmaxtr=0.1/rcontur*istack*3.0;

while (fseek(fp,20,'cof')==0)
    curtrace=curtrace+1;
    if (mod(curtrace,50)==0)
        fprintf('Reading trace no. %d\n',curtrace);
    end
    cdp(curtrace)=fread(fp,1,'int','ieee-be');
    status=fseek(fp,90,'cof');
    ns=fread(fp,1,'int16','ieee-be');
    dt=fread(fp,1,'int16','ieee-be');
    status=fseek(fp,122,'cof');
    traces(1:ns,curtrace)=fread(fp,ns,'float','ieee-be');
    % envelope derived from the analytical trace
    traces(1:ns,curtrace)=rother*traces(1:ns,curtrace);
    envtraces(1:ns,curtrace)=abs(hilbert(traces(:,curtrace)));
    % moving average within trace
    envtraces(:,curtrace)=conv(envtraces(:,curtrace), ones(2*iwin+1,1)/ ...
        (2*iwin+1), 'same');
end
fclose(fp);

fprintf('Stacking...');
% stacks every istack traces into one new trace
no_of_istacks=floor(curtrace/istack);
envtraces_stacked=zeros(ns,no_of_istacks);
for i=1:ns
    b=filter(ones(istack,1),1,envtraces(i,:));
    envtraces_stacked(i,:)=b(istack:istack:end);
end

% lateral moving average, window width 3 samples
envtraces_averaged=zeros(ns,no_of_istacks);
for i=1:ns
    b=filter(ones(3,1),1,envtraces_stacked(i,:));
    envtraces_averaged(i,:)=b;
end

% preserves every nth depth sample
envtraces_averaged=envtraces_averaged(1:ispar:end,:)/rmaxtr;

% decibel (dB) scale used for plotting
% note! original code erroneously uses natural logarithm
fprintf('Calculating 20*log...\n');
nanneja=isnan(envtraces_averaged);
negatiivisia=envtraces_averaged<0;
envtraces_averaged=20*log(envtraces_averaged);
% if sample is NaN or negative, plot as white color
envtraces_averaged(nanneja | negatiivisia)=-999.999;

% make a depth scale in kilometers
depth=linspace(0,(ns-1)*depthsamprate,size(envtraces_averaged,1));
% make a horizontal scale in kilometers
cdp=[1:size(envtraces_averaged,2)]*istack*cdpsamprate;

% finally, plot and whiten the data by 0.8 dB
fprintf('Plotting...\n');
pcolor(cdp,depth,envtraces_averaged-0.8);
set(gca,'Ydir','reverse');

% set reverse horizontal axis, if section should be flipped
switch bflip
    case 'Yes'
        set(gca,'Xdir','reverse');
end

% empirical colormap
colormap gray
colormap(flipud(colormap))
caxis([0 20]);
shading interp
axis equal tight

% axis labeling
xlabel('Distance [km]');
ylabel('Depth [km]');

% fix underscores in figure title
title_filename=strrep(filename,'_','\_');
title(title_filename);

%% Run this section if you wish to export a TIFF figure.
fprintf('Exporting figure...\n');
output_filename=strrep(filename,'.su','.tiff');
print(output_filename,'-dtiff','-r400');

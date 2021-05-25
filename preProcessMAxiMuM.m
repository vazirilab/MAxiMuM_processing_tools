%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preProcessMAxiMuM.m
%
% J.D. 05/19/2020
%
% For each file, import raw .tif files, assemble ROIs, 
% correct scan phase, re-order planes, and save a memory-mapped file.
% Access each plane in each file and concatenate in time
% Motion correct each plane and save to separate .tif files.
%
% Use the 'filePath' input argument to point to a local folder with raw
% .tif files
%
% Use the 'fileNameRoot' input argument to point toward the root for each
% file in the directory. The code will appeand '_00001.tif' to the end, so
% remove this suffix from the root
%
% The 'diagnosticFlag' argument, when set to '1', will report all .tif files in
% the directory specified by 'path'
%
% Each motion-corrected plane is saved as a separate .mat file with the following fields:
% Y: single plane recording data (x,y,T) (single)
% Ym: mean projection image of Y (x,y) (single)
% sizY: array with size of dimension of Y (1,3) 
% volumeRate: volume rate of the recording (1,1) (Hz)
% pixelResolution: size of each pixel in microns (1,1) (um)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function preProcessMAxiMuM(filePath,fileNameRoot,diagnosticFlag,numCores)

if ~exist([filePath fileNameRoot '_00001.tif'],'file')
   error('File does not exist.')
end

if strcmp(diagnosticFlag,'1') % if the diagnostic flag is set to 1, spit out contents of directory specified by 'path'
    dir([filePath,'*.tif'])
else

    tic
    
    clck = clock; % use current time and date to make a log file
    fid = fopen(fullfile(filePath,['matlab_log_' num2str(clck(1)) '_' num2str(clck(2)) '_' num2str(clck(3)) '_' num2str(clck(4)) '_' num2str(clck(5)) '.txt']),'w');
    
    poolobj = gcp('nocreate'); % if a parallel pool is running, kill it and restart it to make sure parameters are correct
    if ~isempty(poolobj)
        disp('Removing existing parallel pool.')
        delete(poolobj)
    end
    
    if str2double(numCores) == 0 || size(str2double(numCores),1) == 0
        numCores = 12;
    else
        numCores = str2double(numCores);
    end
    
    mkdir(fullfile(filePath,'TMP'));
    tmpPath = fullfile(filePath,'TMP');

    %% Set up path dependencies
    addpath(genpath(fullfile('motion_correction')));
    addpath(genpath(fullfile('CaImAn-MATLAB-master','CaImAn-MATLAB-master')));
    addpath(genpath(fullfile('ScanImage')));

    %% Determine how many files there are

    disp('Beginning processing routine.')
    
    date = datetime(now,'ConvertFrom','datenum');
    formatSpec = '%s Beginning processing routine...\n';
    fprintf(fid,formatSpec,date);

    if exist([filePath fileNameRoot '_00001.tif'],'file') > 0
        xx = dir(fullfile(filePath,'*.tif'));
        N = {xx.name};
        X = ~cellfun('isempty',strfind(N,fileNameRoot));
        numFiles = sum(X);
        multiFile = true;
    else
        multiFile = false;
        numFiles = 1;
    end

    tic;

    %% Loop through files

    for ijk = 1:numFiles
        disp(['Loading file ' num2str(ijk) ' of ' num2str(numFiles) '...'])
        date = datetime(now,'ConvertFrom','datenum');
        formatSpec = '%s Beginning plane %u...\n';
        fprintf(fid,formatSpec,date,ijk);

        % Determine current file name
        if multiFile % if multiple files, append '_0000x.tif' or '_000xx.tif'
            if ijk < 10
                fileName = [fileNameRoot '_0000' num2str(ijk) '.tif'];
            else
                fileName = [fileNameRoot '_000'  num2str(ijk) '.tif'];
            end
        else % single file recording, append '.tif' to end
            fileName = [fileNameRoot '.tif'];
        end

        [vol,volumeRate,pixelResolution] = assembleCorrectedROITiff([filePath fileName]);
        tt = toc/3600;
        disp(['Volume loaded and processed. Elapsed time: ' num2str(tt) ' hours. Saving volume to temp...'])

        fullVolumeSize = size(vol);
        numberOfPlanes = fullVolumeSize(3);
        Tfile = fullVolumeSize(4);
        d1file = fullVolumeSize(1);
        d2file = fullVolumeSize(2);

        if numberOfPlanes == 30
            order = [1 5:10 2 11:17 3 18:23 4 24:30];
            order = fliplr(order);
        elseif numberOfPlanes == 15
            order = [1 5 6 2 7:9 3 10:14 4 15];
            order = fliplr(order);
        else
            disp('Number of planes not recognized.')
        end

        vol = vol(:,:,order,:);

        savefast([tmpPath fileName(1:end-3) 'mat'],'vol','volumeRate','pixelResolution','fullVolumeSize');
        tt = toc/3600;
        disp(['Volume loaded, processed, and saved to disk. Elapsed time: ' num2str(tt) ' hours. Processing next volume...'])

        clear vol
        pause(0.5)

    end

    %% (II) Process each plane

    for aaa = 1:numberOfPlanes

        disp(['PROCESSING PLANE ' num2str(aaa) ' OF ' num2str(numberOfPlanes) '...'])

        Y = zeros(d1file,d2file,Tfile*numFiles,'single');
        for ddd = 1:numFiles
            if multiFile
                if ddd < 10
                    matfilename = [tmpPath fileNameRoot '_0000' num2str(ddd) '.mat'];
                else
                    matfilename = [tmpPath fileNameRoot '_000' num2str(ddd) '.mat'];
                end
            else
                matfilename = [tmpPath fileNameRoot '.mat'];
            end

            data = matfile(matfilename,'Writable',true);
            Y(:,:,(ddd-1)*Tfile+1:ddd*Tfile) = single(reshape(data.vol(:,:,aaa,:),d1file,d2file,Tfile));
        end

        tt = toc/3600;
        disp(['Current plane loaded. Time elapsed: ' num2str(tt) ' hours. Beginning motion correction...'])
        
        date = datetime(now,'ConvertFrom','datenum');
        formatSpec = '%s Beginning processing for plane %u...\n';
        fprintf(fid,formatSpec,date,aaa);
        
        sizY = size(Y);
        d1 = sizY(1);
        d2 = sizY(2);
        T = sizY(3);

        Y = Y-min(Y(:));

        %% Motion correction

        disp('Starting the parallel pool...')
        poolobj = parpool('local',numCores);
        tmpDir = tempname();
        mkdir(tmpDir);
        poolobj.Cluster.JobStorageLocation = tmpDir;

        % Rigid motion correction using NoRMCorre algorithm:    
        options_rigid = NoRMCorreSetParms(...
            'd1',d1,...
            'd2',d2,...
            'bin_width',200,...       % Bin width for motion correction
            'max_shift',round(20/pixelResolution),...        % Max shift in px
            'us_fac',20,...
            'init_batch',200,...     % Initial batch size
            'correct_bidir',false... % Correct bidirectional scanning
            );

        [M1,shifts1,~,~] = normcorre_batch(Y,options_rigid);

        disp('Rigid motion correction complete. Beginning non-rigid motion correction...')

        shifts_r = squeeze(cat(3,shifts1(:).shifts));
        shifts_v = movvar(shifts_r,24,1);
    %     [~,minv_idx] = mink(shifts_v,120,1);
        [srt,minv_idx] = sort(shifts_v,120); 
        best_idx = unique(reshape(minv_idx,1,[]));
        template_good = mean(M1(:,:,best_idx),3);

        % No rigid motion correction using the good tamplate from the rigid
        % correction.
          options_nonrigid = NoRMCorreSetParms(...
            'd1',d1,...
            'd2',d2,...
            'bin_width',24,...
            'max_shift',round(20/pixelResolution),...
            'us_fac',20,...
            'init_batch',120,...
            'correct_bidir',false...
            );

        % Data from the motion correction that will be used for the CNMF
        [M2,shifts2,~,~] = normcorre_batch(Y,options_nonrigid,template_good);

        disp('Calculating motion correction metrics...')

        shifts_r = squeeze(cat(3,shifts1(:).shifts));
        shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
        shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
        shifts_x = squeeze(shifts_nr(:,1,:))';
        shifts_y = squeeze(shifts_nr(:,2,:))';

        [cY,~,~] = motion_metrics(Y,10);
        [cM1,~,~] = motion_metrics(M1,10);
        [cM2,~,~] = motion_metrics(M2,10);

        motionCorrectionFigure = figure;

        ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax2 = subplot(312); %plot(shifts_x); hold on; 
        plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
                set(gca,'Xtick',[])
        ax3 = subplot(313); %plot(shifts_y); hold on; 
        plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
                xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax1,ax2,ax3],'x')

        % Figure: Motion correction Metrics
        saveas(motionCorrectionFigure,[tmpPath 'motion_corr_metrics_plane_' num2str(aaa) '.fig']);
        close(motionCorrectionFigure)

        Y = M2;
        clear M2 M1 cM1 cM2 template_good shifts1 shifts2 shifts_nr shifts_r shifts_x shifts_y cY

        tt = toc/3600;
        disp(['Motion correction complete. Time elapsed: ' num2str(tt) ' hours. Saving to disk...'])

        outputFile = [tmpPath fileNameRoot '_plane_' num2str(aaa) '.mat'];

        Ym = mean(Y,3);

        savefast(outputFile,'Y','volumeRate','sizY','pixelResolution','Ym')

        disp('Data saved, beginning next plane...')
    end

    disp('All planes processed. Deleting temporary files...')

    for xyz = 1:numFiles
        if multiFile % if multiple files, append '_0000x.tif' or '_000xx.tif'
                if xyz < 10
                    fileName = [fileNameRoot '_0000' num2str(xyz) '.tif'];
                else
                    fileName = [fileNameRoot '_000'  num2str(xyz) '.tif'];
                end
            else % single file recording, append '.tif' to end
                fileName = [fileNameRoot '.tif'];
        end

        delete([tmpPath fileName(1:end-3) 'mat'])

    end
    t = toc;
    disp(['Routine complete. Total run time ' num2str(t./3600) ' hours.'])
    date = datetime(now,'ConvertFrom','datenum');
    formatSpec = '%s Routine complete.';
    fprintf(fid,formatSpec,date);

end
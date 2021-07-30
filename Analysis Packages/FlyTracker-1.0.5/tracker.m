
% Track (and calibrate) videos.
%
% To run tracker with interface, use:
%
%   tracker
%
% To run tracker without interface, use:
%
%   tracker(videos, [options], [f_calib], [vinfo])
%
%  where [] denotes an optional parameter (default values used if set to []) and:
%
%    videos.            - videos to process through tracking pipeline
%       dir_in          - directory containing input videos
%       dir_out         - directory in which to save results
%       filter          - file filter (eg '*.avi') (default: '*')
%       
%    options.           - cluster processing and output options
%       max_minutes     - maximum number of minutes to process (default: inf)
%       num_cores       - number of workers to process jobs in parallel? (default: 1)
%       granularity     - number of video frames per job (default: 10,000)
%       num_chunks      - number of chunks to process (default: num_frames/granularity)
%       save_JAABA      - write JAABA folders from features? (default: false)
%       save_xls        - save tracks and feats to xls? (default: false)
%       save_seg        - save segmentation from tracking process? (default: false)
%      
%    f_calib            - file containing calibration data (default: [videos.dir_in]/calibration.mat)
%    vinfo              - if specified, ignore videos and use loaded video
%
function tracker(videos, options, f_calib, vinfo)

   % add tracker to path if its not there already
   check = which('is_atomic_detection');
   if isempty(check)
       parentdir = fileparts(mfilename('fullpath'));
       addpath(genpath(parentdir));
   end   

   if nargin == 0
       % If tracker is started with no arguments, load interface
       display_available = feature('ShowFigureWindows');
       if display_available
          tracker_interface();          
       else          
          disp('No display available: run tracker with arguments:') 
          help tracker
          return
       end
   else
       % If tracker is run with arguments, make sure all fields are entered
       if nargin < 2, options = []; end
       if nargin < 3, f_calib = []; end
       if nargin < 4, vinfo = []; end
       run_tracker(videos,options,f_calib,vinfo);
   end
end

function run_tracker(videos, options, f_calib, vinfo)
   % default options
   options_def.granularity = 10000;
   options_def.num_chunks  = [];
   options_def.num_cores   = 1;     
   options_def.max_minutes = inf;      
   options_def.fr_samp     = 100;
   options_def.save_JAABA  = false;
   options_def.save_xls    = false;
   options_def.save_seg    = false;
   
   % set display variables
   display_available = feature('ShowFigureWindows');   
   fs = 72/get(0,'ScreenPixelsPerInch'); % scale fonts to resemble OSX
   dh = [];
   
   % fill in specified options
   if ((nargin < 2) || (isempty(options))), options = options_def; end
   options = set_defaults(options, options_def);
   
   % make sure we don't try to use more workers than available
   n_cores = feature('numCores');
   options.num_cores = min(n_cores,options.num_cores);
   % open parallel pool if not already open
   if options.num_cores > 1
       try
           open_pool = 1;
           if ~isempty(gcp('nocreate'))
               par = gcp;
               n_workers = par.NumWorkers;
               if n_workers == options.num_cores
                   open_pool = 0;
               else
                   delete(gcp);
               end
           end
           if open_pool
               parpool(options.num_cores);
           end
       catch
           options.num_cores = 1;
           str = 'Could not open parallel pool. Using single thread.';
           disp(str);
       end
   end
   
   % collect video information
   if nargin > 3 && ~isempty(vinfo)
       % use given vinfo rather than videos
       [path,filename,ext] = fileparts(vinfo.filename);
       videos.dir_in = path;
       if ~isfield(videos,'dir_out'), videos.dir_out = path; end
       videos.filter = [filename ext];
       vid_files = {vinfo.filename};
   else
       % convert input/output directories to absolute path form
       videos.dir_in  = absolute_path(videos.dir_in);
       videos.dir_out = absolute_path(videos.dir_out);
       % check video list
       if ((~isfield(videos,'filter')) || (isempty(videos.filter)))
          videos.filter = '*';
       end   
       % scan input directory for video files
       vid_files = dir(fullfile(videos.dir_in, videos.filter));
       vid_files([vid_files.isdir]) = [];
       vid_files = { vid_files.name };       
   end
   
   % make sure there are videos to process
   n_vids = numel(vid_files);
   if n_vids < 1
       str = ['No videos to process for: ' fullfile(videos.dir_in, videos.filter)];
       if display_available
          customDialog('warn',str,12*fs);
       else
          disp(str);
       end
       return
   end
   
   % create output directories
   if (~exist(videos.dir_out,'dir')), mkdir(videos.dir_out); end
   for n = 1:n_vids
      [~, name] = fileparts(vid_files{n});
      dir_vid = fullfile(videos.dir_out, name);
      if (~exist(dir_vid,'dir')), mkdir(dir_vid); end
   end
   
   % load calibration file
   if nargin < 3 || isempty(f_calib)
      f_calib = fullfile(videos.dir_in, 'calibration.mat');
   end
   if ~exist(f_calib,'file')
       str = [f_calib ' not found: run calibrator first or input a valid calibration file.'];
       if display_available
          customDialog('warn',str,12*fs);
       else
          disp(str);
       end
       return
   end
   D = load(f_calib); parent_calib = D.calib;

   % compute maximum number of frames to process
   max_frames = options.max_minutes*parent_calib.FPS*60;
   min_chunksize = 100;
   
   % package jobs to be run in sequence
   for n = 1:n_vids
      % get input video file
      [~, name,ext] = fileparts(vid_files{n});
      f_vid = fullfile(videos.dir_in, vid_files{n});
      % display progress
      waitstr = ['Processing ' name ext ...
          '   (movie ' num2str(n) '/' num2str(n_vids) ')'];      
      if display_available 
          if isempty(dh) || ~ishandle(dh)            
            dh = customDialog('wait',waitstr,12*fs);
          else            
            child = get(dh,'Children');
            set(child(1),'string',waitstr);
          end
      else
          disp(['*** ' waitstr])
      end      
      % check whether video has already been tracked
      f_res_final = fullfile(videos.dir_out, name, [name '-track.mat']);
      if exist(f_res_final,'file')
          disp('Movie already tracked')
          % compute features and learning files if specified
          tracker_job('track_features',f_vid,f_res_final,f_calib,options,0);          
          continue;
      end
      % load video
      do_close = 0;
      if nargin < 4 || isempty(vinfo)
        vinfo = video_open(f_vid);
        do_close = 1;
      end
      % get length of video 
      n_frames = min(vinfo.n_frames,max_frames); 
      % output filenames
      f_bg  = fullfile(videos.dir_out, name, [name '-bg.mat']);
      if n_vids > 1 && parent_calib.auto_detect
          % generate new calibration files if multiple videos 
          f_calib = fullfile(videos.dir_out, name, [name '-calibration.mat']);
      end
      % compute background and calibration if needed
      flag = tracker_job('track_calibrate', f_vid, f_bg, f_calib, options, parent_calib, vinfo);
      if check_abort(flag), return; end
      % load calibration 
      D = load(f_calib); calib = D.calib;      
      % compute number of chunks to process
      if ~isempty(options.num_chunks)
          n_chunks = options.num_chunks;
          chunksize = ceil(n_frames/n_chunks);
          options.granularity = max(chunksize,min_chunksize);
      end
      n_chunks = ceil(n_frames./options.granularity);
      % loop through all chambers
      valid = find(calib.valid_chambers);
      n_chambers = numel(valid);      
      % process tracks for all chambers
      % set frame range
      frs = cell(1,n_chunks);
      f_trks = cell(1,n_chunks);
      for c=1:n_chunks
         fr.start = (c-1).*options.granularity;
         fr.step  = 1;
         fr.limit = min(fr.start + options.granularity, n_frames);
         frs{c} = fr;
         f_trks{c} = cell(1,n_chambers);
         for i=1:n_chambers
             if n_chambers == 1
                 chamber_str = '';
             else
                 chamber_str = ['-c' num2str(i)];
             end    
             % determine output filename
             f_trk = fullfile( ...
                videos.dir_out, name, ...
                [name chamber_str '-trk' '-' num2str(frs{c}.start,'%010d') '.mat'] ...
             );     
            f_trks{c}{i} = f_trk;
         end
      end
      % check whether chamber files already exist
      f_res = fullfile(videos.dir_out, name, [name chamber_str '-track.mat']);
      if ~exist(f_res,'file') 
          if options.num_cores > 1
              success = zeros(1,n_chunks);              
              parfor c = 1:n_chunks
                 % store job parameters
                 if strcmp(ext,'.ufmf')
                     % ufmf cannot be processed in parallel with single fid, must create a new one for each chunk
                     flag = tracker_job('track_process', f_vid, f_bg, f_calib, f_trks{c}, frs{c});
                 else
                     flag = tracker_job('track_process', f_vid, f_bg, f_calib, f_trks{c}, frs{c}, vinfo);
                 end
                 success(c) = flag;
              end
              if check_abort(min(success)), return; end
          else
              for c = 1:n_chunks
                 % store job parameters
                 flag = tracker_job('track_process', f_vid, f_bg, f_calib, f_trks{c}, frs{c}, vinfo);
                 if check_abort(flag), return; end
              end
          end      
      end      
      % combine results for all chunks (per chamber)
      f_res_list = cell(1,n_chambers);
      for i=1:n_chambers
          if n_chambers == 1
             chamber_str = '';
          else
             chamber_str = ['-c' num2str(i)];
          end     
          f_trk_list = cell(1,n_chunks);
          for c=1:n_chunks
              f_trk_list{c} = f_trks{c}{i};
          end
          f_res = fullfile(videos.dir_out, name, [name chamber_str '-track.mat']);
          f_res_list{i} = f_res;
          if n_chambers > 1
              chamber_str = ['c' num2str(i) ' - '];
          end
          tracker_job('track_combine', f_res, f_trk_list, f_calib, options, chamber_str);
      end
      % combine results from chambers
      if n_chambers > 1
        tracker_job('track_consolidate', f_res_final, f_res_list, options);
      end
      % compute features and learning files if specified
      tracker_job('track_features', f_vid, f_res_final, f_calib, options, 0);
      
      % close video
      if do_close
        video_close(vinfo);
        vinfo = [];
      end
   end
   
   if ~isempty(dh) && ishandle(dh)      
       child = get(dh,'Children');
       set(child(1),'string','Finished!');
       pause(2)
       delete(dh); 
   end
   
   function do_abort = check_abort(flag)
      do_abort = 0;
      if ~flag
          if ~isempty(dh) && ishandle(dh), 
              child = get(dh,'Children');
              set(child(1),'string','** canceled by user **');
              pause(2)
              delete(dh); 
          end
          do_abort = 1;
      end
   end
end

%% TRACKER INTERFACE
function tracker_interface()
    videos.dir_in = '';
    videos.dir_out = '';
    videos.filter = '*';  
    options = [];
    f_calib = '';
    
    % specify valid extensions
    videoReaderFmts = VideoReader.getFileFormats();
    commonFmts = {'avi','mov','mp4','wmv','m4v','mpg'};
    specialFmts = {'fmf','sbfmf','ufmf','seq','bin'};
    valid_extns = union(commonFmts,specialFmts);
    valid_extns = union({videoReaderFmts.Extension},valid_extns);  
    
    % ----- LOAD INTERFACE -----   
    % MAIN WINDOW
    scrsz = get(0,'ScreenSize');
    fig_width = 620;
    fig_height = 500;
    fig_h = figure('Position',[scrsz(3)/2-fig_width/2 scrsz(4)/2-fig_height/2 fig_width fig_height],...
        'Name','FlyTracker-1.0.5','NumberTitle','off','Color',.94*[1 1 1]);
    set(fig_h,'MenuBar','none')
    set(fig_h,'Resize','off')
    set(fig_h,'CloseRequestFcn',@ui_close)
    figclr = get(fig_h,'color');
    fs = 72/get(0,'ScreenPixelsPerInch'); % scale fonts to resemble the mac
    
    % title
    text_x = 20;
    text_y = fig_height-50;
    uicontrol('Style', 'text', 'String', 'Select files:', ...
    'Position',[text_x text_y fig_width-40 30], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*14, 'FontWeight', 'bold');
    
    % SELECT READ FOLDER
    text_y = text_y - 40;
    uicontrol('Style', 'pushbutton', 'String', 'VIDEO folder', ...
    'Position',[text_x text_y-3 160 35], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*12, ...
    'Callback', @selectReadFolder, ...
    'ToolTipString','Specify folder containing videos to process');
    f_read_h = uicontrol('Style', 'edit', 'String', videos.dir_in, ...
    'Position',[text_x+170 text_y 300 30], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*10);
    extn_h = uicontrol('Style', 'popup', 'String', 'extn', ...
    'Position',[text_x+490 text_y-2 80 30], ...
    'Value',1,'BackgroundColor',figclr,...
    'FontSize',fs*12, ...
    'Callback',@setExtension,...
    'ToolTipString','Only videos with selected extension will be processed');  
    
    % SELECT SAVE FOLDER
    text_y = text_y - 50;
    uicontrol('Style', 'pushbutton', 'String', 'OUTPUT folder', ...
    'Position',[text_x text_y-3 160 35], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*12, ...
    'Callback', @selectSaveFolder, ...
    'ToolTipString','Specify folder to which output will be saved');
    f_save_h = uicontrol('Style', 'edit', 'String', videos.dir_out, ...
    'Position',[text_x+170 text_y 300 30], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*10);   

    % CALIBRATION FILE
    text_y = text_y - 50;
    uicontrol('Style', 'pushbutton', 'String', 'calibration file', ...
    'Position',[text_x text_y-3 160 35], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*12,...
    'Callback', @selectCalibrationFile, ...
    'ToolTipString','Select existing calibration file');     
    f_calib_h = uicontrol('Style', 'edit', 'String',  f_calib, ...
    'Position',[text_x+170 text_y 300 30], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*10);
    uicontrol('Style', 'text', 'String', '/', ...
    'Position',[text_x+470 text_y-15 20 40], ...
    'HorizontalAlignment', 'center', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*12);     
    uicontrol('Style', 'pushbutton', 'String', 'CALIBRATE', ...
    'Position',[text_x+490 text_y-3 80 35], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [.9 .85 .5], ...
    'FontSize',fs*10, ...
    'Callback', @runCalibrator, ...
    'ToolTipString','Create new calibration file');

    % title
    text_y = text_y - 30;
    uicontrol('Style', 'text', 'String', '', ...
    'Position',[text_x text_y fig_width-40 1], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr*.8, ...
    'FontSize',fs*14, 'FontWeight', 'bold');
    text_y = text_y - 50;
    uicontrol('Style', 'text', 'String', 'Options:', ...
    'Position',[text_x text_y fig_width-40 30], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr, ...
    'FontSize',fs*14, 'FontWeight', 'bold');

    % OPTIONS
    text_y = text_y - 40;
    % max minutes to track
    uicontrol('Style', 'text', 'String', 'Process:', ...
        'Position',[text_x text_y 60 25], ...
        'HorizontalAlignment', 'right', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*13);   
    text_x = text_x + 70;
    uicontrol('Style', 'text', 'String', 'max minutes:', ...
        'Position',[text_x text_y 87 25], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12, ...
        'ToolTipString','Upper limit on number of minutes to process');        
    max_h = uicontrol('Style', 'edit', 'String', 'Inf', ...
        'Position',[text_x+87 text_y+2 50 25], ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12);    
    % chunksize 
    chunktype_h = uicontrol('Style', 'popup', ...
        'String', 'chunksize (frames):|number of chunks:', ...
        'Position',[text_x+158 text_y 155 25], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12,...
        'Callback',@setChunkType, ...
        'ToolTipString','Each video is processed in chunks of frames');    
    chunk_h = uicontrol('Style', 'edit', 'String', '10000', ...
        'Position',[text_x+311 text_y+2 54 25], ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12);    
    % use parallel pool     
    n_cores = feature('numCores');
    if n_cores > 1
        popstring = 'use 1 core';
        for i=2:n_cores
            popstring = [popstring '|use ' num2str(i) ' cores'];
        end
        par_h = uicontrol('Style', 'popup', 'String',popstring, ...
            'Position',[text_x+390 text_y 120 25], ...
            'HorizontalAlignment', 'left', ...
            'BackgroundColor', figclr, ...
            'FontSize',fs*12, ...
            'Callback',@setNCores,...
            'ToolTipString','Process chunks in parallel on multiple cores');    
    end
    text_x = text_x - 70;
    % OUTPUT OPTIONS
    text_y = text_y - 45;
    uicontrol('Style', 'text', 'String', 'Extra:', ...
        'Position',[text_x text_y 60 25], ...
        'HorizontalAlignment', 'right', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*13);   
    % output writeJAABA folders
    jab_h = uicontrol('Style', 'checkbox', 'String', 'save JAABA folders', ...
        'Position',[text_x+70 text_y+2 200 30], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12, ...
        'ToolTipString','Save output to JAABA compatible folders');    
    xls_h = uicontrol('Style', 'checkbox', 'String', 'save to .xls', ...
        'Position',[text_x+228 text_y+2 200 30], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12, ...
        'ToolTipString','Save output to .xls');    
    seg_h = uicontrol('Style', 'checkbox', 'String', 'save segmentation', ...
        'Position',[text_x+340 text_y+2 200 30], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', figclr, ...
        'FontSize',fs*12, ...
        'ToolTipString','Save video segmentation (this produces large files)');    

    % CLOSE and TRACK BUTTONS   
    text_y = text_y - 30;
    uicontrol('Style', 'text', 'String', '', ...
    'Position',[text_x text_y fig_width-40 1], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', figclr*.8, ...
    'FontSize',fs*14, 'FontWeight', 'bold');
    text_y = text_y - 70;
    text_x = text_x + 170;
    uicontrol('Style', 'pushbutton', 'String', 'CLOSE', ...
    'Position',[text_x text_y-5 80 40], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [.9 .6 .7], ...
    'FontSize',fs*12, ...
    'Callback', @ui_close, ...
    'ToolTipString','Close without tracking');    
    uicontrol('Style', 'pushbutton', 'String', 'TRACK', ...
    'Position',[text_x+125 text_y-5 140 40], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [.4 .8 .5], ...
    'FontSize',fs*12, ...
    'Callback', @finishAndTrack, ...
    'ToolTipString','Track all videos in READ folder');
    
    % CALLBACK FUNCTIONS
    function setChunkType(hObj,event) %#ok<INUSD>
        value = get(hObj,'Value');
        if value == 1
            set(chunk_h,'String',num2str(10000));
        else
            n_chunks = 10;
            if n_cores > 1
                n_chunks = get(par_h,'Value')*2;                
            end
            set(chunk_h,'String',num2str(n_chunks));
        end
    end
    function setNCores(hObj,event) %#ok<INUSD>
        options.num_cores = get(hObj,'Value');                
        if options.num_cores > 1
            set(chunk_h,'String',num2str(options.num_cores*2));
            set(chunktype_h,'Value',2);
        else
            set(chunk_h,'String',num2str(10000));
            set(chunktype_h,'Value',1);
        end
    end
    function setExtension(hObj,event) %#ok<INUSD>
        string = get(hObj,'String');        
        if ~iscell(string), 
            return; 
        end
        idx = get(hObj,'Value');
        val = string{idx};
        if isempty(val)
            val = '*';
        elseif ~strcmp(val(1),'*')
            val = ['*' val];
        end
        videos.filter = val;
    end
    function updateExtensions()
        directory = videos.dir_in;
        files = dir(fullfile(directory,'*'));
        extns = cell(1,numel(files));
        count_valid = 0;
        for f=1:numel(files)
            [~,~,extn] = fileparts(files(f).name);  
            if any(strcmpi(extn(2:end),valid_extns))
                count_valid = count_valid+1;
                extns{count_valid} = extn;
            end
        end
        extns = extns(1:count_valid);
        extns = unique(extns);        
        if numel(extns)>0
            set(extn_h,'string',extns);  
            videos.filter = ['*' extns{1}];
        end
    end
    function selectReadFolder(hObj,event) %#ok<INUSD>
        directory = uigetdir(videos.dir_out,'Select READ folder');
        if ~directory, directory = ''; end
        videos.dir_in = directory;
        set(f_read_h,'String',directory);
        % update extensions
        updateExtensions;        
        % update save directory
        videos.dir_out = directory;
        set(f_save_h,'String',directory);            
        % update calibration file
        f_calib = fullfile(videos.dir_in,'calibration.mat');
        if ~exist(f_calib,'file')
           f_calib = '';
        else 
           set(f_calib_h,'String',f_calib);
        end
    end 
    function selectSaveFolder(hObj,event) %#ok<INUSD>
        directory = uigetdir(videos.dir_in,'Select SAVE folder');
        if ~directory, directory = ''; end
        videos.dir_out = directory;
        set(f_save_h,'String',directory);
        if numel(get(f_read_h,'String')) == 0
            videos.dir_in = directory;
            set(f_read_h,'String',directory);
        end
    end
    function selectCalibrationFile(hObj,event) %#ok<INUSD>
        [file,path] = uigetfile('*.mat','Select calibration file',videos.dir_in);
        f_calib = fullfile(path,file);
        if ~f_calib, f_calib = ''; end
        set(f_calib_h,'String',f_calib);
    end
    function runCalibrator(hObj,event) %#ok<INUSD>
       vid_files = dir(fullfile(videos.dir_in, videos.filter));
       vid_files([vid_files.isdir]) = [];
       vid_files = { vid_files.name };           
       valid = false(size(vid_files));
       for f=1:numel(vid_files)
          [~,~,extn] = fileparts(vid_files{f});  
          if any(strcmpi(extn(2:end),valid_extns))
              valid(f) = 1;
          end
       end
       vid_files = vid_files(valid);
       if numel(vid_files) == 0
           customDialog('warn','No valid video folder selected',12*fs);
           return;
       end
       f_vid = fullfile(videos.dir_in, vid_files{1});
       f_calib = fullfile(videos.dir_in,'calibration.mat');
       calib_success = calibrator(f_vid,f_calib);
       if ~calib_success, f_calib = ''; end
       set(f_calib_h,'String',f_calib);
    end
    function finishAndTrack(hObj,event) %#ok<INUSD>
        % collect file information
        videos.dir_in = get(f_read_h,'String');
        videos.dir_out = get(f_save_h,'String');
        if isempty(videos.dir_in) || isempty(videos.dir_out) || ...
                isempty(videos.filter) || ~exist(f_calib,'file')
            customDialog('warn','File information incomplete',12*fs);
            return;
        end
        % collect options
        str = get(max_h,'String');
        minutes = str2double(str);
        if ~isempty(minutes)
            options.max_minutes = minutes;
        end        
        str = get(chunk_h,'String');
        value = str2double(str);
        type = get(chunktype_h,'Value');        
        if ~isempty(value)
            if type == 1
                options.granularity = value;
            else
                options.num_chunks = value;
            end    
        end
        value = get(jab_h,'Value');
        options.save_JAABA = value;
        value = get(xls_h,'Value');
        options.save_xls = value;
        value = get(seg_h,'Value');
        options.save_seg = value;
        % track
        delete(fig_h)
        pause(.5)
        run_tracker(videos,options,f_calib);
    end
    function ui_close(hObj,event) %#ok<INUSD>
        delete(fig_h);     
    end    
end

%% Help functions
% Set default values for parameter fields not specified by the user.
%
%    params = set_defaults(params, params_def)
%
% copies any fields specified in params_def but not params into params.  If
% params is empty, then params is set to params_def.
function params = set_defaults(params, params_def)
   if (isempty(params))
      % return default parameters
      params = params_def;
   else
      % set default values for any unspecified parameters
      names = setdiff(fieldnames(params_def),fieldnames(params));
      for n = 1:numel(names)
         params = setfield(params, names{n}, getfield(params_def, names{n}));
      end
   end
end

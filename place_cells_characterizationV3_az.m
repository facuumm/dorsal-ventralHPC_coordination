clear
clc
close all
%% Parameters
path = {'\\Maryjackson\e\Rat127\Ephys\pyr';'\\Maryjackson\e\Rat128\Ephys\in_pyr\ready';'\\Maryjackson\e\Rat103\usable';'\\Maryjackson\e\Rat132\recordings\in_pyr'; '\\Maryjackson\e\Rat165\in_pyr'};%List of folders from the path
% path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr'};%List of folders from the path
% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction
sigma = 2;%round(15/(180/Xedges)); %defined for gauss kernel of 15cm
binSize = 0.001; % bin size for replay events detection
bin_size = 1; % to bin pos ans spks in between/within  
% Behavior
minimal_speed = 2.5;% minimal speed to detect quite periods
minimal_speed_time = 2; % minimal time to detect quite periods

%% Main loop, to iterate across sessions
for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    clear files dirFlags
    
    for t = 1 : length(subFolders)-2
        disp(['-- Initiating analysis of folder #' , num2str(t) , ' from rat #',num2str(tt) , ' --'])
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        %Loading TS of the sessions
        disp('Uploading session time stamps')
        x = dir([cd,'\*.cat.evt']);
        segments = readtable([cd,'\',x.name],'FileType','text');
        clear x
        % TimeStamps of begening and end of the sleep and awake trials
        % Reward and Aversive trials
        aversiveTS = [];
        aversiveTS_run = [];
        rewardTS = [];
        rewardTS_run = [];
        for y = 1 : height(segments)
            % Baseline sleep session TS detection
            if y == 1
                baselineTS(1,1) = segments.Var1(y);
            elseif y ==2
                baselineTS(1,2) = segments.Var1(y);
            end
            % Aversive sleep session TS detection
            if strcmp(segments.Var2{y},'aversive')
                if strcmp(segments.Var3{y},'End')
                    aversiveTS(1,1) = segments.Var1(y+1);
                    aversiveTS(1,2) = segments.Var1(y+2);
                    aversiveTS_run(1,1) = segments.Var1(y-1);
                    aversiveTS_run(1,2) = segments.Var1(y);
                    A = y;
                end
                % Rewarded sleep session TS detection
            elseif strcmp(segments.Var2{y},'reward')
                if strcmp(segments.Var3{y},'End')
                    rewardTS(1,1) = segments.Var1(y+1);
                    rewardTS(1,2) = segments.Var1(y+2);
                    rewardTS_run(1,1) = segments.Var1(y-1);
                    rewardTS_run(1,2) = segments.Var1(y);
                    R = y;
                end
            end
        end
        clear y
        
        %Session output
        dHPC = {};% one cell per pc
        vHPC = {}; 
        
        %% Awake
        disp('Uploading digital imputs')
        % Load digitalin.mat
        load('digitalin.mat')
        %Shocks selection
        Shocks_filt = Restrict(shock,aversiveTS_run ./1000);
        % Keep only the first shock of each TTL (first from 20)
        count = 1;
        deff = [];
        for i = 1:length(Shocks_filt)
            if count == 1
                deff = [deff ; Shocks_filt(i,1)];
            end
            if count ==20
                count = 0;
            end
            count = count + 1;
        end
        Shocks_filt = deff;
        clear count deff shock i
        %Rewards selection
        Rewards_filt = Restrict([leftvalve ; rightvalve],rewardTS_run ./1000);
        disp('Uploading DLC outputs')
        camara = ((camara(:,2)-camara(:,1))/2)+camara(:,1);
        % periods of movment during eacj condition
        if rewardTS_run(1) < aversiveTS_run(1)
            load('laps1.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraR
        else
            load('laps2.mat','posx','posy');
            [camaraR,~] = find((camara(:,1)-rewardTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            %         [camaraR2,~] = find((camara(:,1)-rewardTS_run(2)/1000)<0,1,'last'); %TimeStamp of the ending of aversive
            pos = [camara(camaraR : camaraR+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            behavior.pos.reward = [pos];
            behavior.speed.reward = LinearVelocity(pos,0);
            behavior.speed.reward(:,2) = smoothdata(behavior.speed.reward(:,2),'gaussian',1,'SamplePoints',behavior.speed.reward(:,1));
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , minimal_speed , minimal_speed_time);
            clear pos camaraR posx posy
            load('laps1.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int; %saving corrected pos
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,0);
            behavior.speed.aversive(:,2) = smoothdata(behavior.speed.aversive(:,2),'gaussian',1,'SamplePoints',behavior.speed.aversive(:,1));
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
        % Generation of no-movements periods
        % Reward
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        movement.reward = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
        movement.reward(movement.reward(:,2) - movement.reward(:,1) <1,:)=[];
        clear tmp start stop
        % Aversive
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        movement.aversive = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        movement.aversive(movement.aversive(:,2) - movement.aversive(:,1) <1,:)=[];
        clear tmp start stop

        %% load sleep states
%         disp('Uploading sleep scoring')
%         x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
%         REM.all = ToIntervals(states==5);    NREM.all = ToIntervals(states==3);    WAKE.all = ToIntervals(states==1);
%         clear x states
%         NREM.baseline = Restrict(NREM.all,baselineTS./1000);
%         NREM.aversive = Restrict(NREM.all,aversiveTS./1000);
%         NREM.reward = Restrict(NREM.all,rewardTS./1000);
%         REM.baseline = Restrict(REM.all,baselineTS./1000);
%         REM.aversive = Restrict(REM.all,aversiveTS./1000);
%         REM.reward = Restrict(REM.all,rewardTS./1000);
        
        %% load coordinated ripple bursts
%         load('coordinated_ripple_bursts.mat')
%         ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
%         ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
%         ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
%         ripple_bursts.all = coordinated_ripple_bursts;
%         clear coordinated_ripple_bursts
        
        %% Load ripples
%         ripplesD = table2array(readtable('ripplesD_customized2.csv'));
%         ripplesV = table2array(readtable('ripplesV_customized2.csv'));
%         % coordination
%         coordinated = [];
%         coordinatedV = [];
%         coordinatedV_refined = [];
%         for i = 1:length(ripplesD)
%             r = ripplesD(i,:);
%             tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
%             if tmp>0
%                 z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
%                 coordinatedV = [coordinatedV ; z];
%                 [p,indice] = min(abs(r(2)-z(:,2)));
%                 coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
%                 coordinated = [coordinated ; r];
%                 clear tmp2 tmp1 p indice z
%             end
%             clear r
%         end
%         clear x tmp i
%         
        %% Spikes
        %Load Units
        disp('Uploading Spiking activity')
        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters
        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:7)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
        %Loop to select dorsal or ventral LFP and SU
        % z=1 --> dorsal
        % z=2 --> ventral
        for z = 1:2
            if z == 1
                n_SU_D = [n_SU_D ; length(group_dHPC)];
                spks_dHPC = spks(ismember(spks(:,1),group_dHPC(:,1)),:); %keep spks from good clusters
            else
                n_SU_V = [n_SU_V ; length(group_vHPC)];
                spks_vHPC = spks(ismember(spks(:,1),group_vHPC(:,1)),:); %keep spks from good clusters
            end
        end
        clear z
        spks_vHPC(:,2) = double(spks_vHPC(:,2))./20000;
        spks_dHPC(:,2) = double(spks_dHPC(:,2))./20000;
        % Selection of celltype to analyze
        if criteria_type == 0 %pyr
            cellulartype = [K(:,1) , K(:,3)];
        elseif criteria_type == 1 % int
            cellulartype = [K(:,1) , K(:,4)];
        elseif criteria_type == 2 % all
            cellulartype = [K(:,1) , ones(length(K),1)];
        end
        % Definition of the emotional transition
        if aversiveTS_run(1)<rewardTS_run(1)
            cond = 1; % 1 if the order was Aversive -> Reward
        else
            cond = 2;% 1 if the order was Reward -> Aversive
        end
        
        %% Firing maps calculation
        disp('dHPC Firing rate map calculation')
        maps_dHPC_A = [];        maps_dHPC_R = [];
        maps_dHPC_SA = [];        maps_dHPC_SR = [];
        peak_dHPC_A = [] ;       peak_dHPC_R = [];
        field_dHPC_A = [];       field_dHPC_R = [];
        PC.dHPC.aversive=[];     PC.dHPC.reward=[];
        cluster_dHPC = [];   
      
        
        pc_dHPC_par.within.ave = []; pc_dHPC_par.within.rew = []; 
        pc_dHPC_par.within.ave_tresh = []; pc_dHPC_par.within.rew_tresh = []; 
        pc_dHPC_par.between = []; 
        
        pc_dHPC_par.firingMap.ave = []; 
        pc_dHPC_par.firingMap.rew = [];
        
        pc_dHPC_par.pc_params.ave = []; 
        pc_dHPC_par.pc_params.rew = [];
       
        for ii=1:size(group_dHPC,1)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
            if celltype % check if pyr
                
                spks = spks_dHPC(spks_dHPC(:,1)==cluster,2); % select tspk from cluster
              
                % Select neurons with fr greater than tresh in aversive or reward - PC criteria 1   
                tresh = 0.1; % hz
                spks_tmp = Restrict(spks,movement.reward);
                fr_rew= size(spks_tmp,1)/sum(movement.reward(:,2)-movement.reward(:,1));
                
                spks_tmp = Restrict(spks,movement.aversive);
                fr_ave= size(spks_tmp,1)/sum(movement.aversive(:,2)-movement.aversive(:,1));
                
                if fr_rew >= tresh || fr_ave >= tresh
                    m = 1;
                else 
                    m = nan;
                end 
                
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive); % Restrict pos to movement periods
                    pos_tmp = pos_tmp(and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170),:); % eliminating the extrems of the maze
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
  
                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4, 'minPeak' , 0.2);
                    curveA_dhpc = curveA.rate; 
                    
                    %%%Within-trial pc parameters 
                    [withinA,withinA_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
                    
                    %Store par for between comp
                    spks_ave = spks_tmp;
                    pos_ave = pos_tmp;
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward); % Restrict to movement periods
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.2);
                    curveR_dhpc = curveR.rate;
                    %%%Within-trial pc parameters 
                    [withinR,withinR_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
                   
                    %Store par for between comp
                    spks_rew = spks_tmp;
                    pos_rew = pos_tmp;
                    
                    %%%PC criteria: have a PF bigger or equal to 4 bins 
                    
                    %if there is no PF, assigne 0 to correct assesment in
                    %the next if
                    if isempty(statsA.field)
                        statsA.field=0;
                    end 
                     if isempty(statsR.field)
                        statsR.field=0;
                    end 
                    
                    if or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4) 
                        
                        %Store pc info 
                        n.id = cluster; 
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        dHPC{ii}= n;
                        
                        %Store firing maps
                        pc_dHPC_par.firingMap.ave = [pc_dHPC_par.firingMap.ave;curveA_dhpc]; 
                        pc_dHPC_par.firingMap.rew = [pc_dHPC_par.firingMap.rew;curveR_dhpc];
                        
                       
                        %Sparsity: 
                        sa = sparsity_info(curveA.rate,curveA.time);
                        sr = sparsity_info(curveR.rate,curveR.time);
                        %Store pc parametres 
                        pc_dHPC_par.pc_params.ave = [pc_dHPC_par.pc_params.ave;statsA.specificity,sa,statsA.size(1),statsA.peak(1)]; 
                        pc_dHPC_par.pc_params.rew = [pc_dHPC_par.pc_params.rew;statsR.specificity,sr,statsR.size(1),statsR.peak(1)];
                        
                        % Aversive 
                        [ff,f] = max(curveA.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        fr_ave = nanmean(curveA.rate);
                        PC.dHPC.aversive=[PC.dHPC.aversive ; cluster , ff , fff(f),fr_ave];
                        clear f ff fff
                        
                        maps_dHPC_A = [maps_dHPC_A ; curveA.rate];
%                         maps_dHPC_SA = [maps_dHPC_SA ; curveSA.rate];
                        map_ave = curveA.rate;
                        peak_dHPC_A = [peak_dHPC_A ; statsA.peak(1)];
                        field_dHPC_A = [field_dHPC_A ; statsA.fieldX(1,:)];
                        
                        %Rewarded 
                        [ff,f] = max(curveR.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        fr_rew = nanmean(curveR.rate);
                        PC.dHPC.reward=[PC.dHPC.aversive ; cluster , ff , fff(f),fr_rew];
                        clear f ff fff
                        
                        maps_dHPC_R = [maps_dHPC_R ; curveR.rate];
%                         maps_dHPC_SR = [maps_dHPC_SR ; curveSR.rate];
                        map_rew = curveR.rate;
                        peak_dHPC_R = [peak_dHPC_R ; statsR.peak(1)] ;
                        field_dHPC_R = [field_dHPC_R ; statsR.fieldX(1,:)];
                        cluster_dHPC = [cluster_dHPC ; cluster , cond];
                        
                        %Store within 
                        pc_dHPC_par.within.ave = [pc_dHPC_par.within.ave; withinA];
                        pc_dHPC_par.within.rew = [pc_dHPC_par.within.rew; withinR];
                        
                        pc_dHPC_par.within.ave_tresh = [pc_dHPC_par.within.ave_tresh; withinA_tresh];
                        pc_dHPC_par.within.rew_tresh = [pc_dHPC_par.within.rew_tresh; withinR_tresh];
                        
                        %Between-trial pc parameters (aversive vs. rewarded)
                        [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);

                        %Save between pc parameters 
                        pc_dHPC_par.between = [pc_dHPC_par.between; between];
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
            
        end
        
        disp('vHPC Firing rate map calculation')
        maps_vHPC_A = [];        maps_vHPC_R = [];
        maps_vHPC_SA = [];        maps_vHPC_SR = [];
        peak_vHPC_A = [] ;       peak_vHPC_R = [];
        field_vHPC_A = [];       field_vHPC_R = [];
        PC.vHPC.aversive=[];     PC.vHPC.reward=[];
        cluster_vHPC = [];
        
        pc_vHPC_par.within.ave = []; pc_vHPC_par.within.rew = []; 
        pc_vHPC_par.within.ave_tresh = []; pc_vHPC_par.within.rew_tresh = []; 
        pc_vHPC_par.between = []; 
        pc_vHPC_par.firingMap.ave = []; 
        pc_vHPC_par.firingMap.rew = [];
        
        pc_vHPC_par.pc_params.ave = []; 
        pc_vHPC_par.pc_params.rew = [];
        
    for ii=1:size(group_vHPC,1)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            
               if celltype % check if pyr
                
                spks = spks_vHPC(spks_vHPC(:,1)==cluster,2); % select tspk from cluster
              
                % Select neurons with fr greater than tresh in aversive or reward - PC criteria 1   
                tresh = 0.1; % hz
                spks_tmp = Restrict(spks,movement.reward);
                fr_rew= size(spks_tmp,1)/sum(movement.reward(:,2)-movement.reward(:,1));
                
                spks_tmp = Restrict(spks,movement.aversive);
                fr_ave= size(spks_tmp,1)/sum(movement.aversive(:,2)-movement.aversive(:,1));
                
                if fr_rew >= tresh || fr_ave >= tresh
                    m = 1;
                else 
                    m = nan;
                end 
                
                if ~isnan(m)
                    % --- Aversive ---
                    spks_tmp = Restrict(spks , movement.aversive);
                    pos_tmp = Restrict(behavior.pos.aversive(:,1:2) , movement.aversive); % Restrict pos to movement periods
                    pos_tmp = pos_tmp(and(pos_tmp(:,2)>30 , pos_tmp(:,2)<170),:); % eliminating the extrems of the maze
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveA , statsA] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.1);
                    curveA_vhpc = curveA.rate;
                    %%%Within-trial pc parameters 
                    [withinA,withinA_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
       
                    %Store par for between comp
                    spks_ave = spks_tmp;
                    pos_ave = pos_tmp;
                    
                    % --- Reward ---
                    spks_tmp = Restrict(spks , movement.reward); % Restrict to movement periods
                    pos_tmp = Restrict(behavior.pos.reward(:,1:2) , movement.reward);
                    pos_tmp(:,2) = pos_tmp(:,2)-min(pos_tmp(:,2)); pos_tmp(:,2) = pos_tmp(:,2)/max(pos_tmp(:,2)); %normalization of position
                    
                    %Firing curve construction
                    [curveR , statsR] = FiringCurve(pos_tmp , spks_tmp , 'smooth' , sigma , 'nBins' , Xedges , 'minSize' , 4 , 'minPeak' , 0.1);
                    curveR_vhpc = curveR.rate;
                    %%%Within-trial pc parameters 
                    [withinR,withinR_tresh] = Within_pc(pos_tmp,spks_tmp,1,sigma,Xedges);
                    
                    %Store par for between comp
                    spks_rew = spks_tmp;
                    pos_rew = pos_tmp;
                    
                    %PC cirteria: have a PF of at least 4 bins 
                    %if there is no PF, assigne 0 to correct assesment in
                    %the next if
                    if isempty(statsA.field)
                        statsA.field=0;
                    end 
                     if isempty(statsR.field)
                        statsR.field=0;
                    end 
                    
                    if or(sum(statsA.field(:,:,1))>= 4 , sum(statsR.field(:,:,1))>= 4) 
                        
                        %Save pc info 
                        n.id = cluster; 
                        n.frMap_ave = curveA.rate;
                        n.frMap_rew = curveR.rate;
                        vHPC{ii} = n; 
                        %Store firing maps
                        pc_vHPC_par.firingMap.ave = [pc_vHPC_par.firingMap.ave;curveA_vhpc]; 
                        pc_vHPC_par.firingMap.rew = [pc_vHPC_par.firingMap.rew;curveR_vhpc];
                        
                        %Sparsity: 
                        sa = sparsity_info(curveA.rate,curveA.time);
                        sr = sparsity_info(curveR.rate,curveR.time);
                        %Store pc parametres 
                        pc_vHPC_par.pc_params.ave = [pc_vHPC_par.pc_params.ave;statsA.specificity,sa,statsA.size(1),statsA.peak(1)]; 
                        pc_vHPC_par.pc_params.rew = [pc_vHPC_par.pc_params.rew;statsR.specificity,sr,statsR.size(1),statsR.peak(1)];
                        
                        % Aversive 
                        [ff,f] = max(curveA.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        fr_ave = nanmean(curveA.rate);
                        PC.vHPC.aversive=[PC.vHPC.aversive ; cluster , ff , fff(f),fr_ave];
                        clear f ff fff
                        
                        maps_vHPC_A = [maps_vHPC_A ; curveA.rate];
%                         maps_dHPC_SA = [maps_dHPC_SA ; curveSA.rate];
                        map_ave = curveA.rate;
                        peak_vHPC_A = [peak_vHPC_A ; statsA.peak(1)];
                        field_vHPC_A = [field_vHPC_A ; statsA.fieldX(1,:)];
                        
                        %Rewarded 
                        [ff,f] = max(curveR.rate);
                        fff = [1/Xedges:1/Xedges:1];
                        fr_rew = nanmean(curveR.rate);
                        PC.vHPC.reward=[PC.vHPC.aversive ; cluster , ff , fff(f),fr_rew];
                        clear f ff fff
                        
                        maps_vHPC_R = [maps_vHPC_R ; curveR.rate];
%                         maps_dHPC_SR = [maps_dHPC_SR ; curveSR.rate];
                        map_rew = curveR.rate;
                        peak_vHPC_R = [peak_vHPC_R ; statsR.peak(1)] ;
                        field_vHPC_R = [field_vHPC_R ; statsR.fieldX(1,:)];
                        cluster_dHPC = [cluster_dHPC ; cluster , cond];
                        
                        %Store within 
                        pc_vHPC_par.within.ave = [pc_vHPC_par.within.ave; withinA];
                        pc_vHPC_par.within.rew = [pc_vHPC_par.within.rew; withinR];
                        
                        pc_vHPC_par.within.ave_tresh = [pc_vHPC_par.within.ave_tresh; withinA_tresh];
                        pc_vHPC_par.within.rew_tresh = [pc_vHPC_par.within.rew_tresh; withinR_tresh];
                        
                        %Between-trial pc parameters (aversive vs. rewarded
                        [between] = Between_pc(pos_ave,spks_ave,pos_rew,spks_rew,bin_size,sigma,Xedges);
                        
                        %Save between pc parameters 
                        pc_vHPC_par.between = [pc_vHPC_par.between; between];
                        
                        clear fr_ave fr_rew map_ave map_rew
                        
                        
                    end
                    clear s q curve OccMap Nspikes spks_tmp pos_tmp m curve stats sp curve1
                    clear curveA curveR qA qR curveSA curveSR statsSA statsSR
                end
            end
            clear celltype tmp b cluster
            
    end
        

    % Save output session
    if ~isempty(dHPC)
        dHPC = dHPC(~cellfun('isempty',dHPC));
        save([cd,'\dHPC_pc.mat'],'dHPC'); 
    end
    if ~isempty(vHPC)
        vHPC = vHPC(~cellfun('isempty',vHPC));
        save([cd,'\vHPC_pc.mat'],'vHPC'); 
    end
  
            %% Saving data
            disp('Saving outputs')
            if not(exist('pc_all','var'))
                maps.dHPC.aversive = maps_dHPC_A;
                maps.vHPC.aversive = maps_vHPC_A;
                maps.dHPC.reward   = maps_dHPC_R;
                maps.vHPC.reward   = maps_vHPC_R;
                maps.dHPC.aversiveS = maps_dHPC_SA;
                maps.vHPC.aversiveS = maps_vHPC_SA;
                maps.dHPC.rewardS   = maps_dHPC_SR;
                maps.vHPC.rewardS   = maps_vHPC_SR;
                peak.dHPC.aversive = peak_dHPC_A;
                peak.vHPC.aversive = peak_vHPC_A;
                peak.dHPC.reward   = peak_dHPC_R;
                peak.vHPC.reward   = peak_vHPC_R;
                field.dHPC.aversive = field_dHPC_A;
                field.vHPC.aversive = field_vHPC_A;
                field.dHPC.reward   = field_dHPC_R;
                field.vHPC.reward   = field_vHPC_R;
                clus.dHPC = cluster_dHPC;
                clus.vHPC = cluster_vHPC;
                
                pc_all.dHPC.within.A =  pc_dHPC_par.within.ave;
                pc_all.dHPC.within.R =  pc_dHPC_par.within.rew;
                pc_all.dHPC.within.Atresh = pc_dHPC_par.within.ave_tresh;
                pc_all.dHPC.within.Rtresh = pc_dHPC_par.within.rew_tresh;
                pc_all.dHPC.between =   pc_dHPC_par.between; 
                pc_all.dHPC.firingMap.ave = pc_dHPC_par.firingMap.ave; 
                pc_all.dHPC.firingMap.rew = pc_dHPC_par.firingMap.rew;
                pc_all.dHPC.pc_params.ave=pc_dHPC_par.pc_params.ave; 
                pc_all.dHPC.pc_params.rew=pc_dHPC_par.pc_params.rew;
                
                pc_all.vHPC.within.A =  pc_vHPC_par.within.ave;
                pc_all.vHPC.within.R =  pc_vHPC_par.within.rew;
                pc_all.vHPC.within.Atresh = pc_vHPC_par.within.ave_tresh;
                pc_all.vHPC.within.Rtresh = pc_vHPC_par.within.rew_tresh;
                pc_all.vHPC.between =  pc_vHPC_par.between; 
                pc_all.vHPC.firingMap.ave = pc_vHPC_par.firingMap.ave; 
                pc_all.vHPC.firingMap.rew = pc_vHPC_par.firingMap.rew;
                pc_all.vHPC.pc_params.ave=pc_vHPC_par.pc_params.ave; 
                pc_all.vHPC.pc_params.rew=pc_vHPC_par.pc_params.rew;
                
            else
                maps.dHPC.aversive = [maps.dHPC.aversive ; maps_dHPC_A];
                maps.vHPC.aversive = [maps.vHPC.aversive ; maps_vHPC_A];
                maps.dHPC.reward   = [maps.dHPC.reward ; maps_dHPC_R];
                maps.vHPC.reward   = [maps.vHPC.reward ; maps_vHPC_R];
                maps.dHPC.aversiveS = [maps.dHPC.aversiveS ; maps_dHPC_SA];
                maps.vHPC.aversiveS = [maps.vHPC.aversiveS ; maps_vHPC_SA];
                maps.dHPC.rewardS   = [maps.dHPC.rewardS ; maps_dHPC_SR];
                maps.vHPC.rewardS   = [maps.vHPC.rewardS ; maps_vHPC_SR];                
                peak.dHPC.aversive = [peak.dHPC.aversive ; peak_dHPC_A];
                peak.vHPC.aversive = [peak.vHPC.aversive ; peak_vHPC_A];
                peak.dHPC.reward   = [peak.dHPC.reward ; peak_dHPC_R];
                peak.vHPC.reward   = [peak.vHPC.reward ; peak_vHPC_R];
                field.dHPC.aversive = [field.dHPC.aversive ; field_dHPC_A];
                field.vHPC.aversive = [field.vHPC.aversive ; field_vHPC_A];
                field.dHPC.reward   = [field.dHPC.reward ; field_dHPC_R];
                field.vHPC.reward   = [field.vHPC.reward ; field_vHPC_R];
                clus.dHPC = [clus.dHPC ; cluster_dHPC];
                clus.vHPC = [clus.vHPC ; cluster_vHPC]; 
                
                pc_all.dHPC.within.A = [pc_all.dHPC.within.A ; pc_dHPC_par.within.ave];
                pc_all.dHPC.within.R  = [pc_all.dHPC.within.R ; pc_dHPC_par.within.rew];
                pc_all.dHPC.between = [ pc_all.dHPC.between;  pc_dHPC_par.between]; 
                pc_all.dHPC.within.Atresh = [pc_all.dHPC.within.Atresh; pc_dHPC_par.within.ave_tresh];
                pc_all.dHPC.within.Rtresh = [pc_all.dHPC.within.Rtresh; pc_dHPC_par.within.rew_tresh];
                pc_all.dHPC.firingMap.ave = [ pc_all.dHPC.firingMap.ave;pc_dHPC_par.firingMap.ave]; 
                pc_all.dHPC.firingMap.rew = [pc_all.dHPC.firingMap.rew;pc_dHPC_par.firingMap.rew];
                pc_all.dHPC.pc_params.ave=[pc_all.dHPC.pc_params.ave;pc_dHPC_par.pc_params.ave]; 
                pc_all.dHPC.pc_params.rew=[pc_all.dHPC.pc_params.rew;pc_dHPC_par.pc_params.rew];
                
                pc_all.vHPC.within.A = [pc_all.vHPC.within.A ; pc_vHPC_par.within.ave];
                pc_all.vHPC.within.R  = [pc_all.vHPC.within.R ; pc_vHPC_par.within.rew];
                pc_all.vHPC.between = [ pc_all.vHPC.between;  pc_vHPC_par.between];
                pc_all.vHPC.within.Atresh = [pc_all.vHPC.within.Atresh; pc_vHPC_par.within.ave_tresh];
                pc_all.vHPC.within.Rtresh = [pc_all.vHPC.within.Rtresh; pc_vHPC_par.within.rew_tresh];
                pc_all.vHPC.firingMap.ave = [pc_all.vHPC.firingMap.ave;pc_vHPC_par.firingMap.ave]; 
                pc_all.vHPC.firingMap.rew = [pc_all.vHPC.firingMap.rew;pc_vHPC_par.firingMap.rew];
                pc_all.vHPC.pc_params.ave=[pc_all.vHPC.pc_params.ave;pc_vHPC_par.pc_params.ave]; 
                pc_all.vHPC.pc_params.rew=[pc_all.vHPC.pc_params.rew;pc_vHPC_par.pc_params.rew];
                
            end
        clear specificity_dHPC_A specificity_dHPC_R specificity_vHPC_A specificity_vHPC_R
        clear field_dHPC_A field_dHPC_R field_vHPC_A field_vHPC_R
        clear peak_dHPC_A peak_dHPC_R peak_vHPC_A peak_vHPC_R
        clear maps_dHPC_A maps_dHPC_R maps_vHPC_A maps_vHPC_R
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear behavior camara Cell_type_classification cellulartype cluster
        clear group_dHPC group_vHPC K Kinfo leftvalve rightvalve movement R
        clear Rewards_filt Shocks_filt spks_dHPC spks_vHPC r rZ
        clear cluster cluster_dHPC cluster_vHPC coordinated coordinatedV coordinatedV_refined
        clear camaraA count dX dX_int dY dY_int gc grps i ii MUA segments tmp WAKE REM NREM
        clear ripple_bursts ripplesD ripplesV x w position_shocks posx posy ejeX ejeY PC replay Replay limits
        disp(['-- Finished folder #' , num2str(t) , ' from rat #' , num2str(tt) , ' ---'])
        disp('  ')
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

save('W:\Remapping-analysis-Facu\pc_all_within_between_v3.mat', 'pc_all');

%% Total number of recorded neurons - next time put this inseide the main loop
% c1: rat c2: session c3: #dHPC neurons c4: dHPC #pyr c5: #vHPCneurons  c6:VHPC#pyr 
n_total = nan(2000,6); 

ind_global = 1; 

for tt = 1:length(path)
    %List of folders from the path
    files = dir(path{tt});
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags);
    
    for t = 1 : length(subFolders)-2
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)
        rat = str2num(files(3).name(4:6)); 
        %Output
        n_total(ind_global,1) = rat; % rat
        n_total(ind_global,2) = str2num(session(end-7:end)); %session
        
        %% Spikes

        cd 'Spikesorting'
        spks = double([readNPY('spike_clusters.npy') readNPY('spike_times.npy')]);
        K = tdfread('cluster_group.tsv'); % import clusters ID and groups (MUA,Noise,Good)
        Kinfo = tdfread('cluster_info.tsv'); % import information of clusters
        K = [K.cluster_id(K.group(:,1) == 'g') , Kinfo.ch(K.group(:,1) == 'g'),]; % defining good clusters

        % Load neuronal classification
        load('Cell_type_classification')
        K = [K , Cell_type_classification(:,6:7)];
        group_dHPC = K(K(:,2) > 63,:);
        group_vHPC = K(K(:,2) <= 63,:);
      
        %Output
        n_total(ind_global,3) = nansum([n_total(ind_global,3),size(group_dHPC,1)]); %# neurons per session dhpc
        n_total(ind_global,4) = nansum([n_total(ind_global,4),sum(group_dHPC(:,3))]); %# pyr per session dhpc
        
        n_total(ind_global,5) = nansum([n_total(ind_global,5),size(group_vHPC,1)]); %# neurons per session vhpc
        n_total(ind_global,6) = nansum([n_total(ind_global,6),sum(group_vHPC(:,3))]); %# pyr per session vhpc 
        
       ind_global = ind_global+1;
          
    end
    disp(['-------- Finished rat#' , num2str(tt) , ' --------'])
    disp(' ')
end

n_total(any(isnan(n_total), 2), :) = [];

% Count - change manually
size(unique(n_total(n_total(:,1)==165,2)))

sum(n_total(n_total(:,1)==165,6))
%% Remapping stats 
%dHPC 
% 1 = between 2 = within aversive 3= within reward
% c1 = spatial c2= fr_change c3=  overlap
data = [pc_all.dHPC.between,ones(size(pc_all.dHPC.between,1),1);pc_all.dHPC.within.A,...
ones(size(pc_all.dHPC.within.A,1),1)*2;pc_all.dHPC.within.R, ones(size(pc_all.dHPC.within.R,1),1)*3];

data = [pc_all.vHPC.between,ones(size(pc_all.vHPC.between,1),1);pc_all.vHPC.within.A,...
ones(size(pc_all.vHPC.within.A,1),1)*2;pc_all.vHPC.within.R, ones(size(pc_all.vHPC.within.R,1),1)*3];

%Stats
[P,ANOVATAB,STATS] = kruskalwallis(data(:,2),data(:,4));
c = multcompare(STATS)

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%Boxplot
id= 3; % 1 = spatial 2= fr_change 3=  overlap
data_box = [pc_all.dHPC.between(:,id),pc_all.dHPC.within.A(:,id),pc_all.dHPC.within.R(:,id)];

figure(1);clf;
boxplot(data_box, 'Labels',{'Between', 'Within Ave', 'Within Rew'})
ylabel('Fr change','FontSize',14);
title('dHPC'); 

%Violin plot
figure(1);clf;
label = {'Between', 'Within Ave', 'Within Rew'}; 
[h,L,MX,MED,bw] = violin(data_box, 'xlabel', label,'facecolor',[0 0 0;1 0 0; 0 0 1]);
ylabel('Fr change','FontSize',14);
title('dHPC');

%% PC plot

tempA=[];
tempR=[];

for i =1:size(pc_all.dHPC.firingMap.ave,1)
    A = pc_all.dHPC.firingMap.ave(i,:) - min(pc_all.dHPC.firingMap.ave(i,:));
    A = A ./ max(A);
    tempA=[tempA ; A];
    clear A
    
    R = pc_all.dHPC.firingMap.rew(i,:) - min(pc_all.dHPC.firingMap.rew(i,:));
    R = R ./ max(R);
    tempR=[tempR ; R];
    clear R    
end

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 
figure(1);clf;hold on; 
subplot(1,2,1);imagesc([3:3:180], [1:1:size(pc_all.dHPC.firingMap.ave,1)],tempA(mm,:)), colormap 'jet'
subplot(1,2,2);imagesc([3:3:180], [1:1:size(pc_all.dHPC.firingMap.rew,1)],tempR(mm,:)), colormap 'jet'
sgtitle('dHPC');

tempA=[];
tempR=[];

for i =1:size(pc_all.vHPC.firingMap.ave,1)
    A = pc_all.vHPC.firingMap.ave(i,:) - min(pc_all.vHPC.firingMap.ave(i,:));
    A = A ./ max(A);
    tempA=[tempA ; A];
    clear A
    
    R = pc_all.vHPC.firingMap.rew(i,:) - min(pc_all.vHPC.firingMap.rew(i,:));
    R = R ./ max(R);
    tempR=[tempR ; R];
    clear R    
end

[h idx] = max (tempA, [],2);
[m mm] = sort(idx); 
figure(2);clf;
subplot(1,2,1); imagesc([3:3:180], [1:1:size(pc_all.vHPC.firingMap.ave,1)],tempA(mm,:)), caxis([0 1]),colormap 'jet'
subplot(1,2,2); imagesc([3:3:180], [1:1:size(pc_all.vHPC.firingMap.rew,1)],tempR(mm,:)), caxis([0 1]), colormap 'jet'


%% Aversive  --> Reward  PREGUNTAR  a facu porque cambio de sentido 
criteria = clus.dHPC(:,2) == 1;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
        ca = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.rewardS(i,:));
        cr = corrcoef(maps.dHPC.reward(i,:) , maps.dHPC.aversiveS(i,:));
        SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)/(ca(1,2) + cr(1,2))]; clear c ca cr
    end
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))


criteria = clus.vHPC(:,2) == 1;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        
        c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
        ca = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.rewardS(i,:));
        cr = corrcoef(maps.vHPC.reward(i,:) , maps.vHPC.aversiveS(i,:));
        SpatialCorr.vHPC = [SpatialCorr.dHPC ; c(1,2)/(ca(1,2) + cr(1,2))]; clear c ca cr
    end
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y]),ylim([-3 3])

%% Reward --> Aversive
criteria = clus.dHPC(:,2) == 2;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
        SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)]; clear c
    end
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))


criteria = clus.vHPC(:,2) == 2;
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    if logical(criteria(i))
        ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MA = [MA ; iii];
        tmpA = [tmpA ; ii]; clear ii trash iii
        
        
        ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
        ii = ii ./ max(ii);
        [trash , iii] = max(ii);
        MR = [MR ; iii];
        tmpR = [tmpR ; ii]; clear ii trash iii
        
        c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
        SpatialCorr.vHPC = [SpatialCorr.vHPC ; c(1,2)]; clear c
    end
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y])
[h p] = ttest2(SpatialCorr.dHPC , SpatialCorr.vHPC);

%% All
tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.dHPC = [];
for i = 1 : size(maps.dHPC.aversive,1)
    ii = maps.dHPC.aversive(i,:) - min(maps.dHPC.aversive(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MA = [MA ; iii];
    tmpA = [tmpA ; ii]; clear ii trash iii
    
    
    ii = maps.dHPC.reward(i,:) - min(maps.dHPC.reward(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MR = [MR ; iii];
    tmpR = [tmpR ; ii]; clear ii trash iii  
    c = corrcoef(maps.dHPC.aversive(i,:) , maps.dHPC.reward(i,:));
    SpatialCorr.dHPC = [SpatialCorr.dHPC ; c(1,2)]; clear c
end
clear i
[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

tmpA = [];
MA = [];
tmpR = [];
MR = [];
SpatialCorr.vHPC = [];
for i = 1 : size(maps.vHPC.aversive,1)
    ii = maps.vHPC.aversive(i,:) - min(maps.vHPC.aversive(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MA = [MA ; iii];
    tmpA = [tmpA ; ii]; clear ii trash iii
    
    
    ii = maps.vHPC.reward(i,:) - min(maps.vHPC.reward(i,:));
    ii = ii ./ max(ii);
    [trash , iii] = max(ii);
    MR = [MR ; iii];
    tmpR = [tmpR ; ii]; clear ii trash iii   
    
    c = corrcoef(maps.vHPC.aversive(i,:) , maps.vHPC.reward(i,:));
    SpatialCorr.vHPC = [SpatialCorr.vHPC ; c(1,2)]; clear c
end
clear i

[i ii] = sort(MA);

figure,
subplot(121),imagesc(tmpA(ii,:))
subplot(122),imagesc(tmpR(ii,:))

figure,
x = ones(size(SpatialCorr.dHPC));
y = ones(size(SpatialCorr.vHPC))*2;
boxplot([SpatialCorr.dHPC ; SpatialCorr.vHPC],[x;y])
%%
parfor i=1:10
   disp('Hola') 

end   
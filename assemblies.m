clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\ready';'E:\Rat103\usable'};%List of folders from the path

% for SU
criteria_fr = 0.01; %criteria to include or not a SU into the analysis
criteria_n = 6; % minimal number of neurons from each structure
criteria_type = 0; %criteria for celltype (0:pyr, 1:int, 2:all)
n_SU_V = [];
n_SU_D = [];
Xedges = 60; %number of bins for RateMap construction

binSize = 0.001; % bin size for replay events detection

% Behavior
minimal_speed = 5; % minimal speed to detect active periods
minimal_speed_time = 2; % minimal time to detect active periods

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
        session = [subFolders(t+2).folder,'\',subFolders(t+2).name];
        cd(session)

        %Loading TS of the sessions
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
        %% Awake
        % Load digitalin:mat
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
            behavior.speed.reward = LinearVelocity(pos,15);
            behavior.quiet.reward = QuietPeriods( behavior.speed.reward , 5 , 2);
            clear pos camaraR posx posy
            
            load('laps2.mat','posx','posy');
            [camaraA,~] = find((camara(:,1)-aversiveTS_run(1)/1000)>0,1,'first'); %TimeStamp of the begining of aversive
            pos = [camara(camaraA : camaraA+length(posx)-1),posx,posy];
            
            %interpolation of dropped frames
            ejeX = pos(~isnan(pos(:,2)),1); dX = pos(~isnan(pos(:,2)),2); dX_int = interp1(ejeX , dX , pos(:,1));
            ejeY = pos(~isnan(pos(:,3)),1); dY = pos(~isnan(pos(:,3)),3); dY_int = interp1(ejeY , dY , pos(:,1));
            pos(:,2) =dX_int; pos(:,3) =dY_int;%saving corrected pos
            
            behavior.pos.aversive = [pos];
            behavior.speed.aversive = LinearVelocity(pos,15);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , 5 , 2);
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
            behavior.speed.reward = LinearVelocity(pos,15);
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
            behavior.speed.aversive = LinearVelocity(pos,15);
            behavior.quiet.aversive = QuietPeriods(behavior.speed.aversive , minimal_speed , minimal_speed_time);
            clear pos camaraA posx posy
        end
        
        % Generation of no-movements periods
        start = behavior.speed.reward(1,1);   stop = behavior.speed.reward(end,1);
        tmp = InvertIntervals(behavior.quiet.reward , start , stop); %keep only those higher than criteria
        [movement.reward , ~] =  ExcludeIntervals(tmp ,[Rewards_filt Rewards_filt+4]);
        clear tmp start stop
        
        start = behavior.speed.aversive(1,1);   stop = behavior.speed.aversive(end,1);
        tmp = InvertIntervals(behavior.quiet.aversive , start , stop);%keep only those higher than criteria
        [movement.aversive , ~] =  ExcludeIntervals(tmp ,[Shocks_filt-0.5 Shocks_filt+2]);
        clear tmp start stop
        
        % Definition of the place where the shocks were admninistrated
        position_shocks = [];
        for i = 1 : length(Shocks_filt)
            [~ , ii] = min(abs(behavior.pos.aversive(:,1)-Shocks_filt(i)));
            position_shocks = [position_shocks ; behavior.pos.aversive(ii,2)]; 
            clear ii
        end
        clear i
        position_shocks = position_shocks-min(behavior.pos.aversive(:,2));
        position_shocks = position_shocks./max(behavior.pos.aversive(:,2));

        %% load coordinated ripple bursts
        load('coordinated_ripple_bursts.mat')
        ripple_bursts.baseline = Restrict(coordinated_ripple_bursts,baselineTS./1000);
        ripple_bursts.reward = Restrict(coordinated_ripple_bursts,rewardTS./1000);
        ripple_bursts.aversive = Restrict(coordinated_ripple_bursts,aversiveTS./1000);
        clear coordinated_ripple_bursts
        
        %% Load ripples
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        % coordination
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        %% Spikes
        %Load Units
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

        %% Neuronal Ensembles detection
        % Ensembles during shocks
        % vHPC
        limits = aversiveTS_run./1000;
        events = [];
%         events = movement.aversive;
        [SpksTrains.dHPC.aversive , Bins.aversive , Cluster.dHPC.aversive] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, 0.1, limits, events, true);
        [SpksTrains.vHPC.aversive , Bins.aversive , Cluster.vHPC.aversive] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, 0.1, limits, events, true);
        C = SpksTrains.dHPC.aversive' * SpksTrains.dHPC.aversive ./ size(SpksTrains.dHPC.aversive,1);
%         C = C - diag(diag(C));
%         [C , P]=corr(SpksTrains,SpksTrains);
        
        [EigenValues] = eig(C);
        lambda = marchenko(SpksTrains.dHPC.aversive);
        
        e = EigenValues > lambda(2);% + lambda(2)^(-2/3);
        Coeff = pca(SpksTrains.dHPC.aversive);
%         P = Coeff(e,:);
        reactivation.dHPC = [];
        
       limits = [0 segments.Var1(end)/1000];
        events = [];
%         events = movement.aversive;
        [SpksTrains.dHPC.all , Bins.all , Cluster.dHPC.all] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, 0.1, limits, events, true);
        for i = 1 : length(e)
            if e(i)
                N = (Coeff(:,i)*Coeff(:,i)');%*C(:,i);
                figure,imagesc(N)
                 reactivation.dHPC = [reactivation.dHPC , SpksTrains.dHPC.all .* N .* SpksTrains.dHPC.all'];

            end
        end
        
        reactivation
        
        [SpksTrains , Bins , Cluster] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, 0.1, [0 segments.Var1(end)/1000], [], false);

        SpksTrains = zscore(SpksTrains ./ 0.1);
        
        p = P(1,:).'*P(1,:);
        
        ensembles.vHPC.aversive = [];
        for i = 1 : length(e)
            if logical(e(i))
                tmp = real(Coeff(i,:));
                criteria = lambda(2) + lambda(2)^(-2/3);
                if sum(or(tmp>criteria,tmp<criteria*-1))>1
                    figure,stem(tmp./lambda(2),'filled'),hold on
                    yline(criteria,'--')
                    yline(criteria*-1,'--')
                    c = Cluster(or(tmp>criteria,tmp<criteria*-1));
                    ensembles.vHPC.aversive =[ensembles.vHPC.aversive ; {c}];
                end
            end
        end
        clear i Coeff C lambda EigenValues limits events SpksTrains Bins Cluster
        
        % dHPC
        limits = aversiveTS_run./1000;
        events = [Shocks_filt Shocks_filt+2];
%         events = movement.aversive;
        [SpksTrains , Bins , Cluster] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, 0.1, limits, events, true);
        SpksTrains = SpksTrains';
        C = SpksTrains * SpksTrains' / size(SpksTrains,2);
%         [C , P]=corr(SpksTrains,SpksTrains);
        
        [EigenValues] = eig(C);
        lambda = marchenko(SpksTrains);
        
        e = EigenValues > lambda(2) + lambda(2)^(-2/3);
        
        Coeff = pca(SpksTrains,'eig');
        ensembles.dHPC.aversive = [];
        for i = 1 : length(e)
            if logical(e(i))
                tmp = real(Coeff(i,:));
                criteria = lambda(2) + lambda(2)^(-2/3);
                if sum(or(tmp>criteria,tmp<criteria*-1))>1
                    figure,stem(tmp./lambda(2),'filled'),hold on
                    yline(criteria,'--')
                    yline(criteria*-1,'--')
                    c = Cluster(or(tmp>criteria,tmp<criteria*-1));
                    ensembles.dHPC.aversive =[ensembles.dHPC.aversive ; {c}];
                end
            end
        end
        clear i Coeff C lambda EigenValues limits events SpksTrains Bins Cluster

        
        % Ensembles during valve
        % vHPC
        limits = rewardTS_run./1000;
        events = [Rewards_filt(:,1) Rewards_filt(:,2)+2];
%         events = movement.reward;
        [SpksTrains , Bins , Cluster] = spike_train_construction(spks_vHPC, group_vHPC(:,1), cellulartype, 0.1, limits, events, true);
        SpksTrains = SpksTrains';
        C = SpksTrains * SpksTrains' / size(SpksTrains,2);
%         [C , P]=corr(SpksTrains,SpksTrains);
        
        [EigenValues] = eig(C);
        lambda = marchenko(SpksTrains);
        
        e = EigenValues > lambda(2);
        
        Coeff = pca(SpksTrains,'eig');
        ensembles.vHPC.reward = [];
        for i = 1 : length(e)
            if logical(e(i))
                tmp = real(Coeff(i,:));
                 criteria = lambda(2) + lambda(2)^(-2/3);
                if sum(or(tmp>criteria,tmp<criteria*-1))>1
                    figure,stem(tmp./lambda(2),'filled'),hold on
                    yline(criteria,'--')
                    yline(criteria*-1,'--')
                    c = Cluster(or(tmp>criteria,tmp<criteria*-1));
                    ensembles.vHPC.reward =[ensembles.vHPC.reward ; {c}];
                end
            end
        end
        clear i Coeff C lambda EigenValues limits events SpksTrains Bins Cluster
        
        % dHPC
        limits = rewardTS_run./1000;
        events = [Rewards_filt(:,1) Rewards_filt(:,2)+2];
%         events = movement.reward;
        [SpksTrains , Bins , Cluster] = spike_train_construction(spks_dHPC, group_dHPC(:,1), cellulartype, 0.1, limits, events, true);
        SpksTrains = SpksTrains';
        C = SpksTrains * SpksTrains' / size(SpksTrains,2);
%         [C , P]=corr(SpksTrains,SpksTrains);
        
        [EigenValues] = eig(C);
        lambda = marchenko(SpksTrains);
        
        e = EigenValues > lambda(2) + lambda(2)^(-2/3);
        
        Coeff = pca(SpksTrains,'eig');
        ensembles.dHPC.reward = [];
        for i = 1 : length(e)
            if logical(e(i))
                tmp = real(Coeff(i,:));
                 criteria = lambda(2) + lambda(2)^(-2/3);
                if sum(or(tmp>criteria,tmp<criteria*-1))>1
                    figure,stem(tmp./lambda(2),'filled'),hold on
                    yline(criteria,'--')
                    yline(criteria*-1,'--')
                    c = Cluster(or(tmp>criteria,tmp<criteria*-1));
                    ensembles.dHPC.reward =[ensembles.dHPC.reward ; {c}];
                end
            end
        end
        clear i Coeff C lambda EigenValues limits events SpksTrains Bins Cluster
        
        %% Firing maps calculation
        PHIST.vHPC.aversive = [];        PHIST.vHPC.reward = [];        
        for ii=1:length(ensembles.vHPC.aversive)
            spks = spks_vHPC(ismember(spks_vHPC(:,1) , ensembles.vHPC.aversive{ii}),2);
            % --- Aversive ---
            x = Restrict(spks , [Shocks_filt(:,1)-2 Shocks_filt(:,1)+2]);
            y = [Shocks_filt(:,1)];
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector1] = CCG(times,groups,'binsize',0.1,'duration',4,'smooth',1);
            ccg = ccg(:,1,2)./length(y)./0.05;
            PHIST.vHPC.aversive = [PHIST.vHPC.aversive ; ccg'];
            clear ccg x y celltype tmp b cluster
        end
        
        for ii=1:length(ensembles.vHPC.reward)
            spks = spks_vHPC(ismember(spks_vHPC(:,1) , ensembles.vHPC.reward{ii}),2);
            % --- Reward ---
            x = Restrict(spks , [Rewards_filt(:,1)-2 Rewards_filt(:,1)+2]);
            y = [Rewards_filt(:,1)];
            [times,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,TimeVector1] = CCG(times,groups,'binsize',0.1,'duration',4,'smooth',1);
            ccg = ccg(:,1,2)./length(y)./0.05;
            PHIST.vHPC.reward = [PHIST.vHPC.reward ; ccg'];
            clear ccg x y celltype tmp b cluster
        end        
        
        %% Replay detection
    freq = 1/binSize;
    limits = [0 segments.Var1(end)/1000];
        spks = [];
        %Replay events in dHPC
        for ii=1:length(group_dHPC)
            cluster = group_dHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = [spks;spks_dHPC(spks_dHPC(:,1)==cluster,2)];
            end
            clear celltype
        end
        [MUA.dHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks 
        w = gausswin(30, 5);
        MUA.dHPC = filter(w , 1 ,MUA.dHPC);
        replay.dHPC = MUA.dHPC>mean(MUA.dHPC)+std(MUA.dHPC)*3;
        replay.dHPC = ToIntervals(bins',replay.dHPC);
        replay.dHPC = replay.dHPC(replay.dHPC(:,2)-replay.dHPC(:,1)>0.04,:);    
        
        %filtering replay events by amount of PCs
%         clusterdHPC = unique([PC.dHPC.aversive(:,1);PC.dHPC.reward(:,1)]);
        count = [];
        for i = 1 :  length(cluster_dHPC)
            cluster = cluster_dHPC(i,1);
            spks = spks_dHPC(spks_dHPC(:,1)==cluster,2);
            [status,interval,index] = InIntervals(spks,replay.dHPC);
            interval = unique(interval);
            count = [count ; interval(interval~=0)];
            clear spks cluster status interval index
        end
        [gc,grps] = groupcounts(count);
        replay.dHPC = replay.dHPC(grps(gc>length(cluster_dHPC)*0.15),:);
        
        
        %Replay events in vHPC
        spks = [];
        for ii=1:length(group_vHPC)
            cluster = group_vHPC(ii,1);
            celltype = logical(cellulartype(cellulartype(:,1) == cluster,2));
            if celltype
                spks = [spks;spks_vHPC(spks_vHPC(:,1)==cluster,2)];
            end
            clear celltype
        end
        [MUA.vHPC,bins]=binspikes(sort(spks,'ascend'),freq,limits);clear spks 
        w = gausswin(30, 1);
        MUA.vHPC = filter(w , 1 ,MUA.vHPC); 
        replay.vHPC = MUA.vHPC>mean(MUA.vHPC)+std(MUA.vHPC)*2;
        replay.vHPC = ToIntervals(bins',replay.vHPC);
        replay.vHPC = replay.vHPC(replay.vHPC(:,2)-replay.vHPC(:,1)>0.04,:);

        %filtering replay events by amount of PCs
%         clustervHPC = unique([PC.vHPC.aversive(:,1);PC.vHPC.reward(:,1)]);
        count = [];
        for i = 1 :  length(cluster_vHPC)
            cluster = cluster_vHPC(i,1);
            spks = spks_vHPC(spks_vHPC(:,1)==cluster,2);
            [status,interval,index] = InIntervals(spks,replay.vHPC);
            interval = unique(interval);
            count = [count ; interval(interval~=0)];
            clear spks cluster status interval index
        end
        [gc,grps] = groupcounts(count);
        replay.vHPC = replay.vHPC(grps(gc>length(cluster_vHPC)*0.15),:);
        
        
        %% Bayesian Decoding across conditions
        % Position
        decoded_positions.dHPC.aversive = [];
        decoded_positions.vHPC.aversive = [];
        decoded_positions.dHPC.reward = [];
        decoded_positions.vHPC.reward = [];
        
        % Aversive
      % Aversive
        tmp.vHPC = Restrict(Restrict(replay.vHPC , aversiveTS./1000),[ripple_bursts.aversive(:,1),ripple_bursts.aversive(:,3)]);
        tmp.dHPC = Restrict(Restrict(replay.dHPC , aversiveTS./1000),[ripple_bursts.aversive(:,1),ripple_bursts.aversive(:,3)]);
        
        figure, subplot(221)
        for i = 1 : length(tmp.dHPC)
            % Decoding using dHPC SUs
            start = tmp.dHPC(i,1);
            stop = tmp.dHPC(i,2);
            nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.aversive, start, stop);
            probability = bayesian_replay(PHIST.dHPC.normalizado.aversive, nSpks, start, stop);
            plot([1/60:1/60:1],probability,'r'),hold on
            decoded_positions.dHPC.aversive = [decoded_positions.dHPC.aversive ; probability];
            clear probability nSpks start stop
        end
        title('Decoding - dHPC - Aversive')
        
        subplot(223)
        for i = 1 : length(tmp.vHPC)
            % Decoding using vHPC SUs
            start = tmp.vHPC(i,1);
            stop = tmp.vHPC(i,2);
            nSpks = count_spks(spks_vHPC, clustervHPC.normalizado.aversive, start, stop);
            probability = bayesian_replay(PHIST.vHPC.normalizado.aversive, nSpks, start, stop);
            decoded_positions.vHPC.aversive = [decoded_positions.vHPC.aversive ; probability];
            plot([1/60:1/60:1],probability,'b'),hold on
            clear probability nSpks start stop
        end
        title('Decoding - vHPC - Aversive')

        % Reward
      % Reward
        tmp.vHPC = Restrict(Restrict(replay.vHPC , rewardTS./1000),[ripple_bursts.reward(:,1),ripple_bursts.reward(:,3)]);
        tmp.dHPC = Restrict(Restrict(replay.dHPC , rewardTS./1000),[ripple_bursts.reward(:,1),ripple_bursts.reward(:,3)]);

        subplot(222)
        for i = 1 : length(tmp.dHPC)
            % Decoding using dHPC SUs
            start = tmp.dHPC(i,1);
            stop = tmp.dHPC(i,2);
            nSpks = count_spks(spks_dHPC, clusterdHPC.normalizado.reward, start, stop);
            probability = bayesian_replay(PHIST.dHPC.normalizado.reward, nSpks, start, stop);
            plot([1/60:1/60:1],probability,'r'),hold on
            decoded_positions.dHPC.reward = [decoded_positions.dHPC.reward ; probability];
            clear probability nSpks start stop
        end
        title('Decoding - dHPC - Reward')

        subplot(224)
        for i = 1 : length(tmp.vHPC)
            % Decoding using vHPC SUs
            start = tmp.vHPC(i,1);
            stop = tmp.vHPC(i,2);
            nSpks = count_spks(spks_vHPC, clustervHPC.normalizado.reward, start, stop);
            probability = bayesian_replay(PHIST.vHPC.normalizado.reward, nSpks, start, stop);
            decoded_positions.vHPC.reward = [decoded_positions.vHPC.reward ; probability];
            plot([1/60:1/60:1],probability,'b'),hold on
            clear probability nSpks start stop
        end
        title('Decoding - vHPC - Reward')
        
        
        %% Bayesian Decoding across conditions
        % Events
        decoded_events.dHPC.aversive = [];
        decoded_events.vHPC.aversive = [];
        decoded_events.dHPC.reward = [];
        decoded_events.vHPC.reward = [];
        
        % Aversive
        tmp.vHPC = Restrict(Restrict(replay.vHPC , aversiveTS./1000),[coordinatedV_refined(:,2)-0.2,coordinatedV_refined(:,2)+0.2]);
        tmp.dHPC = Restrict(Restrict(replay.dHPC , aversiveTS./1000),[coordinated(:,2)-0.2,coordinated(:,2)+0.2]);
        figure, subplot(221)
        for i = 1 : length(tmp.dHPC)
            % Decoding using dHPC SUs
            start = tmp.dHPC(i,1);
            stop = tmp.dHPC(i,2);
            nSpks = count_spks(spks_dHPC, clusterdHPC.aversive, start, stop);
            probability = bayesian_replay(PHIST.dHPC.aversive, nSpks, start, stop);
            plot(TimeVector1,probability,'r'),hold on
            decoded_events.dHPC.aversive = [decoded_events.dHPC.aversive ; probability];
            clear probability nSpks start stop
        end
        title('Decoding - dHPC - Aversive')
        
        subplot(223)
        for i = 1 : length(tmp.vHPC)
            % Decoding using vHPC SUs
            start = tmp.vHPC(i,1);
            stop = tmp.vHPC(i,2);
            nSpks = count_spks(spks_vHPC, clustervHPC.aversive, start, stop);
            probability = bayesian_replay(PHIST.vHPC.aversive, nSpks, start, stop);
            decoded_events.vHPC.aversive = [decoded_events.vHPC.aversive ; probability];
            plot(TimeVector1,probability,'b'),hold on
            clear probability nSpks start stop
        end
        title('Decoding - vHPC - Aversive')

        % Reward
        tmp.vHPC = Restrict(Restrict(replay.vHPC , rewardTS./1000),[coordinatedV_refined(:,2)-0.2,coordinatedV_refined(:,2)+0.2]);
        tmp.dHPC = Restrict(Restrict(replay.dHPC , rewardTS./1000),[coordinated(:,2)-0.2,coordinated(:,2)+0.2]);
        subplot(222)
        for i = 1 : length(tmp.dHPC)
            % Decoding using dHPC SUs
            start = tmp.dHPC(i,1);
            stop = tmp.dHPC(i,2);
            nSpks = count_spks(spks_dHPC, clusterdHPC.reward, start, stop);
            probability = bayesian_replay(PHIST.dHPC.reward, nSpks, start, stop);
            plot(TimeVector1,probability,'r'),hold on
            decoded_events.dHPC.reward = [decoded_events.dHPC.reward ; probability];
            clear probability nSpks start stop
        end
        title('Decoding - dHPC - Reward')

        subplot(224)
        for i = 1 : length(tmp.vHPC)
            % Decoding using vHPC SUs
            start = tmp.vHPC(i,1);
            stop = tmp.vHPC(i,2);
            nSpks = count_spks(spks_vHPC, clustervHPC.reward, start, stop);
            probability = bayesian_replay(PHIST.vHPC.reward, nSpks, start, stop);
            decoded_events.vHPC.reward = [decoded_events.vHPC.reward ; probability];
            plot(TimeVector1,probability,'b'),hold on
            clear probability nSpks start stop
        end
        title('Decoding - vHPC - Reward')
         
        
        % implement replay of the position related to the shock or reward
        t
    end
    tt
end
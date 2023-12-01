clear
clc
close all

%% Parameters
path = {'E:\Rat127\Ephys\pyr';'E:\Rat128\Ephys\in_pyr\to_finish_spikesorting\usable';'E:\Rat103\usable';'E:\Rat132\recordings\in_pyr';'E:\Rat165\in_pyr'};%List of folders from the path

%Sleep
time_criteria = 600; %time criteria to define the maximal time of sleep to include

% Ripples
q = 0.25; %quantile to restrict above it ripples according to their peak amplitude
ripples_coordinated_percentage = []; %for storing percentage of coordnated events across conditions
ripples_coordinated_numbers = []; %for storing the number of cooridnated events all conditions pooled
rateV = []; rateD = []; % rate
rateVB = []; rateDB = []; % rate baseline
rateVR = []; rateDR = []; % rate reward
rateVA = []; rateDA = []; % rate aversive

deltaB = []; deltaR = []; deltaA = []; %to store all the time delta beween dorsal-ventral ripples
% for SU
criteria_fr = 0; %criteria to include or not a SU into the analysis
criteria_n = 3; % minimal number of neurons from each structure
pval = 0.001/2; % p value to define if SU are ripple modulated
ss = 2; %smooth level of CCG
n_SU_V = 0;
n_SU_D = 0;
FR_B_V = []; FR_R_V = [];  FR_A_V = []; % Firing Rate during NREM
FR_B_D = []; FR_R_D = [];  FR_A_D = []; % Firing Rate during NREM
poisson_dHPC_split = []; poisson_vHPC_split = []; %poisson results split by conditions


dRipples = []; dRipplesB = []; dRipplesR = []; dRipplesA = [];
dIRI = []; dIRIB = []; dIRIR = []; dIRIA = [];

dBurstIndexB = []; dBurstIndexR = []; dBurstIndexA = [];
vBurstIndexB = []; vBurstIndexR = []; vBurstIndexA = [];


vRipples = []; vRipplesB = []; vRipplesR = []; vRipplesA = [];
vIRI = []; vIRIB = []; vIRIR = []; vIRIA = [];

vRipples_coordinatedB = []; vRipples_coordinatedR = []; vRipples_coordinatedA = []; 
dRipples_coordinatedB = []; dRipples_coordinatedR = []; dRipples_coordinatedA = []; 

vRipples_uncoordinatedB = []; vRipples_uncoordinatedR = []; vRipples_uncoordinatedA = []; 
dRipples_uncoordinatedB = []; dRipples_uncoordinatedR = []; dRipples_uncoordinatedA = [];

durationsD = []; durationsDB = []; durationsDR = []; durationsDA = []; 
durationsV = []; durationsVB = []; durationsVR = []; durationsVA = []; 

amplitudeD = []; amplitudeDB = []; amplitudeDR = []; amplitudeDA = [];
amplitudeV = []; amplitudeVB = []; amplitudeVR = []; amplitudeVA = [];

coordinatedD_IRI_B = []; coordinatedD_IRI_R = []; coordinatedD_IRI_A = [];
coordinatedV_IRI_B = []; coordinatedV_IRI_R = []; coordinatedV_IRI_A = [];

dRipples_coordinated_single = [];   vRipples_coordinated_single = [];
dRipples_uncoordinated_single = []; vRipples_uncoordinated_single = [];
dRipples_coordinated_single_ds = [];   vRipples_coordinated_single_ds = [];
dRipples_uncoordinated_single_ds = []; vRipples_uncoordinated_single_ds = [];

dRipples_coordinated_single_cross_with_all = []; vRipples_coordinated_single_cross_with_all = [];
dRipples_uncoordinated_single_cross_with_all = []; vRipples_uncoordinated_single_cross_with_all = [];

%coordinated against uncoordinated
dRipples_coordinated_single_cross_with_uncoordinated = [];
vRipples_coordinated_single_cross_with_uncoordinated = [];

%Storage of BUrst Index as we discuss with GG
Burst_Index_cooridnated_vRipples = []; Burst_Index_cooridnated_dRipples = [];
Burst_Index_uncooridnated_vRipples = []; Burst_Index_uncooridnated_dRipples = [];

Burst_Index_cooridnatedB_vRipples = []; Burst_Index_cooridnatedB_dRipples = [];
Burst_Index_cooridnatedR_vRipples = []; Burst_Index_cooridnatedR_dRipples = [];
Burst_Index_cooridnatedA_vRipples = []; Burst_Index_cooridnatedA_dRipples = [];

Burst_Index_uncooridnatedB_vRipples = []; Burst_Index_uncooridnatedB_dRipples = [];
Burst_Index_uncooridnatedR_vRipples = []; Burst_Index_uncooridnatedR_dRipples = [];
Burst_Index_uncooridnatedA_vRipples = []; Burst_Index_uncooridnatedA_dRipples = [];

Burst_Index_allB_vRipples = []; Burst_Index_allB_dRipples = [];
Burst_Index_allR_vRipples = []; Burst_Index_allR_dRipples = [];
Burst_Index_allA_vRipples = []; Burst_Index_allA_dRipples = [];

%% Main tloop to iterate across sessions
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
        
        %     load([cd,'\lfp.mat'])
        %     Time = dHPC(:,1);
        
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

        %% Sleep
        x = dir([cd,'\*-states.mat']);    states = load([cd,'\',x.name]);    states = states.states;
        
        REM = ToIntervals(states==5);    NREM = ToIntervals(states==3);    WAKE = ToIntervals(states==1);
        clear x states
        
        %keep only WAKE in HomeCage
%         WAKE = Restrict(WAKE, [aversiveTS ; rewardTS ; baselineTS] ./1000);
        
        % NREM events restriction according conditions
        NREM_B = NREM(NREM(:,2)<baselineTS(1,2)/1000,:);
        NREM_A = NREM(NREM(:,2)>aversiveTS(1,1)/1000 & NREM(:,2)<aversiveTS(1,2)/1000,:);
        NREM_R = NREM(NREM(:,2)>rewardTS(1,1)/1000 & NREM(:,2)<rewardTS(1,2)/1000,:);
        
        
        % REM events restriction according conditions
        REM_B = REM(REM(:,2)<baselineTS(1,2)/1000,:);
        REM_A = REM(REM(:,2)>aversiveTS(1,1)/1000 & REM(:,2)<aversiveTS(1,2)/1000,:);
        REM_R = REM(REM(:,2)>rewardTS(1,1)/1000 & REM(:,2)<rewardTS(1,2)/1000,:);
        %         load('detected_ripples.mat')
        ripplesD = table2array(readtable('ripplesD_customized2.csv'));
        ripplesV = table2array(readtable('ripplesV_customized2.csv'));
        
        %% dRipples
        rateD = [rateD ; length(ripplesD)/sum(NREM(:,2)-NREM(:,1))];
        durationsD = [durationsD ; ripplesD(:,3)-ripplesD(:,1)];
        dRipples = [dRipples ; ripplesD]; %pool of dRipples
        amplitudeD = [amplitudeD ; ripplesD(:,4)];
        
        %IRI from dRipples
        tmp = diff(ripplesD(:,2)); %pool of inter-dRipples-interval
        dIRI = [dIRI ; tmp];
        clear tmp
        
        %Separation of dRipples per condition
        tmp = Restrict(ripplesD,NREM_B);
        tmp1 = Restrict(ripplesD,NREM_R);
        tmp2 = Restrict(ripplesD,NREM_A);
        dRipplesB = [dRipplesB ; tmp]; 
        dRipplesR = [dRipplesR ; tmp1];
        dRipplesA = [dRipplesA ; tmp2];
        %Durations per condition
        durationsDB = [durationsDB ; tmp(:,3)-tmp(:,1)];
        durationsDR = [durationsDR ; tmp1(:,3)-tmp1(:,1)];
        durationsDA = [durationsDA ; tmp2(:,3)-tmp2(:,1)];
        %Amplitude per conditions
        amplitudeDB = [amplitudeDB ; tmp(:,4)];
        amplitudeDR = [amplitudeDR ; tmp1(:,4)];
        amplitudeDA = [amplitudeDA ; tmp2(:,4)];
        rateDB = [rateDB ; length(tmp)/sum(NREM_B(:,2)-NREM_B(:,1))];
        rateDR = [rateDR ; length(tmp1)/sum(NREM_R(:,2)-NREM_R(:,1))];
        rateDA = [rateDA ; length(tmp2)/sum(NREM_A(:,2)-NREM_A(:,1))];

        clear tmp tmp1 tmp2
        
        %dIRI per condtion
        tmp = diff(Restrict(ripplesD(:,2),NREM_B)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        dIRIB = [dIRIB ; tmp1/length(tmp)];
        clear tmp N EDGES ii p iii tmp1
        
        tmp = diff(Restrict(ripplesD(:,2),NREM_R)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.2 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        dIRIR = [dIRIR ; tmp1/length(tmp)];
        clear tmp N EDGES ii p iii tmp1
        
        tmp = diff(Restrict(ripplesD(:,2),NREM_A)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        dIRIA = [dIRIA ; tmp1/length(tmp)];
        clear tmp N EDGES ii p iii tmp1
        
        % Burst-Index dRipples Baseline
        x = Restrict(ripplesD(:,2),NREM_B);
        y = Restrict(ripplesD(:,2),NREM_B);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg'); %ccg calculation
        [i,ii] = min(abs(tttt-0.00)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstD = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1)); %burst index calculation
        dBurstIndexB = [dBurstIndexB ; BurstD];
        clear tttt ccg i ii iii iiii BurstD x y
        
        % Burst-Index dRipples Reward
        x = Restrict(ripplesD(:,2),NREM_R);
        y = Restrict(ripplesD(:,2),NREM_R);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
        [i,ii] = min(abs(tttt-0.00)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstD = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
        dBurstIndexR = [dBurstIndexR ; BurstD];
        clear tttt ccg i ii iii iiii BurstD x y
        
        % Burst-Index dRipples Aversive
        x = Restrict(ripplesD(:,2),NREM_A);
        y = Restrict(ripplesD(:,2),NREM_A);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
        [i,ii] = min(abs(tttt-0.00)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstD = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
        dBurstIndexA = [dBurstIndexA ; BurstD];
        clear tttt ccg i ii iii iiii BurstD x y
        
        %% vRipples
        rateV = [rateV ; length(ripplesV)/sum(NREM(:,2)-NREM(:,1))];
        durationsV = [durationsV ; ripplesV(:,3)-ripplesV(:,1)];
        vRipples = [vRipples ; ripplesV]; %pool of vRipples
        amplitudeV = [amplitudeV ; ripplesV(:,4)];

        %IRI from dRipples
        tmp = diff(ripplesV(:,2)); %pool of inter-dRipples-interval
        vIRI = [vIRI ; tmp];
        clear tmp

        %Separation of vRipples per condition
        tmp = Restrict(ripplesV,NREM_B);
        tmp1 = Restrict(ripplesV,NREM_R);
        tmp2 = Restrict(ripplesV,NREM_A);
        vRipplesB = [vRipplesB ; tmp]; 
        vRipplesR = [vRipplesR ; tmp1];
        vRipplesA = [vRipplesA ; tmp2];
        %Durations per condition
        durationsVB = [durationsVB ; tmp(:,3)-tmp(:,1)];
        durationsVR = [durationsVR ; tmp1(:,3)-tmp1(:,1)];
        durationsVA = [durationsVA ; tmp2(:,3)-tmp2(:,1)];
        %Amplitude per conditions
        amplitudeVB = [amplitudeVB ; tmp(:,4)];
        amplitudeVR = [amplitudeVR ; tmp1(:,4)];
        amplitudeVA = [amplitudeVA ; tmp2(:,4)];
        rateVB = [rateVB ; length(tmp)/sum(NREM_B(:,2)-NREM_B(:,1))];
        rateVR = [rateVR ; length(tmp1)/sum(NREM_R(:,2)-NREM_R(:,1))];
        rateVA = [rateVA ; length(tmp2)/sum(NREM_A(:,2)-NREM_A(:,1))];        
        clear tmp tmp1 tmp2
        
        %vIRI per condtion
        tmp = diff(Restrict(ripplesV(:,2),NREM_B)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        vIRIB = [vIRIB ; tmp1/length(tmp)];
        clear tmp N EDGES ii p iii tmp1

        tmp = diff(Restrict(ripplesV(:,2),NREM_R)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        vIRIR = [vIRIR ; tmp1/length(tmp)];        
        clear tmp N EDGES ii p iii tmp1

        tmp = diff(Restrict(ripplesV(:,2),NREM_A)); %pool of inter-dRipples-interval
        tmp1 = 0;
        for p = 1:length(tmp)
            if and(tmp(p) <= 0.20 , tmp(p)>=0.00)
                tmp1 = tmp1 + 1;
            end
        end
        vIRIA = [vIRIA ; tmp1/length(tmp)];
        clear tmp N EDGES ii p iii tmp1
        
        % Burst-Index vRipples Baseline
        x = Restrict(ripplesV(:,2),NREM_B);
        y = Restrict(ripplesV(:,2),NREM_B);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg'); %ccg calculation
        [i,ii] = min(abs(tttt-0.0)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstV = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1)); %burst index calculation
        vBurstIndexB = [vBurstIndexB ; BurstV];
        clear tttt ccg i ii iii iiii BurstV x y
        
        % Burst-Index vRipples Reward
        x = Restrict(ripplesV(:,2),NREM_R);
        y = Restrict(ripplesV(:,2),NREM_R);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
        [i,ii] = min(abs(tttt-0.0)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstV = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
        vBurstIndexR = [vBurstIndexR ; BurstV];
        clear tttt ccg i ii iii iiii BurstV x y
        
        % Burst-Index vRipples Aversive
        x = Restrict(ripplesV(:,2),NREM_A);
        y = Restrict(ripplesV(:,2),NREM_A);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
        [i,ii] = min(abs(tttt-0.0)); %looking for the index of limits for burst calculation in time
        [i,iii] = min(abs(tttt-0.20));
        [i,iiii] = min(abs(tttt-0));
        BurstV = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
        vBurstIndexA = [vBurstIndexA ; BurstV];
        clear tttt ccg i ii iii iiii BurstV x y
        
        %% Coordinated dHPC ripples
        coordinated = [];
        coordinatedV = [];
        coordinatedV_refined = [];
        tmpB_D = [];        tmpR_D = [];        tmpA_D = [];
        tmpB_V = [];        tmpR_V = [];        tmpA_V = [];

        for i = 1:length(ripplesD)
            r = ripplesD(i,:);
            tmp = sum(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1));
            if tmp>0
                z = ripplesV(and(ripplesV(:,2)>= r(1,2)-0.1, ripplesV(:,2)<= r(1,2)+0.1),:);
                coordinatedV = [coordinatedV ; z];
                [p,indice] = min(abs(r(2)-z(:,2)));
                coordinatedV_refined = [coordinatedV_refined ; z(indice,:)];
                coordinated = [coordinated ; r];
                
                if and(r(1,2)>baselineTS(1)/1000 , r(1,2)<baselineTS(2)/1000)
                    iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                    tmpB_D = [tmpB_D ; iri_t]; clear iri_t
                    
                    iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                    tmpB_V = [tmpB_V ; z(indice,2)]; clear iri_t
                elseif and(r(1,2)>rewardTS(1)/1000 , r(1,2)<rewardTS(2)/1000)
                    iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                    tmpR_D = [tmpR_D ; iri_t]; clear iri_t
                    
                    iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                    tmpR_V = [tmpR_V ; z(indice,2)]; clear iri_t
                elseif and(r(1,2)>aversiveTS(1)/1000 , r(1,2)<aversiveTS(2)/1000)
                    iri_t = ripplesD(and(ripplesD(:,2)>r(1,2)-0.4 , ripplesD(:,2)<r(1,2)+0.4),2);
                    tmpA_D = [tmpA_D ; iri_t]; clear iri_t
                    
                    iri_t = ripplesV(and(ripplesV(:,2)>z(indice,2)-0.4 , ripplesV(:,2)<z(indice,2)+0.4),2);
                    tmpA_V = [tmpA_V ; z(indice,2)]; clear iri_t
                end
                
                clear tmp2 tmp1 p indice z
            end
            clear r
        end
        clear x tmp i
        
        coordinatedB = Restrict(coordinated,NREM_B);    coordinatedA = Restrict(coordinated,NREM_A);    coordinatedR = Restrict(coordinated,NREM_R);
        coordinatedB_V = Restrict(coordinatedV_refined,NREM_B);    coordinatedR_V = Restrict(coordinatedV_refined,NREM_R);    coordinatedA_V = Restrict(coordinatedV_refined,NREM_A);
        coordinatedB_V_non_refined = Restrict(coordinatedV,NREM_B);    coordinatedR_V_non_refined  = Restrict(coordinatedV,NREM_R);    coordinatedA_V_non_refined  = Restrict(coordinatedV,NREM_A);
%         coordinatedB_V = Restrict(coordinatedV,NREM_B);    coordinatedR_V = Restrict(coordinatedV,NREM_R);    coordinatedA_V = Restrict(coordinatedV,NREM_A);

        % Detection of uncoordinated ripples
        uncoordinated = ripplesD(~ismember(ripplesD(:,1),coordinated(:,1)),:);
        uncoordinatedV = ripplesV(~ismember(ripplesV(:,1),coordinatedV_refined(:,1)),:);
        
        %% Detection of ripple burst
        coordinated_ripple_bursts = [];
        for i = 1:length(coordinated)
            [p,ii] = min(abs(coordinatedV_refined(:,2) - coordinated(i,2)));
            
            %Burst detection for dorsal ripples
            tmp = and(ripplesD(:,2) > coordinated(i)-0.2 , ripplesD(:,2) < coordinated(i)+0.2);
            if sum(tmp)>1
                burstD = [coordinated(i,:) ; ripplesD(tmp,:)];
            else
                burstD = [coordinated(i,:)];
            end
            clear tmp
            
            %Burst detection for ventral ripples
            tmp = and(ripplesV(:,2) > coordinatedV_refined(ii,2)-0.2 , ripplesV(:,2) < coordinatedV_refined(ii,2)+0.2);
            if sum(tmp)>1
                burstV = [coordinatedV_refined(ii,:) ; ripplesV(tmp,:)];
            else
                burstV = [coordinatedV_refined(ii,:)];
            end
            clear p ii tmp
            
            if or(size(burstD,1)>1 , size(burstV,1)>1)
%                 %To keep the overlapping area
%                 tmp = [min(burstD(:,1)) , max(burstD(:,3)) ; min(burstV(:,1)) , max(burstV(:,3))];
%                 tmp = [max(tmp(:,1)) , max(tmp(:,1)) + ((min(tmp(:,2))-max(tmp(:,1)))/2) , min(tmp(:,2))];
                % To keep the entier burst
                tmp = [burstD;burstV];
                tmp = [min(tmp(:,1)) , ((max(tmp(:,3))-min(tmp(:,1)))/2)+min(tmp(:,1)), max(tmp(:,3))];
                
                coordinated_ripple_bursts = [coordinated_ripple_bursts ; tmp];
                clear tmp
            end
        end
        save ([cd,'\coordinated_ripple_bursts.mat'],'coordinated_ripple_bursts')   
        
%         tmp = Restrict(ripplesV(:,2),NREM_B);
%         tmp1 = Restrict(ripplesV(:,2),NREM_R);
%         tmp2 = Restrict(ripplesV(:,2),NREM_A);
%         tmp3 = Restrict(ripplesD(:,2),NREM_B);
%         tmp4 = Restrict(ripplesD(:,2),NREM_R);
%         tmp5 = Restrict(ripplesD(:,2),NREM_A);
%         
%         vRipples_coordinatedB = [vRipples_coordinatedB ; (length(coordinatedB_V)*100)/length(tmp)];
%         vRipples_coordinatedR = [vRipples_coordinatedR ; (length(coordinatedR_V)*100)/length(tmp1)];
%         vRipples_coordinatedA = [vRipples_coordinatedA ; (length(coordinatedA_V)*100)/length(tmp2)]; 
%         dRipples_coordinatedB = [dRipples_coordinatedB ; (length(coordinatedB)*100)/length(tmp3)];
%         dRipples_coordinatedR = [dRipples_coordinatedR ; (length(coordinatedR)*100)/length(tmp4)];
%         dRipples_coordinatedA = [dRipples_coordinatedA ; (length(coordinatedA)*100)/length(tmp5)]; 
%         clear tmp tmp1 tmp2 tmp3 tmp4 tmp5
        
        
%         vRipples_coordinatedB = [vRipples_coordinatedB ; coordinatedB_V];
%         vRipples_coordinatedR = [vRipples_coordinatedR ; coordinatedR_V];
%         vRipples_coordinatedA = [vRipples_coordinatedA ; coordinatedA_V]; 
%         dRipples_coordinatedB = [dRipples_coordinatedB ; coordinatedB];
%         dRipples_coordinatedR = [dRipples_coordinatedR ; coordinatedR];
%         dRipples_coordinatedA = [dRipples_coordinatedA ; coordinatedA];         
        
        uncoordinatedB = Restrict(uncoordinated,NREM_B);    uncoordinatedA = Restrict(uncoordinated,NREM_A);    uncoordinatedR = Restrict(uncoordinated,NREM_R);
        uncoordinatedB_V = Restrict(uncoordinatedV,NREM_B);    uncoordinatedR_V = Restrict(uncoordinatedV,NREM_R);    uncoordinatedA_V = Restrict(uncoordinatedV,NREM_A);
 
        vRipples_uncoordinatedB = [vRipples_uncoordinatedB ; uncoordinatedB_V];
        vRipples_uncoordinatedR = [vRipples_uncoordinatedR ; uncoordinatedR_V];
        vRipples_uncoordinatedA = [vRipples_uncoordinatedA ; uncoordinatedA_V]; 
        dRipples_uncoordinatedB = [dRipples_uncoordinatedB ; uncoordinatedB];
        dRipples_uncoordinatedR = [dRipples_uncoordinatedR ; uncoordinatedR];
        dRipples_uncoordinatedA = [dRipples_uncoordinatedA ; uncoordinatedA];       
        
        %% For storing number of coordinated events
        ripples_coordinated_numbers = [ripples_coordinated_numbers ; length(coordinated) , length(ripplesD) , length(coordinatedV) , length(ripplesV)];

        %% Autoccorelogram coordionated and uncooridnated without subsampling
        % Burst Index calculation as we discuss with GG
        % coordinated dRipples
        x = coordinated(:,2);
        y = coordinated(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        dRipples_coordinated_single = [dRipples_coordinated_single , ccg];
        
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnated_dRipples = [Burst_Index_cooridnated_dRipples ; count/length(x)];
        clear ccg time x y count
        
        % coordinated vRipples
        x = unique(coordinatedV(:,2));
        y = unique(coordinatedV(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        vRipples_coordinated_single = [vRipples_coordinated_single , ccg];

        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnated_vRipples = [Burst_Index_cooridnated_vRipples ; count/length(x)];
        clear ccg time x y count
        
        % uncoordinated dRipples
        x = uncoordinated(:,2);
        y = uncoordinated(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        dRipples_uncoordinated_single = [dRipples_uncoordinated_single , ccg];
        
        count = 0;
        y = ripplesD(:,2); %exclude coordinated + associated ripples
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnated_dRipples = [Burst_Index_uncooridnated_dRipples ; count/length(x)];
        clear ccg time x y count
        
        % uncoordinated vRipples
        x = unique(uncoordinatedV(:,2));
        y = unique(uncoordinatedV(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        vRipples_uncoordinated_single = [vRipples_uncoordinated_single , ccg];
        
        count = 0;
        y = ripplesV(:,2);%exclude coordinated + associated ripples
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnated_vRipples = [Burst_Index_uncooridnated_vRipples ; count/length(x)];
        clear ccg time x y count
        
        %% Autoccorelogram coordionated and uncooridnated with subsampling
       % Burst Index calculation as we discuss with GG
        % coordinated dRipples
        x = coordinated(:,2);
        y = coordinated(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        dRipples_coordinated_single_ds = [dRipples_coordinated_single_ds , ccg];
        clear ccg time x y count tmp
        
        % coordinated vRipples
        x = unique(coordinatedV(:,2));
        y = unique(coordinatedV(:,2));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,1);
        vRipples_coordinated_single_ds = [vRipples_coordinated_single_ds , ccg];
        clear ccg time x y count tmp
        
        % uncoordinated dRipples
        tmp = [];
        for i = 1:50
            x = uncoordinated(randperm(length(uncoordinated)),2);
            x = x(1:length(coordinated));
            y = x;
            [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
            [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
            ccg = ccg(:,1,1);
            tmp = [tmp , ccg];
            clear x y ccg time s ids groups
        end
        dRipples_uncoordinated_single_ds = [dRipples_uncoordinated_single_ds , mean(tmp,2)];
        clear ccg time x y count tmp
        
        % uncoordinated vRipples
        if length(coordinatedV) < length(uncoordinatedV)
            tmp = [];
            for i = 1:50
                x = uncoordinatedV(randperm(length(uncoordinatedV)),2);
                x = x(1:length(coordinatedV));
                y = x;
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
                ccg = ccg(:,1,1);
                tmp = [tmp , ccg];
                clear x y ccg time s ids groups
            end
        else
            tmp = [];
            for i = 1:50
                x = uncoordinatedV(randperm(length(uncoordinatedV)),2);
                y = x;
                [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
                [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
                ccg = ccg(:,1,1);
                tmp = [tmp , ccg];
                clear x y ccg time s ids groups
            end
        end
        
        vRipples_uncoordinated_single_ds = [vRipples_uncoordinated_single_ds , mean(tmp,2)];
        clear ccg time x y count tmp
        %% Crosscorrelogram coordionated/uncooridnated with all ripples
        % Burst Index calculation as we discuss with GG
        % coordinated dRipples
        x = coordinated(:,2);
        y = ripplesD(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        dRipples_coordinated_single_cross_with_all = [dRipples_coordinated_single_cross_with_all , ccg];
        clear ccg time x y count tmp
        
        % coordinated vRipples
        x = unique(coordinatedV(:,2));
        y = ripplesV(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        vRipples_coordinated_single_cross_with_all = [vRipples_coordinated_single_cross_with_all , ccg];
        clear ccg time x y count tmp
        
        % uncoordinated dRipples
        x = uncoordinated(:,2);
        y = ripplesD(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        dRipples_uncoordinated_single_cross_with_all = [dRipples_uncoordinated_single_cross_with_all , ccg];
        clear ccg time x y count tmp
        
        % uncoordinated dRipples
        x = uncoordinatedV(:,2);
        y = ripplesV(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        vRipples_uncoordinated_single_cross_with_all = [vRipples_uncoordinated_single_cross_with_all , ccg];
        clear ccg time x y count tmp
        
        %% Crosscorrelogram coordionated vs uncooridnated
        % coordinated dRipples
        x = coordinated(:,2);
        y = uncoordinated(:,2);
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        dRipples_coordinated_single_cross_with_uncoordinated = [dRipples_coordinated_single_cross_with_uncoordinated , ccg];
        clear ccg time x y count tmp
        
        % coordinated vRipples
        x = unique(coordinatedV(:,2));
        y = unique(uncoordinatedV(:,2));
        y = y(~ismember(y,x));
        [s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
        [ccg,time] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',1,'mode','ccg');
        ccg = ccg(:,1,2);
        vRipples_coordinated_single_cross_with_uncoordinated = [vRipples_coordinated_single_cross_with_uncoordinated , ccg];
        clear ccg time x y count tmp
        
        
        %% Burst Index calculation as we discuss with GG
        % per condition
        % coordinated dRipples
        x = coordinatedB(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                break
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedB_dRipples = [Burst_Index_cooridnatedB_dRipples ; count/length(x)];
        clear ccg time x y count
        
        x = coordinatedR(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedR_dRipples = [Burst_Index_cooridnatedR_dRipples ; count/length(x)];
        clear ccg time x y count

        x = coordinatedA(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedA_dRipples = [Burst_Index_cooridnatedA_dRipples ; count/length(x)];
        clear ccg time x y count        
        
        
        x = coordinatedB_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedB_vRipples = [Burst_Index_cooridnatedB_vRipples ; count/length(x)];
        clear ccg time x y count
        
        x = coordinatedR_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedR_vRipples = [Burst_Index_cooridnatedR_vRipples ; count/length(x)];
        clear ccg time x y count

        x = coordinatedA_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_cooridnatedA_vRipples = [Burst_Index_cooridnatedA_vRipples ; count/length(x)];
        clear ccg time x y count               
        
        % uncoordinated
        x = uncoordinatedB(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedB_dRipples = [Burst_Index_uncooridnatedB_dRipples ; count/length(x)];
        clear ccg time x y count
        
        x = uncoordinatedR(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedR_dRipples = [Burst_Index_uncooridnatedR_dRipples ; count/length(x)];
        clear ccg time x y count

        x = uncoordinatedA(:,2);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedA_dRipples = [Burst_Index_uncooridnatedA_dRipples ; count/length(x)];
        clear ccg time x y count        
        
        
        x = uncoordinatedB_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedB_vRipples = [Burst_Index_uncooridnatedB_vRipples ; count/length(x)];
        clear ccg time x y count
        
        x = uncoordinatedR_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedR_vRipples = [Burst_Index_uncooridnatedR_vRipples ; count/length(x)];
        clear ccg time x y count

        x = uncoordinatedA_V(:,2);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_uncooridnatedA_vRipples = [Burst_Index_uncooridnatedA_vRipples ; count/length(x)];
        clear ccg time x y count                 

        % All
        x = Restrict(ripplesD(:,2),NREM_B);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allB_dRipples = [Burst_Index_allB_dRipples ; count/length(x)];
        clear ccg time x y count
        
        x = Restrict(ripplesD(:,2),NREM_R);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allR_dRipples = [Burst_Index_allR_dRipples ; count/length(x)];
        clear ccg time x y count

        x = Restrict(ripplesD(:,2),NREM_A);  
        count = 0;
        y = ripplesD(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allA_dRipples = [Burst_Index_allA_dRipples ; count/length(x)];
        clear ccg time x y count        
        
        
        x = Restrict(ripplesV(:,2),NREM_B);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allB_vRipples = [Burst_Index_allB_vRipples ; count/length(x)];
        clear ccg time x y count
        
        x = Restrict(ripplesV(:,2),NREM_R);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allR_vRipples = [Burst_Index_allR_vRipples ; count/length(x)];
        clear ccg time x y count

        x = Restrict(ripplesV(:,2),NREM_A);  
        count = 0;
        y = ripplesV(:,2);
        for s = 1 : length(x)
            seed = x(s);
            tmp = sum(and(y > seed-0.2 , y < seed+0.2));
            if tmp>1
                count = count + 1;
            end
        end
        Burst_Index_allA_vRipples = [Burst_Index_allA_vRipples ; count/length(x)];
        clear ccg time x y count
        
        t
        clear ripplesD ripplesV
        clear aversiveTS aversiveTS_run rewardTS rewardTS_run baselineTS
        clear coordinated coordinatedA coordinatedB coordinatedR
        clear coordinatedV coordinatedA_V coordinatedB_V coordinatedR_V
        clear coordinatedA_V_non_refined coordinatedB_V_non_refined coordinatedR_V_non_refined
        clear uncoordinated uncoordinatedA uncoordinatedA_V uncoordinatedB uncoordinatedB_V
        clear uncoordinatedR uncoordinatedR_V uncoordinatedV
        clear REM REM_A REM_B REM_R NREM NREM_A NREM_B NREM_R WAKE
        
    end
    tt
    
end


%% ALl
x = ripplesD(:,2);
y = ripplesD(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',0,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstD(1) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
figure,
ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg),hold on
clear tttt ccg i ii iii iiii

x = vRipples(:,2);
y = vRipples(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',0,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstV(1) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg),ylim([0 1.01])%,xlim([-1 1])
xline(0.05,'--'), xline(0.2,'--')
yline(mean(ccg))
clear ccg tttt i ii iii iiii


%% Separated
%dHPC Baseline
x = dRipplesB(:,2);
y = dRipplesB(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg'); %ccg calculation
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstD(2) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1)); %burst index calculation
figure,
ccg = ccg(:,1,1)-min(ccg(:,1,1)); %normalization of ccg to be plotted from 0 to 1 in y
ccg = ccg./max(ccg);
plot(tttt,ccg,'k'),hold on
clear tttt ccg i ii iii iiii

%dHPC Reward
x = dRipplesR(:,2);
y = dRipplesR(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstD(3) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg,'b'),ylim([0 1.01])%,xlim([-1 1])
clear tttt ccg i ii iii iiii

%dHPC Aversive
x = dRipplesA(:,2);
y = dRipplesA(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstD(4) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));
ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg,'r'),ylim([0 1.01])%,xlim([-1 1])
xline(0.05,'--'), xline(0.2,'--')
clear tttt ccg i ii iii iiii

%vHPC Baseline
x = vRipplesB(:,2);
y = vRipplesB(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstV(2) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));

figure,
ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg,'k'),hold on
clear tttt ccg i ii iii iiii

%vHPC Reward
x = vRipplesR(:,2);
y = vRipplesR(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstV(3) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));

ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg,'b'),ylim([0 1.01])%,xlim([-1 1])
clear tttt ccg i ii iii iiii

%vHPC Aversive
x = vRipplesA(:,2);
y = vRipplesA(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',2,'mode','ccg');
[i,ii] = min(abs(tttt-0.05)); %looking for the index of limits for burst calculation in time
[i,iii] = min(abs(tttt-0.15));
[i,iiii] = min(abs(tttt-0));
BurstV(4) = max(ccg(ii:iii,1,1))/mean(ccg(iiii:end,1,1));

ccg = ccg(:,1,1)-min(ccg(:,1,1));
ccg = ccg./max(ccg);
plot(tttt,ccg,'r'),ylim([0 1.01])%,xlim([-1 1])
xline(0.05,'--'), xline(0.2,'--')
clear tttt ccg i ii iii iiii

%% Cross-Correlograms
% dorsal-ventral Ripples per condition
x = dRipples(:,2);
y = vRipples(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',4,'smooth',1,'mode','ccg');
ccg = ccg(:,1,2)%./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral Ripples per condition
x = dRipplesB(:,2);
y = vRipplesB(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

% dorsal-ventral Ripples
x = dRipplesR(:,2);
y = vRipplesR(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'b','LineWidth',2)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;

% dorsal-ventral Ripples
x = dRipplesA(:,2);
y = vRipplesA(:,2);
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'r','LineWidth',2)
ax = gca; % axes handle
ax.YAxis.Exponent = 0;


%% Cross-Correlograms
% dorsal-ventral coordinated and uncoordinated Ripples
data = [dRipples_coordinatedB ; dRipples_coordinatedR ; dRipples_coordinatedA];
data1 = [dRipples_uncoordinatedB ; dRipples_uncoordinatedR ; dRipples_uncoordinatedA];
data2 = [vRipples_coordinatedB ; vRipples_coordinatedR ; vRipples_coordinatedA];
data3 = [vRipples_uncoordinatedB ; vRipples_uncoordinatedR ; vRipples_uncoordinatedA];
x = unique(data2(:,2));
y =  unique(data3(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

x = unique(data2(:,2));
y =  unique(data3(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples
x = unique(data(:,2));
y =  unique(data2(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',2,'smooth',2,'mode','ccg');
ccg = ccg(:,1,2)./sum(ccg(:,1,2));
figure,plot(tttt,ccg,'b','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')


% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedB];
data1 = [dRipples_uncoordinatedB];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
figure,plot(tttt,ccg,'k','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedR];
data1 = [dRipples_uncoordinatedR];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
plot(tttt,ccg,'b','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

% dorsal-ventral coordinated Ripples per condition
data = [dRipples_uncoordinatedA];
data1 = [dRipples_uncoordinatedA];
x = unique(data(:,2));
y =  unique(data1(:,2));
[s,ids,groups] = CCGParameters(y,ones(length(y),1),x,ones(length(x),1)*2);
[ccg,tttt] = CCG(s,ids,'binSize',0.01,'duration',1,'smooth',2,'mode','ccg');
ccg = ccg(:,1,1)./sum(ccg(:,1,1));
plot(tttt,ccg,'r','LineWidth',2),hold on
% ax = gca; % axes handle
% ax.YAxis.Exponent = 0;
xline(0,'--')

%% IRI
tmp1 = dIRI;
tmp2 = vIRI;

[y,x]=histcounts(tmp1,100,'BinLimits',[0 1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'r','LineWidth',2),xlim([0.015 0.1]),hold on
bar(x(2:end),Smooth(y,2),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none'),hold on
clear x y

[y,x]=histcounts(tmp2,100,'BinLimits',[0 1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'b','LineWidth',2),xlim([0.015 0.1])
bar(x(2:end),Smooth(y,2),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none')%,xlim([0.015 0.1])
clear x y tmp1 tmp2


% IRI distributions per conditions
tmp1 = diff(vRipplesB(:,2));
% tmp1(tmp1>2) = [];
tmp2 = diff(vRipplesR(:,2));
% tmp2(tmp2>2) = [];
tmp3 = diff(vRipplesA(:,2));
% tmp3(tmp3>2) = [];

[y,x]=histcounts(tmp1,200,'BinLimits',[0 1],'Normalization','probability')
figure,plot(x(2:end),Smooth(y,1),'k','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','k','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp2,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'b','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','b','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp3,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'r','LineWidth',2)
% bar(x(2:end),y,'FaceColor','r','FaceAlpha',0.5),hold on
clear x y


% IRI distributions per conditions
tmp1 = diff(dRipplesB(:,2));
% tmp1(tmp1>2) = [];
tmp2 = diff(dRipplesR(:,2));
% tmp2(tmp2>2) = [];
tmp3 = diff(dRipplesA(:,2));
% tmp3(tmp3>2) = [];

[y,x]=histcounts(tmp1,200,'BinLimits',[0 1],'Normalization','probability')
figure,plot(x(2:end),Smooth(y,1),'k','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','k','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp2,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'b','LineWidth',2),hold on
% bar(x(2:end),y,'FaceColor','b','FaceAlpha',0.5),hold on
clear x y

[y,x]=histcounts(tmp3,200,'BinLimits',[0 1],'Normalization','probability')
plot(x(2:end),Smooth(y,1),'r','LineWidth',2)
% bar(x(2:end),y,'FaceColor','r','FaceAlpha',0.5),hold on
clear x y

% Plot cummulative IRI graph
figure,boxplot([dIRIB dIRIR dIRIA vIRIB vIRIR vIRIA])
data = [dIRIB dIRIR dIRIA vIRIB vIRIR vIRIA];

%% Durations
tmp1 = dRipples(:,3)-dRipples(:,1);
tmp2 = vRipples(:,3)-vRipples(:,1);

[y,x]=histcounts(tmp1,100,'BinLimits',[0 0.1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'r','LineWidth',2),xlim([0.015 0.1]),hold on
bar(x(2:end),Smooth(y,2),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none'),hold on
clear x y

[y,x]=histcounts(tmp2,100,'BinLimits',[0 0.1],'Normalization','probability')
% plot(x(2:end),Smooth(y,2),'b','LineWidth',2),xlim([0.015 0.1])
bar(x(2:end),Smooth(y,2),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none'),xlim([0.015 0.1])
clear x y tmp1 tmp2

%% Ripples cooridnated
data = [dRipples_coordinatedB dRipples_coordinatedR dRipples_coordinatedA vRipples_coordinatedB vRipples_coordinatedR vRipples_coordinatedA];
figure,boxplot(data)
[p,t,stats] = anova1(data);
[c,m,h,nms] = multcompare(stats);

%% burst
figure,boxplot([dBurstIndexB dBurstIndexR dBurstIndexA vBurstIndexB vBurstIndexR vBurstIndexA])

figure,boxplot([rateD rateV])


%% Burst Index of coordinated dorsal and ventral ripples (if is positive means that dRipples occure first)

data = [Burst_Index_uncooridnated_dRipples , Burst_Index_cooridnated_dRipples Burst_Index_uncooridnated_vRipples , Burst_Index_cooridnated_vRipples];
figure,boxplot(data)

% per condition
% just coordinated
data = [Burst_Index_cooridnatedB_dRipples , Burst_Index_cooridnatedR_dRipples , Burst_Index_cooridnatedA_dRipples Burst_Index_cooridnatedB_vRipples , Burst_Index_cooridnatedR_vRipples , Burst_Index_cooridnatedA_vRipples];
subplot(1,3,1),boxplot(data),ylim([0 0.8])
% anova1(data)

% just uncoordinated
data = [Burst_Index_uncooridnatedB_dRipples , Burst_Index_uncooridnatedR_dRipples , Burst_Index_uncooridnatedA_dRipples Burst_Index_uncooridnatedB_vRipples , Burst_Index_uncooridnatedR_vRipples , Burst_Index_uncooridnatedA_vRipples];
subplot(1,3,2),boxplot(data),ylim([0 0.8])
% anova1(data)


% just coordinated
data = [Burst_Index_allB_dRipples , Burst_Index_allR_dRipples , Burst_Index_allA_dRipples Burst_Index_allB_vRipples , Burst_Index_allR_vRipples , Burst_Index_allA_vRipples];
subplot(1,3,3),boxplot(data)
% anova1(data)

%% Correlograms session by session
% Autocorrelograms without subsampling
x = sum(dRipples_coordinated_single,2)./sum(sum(dRipples_coordinated_single));
y = sum(dRipples_uncoordinated_single,2)./sum(sum(dRipples_uncoordinated_single));
figure,subplot(4,2,1),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single,2)./sum(sum(vRipples_coordinated_single));
y = sum(vRipples_uncoordinated_single,2)./sum(sum(vRipples_uncoordinated_single));
subplot(4,2,2),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

% Autocorrelograms session by session with subsampling
x = sum(dRipples_coordinated_single_ds,2)./sum(sum(dRipples_coordinated_single_ds));
y = sum(dRipples_uncoordinated_single_ds,2)./sum(sum(dRipples_uncoordinated_single_ds));
subplot(4,2,3),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_ds,2)./sum(sum(vRipples_coordinated_single_ds));
y = sum(vRipples_uncoordinated_single_ds,2)./sum(sum(vRipples_uncoordinated_single_ds));
subplot(4,2,4),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

% Correlograms session by session without subsampling
x = sum(dRipples_coordinated_single_cross_with_all,2)./sum(sum(dRipples_coordinated_single_cross_with_all));
y = sum(dRipples_uncoordinated_single_cross_with_all,2)./sum(sum(dRipples_uncoordinated_single_cross_with_all));
subplot(4,2,5),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_cross_with_all,2)./sum(sum(vRipples_coordinated_single_cross_with_all));
y = sum(vRipples_uncoordinated_single_cross_with_all,2)./sum(sum(vRipples_uncoordinated_single_cross_with_all));
subplot(4,2,6),plot([-1:0.01:1],x)
hold on
plot([-1:0.01:1],y),ylim([0 0.02])

%% Correlograms session by session with subsampling
x = sum(dRipples_coordinated_single_cross_with_uncoordinated,2)./sum(sum(dRipples_coordinated_single_cross_with_uncoordinated));
subplot(4,2,7),plot([-1:0.01:1],x), ylim([0 0.02])
clear x y

x = sum(vRipples_coordinated_single_cross_with_uncoordinated,2)./sum(sum(vRipples_coordinated_single_cross_with_uncoordinated));
subplot(4,2,8),plot([-1:0.01:1],x),ylim([0 0.02])
clear x y

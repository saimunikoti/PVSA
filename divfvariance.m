clear all
close all
clc

sim_length = 100000;
actorlist= [22,14,7,17,8];
tic
for countactor = 0:length(actorlist)
    
    tempactorlist =  actorlist;
    
    if countactor ~= 0
       tempactorlist(countactor)=[];
    end
    
    tempactorlist = string(tempactorlist);

    fileext = join(tempactorlist,"_");    

    actornodename = fileext + "_actor";
    powervecname = './powervec_' + fileext + '.mat';
    
    disp('powervecname');
    disp(powervecname);
    
    volt_tensor = zeros(37,3, sim_length);
    voltchange_tensor = zeros(37,3,sim_length);

    temppower = load(powervecname);
    R = temppower.powervec;

    load_data = load('load_data.mat');
    S_org = load_data.S;

    %% basevolt

    S = S_org;
    basevolt = loadflow_function(S);

    %% for loop

    for countrun = 1:sim_length % for different run

        S = S_org; 
        deltas = zeros(37,6);
        deltas(:,1) = R(countrun,1:37);
        deltas(:,2) = R(countrun,38:74);
        deltas(:,3) = R(countrun,75:111);
        deltas(:,4) = R(countrun,112:148);
        deltas(:,5) = R(countrun,149:185);
        deltas(:,6) = R(countrun,186:222);

        S= S + deltas; % add delta load into base load

        try
            V = loadflow_function(S); % call load flow function
        catch
            disp('skip the current loop')
            continue;
        end

         volt_tensor(:,:, countrun) = V;
         voltchange_tensor(:,:, countrun) = V-basevolt;

    %     volt_tensor2(countrun,:) = V(:,2);
    %     volt_tensor3(countrun,:) = V(:,3);
        clear S
        disp(countrun)
    end


    disp(toc)

    % save voltage 
    filename = strcat('v_',actornodename,'.mat');
    save(filename, 'volt_tensor')

    % save base voltage
    filename = strcat('basevolt_',actornodename,'.mat');
    save(filename, 'basevolt')  
    clear S;  

    % save voltage change
    filename = strcat('deltav_',actornodename,'.mat');
    save(filename, 'voltchange_tensor')
    
    disp('iteration of actornodes');
    disp( actornodename);
    
end











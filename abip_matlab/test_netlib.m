clear;
clc; 
close all;

%% Feasible LP instances
% Probname = {'25FV47', '80BAU3B', 'ADLITTLE', 'AFIRO', 'AGG', 'AGG2', 'AGG3', ...
%     'BANDM', 'BEACONFD', 'BLEND', 'BNL1', 'BNL2', 'BOEING1', 'BOEING2', ...
%     'BORE3D', 'BRANDY', 'CAPRI', 'CRE_A', 'CRE_B', 'CRE_C', 'CRE_D', ...
%     'CYCLE', 'CZPROB', 'D2Q06C', 'D6CUBE', 'DEGEN2', 'DEGEN3', 'DFL001', ...
%     'E226', 'ETAMACRO', 'FFFFF800', 'FINNIS', 'FIT1D', 'FIT1P', 'FIT2D', ...
%     'FIT2P', 'FORPLAN', 'GANGES', 'GFRD_PNC', 'GREENBEA', 'GREENBEB', 'GROW7', ...
%     'GROW15', 'GROW22', 'ISRAEL', 'KB2', 'KEN_7', 'KEN_11', 'KEN_13', ...
%     'KEN_18', 'LOTFI', 'MAROS', 'MAROS_R7', 'MODSZK1', 'NESM', 'OSA_07', ...
%     'OSA_14', 'OSA_30', 'OSA_60', 'PDS_02', 'PDS_06', 'PDS_10', 'PDS_20', ...
%     'PEROLD', 'PILOT', 'PILOT4', 'PILOT87', 'PILOT_JA', 'PILOT_WE', 'PILOTNOV', ...
%     'QAP8', 'QAP12', 'QAP15', 'RECIPE', 'SC50A', 'SC50B', 'SC105', ...
%     'SC205', 'SCAGR7', 'SCAGR25', 'SCFXM1', 'SCFXM2', 'SCFXM3', 'SCORPION', ...
%     'SCRS8', 'SCSD1', 'SCSD6', 'SCSD8', 'SCTAP1', 'SCTAP2', 'SCTAP3', ...
%     'SEBA', 'SHARE1B', 'SHARE2B', 'SHELL', 'SHIP04L', 'SHIP04S', 'SHIP08L', ...
%     'SHIP08S', 'SHIP12L', 'SHIP12S', 'SIERRA', 'STAIR', 'STANDATA', 'STANDGUB', ...
%     'STANDMPS', 'STOCFOR1', 'STOCFOR2', 'STOCFOR3', 'TRUSS', 'TUFF', 'VTP_BASE', ...
%     'WOOD1P', 'WOODW'};

%% Infeasible LP instances
% Probname = {'BGDBG1', 'BGETAM', 'BGINDY', 'BGPRTR', 'BOX1', 'CERIA3D', 'CHEMCOM', ...
%     'CPLEX1', 'CPLEX2', 'EX72A', 'EX73A', 'FOREST6', 'GALENET', 'GOSH', ...
%     'GRAN', 'ITEST2', 'ITEST6', 'KLEIN1', 'KLEIN2', 'KLEIN3', 'MONDOU2', ...
%     'PANG', 'PILOT4I', 'QUAL', 'REACTOR', 'REFINERY', 'VOL1', 'WOODINFE'};

% Probname = {'CAPRI', 'CYCLE', 'GREENBEB', 'MODSZK1', 'PEROLD', 'PILOT4', 'PILOT_JA', ...
%     'PILOT_WE', 'STAIR', 'TUFF', 'VTP_BASE'};

Probname = {'ADLITTLE', 'AFIRO', 'AGG', 'AGG2', 'AGG3', 'BLEND', 'BOEING1', 'BOEING2', 'CAPRI', ...
    'DEGEN2', 'DEGEN3', 'GANGES', 'GFRD_PNC', 'KEN_7', 'RECIPE', 'SC50A', 'SC50B', 'SC105', ...
    'SC205', 'SHELL', 'SHIP04L', 'SHIP04S', 'TUFF'}; 

nprob = length(Probname);

Problist = [1:nprob];

abips = 1;

for di = 1:length(Problist) 
    %% load data
    probID = Problist(di);
    name = Probname{probID};
    load(strcat('../example/netlib/feasible/', Probname{Problist(di)},'.mat'));
    A = Problem.A; 
    b = Problem.b; 
    c = Problem.aux.c; 
    lbounds = Problem.aux.lo; 
    ubounds = Problem.aux.hi; 
    [m, n] = size(A);
    
    %% presolve procedure
    if abips
        [A,b,c,info] = presolve(A,b,c,lbounds,ubounds);
        %% sdpt3, scs and abips procedure
        if ~info.feasible
            fprintf('The problem is infeasible!\n'); 
        else
            % abips setting.
            data.A = sparse(A); 
            data.c = full(c); 
            data.b = full(b);
            params_abips = struct('max_iters', 1000000, 'max_outiters', 100);
            
            % abips implementation.
            tic; [x, y, s, info_abips] = abip_matlab_short(data, params_abips); time_abips = toc; 
            [x_abips, objp_abips] = postsolve(x, info);
        end
    end
       
    %% print procedure
    if abips
        fprintf('%10s & %5d & %5d & %3.2e & %3.2e & %3.2e & %3.2e & %5d & %5d & %3.2e\\\\ \\hline \n', ...
            name, m, n, objp_abips, info_abips.resPri, info_abips.resDual, info_abips.relGap, info_abips.outiter, info_abips.iter, time_abips);
    end
end

function gurobiParam = setGurobiParam(param)
%  The params struct contains Gurobi parameters. A full list may be
%  found on the Parameter page of the reference manual:
%  https://www.gurobi.com/documentation/current/refman/parameter_descriptions.html
%  Parameters must be in TimeLimit not timelimit or timeLimit format, see
%  below for full list that are eligible to be passed to the solver.

if isfield(param,'timelimit')
    param.TimeLimit=param.timelimit;
end

if isfield(param,'multiscale') && param.multiscale==1
    param.ScaleFlag=0;
end

if isfield(param,'lifted') && param.lifted==1
    param.Aggregate=1;
    % Presolve
    % Controls the presolve level
    % Type:	int
    % Default value:	-1
    % Minimum value:	-1
    % Maximum value:	2
    % Controls the presolve level. A value of -1 corresponds to an automatic setting.
    % Other options are off (0), conservative (1), or aggressive (2). More aggressive application
    % of presolve takes more time, but can sometimes lead to a significantly tighter model.
    param.Presolve=0;
end

%backward compatibility
if isfield(param,'method')
    if isempty(param.method)
        param = rmfield(param,'method');
    else
        if ~isfield(param,[lower(param.problemType) 'method'])
            param.([lower(param.problemType) 'method'])=param.method;
        end
    end
end

% https://www.gurobi.com/documentation/current/refman/method.html
% params.method gives the method used to solve continuous models
% -1=automatic,
%  0=primal simplex,
%  1=dual simplex,
%  2=barrier,
%  3=concurrent,
%  4=deterministic concurrent
% i.e. params.method     = 1;          % use dual simplex method
if isfield(param,'lpmethod')
    %gurobiAlgorithms = {'AUTOMATIC','PRIMAL','DUAL','BARRIER','CONCURRENT','CONCURRENT_DETERMINISTIC'};
    % -1=automatic,
    % 0=primal simplex,
    % 1=dual simplex,
    % 2=barrier,
    % 3=concurrent,
    % 4=deterministic concurrent
    switch param.lpmethod
        case 'AUTOMATIC'
            param.Method = -1;
        case 'PRIMAL'
            param.Method = 0;
        case 'DUAL'
            param.Method = 1;
        case 'BARRIER'
            param.Method = 2;
        case 'CONCURRENT'
            param.Method = 3;
        case 'DETERMINISTIC_CONCURRENT'
            param.Method = 4;
        otherwise
            %https://www.gurobi.com/documentation/current/refman/method.html
            %Concurrent methods aren't available for QP and QCP.
            warning([param.lpmethod ' is an unrecognised param.qpmethod for gurobi'])
    end
    param = rmfield(param,'lpmethod');
end

% param.qpmethod gives the qpmethod used to solve continuous models
% -1=automatic,
%  0=primal simplex,
%  1=dual simplex,
%  2=barrier,
%  3=concurrent,
%  4=deterministic concurrent
% i.e. param.qpmethod     = 1;          % use dual simplex method
if isfield(param,'qpmethod')
    %gurobi QP algorithms = {'AUTOMATIC','PRIMAL','DUAL','BARRIER'};
    % -1=automatic,
    % 0=primal simplex,
    % 1=dual simplex,
    % 2=barrier,
    switch param.qpmethod
        case 'AUTOMATIC'
            param.Method = -1;
        case 'PRIMAL'
            param.Method = 0;
        case 'DUAL'
            param.Method = 1;
        case 'BARRIER'
            param.Method = 2;
        otherwise
            %https://www.gurobi.com/documentation/current/refman/method.html
            %Concurrent methods aren't available for QP and QCP.
            warning([param.qpmethod ' is an unrecognised param.qpmethod for gurobi'])
    end
    param = rmfield(param,'qpmethod');
end

if ~isfield(param,'OutputFlag')
    switch param.printLevel
        case 0
            param.OutputFlag = 0;
        case 1
            param.OutputFlag = 0;
        otherwise
            % silent
            param.OutputFlag = 0;
    end
end

if ~isfield(param,'DisplayInterval')
    switch param.printLevel
        case 0
            param.DisplayInterval = 1;
        case 1
            param.DisplayInterval = 1;
        otherwise
            % silent
            param.DisplayInterval = 1;
    end
end

if ~isfield(param,'FeasibilityTol')
    % Primal feasibility tolerance
    % Type:	double
    % Default value:	1e-6
    % Minimum value:	1e-9
    % Maximum value:	1e-2
    % All constraints must be satisfied to a tolerance of FeasibilityTol.
    % Tightening this tolerance can produce smaller constraint violations, but for
    % numerically challenging models it can sometimes lead to much larger iteration counts.
     param.FeasibilityTol = param.feasTol;
end
if ~isfield(param,'OptimalityTol')
    % Dual feasibility tolerance
    % Type:	double
    % Default value:	1e-6
    % Minimum value:	1e-9
    % Maximum value:	1e-2
    % Reduced costs must all be smaller than OptimalityTol in the improving direction in order for a model to be declared optimal.
     param.OptimalityTol = param.optTol;
end


% Permitted parameter fields
permittedFields = {...
'AggFill'...
'Aggregate'...
'BarConvTol'...
'BarCorrectors'...
'BarHomogeneous'...
'BarIterLimit'...
'BarOrder'...
'BarQCPConvTol'...
'BestBdStop'...
'BestObjStop'...
'BQPCuts'...
'BranchDir'...
'CliqueCuts'...
'CloudAccessID'...
'CloudHost'...
'CloudSecretKey'...
'CloudPool'...
'ComputeServer'...
'ConcurrentJobs'...
'ConcurrentMethod'...
'ConcurrentMIP'...
'ConcurrentSettings'...
'CoverCuts'...
'Crossover'...
'CrossoverBasis'...
'CSAPIAccessID'...
'CSAPISecret'...
'CSAppName'...
'CSAuthToken'...
'CSBatchMode'...
'CSClientLog'...
'CSGroup'...
'CSIdleTimeout'...
'CSManager'...
'CSPriority'...
'CSQueueTimeout'...
'CSRouter'...
'CSTLSInsecure'...
'CutAggPasses'...
'Cutoff'...
'CutPasses'...
'Cuts'...
'DegenMoves'...
'Disconnected'...
'DisplayInterval'...
'DistributedMIPJobs'...
'DualReductions'...
'FeasibilityTol'...
'FeasRelaxBigM'...
'FlowCoverCuts'...
'FlowPathCuts'...
'FuncPieceError'...
'FuncPieceLength'...
'FuncPieceRatio'...
'FuncPieces'...
'FuncMaxVal'...
'FuncNonlinear'...
'GomoryPasses'...
'GUBCoverCuts'...
'Heuristics'...
'IgnoreNames'...
'IISMethod'...
'ImpliedCuts'...
'ImproveStartGap'...
'ImproveStartNodes'...
'ImproveStartTime'...
'InfProofCuts'...
'InfUnbdInfo'...
'InputFile'...
'IntegralityFocus'...
'IntFeasTol'...
'IterationLimit'...
'JobID'...
'JSONSolDetail'...
'LazyConstraints'...
'LicenseID'...
'LiftProjectCuts'...
'LPWarmStart'...
'LogFile'...
'LogToConsole'...
'MarkowitzTol'...
'MemLimit'...
'Method'...
'MinRelNodes'...
'MIPFocus'...
'MIPGap'...
'MIPGapAbs'...
'MIPSepCuts'...
'MIQCPMethod'...
'MIRCuts'...
'MixingCuts'...
'ModKCuts'...
'MultiObjMethod'...
'MultiObjPre'...
'MultiObjSettings'...
'NetworkAlg'...
'NetworkCuts'...
'NLPHeur'...
'NodefileDir'...
'NodefileStart'...
'NodeLimit'...
'NodeMethod'...
'NonConvex'...
'NoRelHeurTime'...
'NoRelHeurWork'...
'NormAdjust'...
'NumericFocus'...
'OBBT'...
'ObjNumber'...
'ObjScale'...
'OptimalityTol'...
'OutputFlag'...
'PartitionPlace'...
'PerturbValue'...
'PoolGap'...
'PoolGapAbs'...
'PoolSearchMode'...
'PoolSolutions'...
'PreCrush'...
'PreDepRow'...
'PreDual'...
'PreMIQCPForm'...
'PrePasses'...
'PreQLinearize'...
'Presolve'...
'PreSOS1BigM'...
'PreSOS1Encoding'...
'PreSOS2BigM'...
'PreSOS2Encoding'...
'PreSparsify'...
'ProjImpliedCuts'...
'PSDCuts'...
'PSDTol'...
'PumpPasses'...
'QCPDual'...
'Quad'...
'Record'...
'ResultFile'...
'RINS'...
'RelaxLiftCuts'...
'RLTCuts'...
'ScaleFlag'...
'ScenarioNumber'...
'Seed'...
'ServerPassword'...
'ServerTimeout'...
'Sifting'...
'SiftMethod'...
'SimplexPricing'...
'SoftMemLimit'...
'SolutionLimit'...
'SolutionTarget'...
'SolFiles'...
'SolutionNumber'...
'StartNodeLimit'...
'StartNumber'...
'StrongCGCuts'...
'SubMIPCuts'...
'SubMIPNodes'...
'Symmetry'...
'Threads'...
'TimeLimit'...
'TokenServer'...
'TSPort'...
'TuneBaseSettings'...
'TuneCleanup'...
'TuneCriterion'...
'TuneDynamicJobs'...
'TuneJobs'...
'TuneMetric'...
'TuneOutput'...
'TuneResults'...
'TuneTargetMIPGap'...
'TuneTargetTime'...
'TuneTimeLimit'...
'TuneTrials'...
'TuneUseFilename'...
'UpdateMode'...
'UserName'...
'VarBranch'...
'WLSAccessID'...
'WLSSecret'...
'WLSToken'...
'WLSTokenDuration'...
'WLSTokenRefresh'...
'WorkerPassword'...
'WorkerPool'...
'WorkLimit'...
'ZeroHalfCuts'...
'ZeroObjNodes'};

fields = fieldnames(param);
fieldsToRemove = setdiff(fields, permittedFields);
%remove fields that should not be present
gurobiParam = rmfield(param, fieldsToRemove);

if isfield(gurobiParam,'logFile') && isempty(gurobiParam.logFile)
    gurobiParam = rmfield(gurobiParam,'logFile');
end
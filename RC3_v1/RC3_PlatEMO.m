classdef RC3_PlatEMO < PROBLEM
%RC3_PLATEMO  PlatEMO problem class of the distribution-network model built by
%RC3_setup (decision vector, MP-Opt-Model costs and constraints).
%
%   P = RC3_setup(24, opts);
%   [Dec, Obj, Con] = platemo('problem', {@RC3_PlatEMO, P}, ...
%                             'algorithm', @NSGAII, 'N', 100, 'maxFE', 30000);
%
%   Requires PlatEMO 4.x (mixed encoding: 1 real, 2 integer, 4 binary).
%   M = 1 for the weighted single objective (opts.m = 1), otherwise the
%   number of user objectives.
%
%   The evaluation is identical to the built-in solver (RC3_evaluate):
%     Initialization  rng(P.seed) once (platemo() calls rng('shuffle')),
%                     then the common initial population RC3_init_pop
%     CalDec          RC3_repair (stored: Lamarckian repair)
%     CalObj, CalCon  RC3_fitness (MP-Opt-Model costs and constraints)
    properties
        Data            % problem data from RC3_setup
        nInit = 0       % number of calls of Initialization
    end
    methods
        function Setting(obj)
            Pd = obj.ParameterSet(struct());
            if ~isfield(Pd, 'D')
                error('RC3_PlatEMO: pass the struct from RC3_setup, e.g. {@RC3_PlatEMO, P}.');
            end
            obj.Data     = Pd;
            obj.M        = Pd.nObj;
            obj.D        = Pd.D;
            obj.lower    = Pd.lower;
            obj.upper    = Pd.upper;
            obj.encoding = Pd.enc;
        end
        function Population = Initialization(obj, N)
            if nargin < 2, N = obj.N; end
            obj.nInit = obj.nInit + 1;
            if obj.nInit == 1
                rng(obj.Data.seed, 'twister');      % same random stream as RC3_nsga2
            end
            Population = obj.Evaluation(RC3_init_pop(obj.Data, N));
        end
        function PopDec = CalDec(obj, PopDec)
            PopDec = RC3_repair(PopDec, obj.Data);
        end
        function PopObj = CalObj(obj, PopDec)
            PopObj = zeros(size(PopDec, 1), obj.M);
            for i = 1:size(PopDec, 1)
                PopObj(i, :) = RC3_fitness(PopDec(i, :), obj.Data);
            end
        end
        function PopCon = CalCon(obj, PopDec)
            PopCon = zeros(size(PopDec, 1), obj.Data.nCon);
            for i = 1:size(PopDec, 1)               % objectives come from the cache
                [~, PopCon(i, :)] = RC3_fitness(PopDec(i, :), obj.Data);
            end
        end
    end
end

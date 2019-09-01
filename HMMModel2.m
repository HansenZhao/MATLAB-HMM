classdef HMMModel2 < handle
    
    properties(SetAccess=private)
        transProb,
        emitProb,
        initProb,
        cachedAlpha,
        cachedBeta,
        isCache,
    end
    
    properties(Dependent)
        nStates,
        nObserve,
    end
    
    properties(Constant)
        ITER_TOLERENCE = 1e-4;
        DEBUG = 1;
    end
    
    methods
        function obj = HMMModel2(nstate,nobserve)
            f = @(x)x./sum(x,2);
            obj.transProb = f(rand(nstate));
            obj.emitProb = f(rand(nstate,nobserve));
            obj.initProb = f(rand(1,nstate));
            obj.isCache = 0;
        end
        
        function ns = get.nStates(obj)
            ns = size(obj.transProb,1);
        end
        
        function no = get.nObserve(obj)
            no = size(obj.emitProb,2);
        end
        
        function setModel(obj,option,mat)
            if any(abs(sum(mat,2)-1)>1e-3)
                error('HMMModel: sum of probability should close to 1');
            end
            % option: trans,emit,init
            switch option
                case 'trans'
                    if all(size(mat)==[obj.nStates,obj.nStates])
                        obj.transProb = mat;
                    else
                        error('HMMModel: trans probability should be in size of %dx%d, given %dx%d',...
                            obj.nStates,obj.nStates,size(mat,1),size(mat,2));
                    end
                case 'emit'
                    if all(size(mat)==[obj.nStates,obj.nObserve])
                        obj.emitProb = mat;
                    else
                        error('HMMModel: emit probability should be in size of %dx%d, given %dx%d',...
                            obj.nStates,obj.nObserve,size(mat,1),size(mat,2));
                    end
                case 'init'
                    if all(size(mat)==[1,obj.nStates])
                        obj.initProb = mat;
                    else
                        error('HMMModel: init probability should be in size of 1x%d, given %dx%d',...
                            obj.nStates,size(mat,1),size(mat,2));
                    end
                otherwise
                    error('HMMModel: setModel by trans, emit or init option, given: %s',option);
            end
        end
        
        function p = calForward(obj,obsVec,t,tState)
            if obj.isCache && (~isempty(obj.cachedAlpha{tState,t}))
                p = obj.cachedAlpha{tState,t};
            else
                if t>1
                    prev = cell(1,obj.nStates);
                    for m = 1:obj.nStates
                        prev{m} = obj.calForward(obsVec,t-1,m);
                        prev{m} = prev{m}.multiValue(obj.transProb(m,tState)).multiValue(obj.emitProb(tState,obsVec(t)));
                    end
                    p = ExtNum.ExtSum(prev);
                else
                    p = ExtNum(obj.initProb(tState) * obj.emitProb(tState,obsVec(1)));
                end
                if obj.isCache
                    obj.cachedAlpha{tState,t} = p;
                end
            end           
        end
        
        function p = calBackward(obj,obsVec,t,tState)
            if obj.isCache && (~isempty(obj.cachedBeta{tState,t}))
                p = obj.cachedBeta{tState,t};
            else
                L = length(obsVec);
                if t<L
                    recurFactor = cell(1,obj.nStates);
                    for m = 1:obj.nStates
                        recurFactor{m} = obj.calBackward(obsVec,t+1,m);
                        recurFactor{m} = recurFactor{m}.multiValue(obj.emitProb(m,obsVec(t+1))).multiValue(obj.transProb(tState,m));
                    end
                    p = ExtNum.ExtSum(recurFactor);
                else
                    p = ExtNum(1);
                end
                if obj.isCache
                    obj.cachedBeta{tState,t} = p;
                end
            end
        end
        
        function p = calObserveProb(obj,obsVec,method)
            % calculate P(O|lameda), method: forward, backward
            if ~exist('method','var')
                method = 'forward';
            end
            L = length(obsVec);
            switch method
                case 'forward'
                    p = cell(obj.nStates,1);
                    for m = 1:obj.nStates
                        p{m} = obj.calForward(obsVec,L,m);
                    end
                    p = ExtNum.ExtSum(p);
                case 'backward'
                    p = cell(1,obj.nStates);
                    for m = 1:obj.nStates
                        p{m} = obj.calBackward(obsVec,1,m);
                        p{m} = p{m}.multiValue(obj.initProb(m)).multiValue(obj.emitProb(m,obsVec(1)));
                    end
                    p = ExtNum.ExtSum(p);
                otherwise
                    error('HMMModel: cal observe probability only accept option: forward,backward. Given: %s',method);
            end
        end
        
        function p = calGammaProb(obj,obsVec,t,tState)
            forwardVec = cell(1,obj.nStates);
            backwardVec = cell(obj.nStates,1);
            for m = 1:obj.nStates
                forwardVec{m} = obj.calForward(obsVec,t,m);
                backwardVec{m} = obj.calBackward(obsVec,t,m);
            end
            p = (forwardVec{tState}*backwardVec{tState})/ExtNum.cellVecTimes(forwardVec,backwardVec);
        end
        
        function p = calKasaiProb(obj,obsVec,t,curState,nextState)
            probMat = cell(obj.nStates);
            alpha_t = cell(obj.nStates,1);
            beta_t_1 = cell(obj.nStates,1);
            for m = 1:obj.nStates
                alpha_t{m} = obj.calForward(obsVec,t,m);
                beta_t_1{m} = obj.calBackward(obsVec,t+1,m);
            end
            for m = 1:obj.nStates
                for n = 1:obj.nStates
                    p = alpha_t{m}.multiValue(obj.transProb(m,n)).multiValue(obj.emitProb(n,obsVec(t+1)));
                    probMat{m,n} = p*beta_t_1{n};
                end
            end
            p = probMat{curState,nextState}/ExtNum.ExtSum(probMat(:));
        end
        
        function BWFit(obj,obsVec)
            L = length(obsVec);
            f = @(x)x./sum(x,2);
            tmpTransProb = f(rand(obj.nStates));
            tmpEmitProb = f(rand(obj.nStates,obj.nObserve));
            tmpInitProb = f(rand(1,obj.nStates));
            obj.enableCache(obsVec);
            
            while obj.getModelDiff(tmpTransProb,tmpEmitProb,tmpInitProb) > HMMModel.ITER_TOLERENCE
                obj.transProb = tmpTransProb;
                obj.emitProb = tmpEmitProb;
                obj.initProb = tmpInitProb;
                
                gammaCache = cell(obj.nStates,L);
                kasaiCache = cell(obj.nStates,obj.nStates,L-1);
                
                obj.refreshCache(obsVec);
                for t = 1:L
                    for s = 1:obj.nStates
                        gammaCache{s,t} = obj.calGammaProb(obsVec,t,s);
                        if t < L
                            for m = 1:obj.nStates
                                kasaiCache{s,m,t} = obj.calKasaiProb(obsVec,t,s,m);
                            end
                        end
                    end
                    if HMMModel.DEBUG && mod(t,500)==0
                        fprintf(1,'cache parameters: %d/%d\n',t,L);
                    end
                end
                [tmpTransProb,tmpEmitProb,tmpInitProb] = obj.calNewModelParam(obsVec,gammaCache,kasaiCache);
            end
            
            obj.transProb = tmpTransProb;
            obj.emitProb = tmpEmitProb;
            obj.initProb = tmpInitProb;
            obj.disableCache();
            obj.refreshCache(obsVec);
        end
        
        function [obs,states] = simulateObserve(obj,steps)
            [obs,states] = deal(zeros(steps,1));
            states(1) = genIndexByProb(obj.initProb);
            obs(1) = genIndexByProb(obj.emitProb(states(1),:));
            for m = 2:steps
                states(m) = genIndexByProb(obj.transProb(states(m-1),:));
                obs(m) = genIndexByProb(obj.emitProb(states(m),:));
            end
        end
        
        function disableCache(obj)
            obj.isCache = 0;
        end
        
        function enableCache(obj,obs)
            obj.isCache = 1;
            obj.refreshCache(obs);
        end
        
        function I = vPredict(obj,obsVec)
            L = length(obsVec);
            delta = cell(obj.nStates,L);
            phi = zeros(obj.nStates,L);
            for m = 1:obj.nStates
                delta{m,1} = ExtNum(obj.initProb(m)*obj.emitProb(m,obsVec(1)));
            end
            for t = 2:L
                for m = 1:obj.nStates
                    tmp = cell(1,obj.nStates);
                    for n = 1:obj.nStates
                        tmp{n} = delta{n,t-1}.multiValue(obj.transProb(n,m));
                    end
                    [v,I] = ExtNum.ExtMax(tmp);
                    delta{m,t} = v.multiValue(obj.emitProb(m,obsVec(t)));
                    phi(m,t) = I;
                end
            end
            I = zeros(L,1);
            [~,tmp] = ExtNum.ExtMax(delta(:,t));
            I(L) = tmp;
            for t = fliplr(1:(L-1))
                I(t) = phi(I(t+1),t+1);
            end
        end
    end
    
    methods(Access=private)
        function d = getModelDiff(obj,trans,emit,init)
            diffMat = abs([trans;emit';init]-[obj.transProb;obj.emitProb';obj.initProb]);
            d = max(diffMat(:));
            if HMMModel.DEBUG
                fprintf(1,'model differences: %.4f\n',d);
            end
        end
        
        function [trans,emit,init] = calNewModelParam(obj,obsVec,gammaCache,kasaiCache)
            init = zeros(1,obj.nStates);
            trans = zeros(obj.nStates);
            emit = zeros(obj.nStates,obj.nObserve);
            for m = 1:obj.nStates
                init(m) = gammaCache{m,1}.value;
                for n = 1:obj.nStates
                    p = ExtNum.ExtSum(kasaiCache(m,n,:))/ExtNum.ExtSum(gammaCache(m,1:(end-1)));
                    trans(m,n) = p.value;
                end
                for k = 1:obj.nObserve
                    p = ExtNum.ExtSum(gammaCache(m,obsVec==k))/ExtNum.ExtSum(gammaCache(m,:));
                    emit(m,k) = p.value;
                end
            end
        end
        
        function refreshCache(obj,obsVec)
            L = length(obsVec);
            obj.cachedAlpha = cell(obj.nStates,L);
            obj.cachedBeta = cell(obj.nStates,L);
        end
    end
end


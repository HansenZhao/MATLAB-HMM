classdef HMMModel < handle
    
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
        function obj = HMMModel(nstate,nobserve)
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
            if obj.isCache && obj.cachedAlpha(tState,t)>=0
                p = obj.cachedAlpha(tState,t);
            else
                if t>1
                    prev = zeros(1,obj.nStates);
                    for m = 1:obj.nStates
                        prev(m) = obj.calForward(obsVec,t-1,m);
                    end
                    p = prev*obj.transProb(:,tState)*obj.emitProb(tState,obsVec(t));
                else
                    p = obj.initProb(tState) * obj.emitProb(tState,obsVec(1));
                end
                if obj.isCache
                    obj.cachedAlpha(tState,t) = p;
                end
            end           
        end
        
        function p = calBackward(obj,obsVec,t,tState)
            if obj.isCache && obj.cachedBeta(tState,t)>=0
                p = obj.cachedBeta(tState,t);
            else
                L = length(obsVec);
                if t<L
                    transFactor = obj.transProb(tState,:);
                    emitFactor = obj.emitProb(:,obsVec(t+1))';
                    recurFactor = zeros(1,obj.nStates);
                    for m = 1:obj.nStates
                        recurFactor(m) = obj.calBackward(obsVec,t+1,m);
                    end
                    p = sum(transFactor.*emitFactor.*recurFactor);
                else
                    p = 1;
                end
                if obj.isCache
                    obj.cachedBeta(tState,t) = p;
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
                    p = zeros(obj.nStates,1);
                    for m = 1:obj.nStates
                        p(m) = obj.calForward(obsVec,L,m);
                    end
                    p = sum(p);
                case 'backward'
                    p = zeros(1,obj.nStates);
                    for m = 1:obj.nStates
                        p(m) = obj.calBackward(obsVec,1,m);
                    end
                    p = sum(obj.initProb.*obj.emitProb(:,obsVec(1))'.*p);
                otherwise
                    error('HMMModel: cal observe probability only accept option: forward,backward. Given: %s',method);
            end
        end
        
        function p = calGammaProb(obj,obsVec,t,tState)
            forwardVec = zeros(1,obj.nStates);
            backwardVec = zeros(obj.nStates,1);
            for m = 1:obj.nStates
                forwardVec(m) = obj.calForward(obsVec,t,m);
                backwardVec(m) = obj.calBackward(obsVec,t,m);
            end
            p = (forwardVec(tState)*backwardVec(tState))/(forwardVec*backwardVec);
        end
        
        function p = calKasaiProb(obj,obsVec,t,curState,nextState)
            probMat = zeros(obj.nStates);
            for m = 1:obj.nStates
                alpha_t = obj.calForward(obsVec,t,m);
                transFactor = obj.transProb(m,:); %nstate x 1
                emitFactor = obj.emitProb(:,obsVec(t+1))'; %nstate x 1
                backwardFactor = zeros(1,obj.nStates);
                for n = 1:obj.nStates
                    backwardFactor(n) = obj.calBackward(obsVec,t+1,n);
                end
                probMat(m,:) = alpha_t * (transFactor.*emitFactor.*backwardFactor);
            end
            p = probMat(curState,nextState)/sum(probMat(:));
        end
        
        function BWFit(obj,obsVec)
            L = length(obsVec);
            f = @(x)x./sum(x,2);
            tmpTransProb = f(rand(obj.nStates));
            tmpEmitProb = f(rand(obj.nStates,obj.nObserve));
            tmpInitProb = f(rand(1,obj.nStates));
            obj.isCache = 1;
            
            while obj.getModelDiff(tmpTransProb,tmpEmitProb,tmpInitProb) > HMMModel.ITER_TOLERENCE
                obj.transProb = tmpTransProb;
                obj.emitProb = tmpEmitProb;
                obj.initProb = tmpInitProb;
                
                gammaCache = zeros(obj.nStates,L);
                kasaiCache = zeros(obj.nStates,obj.nStates,L-1);
                
                obj.refreshCache(obsVec);
                for t = 1:L
                    for s = 1:obj.nStates
                        gammaCache(s,t) = obj.calGammaProb(obsVec,t,s);
                        if t < L
                            for m = 1:obj.nStates
                                kasaiCache(s,m,t) = obj.calKasaiProb(obsVec,t,s,m);
                            end
                        end
                    end
                    if HMMModel.DEBUG && mod(t,20)==0
                        fprintf(1,'cache parameters: %d/%d\n',t,L);
                    end
                end
                [tmpTransProb,tmpEmitProb,tmpInitProb] = obj.calNewModelParam(obsVec,gammaCache,kasaiCache);
            end
            
            obj.transProb = tmpTransProb;
            obj.emitProb = tmpEmitProb;
            obj.initProb = tmpInitProb;
            obj.isCache = 0;
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
            init = gammaCache(:,1)';
            trans = zeros(obj.nStates);
            emit = zeros(obj.nStates,obj.nObserve);
            for m = 1:obj.nStates
                for n = 1:obj.nStates
                    trans(m,n) = sum(kasaiCache(m,n,:))/sum(gammaCache(m,1:(end-1)));
                end
                for k = 1:obj.nObserve
                    emit(m,k) = sum(gammaCache(m,obsVec==k))/sum(gammaCache(m,:));
                end
            end
        end
        
        function refreshCache(obj,obsVec)
            L = length(obsVec);
            obj.cachedAlpha = -1*ones(obj.nStates,L);
            obj.cachedBeta = -1*ones(obj.nStates,L);
        end
    end
end


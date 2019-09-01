function [ I ] = genIndexByProb( prob,indexVec,num )
    % tic; I = genIndexByProb([0.35,0.45,0.2],[1,3,4],10000); toc;
    % histogram(I,'Normalization','Probability')
    L = length(prob);
    if ~exist('indexVec','var')
        indexVec = 1:L;
    end
    if ~exist('num','var')
        num = 1;
    end
    tmp = rand(num,1);
    prob = cumsum(prob(:)'/sum(prob));
    I = sum(tmp>prob,2)+1;
    I = indexVec(I);
end


function [pathways, scores] = prepare_FVA(pathways, scores, nRequested)

    [scores, idx] = sort(scores,'ascend');
    pathways = pathways(idx);

    n = min(nRequested, numel(scores));

    pathways = pathways(1:n);
    scores = scores(1:n);

end
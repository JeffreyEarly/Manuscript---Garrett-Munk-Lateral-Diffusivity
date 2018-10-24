function [xMean,yMean,yStdErr] = BinDataWithErrorBars(x,y,xBins)
%BinDataWithErrorBars 
    [~,edges,bin] = histcounts(x,xBins);
    yMean = zeros(length(edges)-1,1);
    yStdErr = zeros(length(edges)-1,1);
    for i=1:max(bin)
        yMean(i) = mean(y(bin==i));
        yStdErr(i) = std(y(bin==i))/sqrt(sum(bin==i));
    end
    xMean = ((edges(1:end-1)+edges(2:end))/2)';
end


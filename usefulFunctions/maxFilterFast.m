function maxdata=maxFilterFast(B,win)

startPadding=max(B(1:win))*ones(1,floor(win/2));
endPadding=max(B(end-win+1:end))*ones(1,win-floor(win/2)/1);

maxdata = zeros(1,length(B)-win+1);
    m = max(B(1:win));
    for k=win+1:length(B)
        maxdata(k-win) = m;
        if B(k-win) < m
            m = max(m, B(k));
        else
            % Matt Fig's improvement
            m = B(1+k-win);
            for ii = 1+k-win+1:k
                if B(ii)>m
                    m = B(ii);
                end
            end
            % m = max(B(1+k-win:k));
        end
    end
    maxdata(end) = m;
    maxdata=[startPadding maxdata  endPadding];
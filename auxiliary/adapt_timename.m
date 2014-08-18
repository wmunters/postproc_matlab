function timename = adapt_timename(i, decimals)

    timename = num2str(i);

    if(length(timename)==1) 
        % No decimals
        timename = strcat(timename,'.');
    end

    while(length(timename)<decimals+2)
        timename = strcat(timename,'0');
    end
    

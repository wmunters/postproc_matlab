function Ihat = cyclicshift(fhat)

    N = size(fhat,1);


    for l=N/2:N-1
        for m=N/2:N-1
            lindex = l+1;
            mindex = m+1;
            Ihat(lindex,mindex) = fhat(lindex-N/2,mindex-N/2);
        end
    end

    for l=0:N/2-1
        for m=N/2:N-1
            lindex = l+1;
            mindex = m+1;
            Ihat(lindex,mindex) = fhat(lindex+N/2,mindex-N/2);
        end
    end

    for l=N/2:N-1
        for m=0:N/2-1
            lindex = l+1;
            mindex = m+1;
            Ihat(lindex,mindex) = fhat(lindex-N/2,mindex+N/2);
        end
    end

    for l=0:N/2-1
        for m=0:N/2-1
            lindex = l+1;
            mindex = m+1;
            Ihat(lindex,mindex) = fhat(lindex+N/2,mindex+N/2);
        end
    end

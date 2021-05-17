function mat2dss(h, A)
    fprintf(h, '(');
    np = size(A,1);

    for r = 1:np
        if r > 1
            fprintf(h, ' | ');
        end

        for c = 1:r
            if c > 1
                fprintf(h, ' ');
            end

            fprintf(h, '%0.6f', A(r,c));
        end
    end

    fprintf(h, ')');
 end
function ps = phases2dss(pp)
    ps = '';

    if pp(1) == 1
        ps = '1';
    end

    if pp(2) == 1
        if length(ps) >= 1
            ps = '1.2';
        else
            ps = '2';
        end
    end

    if pp(3) == 1
        if length(ps) >= 1
            ps = [ps '.3'];
        else
            ps = '3';
        end
    end
end
function stype = parse_type(ttype)
    switch ttype
        case 'double'
            stype = 8;
        case 'single'
            stype = 4;
        case 'uint32'
            stype = 4;
        case 'uint16'
            stype = 2;
        otherwise
            stype = 1;
    end
end
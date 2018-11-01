function flag = judge_file(filename, msg)
    flag = false;
    if exist(filename, 'file')
        w = input(msg, 's');
        switch w
            case 'y'
                flag = true;
        end
    else
        flag = true;
    end
end
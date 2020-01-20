function nbatch = batch_compute(nsize)
    if ispc
        memo = pc_mem / 2;
    elseif ismac
        memo = mac_mem / 2;
    elseif isunix
        memo = unix_mem / 2;
    end
    nbatch = ceil(nsize / memo);
end

function mem = pc_mem()
    [~, mem] = memory;
    mem = mem.PhysicalMemory.Available;
end

function mem = unix_mem()
    [~, out] = system('cat /proc/meminfo |grep MemAvailable');
    mem = sscanf(out, 'MemAvailable: %f'); % in KiB
    mem = mem * 1024; % in bytes
end

function mem = mac_mem()
    [~, m] = unix('vm_stat | grep free');
    spaces = strfind(m, ' ');
    mem = str2num(m(spaces(end):end)) * 4096;
end

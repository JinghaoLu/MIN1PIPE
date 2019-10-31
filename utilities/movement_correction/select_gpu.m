function select_gpu
    n = gpuDeviceCount;
    if n>0
        try
            s = zeros(1, n);
            for i = 1: n
                d = gpuDevice(i);
                s(i) = d.AvailableMemory;
            end
            [~, id] = max(s);
            gpuDevice(id);
            disp(['Using GPU Device #', num2str(id), ' Total Memory: ', num2str(s(id))])
        catch
            disp('GPU Device available but outdated.')
        end
    else
        disp('No GPU Device available.')
    end
end


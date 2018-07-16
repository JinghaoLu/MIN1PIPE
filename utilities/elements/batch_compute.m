function nbatch = batch_compute(nsize)
    [~, memo] = memory;
    memo = memo.PhysicalMemory.Available / 2;
    nbatch = ceil(nsize / memo);
end
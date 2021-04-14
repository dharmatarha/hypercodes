%timingInfo_runs34_wrapper -- runs the parseEPI functions for each pair

for PAIR = 2:9
    
    % listening task
    parseEPI_listeningTask(PAIR);
    
    % reading task
    parseEPI_readingTask(PAIR);
    
end
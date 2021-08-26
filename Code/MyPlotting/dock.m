function dock()

    figHandles = findall(0,'Type','figure');
    for n=1:length(figHandles),
        set(figHandles(n), 'WindowStyle', 'docked');
    end
    
end
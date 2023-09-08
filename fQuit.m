function [value,isTerminal,direction] = fQuit(t,y)
    try value = y(1) - pi;
    catch
        keyboard
    end
    isTerminal = true;
    direction = 1;
end
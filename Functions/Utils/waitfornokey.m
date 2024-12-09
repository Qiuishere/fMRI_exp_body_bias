function [] = waitfornokey
 % wait for key release, used after a key press collection to make sure key
 % press is released before the next scripts are executed
while 1
    if IsOSX
        [a,~,~] = KbCheck(-1);
    else
        [a,~,~] = KbCheck;
    end
    if a ==0
        break
    end
end
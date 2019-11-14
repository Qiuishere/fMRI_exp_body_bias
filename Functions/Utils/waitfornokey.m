function [] = waitfornokey

while 1
    if IsOSX
        [a,~,~] = KbCheck(-1);
    else
        [a,~,~] = KbCheck;
    end
    if a ==0;
        break
    end
end
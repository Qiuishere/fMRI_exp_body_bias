clear
filedir = listdir('./Vanrie_mat/*.mat');
nPoint = 13;

for themov = 1:length(filedir)

    MOV_resamp = zeros(3,nPoint*60);

load(['./Vanrie_mat/' filedir{themov}])
Nframe = size(MOV,2)/nPoint;
for j = 1:3
    for i = 1:nPoint
        v = [MOV(j,i:13:end) MOV(j,i)];
        xq = 1:0.0002: (Nframe+1);% Interpolate at the query points, and specify cubic interpolation.
        Trajectory_long  = round(interp1(1:Nframe+1,v ,xq,'spline'));
        MOV_resamp(j, i:13:end) = Trajectory_long(2:2500:end);
    end
end
    for i = 1:size(MOV,1)
        MOV(i,:) = MOV(i,:) - mean(MOV(i,:));
    end
MOV = RotationMatrix(pi/2, [1 0 0]);
    
newName = ['Vanrie_60fr/' filedir{themov}]
    save(newName, 'MOV', 'Nframe')
    cd ..
    clear MOV
end

fclose('all');

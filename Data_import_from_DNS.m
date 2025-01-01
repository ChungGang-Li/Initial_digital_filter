clear all
clc

% https://turbulence.oden.utexas.edu/  

Re_tau=180

% 讀取數據
data = load("Turbulence_intensity.dat"); % 如果文件是純文本數據，使用 load

% 提取 Y 座標和其他列數據
Y_original = data(:, 1); % 假設第一列是 Y 座標
Data = data(:, 2:end); % 其他列為數據

% 定義等間隔數目
num_points = 40; % 修改此處設置等間隔點的數目

% 定義新的等間隔 Y
Y_equal = linspace(min(Y_original), max(Y_original), num_points);

% 內插其他列的數據
DI = zeros(size(Data, 2), num_points); % 預分配空間
for i = 1:size(Data, 2)
    DI(i, :) = interp1(Y_original, Data(:, i), Y_equal, 'linear');
end

Y_plus = DI(1, :);

for n = 1:num_points
    
    if (Y_plus(n) < 10)
        
        U(n) = Y_plus(n);
        
    else
        
        U(n) = (1/0.41*log(Y_plus(n))+5.5);
        
    end
    
end



% 合併結果 R11 R21 R22 R31 R32 R33
Half = [DI(1,:)'/(2*Re_tau), U', DI(2,:)', DI(5,:)', DI(3,:)', DI(6,:)', DI(7,:)', DI(4,:)'];

additional_rows = Half(end-1:-1:1, :);

Results = [Half; additional_rows];

Results(1:2*num_points-1,1) = linspace(0, 1, 2*num_points-1);

% Results(:,3) = sqrt(Results(:,3))

Results(num_points+1:2*num_points,4) = -Results(num_points:-1:1,4);









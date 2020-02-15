%                 *    (      (   (       )                         
%           (   (  `   )\ )   )\ ))\ ) ( /(             (    *   )  
%         ( )\  )\))( (()/(  (()/(()/( )\())   (  (     )\ ` )  /(  
%         )((_)((_)()\ /(_))  /(_))(_)|(_)\    )\ )\  (((_) ( )(_)) 
%        ((_)_ (_()((_|_))_| (_))(_))   ((_)  ((_|(_) )\___(_(_())  
%         / _ \|  \/  | |_   | _ \ _ \ / _ \ _ | | __((/ __|_   _|  
%        | (_) | |\/| | __|  |  _/   /| (_) | || | _| | (__  | |    
%         \__\_\_|  |_|_|    |_| |_|_\ \___/ \__/|___| \___| |_|                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        scrity by Francesco Morabito                 %%%
%%%                         last edit 29/11                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% importing data %%% 
%%%%%%%%%%%%%%%%%%%%%%

clear; clf;

portfolios = importdata('C:\Users\franc\Documents\portfolios.csv',',');
factors = importdata('C:\Users\franc\Documents\factors.CSV',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ex 1 %%%
%%%%%%%%%%%%

% summary statistics
rf = factors.data(:,5);

summary1 = [summary_func(portfolios.data(:, 2) - rf);...
            summary_func(portfolios.data(:, 6) - rf);... 
            summary_func(portfolios.data(:,22) - rf);... 
            summary_func(portfolios.data(:,26) - rf)];

% plots
portfolios_data = summary_func(portfolios.data(:,2:26) - rf);

subplot(2,2,1)
stem(portfolios_data(1,1:25), "filled")
title("Mean")
subplot(2,2,2)
stem(portfolios_data(1,26:50), "filled")
title("Variance")
subplot(2,2,3)
stem(portfolios_data(1,51:75), "filled")
title("Skewness")
subplot(2,2,4)
stem(portfolios_data(1,76:100), "filled")
title("Kurtosis")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ex 2 %%%
%%%%%%%%%%%%

R_port = portfolios.data(:,2:26);

Y = R_port - rf;
X = [ones(length(portfolios.data),1) factors.data(:,2:4)];

% OLS
B = (inv(X' * X) * X' * Y);
B_prime = B';

index_row = 0;
for index = [1 5 21 25] 
    index_row = index_row + 1;
    summary_beta(index_row,:) = [B_prime(index,1) B_prime(index,2)...
                                 B_prime(index,3) B_prime(index,4)]; 
end

% ploting B'
figure;
hold on
for index = [1 2 3 4]
    plot(B_prime(:,index))
end
hold off
legend('alpha','Mkt-R',' SMB',' HML','NumColumns',2)
legend('boxoff')

% residuals
index_row = 0;
for index = (1:25)
    index_row = index_row + 1;
    E = Y(:,index) - X * B(:,index);  % cell2mat(E_matrix(1,2))
    E_matrix(index_row,:) = [portfolios.textdata(1,index+1) E]; 
end

% variance
betas = 3;
index_row = 0;
for index = 1:length(B_prime)   
    index_row = index_row + 1;
    SSR = sum(cell2mat(E_matrix(index,2)) .^ 2);
    S2_u = 1 / (length(Y) - (betas + 1)) * SSR;
    Var_beta = S2_u * inv(X' * X);
    error = sqrt(diag(Var_beta)');
    e_matrix(index_row,:) = [error(1,1) error(1,2) error(1,3) error(1,4)];
end
summaryE = [e_matrix(1,1)  e_matrix(1,2)  e_matrix(1,3)  e_matrix(1,4);... 
            e_matrix(5,1)  e_matrix(5,2)  e_matrix(5,3)  e_matrix(5,4);...
            e_matrix(21,1) e_matrix(21,2) e_matrix(21,3) e_matrix(21,4);... 
            e_matrix(25,1) e_matrix(25,2) e_matrix(25,3) e_matrix(25,4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ex 3 %%%
%%%%%%%%%%%%

index_col = 0;
for index = 1:length(portfolios.data)
    index_col = index_col + 1;
    Y3 = portfolios.data(index,2:26);
    Y3 = Y3';
    B_prime_new = [ones(length(B_prime),1) B_prime(:,2:4)];
    BX = (inv(B_prime' * B_prime) * B_prime' * Y3);
    RP_matrix(:,index_col) = [BX(1,1) BX(2,1) BX(3,1) BX(4,1)];
end

% summary statistics
for index = 1:4
    if index == 1
        name = "alpha";
    end
    if index == 2
         name = "market";
    end
    if index == 3
         name ="size";
    end
    if index == 4
         name = "value";
    end    
    mean_port = mean(RP_matrix(index,:));
    var_port = var(RP_matrix(index,:));
    ske = skewness(RP_matrix(index,:));
    kurt = kurtosis(RP_matrix(index,:));
    summary3(index,:) = [mean_port var_port ske kurt];
end

% plot
for index = 1:4
    figure;
    plot(RP_matrix(index,:))
    if index == 1
        title("alpha")
       % txt = ["mean" summary3(index,1),...
       %        "variance" summary3(index,2), ...
       %        "skewness" summary3(index,3), ...
       %        "kurtosis" summary3(index,4)];
       % text(4,0.5,txt)
    end
    if index == 2
        title("market")
       % txt = ["mean" summary3(index,1),...
       %        "variance" summary3(index,2), ...
       %        "skewness" summary3(index,3), ...
       %        "kurtosis" summary3(index,4)];
       % text(4,0.5,txt)
    end
    if index == 3
        title("size")
       % txt = ["mean" summary3(index,1),...
       %        "variance" summary3(index,2), ...
       %        "skewness" summary3(index,3), ...
       %        "kurtosis" summary3(index,4)];
       % text(4,0.5,txt)
    end
    if index == 4
        title("value")
       % txt = ["mean" summary3(index,1),...
       %        "variance" summary3(index,2), ...
       %        "skewness" summary3(index,3), ...
       %        "kurtosis" summary3(index,4)];
       %text(4,0.5,txt)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ex 4 %%%
%%%%%%%%%%%%

% creating COV and COVS matrix
start_time1 = sequencer(200801); end_time1 = sequencer(201312);

ret_port = portfolios.data(start_time1:end_time1,2:26)...
            - factors.data(start_time1:end_time1,5);

COV = cov(ret_port);
COVS = (1 - 0.5) * mean(diag(COV)) * eye(25) + 0.5 * COV;

% creating portfolios
w = port_assembler(COV);
ws = port_assembler(COVS);
wEW = 1 / length(COV) * ones(length(COV),1);

% backtest
start_time = sequencer(201401); end_time = sequencer(201808);

ret_portb = portfolios.data(start_time:end_time,2:26)... 
            - rf(start_time:end_time,:);

r_w = w_portfolios(w, ret_portb);
r_ws = w_portfolios(ws, ret_portb);
r_wEW = w_portfolios(wEW, ret_portb);

% cumulative
cumu_r_w = cumulative_builder(r_w);
cumu_r_ws = cumulative_builder(r_ws);
cumu_r_wEW = cumulative_builder(r_wEW);
cumu_benckmark = cumulative_builder(rf(start_time:end_time,:));

% plot
figure;
plot(cumu_r_w)
title("w")
figure;
plot(cumu_r_ws)
title("ws")
figure;
plot(cumu_r_wEW) 
title("wEW")
figure;
plot(cumu_benckmark); 
title("benchmark")

% summary 
summaryRW = [summary_func(r_w); summary_func(r_ws);...
             summary_func(r_wEW); summary_func(rf(start_time:end_time,:))];
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ex 5 %%%
%%%%%%%%%%%%

for i = 1:length(e_matrix)
   test(i,:) = sum(e_matrix(i,1:4)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%  RESULTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("EX 1") 
summary1
fprintf("EX 2")
B_prime([1 5 21 25],:)
summaryE
fprintf("EX 3")
summary3
fprintf("EX 4")
summaryRW

%%%%%%%%%%%%%%%%%%%%%%%%%%  FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cumulative_vector = cumulative_builder(returns)
cummulative = 1;
for i = 1:length(returns)
    cummulative = cummulative * 1 + returns(i,1);
    if i == 1
        cummulative = 1;
    end
    cumu = 1 + returns(i, 1) * cummulative;
    cumulative_vector(i,:) = cumu;
end
end

function summa = summary_func(values)
mean_w_ret = mean(values);
var_w_ret = var(values);
skewness_w_ret = skewness(values);
kurtosis_w_ret = kurtosis(values);
summa = [mean_w_ret var_w_ret skewness_w_ret kurtosis_w_ret];
end

function Rport_vect = w_portfolios(weights, returns)
row_index = 0;
for i = 1:length(returns)
    row_index = row_index + 1;
    R_Port = sum(weights' .* returns(i,:) / 100);
    Rport_vect(row_index,:) = R_Port;
end
end

function w_port = port_assembler(varcov_matrix)
len_matrix = length(varcov_matrix);
w_port = 1 / (ones(len_matrix,1)' * inv(varcov_matrix) * ones(len_matrix,1))...
     * inv(varcov_matrix) * ones(len_matrix,1);
end

function time1 = sequencer(time_cut)
time1 = 0; 
portfolios = importdata('C:\Users\franc\Documents\portfolios.csv',',');
for index = 1:length(portfolios.data(:,1))
    if portfolios.data(index,1) == time_cut
        time1 = time1 + 1;
        break
    end
    if portfolios.data(index,1) ~= time_cut
        time1 = time1 + 1;
    end
end
end

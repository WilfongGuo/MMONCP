clc
clear

                 
%% ***************************Input*************************************

path='E:\Desktop\LSCV_MCEA_code\';        % The path of 'Main,m' on the user's computer and '\' need reserve.
unzip('BRCA_tumor.zip');              % Take BRCA as an example
unzip('BRCA_normal.zip');
 
expression_tumor_fileName = strcat('BRCA_tumor.txt');
expression_normal_fileName = strcat('BRCA_normal.txt');

%% *******************The function of DMOP_LSCV****************************

%% default parameters
popnum=300;
Max_CalNum=200000;
Experiment_num=30;

gen_name=model(expression_tumor_fileName,expression_normal_fileName,path,popnum,Max_CalNum,Experiment_num);


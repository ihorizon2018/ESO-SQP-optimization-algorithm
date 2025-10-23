% This function draw the benchmark functions

function func_plot(func_name)
[lb,ub,dim,fun]=get_benchmark(func_name);

switch func_name 
    case 'sphere' 
        x=-100:2:100; y=x; %[-100,100]
    case 'sumproduct' 
        x=-10:1:10; y=x; %[-10,10] 
    case 'rosenbrock' 
        x=-200:2:200; y=x; %[-5,5]
    case 'schwefel' 
        x=-500:10:500;y=x; %[-500,500]
    case 'rastrigin' 
        x=-5:0.1:5;   y=x; %[-5,5]    
    case 'griewank' 
        x=-500:10:500; y=x;%[-0.5,0.5]
end    

    

L=length(x);
f=[];

for i=1:L
    for j=1:L
        if strcmp(func_name,'F15')==0 && strcmp(func_name,'F19')==0 && strcmp(func_name,'F20')==0 && strcmp(func_name,'F21')==0 && strcmp(func_name,'F22')==0 && strcmp(func_name,'F23')==0
            f(i,j)=fun([x(i),y(j)]);
        end
        if strcmp(func_name,'F15')==1
            f(i,j)=fun([x(i),y(j),0,0]);
        end
        if strcmp(func_name,'F19')==1
            f(i,j)=fun([x(i),y(j),0]);
        end
        if strcmp(func_name,'F20')==1
            f(i,j)=fun([x(i),y(j),0,0,0,0]);
        end       
        if strcmp(func_name,'F21')==1 || strcmp(func_name,'F22')==1 ||strcmp(func_name,'F23')==1
            f(i,j)=fun([x(i),y(j),0,0]);
        end          
    end
end

surfc(x,y,f,'LineStyle','none');

end

func_plot('sphere')
% func_plot('sumproduct')

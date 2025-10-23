% N,T,lb,ub, dim,fobj,sita
function [Xfood, fval,gbest_t,count,X0] = myESO(N,T,lb,ub,dim,fobj)
count=0;
noImprovementCount=0;
vec_flag=[1,-1];
Threshold=0.25;
Thresold2= 0.6;
C3=1;
ub1 = ub.*ones(1,dim);
lb1 = lb.*ones(1,dim);

if(max(size(ub)) == 1)
    ub2 = ub.*ones(1,dim);
    lb2 = lb.*ones(1,dim);
end
%% X=initializationNew(N,dim,ub,lb,fobj);
% 
if isscalar(lb), lbv = lb*ones(1,dim); else, lbv = lb(:)'; end
if isscalar(ub), ubv = ub*ones(1,dim); else, ubv = ub(:)'; end
X = lbv + rand(N,dim).*(ubv - lbv);
fitness = zeros(1, N);
for i = 1:N
    fitness(i) = feval(fobj, X(i,:));
    count = count + 1;
end

t1=zeros(1,T);
t2=zeros(1,T);
a=0.05;
for l= 1:T
    t1(1,l)=1;
    t2(1,l)=1;
end

[GYbest, gbest] = min(fitness);
bestFitness=GYbest;

Xfood = X(gbest,:);
%% Diving the swarm into two equal groups males and females
Nm=round(N/2);%eq.(2&3)
Nf=N-Nm;
Xm=X(1:Nm,:);
Xf=X(Nm+1:N,:);
fitness_m=fitness(1:Nm);
fitness_f=fitness(Nm+1:N);
[fitnessBest_m, gbest1] = min(fitness_m);
BestIndex1=gbest1;
Xbest_m = Xm(gbest1,:);
[fitnessBest_f, gbest2] = min(fitness_f);
BestIndex2=gbest2;
Xbest_f = Xf(gbest2,:);

gbest_t(1)=min([fitnessBest_f fitnessBest_m]);
%% Main loop
for t = 2:T
    Positions=[Xm;Xf];

    C1=0.5;%eq.(14)
    C2=0.05;%eq.(15)

    Temperature=exp(-((t)/T));  %eq.(4)
    Q=C1*exp(((t-T)/(T)));%eq.(5)
    if Q>1        Q=1;    end
    %% Exploration Phase (No Food)
    if Q<Threshold
        for i=1:Nm
            for j=1:1:dim
                rand_leader_index = floor(Nm*rand()+1);
                X_randm = Xm(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Am=exp(-fitness_m(rand_leader_index)/(fitness_m(i)+eps));
                Xnewm(i,j)=X_randm(j)+Flag*C2*Am*((ub-lb)*rand+lb);%eq.(5)
            end
        end
        for i=1:Nf
            for j=1:1:dim
                rand_leader_index = floor(Nf*rand()+1);
                X_randf = Xf(rand_leader_index, :);
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                Af=exp(-fitness_f(rand_leader_index)/(fitness_f(i)+eps));
                Xnewf(i,j)=X_randf(j)+Flag*C2*Af*((ub-lb)*rand+lb);%eq.(6)
            end
        end
        %% Exploitation Phase (Food Exists)
    else
        if Temperature>Thresold2  %hot
            for i=1:Nm
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewm(i,j)=Xfood(j)+C3*Flag*Temperature*rand*(Xfood(j)-Xm(i,j));%eq.(7)
                end
            end
            for i=1:Nf
                flag_index = floor(2*rand()+1);
                Flag=vec_flag(flag_index);
                for j=1:1:dim
                    Xnewf(i,j)=Xfood(j)+Flag*C3*Temperature*rand*(Xfood(j)-Xf(i,j));%eq.(7)
                end
            end
        else %cold
            if rand>0.6 %fight
                for i=1:Nm
                    for j=1:1:dim
                        FM=exp(-(fitnessBest_f)/(fitness_m(i)+eps));
                        Xnewm(i,j)=t1(1,t)*Xm(i,j) +C3*FM*rand*(Q*Xbest_f(j)-Xm(i,j));%eq.(8)

                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        FF=exp(-(fitnessBest_m)/(fitness_f(i)+eps));
                        Xnewf(i,j)=t2(1,t)*Xf(i,j)+C3*FF*rand*(Q*Xbest_m(j)-Xf(i,j));%eq.(9)
                    end
                end
            else%mating
                for i=1:Nm
                    for j=1:1:dim
                        Mm=exp(-fitness_f(i)/(fitness_m(i)+eps));
                        Xnewm(i,j)=Xm(i,j) +C3*rand*Mm*(Q*Xf(i,j)-Xm(i,j));%eq.(10)
                    end
                end
                for i=1:Nf
                    for j=1:1:dim
                        Mf=exp(-fitness_m(i)/(fitness_f(i)+eps));
                        Xnewf(i,j)=Xf(i,j) +C3*rand*Mf*(Q*Xm(i,j)-Xf(i,j));%eq.(11)
                    end
                end
                [~, index]=sort(fitness_m);
                [~, index1]= sort(fitness_f);%排序
                Xnewm(index(end-3),:)=lb+rand*(ub-lb);
                Xnewm(index(end-2),:)=lb+rand*(ub-lb);
                Xnewm(index(end-1),:)=lb+rand*(ub-lb);
                Xnewm(index(end),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-3),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-1),:)=lb+rand*(ub-lb);
                Xnewf(index1(end-2),:)=lb+rand*(ub-lb);
                Xnewf(index1(end),:)=lb+rand*(ub-lb);
            end
        end
    end
    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nm
        Flag4ub=Xnewm(j,:)>ub;
        Flag4lb=Xnewm(j,:)<lb;
        Xnewm(j,:)=(Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewm(j,:));
        count=count+1;%苏晓宇的计数器
        if y<fitness_m(j)
            fitness_m(j)=y;
            Xm(j,:)= Xnewm(j,:);
        end
    end

    %% Return back the search agents that go beyond the boundaries of the search space
    for j=1:Nf
        Flag4ub=Xnewf(j,:)>ub;
        Flag4lb=Xnewf(j,:)<lb;
        Xnewf(j,:)=(Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        y = feval(fobj,Xnewf(j,:));
        count=count+1;%
        if y<fitness_f(j)
            fitness_f(j)=y;
            Xf(j,:)= Xnewf(j,:);

        end
    end
     SSF=0.45;pc=0.8;
    for j = 1:Nm
        if rand>0.8 %rand>0.5
            %% Cauchy mutation
            Temp = Xm(j,:).*(1 + tan(pi*(rand-0.5)));%eq.(28)
            %% Return back the search agents that go beyond the boundaries of the search space
            Temp(Temp>ub) = ub2(Temp>ub);
            Temp(Temp<lb) = lb2(Temp<lb);
            ftemp = fobj(Temp);
            if(ftemp<fitness_m(j))
                fitness_m(j)= ftemp;
                Xm(j,:) = Temp;
            end
        else
            %% Tent-chaos
           % 随机选择三个不同个体进行差分变异
            inx = randperm(Nm, 3);
            while any(inx == j) % 确保不包含自身
                inx = randperm(Nm, 3);
            end
            r1 = inx(1);
            r2 = inx(2);
            r3 = inx(3);

            % 计算差分变异向量
            TempDE = Xm(r1, :) + SSF * (Xm(r2, :) - Xm(r3, :));

            %% Crossover
            del = randi(dim,1);           % genertaing the random variable delta
            for i = 1:dim

                if (rand <= pc) || del == i
                    Temp(i) = TempDE(i);      % check for donor vector or target vector
                else
                    Temp(i) = Xm(j,i);    % check for donor vector or target vector
                end
            end

            ftemp = fobj(Temp);
            %% Return back the search agents that go beyond the boundaries of the search space
            Temp(Temp>ub) = ub2(Temp>ub);
            Temp(Temp<lb) = lb2(Temp<lb);
            if(ftemp<fitness_m(j))
                fitness_m(j) = ftemp;
                Xm(j,:) = Temp;
            end
        end
    end

   avgF = mean(fitness_f);
    for j = 1:Nf
        if rand>0.8 %rand>0.5
            %% Cauchy mutation
            Temp = Xf(j,:).*(1 + tan(pi*(rand-0.5))); %eq.(28)
            %% Return back the search agents that go beyond the boundaries of the search space
            Temp(Temp>ub) = ub2(Temp>ub);
            Temp(Temp<lb) = lb2(Temp<lb);
            ftemp = fobj(Temp);
            if(ftemp<fitness_f(j))
                fitness_f(j) = ftemp;
                Xf(j,:) = Temp;
            end
        else
            % 随机选择三个不同个体进行差分变异
            inx = randperm(Nf, 3);
            while any(inx == j) % 确保不包含自身
                inx = randperm(Nf, 3);
            end
            r1 = inx(1);
            r2 = inx(2);
            r3 = inx(3);

            % 计算差分变异向量
            TempDE = Xf(r1, :) + SSF * (Xf(r2, :) - Xf(r3, :));

            %% Crossover
            del = randi(dim,1);           % genertaing the random variable delta
            for i = 1:dim

                if (rand <= pc) || del == i
                    Temp(i) = TempDE(i);      % check for donor vector or target vector
                else
                    Temp(i) = Xf(j,i);    % check for donor vector or target vector
                end
            end

            ftemp = fobj(Temp);
            %% Return back the search agents that go beyond the boundaries of the search space
            Temp(Temp>ub) = ub2(Temp>ub);
            Temp(Temp<lb) = lb2(Temp<lb);
            if(ftemp<fitness_f(j))
                fitness_f(j) = ftemp;
                Xf(j,:) = Temp;
            end
        end
    end


    [Ybest1,gbest1] = min(fitness_m);
    BestIndex1=gbest1;
    [Ybest2,gbest2] = min(fitness_f);
    BestIndex2=gbest2;
    if Ybest1<fitnessBest_m
        Xbest_m = Xm(gbest1,:);
        fitnessBest_m=Ybest1;
    end
    if Ybest2<fitnessBest_f
        Xbest_f = Xf(gbest2,:);
        fitnessBest_f=Ybest2;

    end
    if Ybest1<Ybest2
        gbest_t(t)=min(Ybest1);
    else
        gbest_t(t)=min(Ybest2);

    end
    if fitnessBest_m<fitnessBest_f
        GYbest=fitnessBest_m;
        Xfood=Xbest_m;
    else
        GYbest=fitnessBest_f;
        Xfood=Xbest_f;
    end


end
fval = GYbest;
X0=[Xm;Xf];



end







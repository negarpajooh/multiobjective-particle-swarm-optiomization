clc;
clear;
close all;
global teta L E v n Kc Kb Ks M Sigmay y z x1 x2 x3 Sai Pi
%% Constant
E=70e9;
n=1;
Ks=8.98;
Kb=41.5;
L=1;
teta=54.7*(pi/180);
y=0;
M=L*v;
Sigmay=.007*E;
%% Problem Definition
Function=questdlg('Plz Select Model?','Selection','Longitudinal','Flexural','Longitudinal')
switch Function
    case 'Longitudinal'
        CostFunction=@(x) Objective1(x);
        v=17920;
        Kc=6.97;
        nPop=500;            % Population Size
        nRep=100;            % Repository Size
    case 'Flexural'
        CostFunction=@(x) Objective2(x);
        v=1792;
        nPop=100;            % Population Size
        nRep=50;            % Repository Size
end

nVar=3;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
VarMin=0.01;          % Lower Bound of Variables
VarMax=0.2;          % Upper Bound of Variables
%% MOPSO Parameters
MaxIt=50;          % Maximum Number of Iterations
% w=0.7298;                % Inertia Weight
% wdamp=1;         % Intertia Weight Damping Rate
% c1=1.4962;               % Personal Learning Coefficient
% c2=c1;  
%% Constriction Coefficient
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt((phi^2)-4*phi));
w=chi;
wdamp=1;
c1=chi*phi1;
c2=chi*phi2;
%%
nGrid=5;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.2;             % Mutation Rate

%% Initialization
empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.X=[];
empty_particle.OutPut=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];
empty_particle.Best.X=[];
empty_particle.Best.OutPut=[];
pop=repmat(empty_particle,nPop,1);

for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Velocity=zeros(VarSize);
    pop(i).particle.X=[abs(x1),abs(x2),abs(x3)];
    pop(i).particle.OutPut=[abs(Sai),abs(Pi)];
    pop(i).Cost=CostFunction(pop(i).Position);
    % Update Personal Best
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost; 
    pop(i).particle.Best.X=pop(i).particle.X;
    pop(i).particle.Best.OutPut=pop(i).particle.OutPut;
end

% Determine Domination
pop=DetermineDomination(pop);

rep=pop(~[pop.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end
%% MOPSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        leader=SelectLeader(rep,beta);
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Cost = CostFunction(pop(i).Position);
        pop(i).particle.X=[abs(x1),abs(x2),abs(x3)];
        pop(i).particle.OutPut=[abs(Sai),abs(Pi)];
        % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
        NewSol.Cost=CostFunction(NewSol.Position);
        if Dominates(NewSol,pop(i))
            pop(i).Position=NewSol.Position;
            pop(i).Cost=NewSol.Cost;
      elseif Dominates(pop(i),NewSol)
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Position=NewSol.Position;
                pop(i).Cost=NewSol.Cost;
            end
        end
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            pop(i).particle.Best.X=pop(i).particle.X;
            pop(i).particle.Best.OutPut=pop(i).particle.OutPut;
  
        elseif Dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
                pop(i).particle.Best.X=pop(i).particle.X;
                pop(i).particle.Best.OutPut=pop(i).particle.OutPut;
            end
        end
        
    end
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep; pop(~[pop.IsDominated])];
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    %%

  
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
    figure(1);
    PlotCosts(pop,rep);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    % Damping Inertia Weight
    w=w*wdamp;

end

%% Resluts
    switch Function
        case 'Longitudinal'
            save('Longitudinal','rep')
        case 'Flexural'
            save('Flexural','rep')
    end



















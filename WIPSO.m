% ***************************************************************
% *** Matlab function for Weight improved Particle Swarm Optimization
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Dr. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
% ***       Mr. M. Prasad (email:prasadgoud333@gmail.com)
% ***       Crustal Processes Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************
function [best_var, best_cost,iter_count,error_energy,all_var,all_cost]  = WIPSO(CostFunction,nVar,MaxIt,nPoP,c1,c2)
	% WIPSO calculates the optimized parameters (best_var) for a given objective function (CostFunction)
	% using Particle Swarm Optimization.
%Inputs
%	Costfunction = Objective function of the optimization problem
%	nVar		 = Number of parameters of the optimization problem
%   Maxit		 = Maximum Generations for PSO 	
%	nPoP		 = Number of particles of the swarm in PSO
%	c1           = cognitive component of PSO
%	c2           = social component of PSO
%   all_var      = all variables at each iterations
%   all_cost     = cost in each iterations

%Outputs 
%best_var		= optimized parameters for PSO
%best_cost		= value of objective function for optimized parameters
%iter_count		= required number of generations for optimization
%error_energy   = Error energy after each generations. 
	
    VarSize = [1 nVar];     %Matrix size of Decision variables
    VarMin= -ones(1,nVar);          %Lower Bound of Unknown variable
    VarMax= ones(1,nVar);           %Upper Bound of Unknown variable

    %% Parameters of PSO
    w_min=0.1; w_max=0.9;                     %inertia coefficient
         
    %% Initialization
    Empty.Particle.Position =[];
    Empty.Particle.Velocity =[];
    Empty.Particle.Cost     =[];

    Empty.Particle.Best.Position =[];
    Empty.Particle.Best.Cost     =[];

    Particle=repmat(Empty.Particle, nPoP,1);

    %initial global best
    GlobalBest.Cost= Inf;
    for i=1:nPoP
        %initialize position with random number from VarMin and VarMax
        for j=1:nVar
            Particle(i).Position(j) =(VarMax(j)-VarMin(j))*rand(1) + VarMin(j);
        end
        %Initialize Velocity
        Particle(i).Velocity =zeros(VarSize);
        %checking cost function value
        Particle(i).Cost = CostFunction(Particle(i).Position);
        %update personal best
        Particle(i).Best.Position =Particle(i).Position;
        Particle(i).Best.Cost =Particle(i).Cost;
        %Update global best
        if Particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest= Particle(i).Best;
        end
    end

    %Best cost value in each iterations
    BestCost=zeros(MaxIt,1);

    %%  Main loop of PSO
%loop for each time steps 
    for it=1: MaxIt
        w=w_max-((w_max-w_min)/MaxIt)*it;
        %loop for all particles of the swarm
        for i=1:nPoP
            %Update Velocity
            Particle(i).Velocity= (w).*Particle(i).Velocity+ ...
                c1*rand(VarSize).*(Particle(i).Best.Position-Particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position-Particle(i).Position);
            %Update Position
            Particle(i).Position=Particle(i).Position+Particle(i).Velocity;
            %cost for this iteration
            Particle(i).Cost = CostFunction(Particle(i).Position);
            %Update Personal Best
            if Particle(i).Cost < Particle(i).Best.Cost
                Particle(i).Best.Position =Particle(i).Position;
                Particle(i).Best.Cost =Particle(i).Cost;
                %Update Global Best
                if Particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest= Particle(i).Best;
                end
            end
        end
        %Display iteration information
        %Store Best Cost value
        BestCost(it)=GlobalBest.Cost;
        BestVar(it,:)=GlobalBest.Position;
        %break the process if misfit achieved 0.01
        if BestCost(it)<0.01
            break
        end
        %printing the cost after each iterations 
        fprintf('After %d iterations Best Cost Value= %.7f\n',it,BestCost(it))  
        %value for error energy plot
        error_energy(it)=(BestCost(it)).^2;
    end
    %best parameter values for the optimization
        best_var= (GlobalBest.Position)'; 
        best_cost= BestCost(it);
        iter_count= it; 
        all_var=BestVar';
        all_cost=BestCost';
end
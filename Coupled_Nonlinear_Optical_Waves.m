classdef Coupled_Nonlinear_Optical_Waves < handle
    properties (Access = public)
        A;
        z;
        convergence_condition = 1e-3;
        max_iterations = 100;
        tol_freq_matching = 1e-3
    end
    properties (Access = private)
        k;
        omega;
        d;
        N;
        delta_p;
        delta_m;
    end
    methods
        function obj = Coupled_Nonlinear_Optical_Waves(d,z,A0,k,omega)
            obj.N = length(A0);
            obj.k = k;
            obj.omega = omega;
            obj.d = d;
            obj.z = z;
            obj.A = NaN(obj.N,length(z));
            obj.A(:,1) = A0;
            
            if obj.N ~= length(obj.omega)
                obj.error('Mismatch between number of initial wave amplitudes and number of frequencies');
            end
            
            if ~isa(obj.k,'function_handle')
                obj.error('k is not a function handle');
            end
            
            obj.delta_p = zeros(length(obj.omega),length(obj.omega),length(obj.omega));
            obj.delta_m = zeros(length(obj.omega),length(obj.omega),length(obj.omega));
            tol = min(gradient(omega))*obj.tol_freq_matching;
            for p = 1:length(obj.omega)
                for m = 1:length(obj.omega)
                    for n = 1:length(obj.omega)
                        diffp = obj.omega(m)+obj.omega(n)-obj.omega(p);
                        obj.delta_p(p,m,n) = (abs(diffp) < tol);
                        diffm = obj.omega(m)-obj.omega(n)-obj.omega(p);
                        obj.delta_m(p,m,n) = (abs(diffm) < tol);
                    end
                end
            end
        end
        
        function solve(obj)
            obj.A(:,2:end)=nan;
            
            multiWaitbar('Solve','Close');
            multiWaitbar('Solve',0,'Color','g');
            for i = 1:length(obj.z)-1
                A1       = obj.A(:,i);
                A_approx = A1;
                
                converge = 1;
                counter = 0;
                while converge > obj.convergence_condition
                    if counter > obj.max_iterations
                        obj.error('Maximum number of iterations reached');
                    end
                    F = obj.eval_F(i+1,A_approx,A1);
                    J = obj.eval_J(i+1,A_approx);
                    dA = J\-F;
                    A_approx = A_approx + dA;
                    converge = norm(dA);
                    counter=counter+1;
                end
                obj.A(:,i+1) = A_approx;
                multiWaitbar('Solve',(i+1)/length(obj.z));
            end
            multiWaitbar('Solve','Close');
        end
    end
    
    methods (Access = private)
        function F = eval_F(obj,i,Ai,A1)
            z1 = obj.z(i-1);
            zi = obj.z(i);
            dz = obj.z(2)-obj.z(1);
            
            F = zeros(obj.N,1);
            
            for p = 1:obj.N
                inds = 1:obj.N;
                inds(p) = [];
                
                thetap = 2*1j*obj.d/299792458^2*obj.omega(p)^2/obj.k(p,zi);
                
                for m = inds
                    for n  = inds
                        Am  = Ai(m);
                        An  = Ai(n);
                        Am1 = A1(m);
                        An1 = A1(n);
                        
                        if obj.delta_p(p,m,n)
                            F(p) = F(p)-dz*thetap/2*( 1/2*Am1*An1*exp(1j*obj.dk_p(p,m,n,z1)*z1) + 1/2*Am *An *exp(1j*obj.dk_p(p,m,n,zi)*zi) );
                        end
                        if obj.delta_m(p,m,n)
                            F(p) = F(p)-dz*thetap/2*( Am1*conj(An1)*exp(1j*obj.dk_m(p,m,n,z1)*z1)+ Am *conj(An )*exp(1j*obj.dk_m(p,m,n,zi)*zi) );
                        end
                    end
                end
                F(p)=F(p)+Ai(p)-A1(p);
            end
        end
        
        function J = eval_J(obj,i,Ai)
            zi = obj.z(i);
            dz = obj.z(2)-obj.z(1);
            
            J = zeros(obj.N,obj.N);
            
            for p = 1:obj.N
                thetap = 2*1j*obj.d/299792458^2*obj.omega(p)^2/obj.k(p,zi);
                for q = 1:obj.N
                    if p == q
                        J(p,q) = 1;
                    else
                        inds = 1:obj.N;
                        inds(p) = [];
                        for m = inds
                            for n  = inds
                                Am  = Ai(m);
                                An  = Ai(n);
                                if obj.delta_p(p,m,n)
                                    J(p,q) = J(p,q) - dz*thetap/2*1/2*(Am*(n==q)+     An*(m==q)) *exp(1j*obj.dk_p(p,m,n,zi)*zi);
                                end
                                if obj.delta_m(p,m,n)
                                    J(p,q) = J(p,q) - dz*thetap/2*    (Am*(n==q)+conj(An)*(m==q))*exp(1j*obj.dk_m(p,m,n,zi)*zi);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function dk = dk_p(obj,p,m,n,zi)
            dk = obj.k(m,zi)+obj.k(n,zi)-obj.k(p,zi);
        end
        
        function dk = dk_m(obj,p,m,n,zi)
            dk = obj.k(m,zi)-obj.k(n,zi)-obj.k(p,zi);
        end
        
        function error(obj,str)
            multiWaitbar('CloseAll');
            error(str);
        end
    end
end

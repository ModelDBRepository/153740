%This programs solves the diffusion-reaction-advection equations related in the transport and conversion of JNK to activated JNK along microtubules to neuron ends
%A refers to JNK; AA refers to activated JNK; F refers to scaffold JIP1; E refers to enzyme MKK7; M refers to motor protein KIF1; E2 refers to phosphatase M3/6 
%Combination of symbols refers to complexes composed of the corresponding compounds. 
clear all

%x is the length of the neuron
x_length = 10; %Each unit is 10um.

%discretizing x
x_step_precount = 200;
x_step_count = x_step_precount+1;
x_step = x_length/(x_step_precount);
x_pos = 0:x_step:x_length;

%setting interval time for calculations
dt = 0.01;
timesteps = 100000;

%recording values at specific time intervals
record = 1000;
record_blocks = timesteps/record;

%s is the variable assigned for recording of values
s_t = zeros(record_blocks,x_step_count,15);

%initializing
countj=0;
%initializaing initial gausian mean concentration of JNK(A), KIF5(M) and MKK7(Enzy)
A = 10;
M = 20;
Enzy = 10/6.266;

for F = 7.5 %set F = 0 to simulate a case where f = 0 and consequently p = 0
    
    %S refers to speed of motor protein translocation
    for S = 0.025

        count_step = 0;
        
        s_0 = zeros(16, x_step_count);
        s_0(1,:) = A*exp(-(x_pos-0).^2/0.5);%A
        
        s_0(3,:) = F*exp(-(x_pos-0).^2/0.5); %F
        
        s_0(4,:) = Enzy*exp(-(x_pos-0).^2/0.5); %E
        s_0(5,:) = 1; %E2
        
        s_0(11,:) = M.*exp(-(x_pos-0).^2/0.05);%M
        
        s_1 = s_0;
        
        %setting diffusion constants
        D = zeros(1,15);
        D(1:5) = 0.1;
        D(6:8) = 0.1/sqrt(2);
        D(9) = 0.1/sqrt(3);
        D(10) = 0.1/sqrt(2);
        D(11:15) = 5e-5;
        
        %setting reaction constants
        %rate constants
        k1 = 0.1;
        k2 = 0.4;
        k4 = 0.1;
        k5 = k2;
        %association constants
        kf1 = 1;
        kf2 = 1;
        kf3 = 1;
        kf4 = 1;
        kf5 = 1;
        kf6 = 1;
        %dissociation constants
        kr1 = 1;
        kr2 = 1;
        kr3 = 1;
        kr4 = 1;
        kr5 = 1;
        kr6 = 1;
        %binding constants
        b1 = 0.1;
        b2 = 0.1;
        b3 = 0.5;
        b4 = 0.5;
        m1 = 0.5;
        m2 = 0.5;
        m3 = 0.5;
        %unbinding constants
        u1 = 0.1;
        u2 = 0.1;
        u3 = 0.1;
        u4 = 0.1;
        n1 = 0.1;
        n2 = 0.1;
        n3 = 0.1;
        
        for t = 1: timesteps
            
            %t
            
            s = zeros(15,x_step_count);
            
            A_0 = s_0(1,:);
            AA_0 = s_0(2,:);
            F_0 = s_0(3,:);
            E_0 = s_0(4,:);
            E2_0 = s_0(5,:);
            FA_0 = s_0(6,:);
            FE_0 = s_0(7,:);
            AE_0 = s_0(8,:);
            FAE_0 = s_0(9,:);
            AAE2_0 = s_0(10,:);
            M_0 = s_0(11,:);
            MF_0 = s_0(12,:);
            MFA_0 = s_0(13,:);
            MFE_0 = s_0(14,:);
            MFAE_0 = s_0(15,:);
            
            %characterizing reactions
            %forward reactions
            v1 = k1.*AE_0;
            v2 = k2.*FAE_0;
            v5 = k5.*MFAE_0;
            %backward reactions
            v4 = k4.*AAE2_0;
            %association dissociation reactions
            ad1 = kf1.*A_0.*E_0 - kr1.*AE_0;
            ad2 = kf2.*FA_0.*E_0 - kr2.*FAE_0;
            ad3 = kf3.*FE_0.*A_0 - kr3.*FAE_0;
            ad4 = kf4.*AA_0.*E2_0 - kr4.*AAE2_0;
            ad5 = kf5.*MFA_0.*E_0 - kr5.*MFAE_0;
            ad6 = kf6.*MFE_0.*A_0 - kr6.*MFAE_0;
            %binding/unbinding reactions
            BU1 = b1.*A_0.*F_0 - u1.*FA_0;
            BU2 = b2.*E_0.*F_0 - u2.*FE_0;
            BU3 = b3.*MF_0.*A_0 - u3.*MFA_0;
            BU4 = b4.*MF_0.*E_0 - u4.*MFE_0;
            M1 = m1.*FA_0.*M_0 - n1.*MFA_0;
            M2 = m2.*F_0.*M_0 - n2.*MF_0;
            M3 = m3.*FE_0.*M_0 - n3.*MFE_0;
           
            %setting reactions for each component in the simulation
            rxns = zeros(15,x_step_count);
            rxns(1,:) = -ad1 -ad3 +v4 - BU1 - BU3 - ad6; %A
            rxns(2,:) = v1 + v2 -ad4 + v5; %AA
            rxns(3,:) = v2-BU1 -BU2 - M2; %F
            rxns(4,:) = -ad1+ v1 -ad2 + v2 -BU2 - BU4 -ad5 +v5; %E
            rxns(5,:) = -ad4 + v4; %E2
            rxns(6,:) = -ad2 + BU1 - M1; %FA
            rxns(7,:) = -ad3 + BU2 - M3; %FE
            rxns(8,:) = ad1- v1; %AE
            rxns(9,:) = ad2 - v2 + ad3;  %FAE
            rxns(10,:) = ad4 - v4; %AAE2
            rxns(11,:) = -M1 - M2 -M3; %M
            rxns(12,:) = -BU3 - BU4 + M2 + v5; %MF
            rxns(13,:) = BU3 +M1 -ad5; %MFA
            rxns(14,:) = BU4 + M3 - ad6; %MFE
            rxns(15,:) = ad5 - v5 + ad6; %MFEA
            
            %initializing boundary conditions
            s_L = zeros(15,1);
            s_R = zeros(15,1);
            
            %updating species,j concentrations due to reaction and diffusion
            
            %non motor protein species
            %Neumann boundary condition for non motor protein species
            s_L(1:10) = s_0(1:10,1);
            s_R(1:10) = s_0(1:10,x_step_count);    
            
            %Forward-Time Central-Space(FTCS)scheme was implemented
            for j = 1:10                
                for x = 1;
                    s(j,x) = ( D(j)*(s_0(j,x+1)-2*(s_0(j,x))+s_L(j))/x_step/x_step + rxns(j,x) )*dt + s_0(j,x);
                end
                
                for x = 2:x_step_count-1;
                    s(j,x) = ( D(j)*(s_0(j,x+1)-2*(s_0(j,x))+s_0(j,x-1))/x_step/x_step + rxns(j,x) )*dt + s_0(j,x);
                end
                
                for x = x_step_count
                    s(j,x) = ( D(j)*(s_R(j)-2*(s_0(j,x))+s_0(j,x-1))/x_step/x_step + rxns(j,x) )*dt + s_0(j,x);
                end                
            end

            %motor protein
            for j = 11
                s(j,:) = s_0(j,:) + rxns(j,:)*dt;
                
            end            
            
            d = D.*dt/x_step/x_step;
            c = S*dt/x_step;
            %motor protein species
            %Dirichlet boundary condition for motor protein species 
            for j = 12:15                
                c = S*dt/x_step;                
                s_L(j) = -s_0(j,1)*(-(c^2)/2 + c/2 -D(j)/x_step/x_step)/((c^2)/2 + c/2 +D(j)/x_step/x_step);   
                
                %second order Lax-Wendroff scheme was implemented
                for x = 1
                    s(j,x) = s_0(j,x) - (c/2)*(s_0(j,x+1)-s_L(j)) + ((c^2)/2)*(s_0(j,x+1)-2*s_0(j,x)+s_L(j))+ D(j)*(s_0(j,x+1)-2*(s_0(j,x))+s_L(j))/x_step/x_step + rxns(j,x)*dt;
                end                
                
                for x = 2:x_step_count-1;
                    s(j,x) = s_0(j,x) - (c/2)*(s_0(j,x+1)-s_0(j,x-1)) + ((c^2)/2)*(s_0(j,x+1)-2*s_0(j,x)+s_0(j,x-1)) + D(j)*(s_0(j,x+1)-2*(s_0(j,x))+s_0(j,x-1))/x_step/x_step + rxns(j,x)*dt;
                end
                
                for x =  x_step_count
                    s(j,x) = (-S*(s_0(j,x)-s_0(j,x-1))/x_step)*dt + s_0(j,x)+ rxns(j,x)*dt;
                end               
            end
            
            %recording of values into s_t
            if rem(t,record) == 0                
                count_step = count_step + 1                
                for j = 1:15
                    s_t(count_step,:,j)= s(j,:);
                end                
            end
            
            s_0 = s;
            
            
        end
        
        %saving variables into file
        filetitle = strcat('S',num2str(S*1000),'_M',num2str(M),'_F',num2str(F*10));
        save(filetitle,'s_t');
    end
end

%display kymograph of JNK
JNK = s_t(:,:,2);
timegaps = (0:record:timesteps).*dt;
spacegaps = x_pos.*10;
imagesc(spacegaps,timegaps,JNK);caxis([0 0.06]); colormap('autumn');colorbar
axis xy
xlabel('Time (s)');
ylabel('x (µm)');
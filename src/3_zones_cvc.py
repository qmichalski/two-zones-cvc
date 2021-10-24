# Cantera sliced constant volume combustion Q. Michalski tests


% Import of initial properties from experimental 0D cycle
load('n0072_0D_essai_n0169_analysis','mech','Ymix','Tmix','p_r','time','CVCC');
nCycle=7;
ignPos=find(time>=CVCC.allumage.instant,1,'first')-1;
Texp=Tmix(ignPos,nCycle);
pExp=p_r(ignPos,nCycle,1);
pExp_r=p_r;
Yexp=squeeze(Ymix(:,ignPos,nCycle,1));

% Declaration of Cantera environnement
% Choice of Cantera mechanism (only species and thermodynamic data are
% relevant)
mech='gri30_highT_mechII.xml';
% Declaration of 3 cantera gas pointers required for 2-zones combustion :  
cvc_gas{1}=Solution(mech); % description of fresh gases
cvc_gas{2}=Solution(mech); % description of mixed burnt gases
cvc_gas{3}=Solution(mech); % description of current burning slice of fresh gas
UVref=Solution(mech); % gas-set used to compute the reference internal energy conservation associated properties
gas_mix=Solution(mech); % gas-set used to compute an averaged mixture (fresh+burnt)

% Initialisation
% Choice of fuel
fuel='C3H8';
switch fuel % set the stoechiometric chemical equation
    case 'C3H8'
        nO2=5;
    case 'NC10H22'
        nO2=15.5;
end
% Choice of initial conditions (T,P,X)
% Construction of mol fraction vector X based on equivalence ratio ER
ER=1;
X.(fuel)=ER/(ER+nO2*(1+0.79/0.21)) ;
X.O2=nO2/(ER+nO2*(1+0.79/0.21)) ;
X.N2=nO2*0.79/0.21/(ER+nO2*(1+0.79/0.21)) ;
% Choice of initial pressure p0
p0=pExp; % 3.5e5;
% Choice of initial temperature T0
T0=Texp; % 400;
% Initializing the sets of gas
set(cvc_gas{1},'X',[fuel ':' num2str(X.(fuel)) ' ; O2:' num2str(X.O2) ' ; N2:' num2str(X.N2) ],'P',p0,'T',T0) ;
set(UVref,'X',[fuel ':' num2str(X.(fuel)) ' ; O2:' num2str(X.O2) ' ; N2:' num2str(X.N2) ],'P',p0,'T',T0) ;
set(cvc_gas{2},'X',[fuel ':' num2str(X.(fuel)) ' ; O2:' num2str(X.O2) ' ; N2:' num2str(X.N2) ],'P',p0,'T',T0) ;
% equilibrate(cvc_gas{2},'HP');
set(gas_mix,'X',[fuel ':' num2str(X.(fuel)) ' ; O2:' num2str(X.O2) ' ; N2:' num2str(X.N2) ],'P',p0,'T',T0) ;
% Get the final properties of the CVC
equilibrate(UVref,'UV'); % equilibrate the mixture
Tref=temperature(UVref); % get the reference final temperature
Pref=pressure(UVref); % get the reference final pressuer
% Declaration of chamber dimensions values (required for heat exchanges)
L=100e-3; % Length of chamber [m]
l=50e-3; % Square width of chamber [m]
Vt=L*l^2; % Volume of chamber [m]
Twall=T0; % Wall temperature [K]

% Initilization
t=0;
% Choice of time step dt [s]
dt=5e-5; % time step [s]
y=0; % extent of reaction [-]
% Volume vector [m3]
V(1)=Vt; % volume de gaz [m3]
V(2)=0; % volume of mixed burnt gases [m3]
V(3)=0; % volume of current burnt slice [m3]
rhoRef=density(cvc_gas{1}); % density conserved through calculation [kg/m3]
u0=intEnergy_mass(cvc_gas{1}); % internal energy conserved through calculation [J/kg]
m=V*rhoRef; % Mass vector [kg]

% Export variables
ii=1; % current step
timer(ii)=t; % time record
p_r(ii)=pressure(cvc_gas{1}); % pressure record
T_r(ii,1)=temperature(cvc_gas{1}); % fresh gases temperature record 
T_r(ii,2)=temperature(cvc_gas{2}); % mixed burnt gases temperature record
T_mix_r(ii)=temperature(gas_mix); % averaged mixture temperature record
s_mix_r(ii)=entropy_mass(gas_mix); % averaged mixture entropy record
h_mix_r(ii)=enthalpy_mass(gas_mix); % averaged mixture enthalpy record
Y_mix(:,ii)=massFractions(gas_mix); % averaged mixture mass fractions record
gm(ii)=cp_mass(cvc_gas{1})/cv_mass(cvc_gas{1}); % fresh gases gamma record
V_r1(ii,1:numel(V))=V;
m_r(1,1:numel(V))=m; % mass record

% Definition of the flame surface model used for calculation
load('spherical_flame_case_9','rset','St','Stt');
r=0; % flame surface radius [m] 
St=St/(10^6); % flame surface [m²]
Stt=Stt/(10^6); % fame surface + burnt gas surface in contact with chamber [m²]
Sb=Stt-St; % burnt gas surface in contact with chamber [m²]
% Construction of interpolated functions
rset=rset/1000; % radius for construction of interpolated functions
pp=interp1(rset,St,'linear','pp');
St=St*Vt/max(cumtrapz(rset,St));
Sf=@(r) ppval(pp,r); % interpolated flame surface
pp=interp1(rset,cumtrapz(rset,St),'linear','pp');
Vf=@(r) ppval(pp,r); % interpolated burnt volume
pp=interp1(rset,max(Sb)-Sb,'linear','pp');
Sex{1}=@(r) ppval(pp,r); % interpolated contact surface with fresh gases
pp=interp1(rset,Sb,'linear','pp');
Sex{2}=@(r) ppval(pp,r); % interpolated contact surface with burnt gases

flag=true; % stop flag
VUnb=1e-6; % unburnt volume
flux=1; % heat exchange on (flux=1) or off (flux=0)
% Constant volume sliced calculation
while flag
    xi=8; % flame corrugation factor [-]
    vF=flame_velocity('C3H8',temperature(cvc_gas{1}),pressure(cvc_gas{1}),ER); % fundamental flame velocity [m/s]
%     vF=1;% m/s
    vF=vF*xi; % turbulent flame velocity [m/s]
    S=Sf(r+vF*dt); % projected flame surface [m²]
    dv=vF*dt*S; % volume burnt [m3]
    if (V(1)-dv)<=VUnb
        newdt=(V(1)-VUnb)/(vF*S);
        flag=false;
        mUnb=VUnb*density(cvc_gas{1});
    else
        newdt=dt;
        flag=true;
    end
    r=r+vF*newdt; % new flame radius [m]
    S=Sf(r); % effective flame surface [m²]
    dv=vF*newdt*S;
    vF_r(ii)=vF;
    t=t+dt;
    dispstat(['Combustion at ... ' num2str(y*100,'%3.0f') '%']) ;
    %%Burnt the fresh gas slice the flame just moved on
    V(1)=V(1)-dv;
    set(cvc_gas{3},'X',moleFractions(cvc_gas{1}),'P',pressure(cvc_gas{1}),'T',temperature(cvc_gas{1})) ;
    equilibrate(cvc_gas{3},'HP');
    rhotrf(ii)=density(cvc_gas{1});
    rhotrb(ii)=density(cvc_gas{3});
    dv=rhotrf(ii)/rhotrb(ii);
    V(3)=vF*newdt*S*dv;
    dv_s(ii)=dv;
    %% Calculate new pressure
    p0=pressure(cvc_gas{1});
    % density prior to compression
    for i=1:numel(V)
        rho0(i)=density(cvc_gas{i});
    end
    fun=@(p) isentropic_compression_gas_set_rev2(p,V/Vt,rho0,cvc_gas,[0 1 1]);
    p=fzero(fun,p0);
    % density after the compression
    for i=1:numel(V)
        rho(i)=density(cvc_gas{i});
    end
    V=V.*rho0./rho;
    %%Mix burnt slice into burnt gases
    % Overal mass fraction based on overall volume fraction and density
    m=V.*rho;
    % Species mass fraction calculation
    Y=(m(2)*massFractions(cvc_gas{2})+m(3)*massFractions(cvc_gas{3}))/(m(2)+m(3));
    H=(m(2)*enthalpy_mass(cvc_gas{2})+m(3)*enthalpy_mass(cvc_gas{3}))/(m(2)+m(3));
    set(cvc_gas{2},'Y',Y,'H',H,'P',p);
    V(2)=V(2)+V(3);
    y=V(2)/Vt;
    V(3)=0;
    m=V.*rho;
    if r>0
        rRec=@(r) Vf(r)-V(2);
        r=fzero(rRec,r);
    end
    %%Echanges aux parois
    if flux
        p=pressure(cvc_gas{1});
        eps=0;
        Vlg=vF*(dv-1);
        for i=1:2
            if m(i)/sum(m)>eps
                alpha=3;
                rho0=density(cvc_gas{i});
                h=alpha*130*sum(V(1,:))^(-0.06)*(p/10^5)^0.8*temperature(cvc_gas{1,i})^(-0.4)*(Vlg+1.4)^0.8;
                dH=Sex{i}(r)*h*(Twall-temperature(cvc_gas{1,i}))*newdt;
                dH_r(i,ii)=dH;
                h=(enthalpy_mass(cvc_gas{1,i})*m(i)+dH)/m(i);
                set(cvc_gas{1,i},'P',p,'H',h);
                dv=rho0/density(cvc_gas{1,i});
                V(i)=V(i)*dv;
            end
        end
        p0=p;
        % density prior to compression
        for i=1:numel(V)
            rho0(i)=density(cvc_gas{i});
        end
        fun=@(p) isentropic_compression_gas_set_rev2(p,V/Vt,rho0,cvc_gas,[0 0 0]);
        p=fzero(fun,p0);
        % density after the compression
        for i=1:numel(V)
            rho(i)=density(cvc_gas{i});
        end
        V=V.*rho0./rho;
    end
    if r>0
        rRec=@(r) Vf(r)-V(2);
        r=fzero(rRec,r);
    end
    %%Export
    ii=ii+1;
    timer(ii)=t;
    m_r(ii,1:numel(V))=m;
    T_r(ii,1)=temperature(cvc_gas{1});
    T_r(ii,2)=temperature(cvc_gas{2});
    rho0_r(ii)=density(cvc_gas{1});
    Vlg_r(ii)=Vlg;
    p_r(ii)=p;
    V_r1(ii,1:numel(V))=V;
    resV(ii)=sum(V)/Vt-1;
    m=V.*rho/rhoRef;
    Y=(m(1)*massFractions(cvc_gas{1})+m(2)*massFractions(cvc_gas{2}))/(m(1)+m(2));
    H=(m(1)*enthalpy_mass(cvc_gas{1})+m(2)*enthalpy_mass(cvc_gas{2}))/(m(1)+m(2));
    set(gas_mix,'Y',Y,'H',H,'P',p);
    cpm=m(1)*cp_mass(cvc_gas{1})+m(2)*cp_mass(cvc_gas{2});
    cvm=m(1)*cv_mass(cvc_gas{1})+m(2)*cv_mass(cvc_gas{2});
    gm(ii)=cpm/cvm;
    T_mix_r(ii) = temperature(gas_mix);
    s_mix_r(ii) = entropy_mass(gas_mix);
    h_mix_r(ii) = enthalpy_mass(gas_mix);
    Y_mix(:,ii) = massFractions(gas_mix);
end

disp(['Conservation densite : ' num2str((rhoRef-density(cvc_gas{2}))/rhoRef*100) '%'])
disp(['Ecart pression : ' num2str((Pref-pressure(cvc_gas{2}))/Pref*100) '%'])
disp(['Ecart temperature : ' num2str((Tref-temperature(cvc_gas{2}))/Tref*100) '%'])
disp(['Conservation energie interne : ' num2str((u0-intEnergy_mass(cvc_gas{2}))/u0*100) '%'])

figure
plot(time(ignPos:end,nCycle)*1000,pExp_r(ignPos:end,nCycle,1)/10^5,'-x')
hold on
plot(timer*1000,p_r);
xlim([0 timer(end)*1000])
% plot(timer(1:end-1)*1000,(gm(1:end-1)-1).*diff(p_r)./diff(timer),'-x')
% figure
% plot(timer(1:end-1),(gm(1:end-1)-1).*diff(p_r)./diff(timer)./rho0_r(1:end-1),'-x')
% figure
% plot(timer,V_r1,'-x')
% figure
% plot(pos,m_r,'-x')
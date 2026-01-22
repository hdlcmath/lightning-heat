clear; clc;

% Heat solver

% Important Parameters
m=50; % Newman poles per corner
t=0.3; % Time of heat solution
M=9; % Evals/points in inverse laplace transform quadrature

vv=[1+1i,2i,-1+1i,-1-1i,1-1i]; % Vertices
r_pols=0.25i;
gamma=cell(1,length(vv));

% Create parametrization of each side
for j=1:length(vv)
	a=vv(j);
	jplus=mod(j,length(vv))+1;
	b=vv(jplus);
	gamma{j}=@(t)a+(b-a)*t;
end
col_factor=5; % This is admittedly aggressively high; optimize with a lower value
runge_factor=3.5;
N_r=ceil(runge_factor*sqrt(m)); % Set Runge order to be O(sqrt(m))
N_nm=length(vv)*m;
N_terms=2*(N_nm+N_r*length(r_pols))+1;
N_dofs=2*N_terms;
N_c=col_factor*N_dofs;
r_orders=N_r*ones(1,length(r_pols));

% Get the Newman poles
edge_lengths=abs(vv-vv([2:length(vv),1]));
r=0.7071*min(edge_lengths);
nm_pols=[];
for j=1:length(vv)
	jplus=mod(j,length(gamma))+1;
	epsilon=1e-5;

	vbwd=gamma{j}(1-epsilon)-vv(jplus);
	vfwd=gamma{jplus}(epsilon)-vv(jplus);
	vang=-vfwd/vbwd; 


	outward=1i*vbwd*sqrt(vang);
	outward=r*outward/abs(outward);
	beta=abs(angle(vang))/pi+1;

	b=vv(jplus);
	a=vv(jplus)-outward;

	tmp_pols=nm_cluster_to(a,b,m,beta); 
	inds=find(abs(tmp_pols-b)<2e-16); tmp_pols(inds)=0;
	nm_pols=[nm_pols,tmp_pols];
end

% Get collocation points
N_c_halfedge=ceil(N_c*edge_lengths/sum(edge_lengths)/2);
cols=[];
for j=1:length(vv)
	N=N_c_halfedge(j);
	c1=col_cluster_to(.5,0,N,m);
	c2=col_cluster_to(.5,1,N,m);
	tt=[flip(c1),c2];
	cols=[cols,gamma{j}(tt)];
end

% Plot geometry, collocation and Newman poles
figure(1); clf;
fill(real(vv),imag(vv),'w','handlevisibility','off'); hold on; grid on;
axis padded equal;
plot(nm_pols,'.','displayname','Newman Poles');
plot(real(r_pols),imag(r_pols),'*','displayname','Runge Pole(s)');
plot(cols,'r.','displayname','Collocation Points');
legend

% Setting up the problem
source=3+2i;
up=@(z,s)psi(0,z,source,s);

uh_mhe=cell(1,M);
[ww,ss]=ilt_quadrature(t,M);
hsolve=tic;
for j=1:M
	F=@(z)-up(z,ss(j));
	tic
	[uh_tmp,residuals]=solve_mhe(cols,F,r_pols,r_orders,nm_pols,ss(j));
	fprintf('Modified Helmholtz problem %.f/%.f solved in %.2f\n',j,M,toc)
	uh_mhe{j}=uh_tmp;
end	
fprintf('Heat problem solved in %.2f\n',toc(hsolve))

% Plot particular/homogenous part and solution
bds=[-4,4];
gridres=200;
zz=linspace(bds(1),bds(2),gridres);
Z=zz+1i*zz.';

% Get helmholtz solutions
Up_mhe=cell(1,M);
Uh_mhe=cell(1,M);
U_mhe=cell(1,M);
heval=tic;
for j=1:M
	tic
	Up_mhe{j}=up(Z,ss(j));
	Uh_mhe{j}=uh_mhe{j}(Z);
	U_mhe{j}=Up_mhe{j}+Uh_mhe{j};
	fprintf('Modified Helmholtz problem %.f/%.f evaluated on grid in %.2f\n',j,M,toc)
end
fprintf('Heat problem evaluted on grid in %.2f\n',toc(heval))

figure(2); clf
tiledlayout('flow')

for j=1:M
	nexttile
	I=pcolor(real(Z),imag(Z),abs(real(U_mhe{j}))); set(I,'EdgeColor','none')
	colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
	axis equal; axis([bds,bds]); 
	title(['$\tilde{u}_p$ for $s=$',num2str(ss(j),'%.2f')],'interpreter','latex');
	hold on; fill(real(vv),imag(vv),'w');
end

% Sum MHE solutions to get heat solution
Up=zeros(size(Up_mhe{1}));
Uh=Up; U=Up;
for j=1:M
	Up=Up+ww(j)*Up_mhe{j};
	Uh=Uh+ww(j)*Uh_mhe{j};
	U=Up+Uh;
end
Up=real(Up/(M*1i)); Uh=real(Uh/(M*1i)); U=real(U/(M*1i));

% The above as a (vectorized!) anonymous function
u_heat=@(z)eval_heat(z,uh_mhe,up,ww,ss);

figure(3); clf
tiledlayout(1,3)

nexttile
I=pcolor(real(Z),imag(Z),abs(real(Up))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$u_h$','interpreter','latex')
hold on; fill(real(vv),imag(vv),'w')

nexttile
I=pcolor(real(Z),imag(Z),abs(real(Uh))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$u_h$','interpreter','latex')
hold on; fill(real(vv),imag(vv),'w')

nexttile
I=pcolor(real(Z),imag(Z),abs(real(U))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$u=u_p+u_h$','interpreter','latex')
hold on; fill(real(vv),imag(vv),'w')

% Get oversampled grid on boundary and evaluate
N_s_halfedge=3*N_c_halfedge;
samps=[];
bdp=[]; % Boundary parametrization
for j=1:length(vv)
	N=N_s_halfedge(j);
	c1=samp_cluster_to(.5,0,N,m);
	c2=samp_cluster_to(.5,1,N,m);
	tt=[flip(c1),c2];
	samps=[samps,gamma{j}(tt)];
	bdp=[bdp,tt+j-1];
end
U_oversampled=u_heat(samps);

figure(4); clf;
semilogy(bdp,abs(U_oversampled),'k.-')
hold on
grid on
yline(max(abs(U_oversampled)),'k--','Einf')
title('Relative boundary error')
fprintf('Einf error: %.6e\n',max(abs(U_oversampled)))


%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=eval_heat(z,uh_mhe,up,ww,ss)
	zz=z(:);
	heatsol=0;
	M=length(ss);
	for j=1:M
		uuh=uh_mhe{j}(zz);
		uup=up(zz,ss(j));
		uu_mhe=uuh+uup;
		heatsol=heatsol+ww(j)*uu_mhe;
	end
	heatsol=real(heatsol/(M*1i));
	out=reshape(heatsol,size(z));
end

function out=psi(n,z,xi,s)
	dd=z-xi;
	arg=sqrt(s)*abs(dd);
	out=besselk(abs(n),arg).*dd.^n./abs(dd).^n/(2*pi);
end

function [qW,qP] = ilt_quadrature(t,M)
% Given a value of t, output a vector of weights and M points for inverse laplace transform

th = (.5:M-.5)*pi/M; % Vector of theta to go in parametrization
sg = -0.6122; mu = 0.5017; bt = 0.6407; nu = 0.2645;
z=@(th,t) 2*M./t.*(sg + mu*th.*cot(bt*th)+nu*1i*th);
dz=@(th,t) 2*M./t.*(mu*cot(bt*th)...
		-mu*bt*th.*csc(bt*th).^2+nu*1i);
	
qP = z(th,t);
qW = exp(z(th,t)*t).*dz(th,t);

end


function [uh,residuals]=solve_mhe(col,F,r_pols,r_order,nm_pols,s)
	cc=col(:);
	col_val=F(cc);
	A=get_problem_matrix(cc,r_pols,r_order,nm_pols,s);
	normA=vecnorm(A); scaledA=A./normA;
	coefs=scaledA\col_val;
	
	uh=@(z)eval_solution(z,coefs./normA(:));
	residuals=0;

	function out=eval_solution(z,coefs)
		zz=z(:);
		A=get_problem_matrix(zz,r_pols,r_order,nm_pols,s);
		U=A*coefs;
		out=reshape(U,size(z));
	end

end

function out=get_problem_matrix(in_pts,r_pols,r_order,nm_pols,s)
	zz=in_pts;
	nm_part=get_newman_part(zz,nm_pols,s);
	rg_part=get_runge_part(zz,r_pols,r_order,s);
	out=[nm_part,rg_part];
end

function out=get_runge_part(in_pts,r_pols,r_order,s)
	zz=in_pts(:);
	
	out=[];
	for NR=length(r_pols)
		N=r_order(NR);
		Nvec=-N:N;
		xi=r_pols(NR);
		tmp=zeros(length(in_pts),length(Nvec));
		for j=1:length(Nvec);
			tmp(:,j)=psi(Nvec(j),zz,xi,s);
		end
		out=[out,tmp];
	end
end


function out=get_newman_part(in_pts,nm_pols,s)
	zz=in_pts(:);
	xi=nm_pols;

	o1=psi(1,zz,xi,s);
	oneg1=psi(-1,zz,xi,s);

	out=[oneg1,o1];
end

function out = nm_cluster_to(a,b,N,beta)
	sg=sqrt(2*(2-beta)*beta)*pi; % from Herremans paper

	% Tapered Lightning
	nn=1:N;
	rate=exp(-sg*(sqrt(N+1)-sqrt(nn)));

	rate(rate==0)=[];
	tmp=(b-a)*(1-rate)+a;
	out=tmp;
end

function out = col_cluster_to(a,b,N,m)
	sg=sqrt(2*(m+1))*pi;
	v=linspace(0,1,N);
	rate=1-exp(-sg*v);
	tmp = (b-a)*rate+a;
	out = tmp;
end

function out = samp_cluster_to(a,b,N,m)
	sg=sqrt(2*(m+1))*pi+log(10);
	v=linspace(0,1,N);
	rate=1-exp(-sg*v);
	tmp = (b-a)*rate+a;
	out = tmp;
end

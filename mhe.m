clear; clc;

% Modified Helmholtz Equation solver

% Important Parameters
m=50; % Newman poles per corner

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
	inds=find(abs(tmp_pols-b)<2*eps); tmp_pols(inds)=0;
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

s=-10+8i;
F=@(z)-up(z,s);
tic
[uh,residuals]=solve_mhe(cols,F,r_pols,r_orders,nm_pols,s);
u=@(z)uh(z)+up(z,s); % Function handle of the solution
fprintf('Modified Helmholtz problem solved in %.2f\n',toc)

% Plot particular/homogenous part and solution
bds=[-4,4];
gridres=200;
zz=linspace(bds(1),bds(2),gridres);
Z=zz+1i*zz.';

tic
Up=up(Z,s);
Uh=uh(Z);
U=Up+Uh;
fprintf('Modified Helmholtz problem evaluted on grid in %.2f\n',toc)

figure(2); clf
tiledlayout(1,3)

nexttile
I=pcolor(real(Z),imag(Z),abs(real(Up))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$\tilde{u}_p$','interpreter','latex')
hold on; fill(real(vv),imag(vv),'w');

nexttile
I=pcolor(real(Z),imag(Z),abs(real(Uh))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$\tilde{u}_h$','interpreter','latex')
hold on; fill(real(vv),imag(vv),'w')

nexttile
I=pcolor(real(Z),imag(Z),abs(real(U))); set(I,'EdgeColor','none')
colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
axis equal; axis([bds,bds]); title('$\tilde{u}=\tilde{u}_p+\tilde{u}_h$','interpreter','latex')
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
U_oversampled=u(samps);

figure(3); clf;
semilogy(bdp,abs(U_oversampled),'k.','displayname','abs')
hold on
semilogy(bdp,abs(real(U_oversampled)),'r-','displayname','real')
semilogy(bdp,abs(imag(U_oversampled)),'b-','displayname','imag')
grid on
yline(max(abs(U_oversampled)),'k--','Einf')
title('Relative boundary error')
fprintf('Einf error: %.6e\n',max(abs(U_oversampled)))

%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function out=psi(n,z,xi,s)
	dd=z-xi;
	arg=sqrt(s)*abs(dd);
	out=besselk(abs(n),arg).*dd.^n./abs(dd).^n/(2*pi);
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

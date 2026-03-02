clear; clc;

% Heat solver

% Important Parameters
m=20; % Newman poles per corner
TT=logspace(-2,4,20);
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

% Get oversampled grid and normals for flux calculation
N_s_halfedge=3*N_c_halfedge;
samps=[];
bdp=[]; % Boundary parametrization
nn_s=[];
cum_side=0;
for j=1:length(vv)
	N=N_s_halfedge(j);
	c1=samp_cluster_to(.5,0,N,m);
	c2=samp_cluster_to(.5,1,N,m);
	tt=[flip(c1),c2];
	samps=[samps,gamma{j}(tt)];
	n_direction=gamma{j}(1)-gamma{j}(0);
	side_length=abs(n_direction);
	n_direction=-1i*n_direction;
	n_direction=n_direction./abs(n_direction);
	nn_s=[nn_s,ones(size(tt))*n_direction];
	bdp=[bdp,cum_side+tt*side_length];
	cum_side=cum_side+side_length;
end

% Plot geometry, collocation and Newman poles
figure(1); clf;
fill(real(vv),imag(vv),'w','handlevisibility','off'); hold on; grid on;
axis padded equal;
plot(nm_pols,'.','displayname','Newman Poles');
plot(real(r_pols),imag(r_pols),'*','displayname','Runge Pole(s)');
plot(cols,'r.','displayname','Collocation Points');
legend

figure(2); clf


cumulative_fluxes=zeros(size(TT));
subindices=find(mod((1:length(TT)),3)==0);
tiledlayout(2,length(subindices))

for I_t = 1:length(TT)
	t=TT(I_t); % Time of heat solution

	% Setting up the problem
	source=2+2i;
	up=@(z,s)psi(0,z,source,s);
	
	uh_mhe=cell(1,M);
	D_uh_mhe=cell(1,M);
	[ww,ss]=ilt_quadrature(t,M);
	hsolve=tic;
	for j=1:M
		F=@(z)-up(z,ss(j));
		tic
		[uh_tmp,D_uh_tmp,residuals]=solve_mhe(cols,F,r_pols,r_orders,nm_pols,ss(j));
		fprintf('Modified Helmholtz problem %.f/%.f solved in %.2f\n',j,M,toc)
		uh_mhe{j}=uh_tmp;
		D_uh_mhe{j}=D_uh_tmp;
	end	
	fprintf('Heat problem solved in %.2f\n',toc(hsolve))
	
	% Plot particular/homogenous part and solution
	bds=[-4,4];
	%gridres=200;
	gridres=64;
	zz=linspace(bds(1),bds(2),gridres);
	Z=zz+1i*zz.';
	
	% Get helmholtz solutions
	U_mhe=cell(1,M); U_mhe_normals=cell(1,M); 
	heval=tic;
	for j=1:M
		tic
		U_mhe{j}=up(Z,ss(j))+uh_mhe{j}(Z);
		Up_mhe_normals = -real(nn_s).*Dx_psi(0,samps,source,ss(j)) ...
			- imag(nn_s).*Dy_psi(0,samps,source,ss(j));
		Uh_mhe_normals = D_uh_mhe{j}(samps,nn_s);
		
		% For the fluxes, we divide our normals for each MHE problem
		% by s, which gives us cumulative fluxes and is more numerically stable
		U_mhe_normals{j}=(Up_mhe_normals+Uh_mhe_normals)/ss(j); 
		fprintf('Modified Helmholtz problem %.f/%.f evaluated on grid in %.2f\n',j,M,toc)
	end
	fprintf('Heat problem evaluted on grid in %.2f\n',toc(heval))
	
	% Sum MHE solutions to get heat solution
	U=zeros(size(U_mhe{1})); U_normals=zeros(size(U_mhe_normals{1}));
	for j=1:M
		U=U+ww(j)*U_mhe{j};
		U_normals=U_normals+ww(j)*U_mhe_normals{j};
	end
	U=real(U/(M*1i));
	U_normals=real(U_normals/(M*1i));
	
	% The above as a (vectorized!) anonymous function
	u_heat=@(z)eval_heat(z,uh_mhe,up,ww,ss);
	
	if any(I_t==subindices)
		nexttile(I_t/3)
		I=pcolor(real(Z),imag(Z),abs(real(U))); set(I,'EdgeColor','none')
		colorbar; colormap jet; set(gca,'colorscale','log'); clim([1e-10,1])
		axis equal; axis([bds,bds]); title('$u=u_p+u_h$','interpreter','latex')
		hold on; fill(real(vv),imag(vv),'w')

		nexttile(I_t/3+length(subindices))
		semilogy(bdp,U_normals,'k-')
		grid on
		ylim([1e-5,1])
		xlim([min(bdp),max(bdp)])
	end

	cumulative_fluxes(I_t)=trapz(bdp,U_normals);
end

figure(3); clf;

semilogx(TT,cumulative_fluxes,'*-')
grid on

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

function out=Dx_psi(n,z,xi,s)
	dd=z-xi;
	arg=sqrt(s)*abs(dd);
	out=dd.^n ./ abs(dd).^n .* (...
		sqrt(s)*dbesselk(n,arg).* real(dd) ./ abs(dd) +...
		besselk(n,arg).*(n./dd - n*real(dd)./abs(dd).^2))/(2*pi);
end

function out=Dy_psi(n,z,xi,s)
	dd=z-xi;
	arg=sqrt(s)*abs(dd);
	out=dd.^n ./ abs(dd).^n .* (...
		sqrt(s)*dbesselk(n,arg).* imag(dd) ./ abs(dd) +...
		besselk(n,arg).*(1i*n./dd - n*imag(dd)./abs(dd).^2))/(2*pi);
end



function out = dbesselk(N,z)
	if N == 0
		out = -besselk(1,z);
	else
		out = -.5*(besselk(N-1,z)+besselk(N+1,z));
	end
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


function [uh,D_uh,residuals]=solve_mhe(col,F,r_pols,r_order,nm_pols,s)
	cc=col(:);
	col_val=F(cc);
	A=get_problem_matrix(cc,r_pols,r_order,nm_pols,s);
	normA=vecnorm(A); scaledA=A./normA;
	coefs=scaledA\col_val;
	
	uh=@(z)eval_solution(z,coefs./normA(:));
	D_uh=@(z,n)D_eval_solution(z,n,coefs./normA(:));
	residuals=scaledA*coefs-col_val;

	function out=eval_solution(z,coefs)
		zz=z(:);
		A=get_problem_matrix(zz,r_pols,r_order,nm_pols,s);
		U=A*coefs;
		out=reshape(U,size(z));
	end

	function out=D_eval_solution(z,n,coefs)
		zz=z(:);
		nn=n(:);
		A=get_D_problem_matrix(zz,nn,r_pols,r_order,nm_pols,s);
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

function out=get_D_problem_matrix(in_pts,nn,r_pols,r_order,nm_pols,s)
	zz=in_pts;
	nn=nn(:);
	nm_part=get_D_newman_part(zz,nn,nm_pols,s);
	rg_part=get_D_runge_part(zz,nn,r_pols,r_order,s);
	out=[nm_part,rg_part];
end

function out=get_D_runge_part(in_pts,nn,r_pols,r_order,s)
	zz=in_pts(:);
	
	ox=[];
	oy=[];
	for NR=length(r_pols)
		N=r_order(NR);
		Nvec=-N:N;
		xi=r_pols(NR);
		tmpx=zeros(length(in_pts),length(Nvec));
		tmpy=zeros(length(in_pts),length(Nvec));
		for j=1:length(Nvec);
			tmpx(:,j)=Dx_psi(Nvec(j),zz,xi,s);
			tmpy(:,j)=Dy_psi(Nvec(j),zz,xi,s);
		end
		ox=[ox,tmpx];
		oy=[oy,tmpy];
	end
	out=[real(nn).*ox+imag(nn).*oy];
end


function out=get_D_newman_part(in_pts,nn,nm_pols,s)
	zz=in_pts(:);
	xi=nm_pols;

	ox=[Dx_psi(-1,zz,xi,s),Dx_psi(1,zz,xi,s)];
	oy=[Dy_psi(-1,zz,xi,s),Dy_psi(1,zz,xi,s)];

	out=[real(nn).*ox+imag(nn).*oy];
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

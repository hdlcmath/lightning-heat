clear; clc;

m=50; % Newman poles per corner

vv=[1+1i,2i,-1+1i,-1-1i,1-1i]; % Vertices
r_pols=0;
gamma=cell(1,length(vv));

% Create parametrization of each side
for j=1:length(vv)
	a=vv(j);
	jplus=mod(j,length(vv))+1;
	b=vv(jplus);
	gamma{j}=@(t)a+(b-a)*t;
end
col_factor=5;
runge_factor=3.5;
N_r=ceil(runge_factor*sqrt(m)); % Set Runge order to be O(sqrt(m))
N_nm=length(vv)*m;
N_terms=2*(N_nm+N_r*length(r_pols))+1;
N_dofs=2*N_terms;
N_c=col_factor*N_dofs;
runge_orders=N_r*ones(1,length(vv));

% Get the Newman poles
edge_lengths=abs(vv-vv([2:length(vv),1]));
r=0.7071*min(edge_lengths);
nm_pols=[];
for j=1:length(vv)
	jplus=mod(j,length(gamma))+1;
	epsilon=1e-5;

	vbwd=gamma{j}(1-epsilon)-vv(jplus);
	vfwd=gamma{jplus}(epsilon)-vv(jplus);

	outward=1i*vbwd*sqrt(-vfwd/vbwd);
	outward=r*outward./abs(outward);
	vang=vfwd/vbwd; beta=abs(angle(vang))/pi+1;

	b=vv(jplus);
	a=vv(jplus)-outward;

	tmp_pols=nm_cluster_to(a,b,m,beta); 
	nm_pols=[nm_pols,tmp_pols];
end

% Get collocation points
N_c_halfedge=ceil(N_c*edge_lengths/max(edge_lengths)/2);
cols=[];
for j=1:length(vv)
	N=N_c_halfedge(j);
	c1=col_cluster_to(.5,0,N,m);
	c2=col_cluster_to(.5,1,N,m);
	tt=[flip(c1),c2];
	cols=[cols,gamma{j}(tt)];
end

% Setting up the problem

source=2+1i;
up=@(z,s)psi(0,z,source,s)


t=0.3;
[ww,ss]=ilp_getquadrature


%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function out=psi(n,z,xi,s)
	dd=z-xi;
	arg=alpha*abs(dd);
	out=besselk(abs(n),arg).*dd.^n./abs(dd).^n/(2*pi);
end

function [solution,residuals]=solve(col,F,r_pols,r_order,nm_pols,s)
	cc=col(:);
	col_val=F(col(:)).';
	A=get_problem_matrix(cc,r_pols,r_order,nm_pols,s);
	keyboard

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
	alpha=sqrt(s);
	
	out=[];
	for NR=length(r_pols)
		N=r_order(NR);
		Nvec=-N:N;
		xi=r_pols(NR);
		tmp=zeros(length(in_pts),length(Nvec));
		for j=1:length(Nvec);
			tmp(:,j)=psi(Nvec(j),zz,xi,alpha);
		end
		out=[out,tmp];
	end
end


function out=get_newman_part(in_pts,nm_pols,s)
	zz=in_pts(:);
	xi=nm_pols;
	alpha=sqrt(s);

	o1=psi(1,zz,xi,alpha);
	oneg1=psi(-1,zz,xi,alpha);

	out=[oneg1,o1];
end

function out = nm_cluster_to(a,b,N,beta)
	sg=sqrt(2*(2-beta)*beta)*pi; % from Herremans paper

	% Tapered Lightning
	nn=1:N;
	rate=exp(-sg*(sqrt(N+1)-sqrt(nn)));

	rate(rate==0)=[];
	tmp=(b-a)*(1-rate)+a;
	% inds = abs(b - tmp) > 1e-9;
	out=tmp;
end


function out = col_cluster_to(a,b,N,m)
	sg=sqrt(2*(m+1))*pi;
	v=linspace(0,1,N);
	rate=1-exp(-sg*v);
	tmp = (b-a)*rate+a;
	out = tmp;
end

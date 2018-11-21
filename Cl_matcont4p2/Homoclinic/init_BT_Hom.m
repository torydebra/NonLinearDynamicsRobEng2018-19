function [x,v] = init_BT_Hom(odefile, x, s, p, ap, ntst, ncol, eps, extravec)
global homds cds

% check input
n_par = size(ap,2);
if n_par ~= 2
    error('2 free system parameters are needed');
end

if extravec==zeros(1,3)
    error('at least 1 free homoclinic parameter is needed');
end

if extravec==ones(1,3)
    error('at most 2 free homoclinic parameters are needed');
end

eps = eps^2;

for i=1:length(s)
    if ~isempty(s.data) && isfield(s(i).data,'evec')
        nph = size(s(i).data.evec,2);
    elseif isfield(s(i).data,'x')
        nph = size(s(i).data.x,2);
    else
        
    end
end

if isempty(cds) || ~isfield(cds,'options')
    cds.options = contset();
end
cds.curve = @homoclinic;
curvehandles = feval(cds.curve);
cds.curve_func = curvehandles{1};
cds.curve_jacobian = curvehandles{4};
cds.curve_hessians = curvehandles{5};

homds = [];
% initialize homds
eps0=1e-2;
eps1=1e-2;
T = 100;
init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph);

func_handles = feval(odefile);
symord = 0; 
symordp = 0;

if     ~isempty(func_handles{9}),   symord = 5; 
elseif ~isempty(func_handles{8}),   symord = 4; 
elseif ~isempty(func_handles{7}),   symord = 3; 
elseif ~isempty(func_handles{5}),   symord = 2; 
elseif ~isempty(func_handles{3}),   symord = 1; 
end
if     ~isempty(func_handles{6}),   symordp = 2; 
elseif ~isempty(func_handles{4}),   symordp = 1; 
end
cds.options = contset(cds.options, 'SymDerivative', symord);
cds.options = contset(cds.options, 'SymDerivativeP', symordp);
cds.symjac = 1;
cds.symhess = 0;


homds.odefile = odefile;
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3 = func_handles{7};
homds.Der4 = func_handles{8};
homds.Der5 = func_handles{9};

cds.oldJac = [];
cds.oldJacX = [];
xp = [x;p(ap)];
cds.ndim = length(xp);

homds.x0 = x(1:homds.nphase);
    
pcell = num2cell(p);
% Compute A, A1
A = cjac(homds.func,homds.Jacobian,homds.x0,pcell,homds.ActiveParams);
A1 = cjacp(homds.func,homds.JacobianP,homds.x0,pcell,ap);
    
% Compute B, B1 and B2
B = chess(homds.func,homds.Jacobian,homds.Hessians,homds.x0,pcell,homds.ActiveParams);
B1 = chessp(homds.func,homds.Jacobian,homds.HessiansP,homds.x0,pcell,homds.ActiveParams);
if size(B1,3) > length(ap)
    B1 = B1(:,:,ap);
end

for i=homds.ActiveParams
    p1 = pcell; p1{i} = p1{i}-cds.options.Increment;
    p2 = pcell; p2{i} = p2{i}+cds.options.Increment;
    Hjp2 = cjacp(homds.func,homds.JacobianP,homds.x0,p2,ap); 
    Hjp1 = cjacp(homds.func,homds.JacobianP,homds.x0,p1,ap); 
    tmph(:,:,i) = Hjp2 - Hjp1;
end
tmph = tmph/(2*cds.options.Increment);
B2 = tmph;
if size(B2,3) > length(ap)
    B2 = B2(:,:,ap);
end
    
  [X,D] = eig(A);
  index1 = find(abs(diag(D)) < 1e-3);%If ok, index1 is 1x2 array otherwise
  vext = real(X(:,index1(1)));
  [X,D] = eig(A');
  index1 = find(abs(diag(D)) < 1e-3);
  wext = real(X(:,index1(1)));
  Bord = [ A wext; vext' 0];
  bunit=[zeros(homds.nphase,1);1];
  q0=Bord\bunit; 
  q0=q0(1:homds.nphase);          % A q0 = 0, <vext,q0> = 1
  p1=Bord'\bunit;
  p1=p1(1:homds.nphase);          % A'p1 = 0, <wext,p1> = 1
  Bord = [ A p1; q0' 0];
  q1 = Bord\[q0; 0];		
  q1 = q1(1:homds.nphase);		% A q1 = q0, <q0,q1> = 0
  p0 = Bord'\[p1; 0];
  p0 = p0(1:homds.nphase);		% A'p0 = p1, <p0,p1> = 0

% normalize so that <p0,q0>=<p1,q1>=1, <p0,q1>=<p1,q0>=0 and <q0,q0>=1, <q0,q1>=0
  mu = sqrt(abs(q0'*q0));
  q0 = (1/mu)*q0;
  q1 = (1/mu)*q1;
  q1 = q1 - (q0'*q1)*q0;
  nu = q0'*p0;
  p1 = (1/nu)*p1;
  p0 = p0 - (p0'*q1)*p1;
  p0 = (1/nu)*p0;
  
% Compute a
a = zeros(1,homds.nphase);
for i=1:homds.nphase
    a(i) = p1' * B(:,:,i) * q0;
end    
a = 1/2 * a * q0;
    
% K1, H01 and K2
gamma = (p1' * A1);
K1 = 1/(gamma(1)^2 + gamma(2)^2) * [gamma(1)  -gamma(2);  gamma(2)   gamma(1)];
tmpH = [q1 zeros(size(q1))] - A1 * K1;
H01 = Bord \ [tmpH; zeros(1,size(tmpH,2))];
H01 = H01(1:end-1,:);
for i=1:homds.nphase
    p1B(:,i) = p1' * B(:,:,i) * H01(:,2);
end
for i=1:length(ap)
    p1B1(:,i) = p1' * B1(:,:,i) * H01(:,2);
    p1B2(:,i) = p1' * B2(:,:,i) * K1(:,2);
end
p1B = p1B * H01(:,2);
p1B1 = p1B1 * K1(:,2);
p1B2 = p1B2 * K1(:,2);
K2 = -(p1B + 2 * p1B1 + p1B2) * K1(:,1);
    
% Actual approximation: alpha 
alpha = -5/7 * eps^2 * K1(:,2) + eps^4/49 * (-6 * K1(:,1) + 12.5 * K2);
p(ap) = p(ap) + alpha; %   tmpfreep = tmpfreep + alpha;
homds.P0 = p;

homds.x0 = homds.x0 + eps^2 * (-5/7*H01(:,2) + 6/(7*a) * q0);
x0 = homds.x0;
tmpmiddle = eps^2 * (-5/7*H01(:,2) + 1/(4*a) * ...
        (2 - 6 + 10/7) * q0) + homds.x0;

dist = norm(tmpmiddle - x0);

delta0 = dist / 200;
myarg = (delta0 * norm(2 * a)) / (3 * eps^2);
myarg = sqrt(myarg);
myargn = log(1/myarg + ((1/myarg)^2 - 1));
if ~isempty(imag(myargn)) && (imag(myargn) ~= 0)
    myargn = log(1/myarg - ((1/myarg)^2 - 1));
end
homds.T = 2 * myargn / eps;


% Actual approximation: cycle
for i=1:length(homds.finemsh)
    % conversion from [0,1] -> [-T,T]
    t = (2*homds.finemsh(i) - 1) * (homds.T);
        
    ups(:,i) = eps^2 * (-5/7*H01(:,2) + 1/(4*a) * ...
        (2 - 6 * (sech(t * eps/2))^2 + 10/7) * q0) + eps^3 / (8*a) *...
        (12 * (sech(t * eps/2))^2 * tanh(t * eps/2)) * q1;
    ups(:,i) = ups(:,i) + homds.x0;
end
% Actual approximation: equilibrium
    
Hom_calc_weights;

A = cjac(homds.func,homds.Jacobian,x0,num2cell(p),homds.ActiveParams);
D = eig(A);
% nneg = dimension of stable subspace
homds.nneg = sum(real(D) < 0);
homds.npos = homds.nphase-homds.nneg;
homds.Ysize = homds.nneg*homds.npos;

homds.eps0 = norm(ups(:,1) - homds.x0);
homds.eps1 = norm(ups(:,end) - homds.x0);

% COMPOSE X0
% ----------

% 1. cycle 
x1 = reshape(ups,size(ups,1)*size(ups,2),1);
v = []; 
[x1,v]=Hom_new_mesh(x1,v,ntst,ncol);

% 2. equilibrium coordinates
x1 = [x1; x0];
% 3. (two) free parameters
x1 = [x1; homds.P0(homds.ActiveParams)];
% 4. extra free parameters
extravec = [homds.T; homds.eps0; homds.eps1];
x1 = [x1; extravec(find(homds.extravec))];
% 5. YS and YU, initialized to 0
for i=1:homds.nneg
    x1 = [x1; zeros(homds.npos,1)];
end
for i=1:homds.npos
    x1 = [x1; zeros(homds.nneg,1)];
end

% ASSIGN SOME VALUES TO HOMOCLINIC FIELDS
% ---------------------------------------

homds.YS = zeros(homds.npos,homds.nneg);
homds.YU = zeros(homds.nneg,homds.npos);

% Third parameter = unstable_flag, 
% 1 if we want the unstable space, 0 if we want the stable one
[QU, se] = computeBase(A,0,homds.npos);
[QS, se] = computeBase(A,1,homds.nneg);

homds.oldStableQ = QS;
homds.oldUnstableQ = QU;
homds.ups = [];
homds.ndim = length(x1);
cd.ndim = homds.ndim;

x = x1;
v=[];

%-----------------------------------------------------------------
function init_homds(odefile,x,p,ap,ntst,ncol,extravec,T,eps0,eps1,nph)
global homds 
homds.odefile = odefile;
func_handles = feval(homds.odefile);
homds.func = func_handles{2};
homds.Jacobian  = func_handles{3};
homds.JacobianP = func_handles{4};
homds.Hessians  = func_handles{5};
homds.HessiansP = func_handles{6};
homds.Der3=[];
siz = size(func_handles,2);
if siz > 9
    j=1;
    for k=10:siz
        homds.user{j}= func_handles{k};
        j=j+1;
    end
else homds.user=[];
end
homds.nphase = nph;
homds.ActiveParams = ap;
homds.P0 = p;
homds.extravec = extravec;
Hom_set_ntst_ncol(ntst,ncol,(0:ntst)/ntst);
homds.T = T;
homds.eps0 = eps0;
homds.eps1 = eps1;
homds.cols_p1 = 1:(homds.ncol+1);
homds.cols_p1_coords = 1:(homds.ncol+1)*homds.nphase;
homds.ncol_coord = homds.ncol*homds.nphase;
homds.col_coords = 1:homds.ncol*homds.nphase;
homds.pars = homds.ncoords+(1:3);
homds.phases = 1:homds.nphase;
homds.ntstcol = homds.ntst*homds.ncol;
homds.wp = kron(homds.wpvec',eye(homds.nphase));
homds.pwwt = kron(homds.wt',eye(homds.nphase));
homds.pwi = homds.wi(ones(1,homds.nphase),:);

homds.bialt_M1 = [];
homds.bialt_M2 = [];
homds.bialt_M3 = [];
homds.bialt_M4 = [];
homds.multipliers = nan;
homds.monodromy = [];
homds.multi_r1 = [];
homds.multi_r2 = [];
homds.ups = [];
homds.vps = [];
homds.tsts = 1:homds.ntst;
homds.cols = 1:homds.ncol;

homds.HTPstep = 0;
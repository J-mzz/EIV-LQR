function [p_psatz, cons_psatz, Gram, Coef] = psatz_xnoise_putinar(p, C, d, vars, A, n, T, B)
% PSATZ_XNOISE_PUTINAR  Putinar Positivstellensatz certificate for the
% scalarized condition z'*F(A,B)*z >= 0 on the consistency set times the
% sphere z'*z = 1, when only the state measurements are noisy (used by
% EIVLQR_xnoise_putinar).

% nAB = length(vars) - 2*n;
% vars_AB = vars(1:nAB);
vars_AB = recover(depends([A(:);B(:)]));
vars_AB = vars_AB(:).';
nAB = length(vars_AB);

n_g = length(C.ineq);
n_h = length(C.eq);

n_var = length(vars);

d_g = 0; %degree(C.ineq);
d_h = 1; %degree(C.eq);

p_psatz = 0;

Gram = {};
Coef = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sigma0 SOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pow0,~] = momentPowers(0,n_var,d);
v0 = recovermonoms(pow0,vars);
l0 = length(v0);

Gram{1} = sdpvar(l0);

sigma0 = v0'*Gram{1}*v0;

cons_psatz = (Gram{1} >= 0):'Gram1';

p_psatz = p_psatz + sigma0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lambda_z (sphere equality multiplier)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hz = C.zeq;

[pow_lz,~] = momentPowers(0,n_var,max(d-1,0));
vlz = recovermonoms(pow_lz,vars);

Coef_lambda = sdpvar(length(vlz),1,'full');

lambda_z = vlz' * Coef_lambda;

p_psatz = p_psatz + lambda_z * hz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta multipliers for inequalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(d_g,2)==0
    [pow_zeta,~] = momentPowers(0,nAB,floor((2*d-d_g)/2));
else
    [pow_zeta,~] = momentPowers(0,nAB,floor((2*d-d_g-1)/2));
end

vzet = recovermonoms(pow_zeta,vars_AB);
lzet = length(vzet);

zeta_p = cell(n_g/2,1);
zeta_n = cell(n_g/2,1);

for i = 1:n_g/2

    Gram{end+1} = sdpvar(lzet);

    zeta_p{i} = vzet'*Gram{end}*vzet;

    cons_psatz = [cons_psatz; Gram{end} >= 0];

    p_psatz = p_psatz + zeta_p{i} * C.ineq(i);

end

zeta_p = reshape(zeta_p,n,T);

for i = 1:n_g/2

    Gram{end+1} = sdpvar(lzet);

    zeta_n{i} = vzet'*Gram{end}*vzet;

    cons_psatz = [cons_psatz; Gram{end} >= 0];

    p_psatz = p_psatz + zeta_n{i} * C.ineq(i+n*T);

end

zeta_n = reshape(zeta_n,n,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% equality multipliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pow_mu,~] = momentPowers(0,nAB,max(2*d-d_h,0));

vmu = recovermonoms(pow_mu,vars_AB);

lmu = length(vmu);

mat_mu = cell(n_h,1);

for i = 1:n_h

    Coef{i} = sdpvar(lmu,1,'full');

    mat_mu{i} = vmu'*Coef{i};

    p_psatz = p_psatz + mat_mu{i} * C.eq(i);

end

mat_mu = reshape(mat_mu,n,T-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefficient matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cons_psatz = [cons_psatz; coefficients(p - p_psatz, vars) == 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% same coupling constraints as Scherer version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n

    Amu = 0;

    for j = 1:n
        Amu = Amu - mat_mu{j,1}*A(j,i);
    end

    cons_psatz = [cons_psatz;
        coefficients(zeta_p{i,1} - zeta_n{i,1} + Amu,vars) == 0];

end

for k = 2:T-1

    for i = 1:n

        Amu = 0;

        for j = 1:n
            Amu = Amu - mat_mu{j,k}*A(j,i);
        end

        Amu = Amu + mat_mu{i,k-1};

        cons_psatz = [cons_psatz;
            coefficients(zeta_p{i,k} - zeta_n{i,k} + Amu,vars) == 0];

    end

end

for i = 1:n

    cons_psatz = [cons_psatz;
        coefficients(zeta_p{i,end} - zeta_n{i,end} + mat_mu{i,end},vars) == 0];

end

end

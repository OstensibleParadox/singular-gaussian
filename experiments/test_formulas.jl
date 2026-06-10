using LinearAlgebra

N = 128
kappa = 1.0
mu = [pi^2 * j^2 + kappa^2 for j in 1:N]
A_prior = Diagonal(-1 ./ mu)
D_prior = Diagonal(fill(1.0, N))

c = zeros(N)
for j in 1:N
    if j % 2 != 0
        c[j] = 4 * sqrt(2.0) / (pi^3 * j^3)
    end
end
v = c ./ norm(c)
P_E = v * v'
P_I = I - P_E

A_II = P_I * A_prior * P_I
D_II = P_I * D_prior * P_I

# Q0 solves A_II Q0 + Q0 A_II + 2 D_II = 0
# Since A_II is symmetric, Q0 and A_II share eigenvectors? Not necessarily because D_II is P_I.
# Wait, A_II and D_II share the same nullspace (v).
# On the subspace I, D_II = I.
# So Q0 = (-A_II)^-1 on the subspace I!
# Let's verify this!
F = eigen(Symmetric(A_II))
lam = F.values
U = F.vectors

Q0 = similar(A_II)
C_tilde = 2 * (U' * D_II * U)
for i in 1:N
    for j in 1:N
        if abs(lam[i]) < 1e-10 && abs(lam[j]) < 1e-10
            Q0[i, j] = 0.0
        else
            Q0[i, j] = - C_tilde[i, j] / (lam[i] + lam[j])
        end
    end
end
Q0 = U * Q0 * U'

F2 = eigen(Symmetric(Q0))
inv_lam2 = [abs(l) > 1e-10 ? 1/l : 0.0 for l in F2.values]
Q0_inv = F2.vectors * Diagonal(inv_lam2) * F2.vectors'

A_IE_v = P_I * A_prior * v
limit_val = dot(A_IE_v, Q0_inv * A_IE_v)

println("Empirical limit (A_IE^T Q0^-1 A_IE): ", limit_val)

# User's formula: sum( c_j^2 / mu_j^3 )
limit_user = sum(v.^2 ./ mu.^3)
println("User formula limit: ", limit_user)

# My formula: sum v_j^2 Lambda_j^3 - ...
Lambda = 1 ./ mu
bar_Lambda = dot(v, Lambda .* v)
limit_mine = sum(v.^2 .* Lambda.^3) - 2*bar_Lambda * sum(v.^2 .* Lambda.^2) + bar_Lambda^2 * sum(v.^2 .* Lambda)
println("My analytical formula limit: ", limit_mine)

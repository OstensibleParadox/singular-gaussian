# ==============================================================================
# 🧪 Experiment 6: Asymptotic Coupling and Fredholm Finite Part
# ==============================================================================
import Pkg
Pkg.add("GenericLinearAlgebra")

using LinearAlgebra
using GenericLinearAlgebra
using Printf

function run_experiment()
    setprecision(BigFloat, 256)
    
    N = 128
    kappa = BigFloat(1.0)
    omega = BigFloat(1.0)
    d_bg = BigFloat(1.0)
    Delta_E = BigFloat(1.0)
    
    # 1. 谱基底与特征值
    mu = [BigFloat(pi)^2 * j^2 + kappa^2 for j in 1:N]
    
    # 构建漂移与扩散算子
    A_prior = Diagonal(-1 ./ mu)
    D_prior = Diagonal(fill(d_bg, N))
    
    # 2. 不对齐目标子空间 E
    c = zeros(BigFloat, N)
    for j in 1:N
        if j % 2 != 0
            c[j] = 4 * sqrt(BigFloat(2)) / (BigFloat(pi)^3 * j^3)
        end
    end
    
    v = c ./ norm(c)
    
    P_E = v * v'
    P_I = I - P_E
    
    # 3. 理论预测值计算 (包含 P_I 子空间正交化带来的精确相消)
    Lambda = 1 ./ mu
    bar_Lambda = dot(v, Lambda .* v)
    exact_limit_val = sum((v[j]^2 * Lambda[j]^3) for j in 1:N) - 2 * bar_Lambda * sum((v[j]^2 * Lambda[j]^2) for j in 1:N) + bar_Lambda^2 * sum((v[j]^2 * Lambda[j]) for j in 1:N)
    C_theory = (omega / d_bg) * Delta_E * exact_limit_val
    
    println("======================================================")
    println("🧪 实验设计：无穷维谱解析下的迹反常级数验证")
    println("======================================================")
    @printf("截断维度 N = %d (BigFloat 256-bit)\n", N)
    @printf("理论极限 lim λ^4 tr(ρ_E): %.10e\n\n", Float64(C_theory))
    
    # 预计算固定矩阵
    A_II = P_I * A_prior * P_I
    D_II = P_I * D_prior * P_I
    A_IE_v = P_I * A_prior * v
    
    println(">>> 正在对 A_II 进行特征值分解 (这在 N=2048 时可能需要数分钟)...")
    flush(stdout)
    
    time_start = time()
    F = eigen(Symmetric(A_II))
    lam = F.values
    U = F.vectors
    println(">>> 特征值分解完成。耗时: ", round(time() - time_start, digits=2), " 秒。开始迭代 λ...\n")
    flush(stdout)
    
    lambdas = [BigFloat(10)^k for k in 1:5]
    
    println(" λ            | λ^4 tr(ρ_E) 观测值 | 相对误差 ")
    println("------------------------------------------------------")
    
    # Pre-allocate for optimization
    u_D_II_u = U' * D_II * U
    u_A_IE_v = U' * A_IE_v
    
    for lambda in lambdas
        Sigma_EE = Delta_E / lambda^2
        
        # 交叉协方差 Sigma_IE
        inv_lam = 1 ./ (lambda .- lam)
        Sigma_IE = U * (inv_lam .* (U' * (A_IE_v * Sigma_EE)))
        
        # 为了高效解 Sigma_II，利用预计算的 U' * D_II * U
        u_Sigma_IE = U' * Sigma_IE
        C_tilde = 2 * u_D_II_u + u_A_IE_v * u_Sigma_IE' + u_Sigma_IE * u_A_IE_v'
        
        # 直接在特征基下求解 X_tilde
        X_tilde = similar(C_tilde)
        for i in 1:N
            for j in 1:N
                if abs(lam[i]) < 1e-10 && abs(lam[j]) < 1e-10
                    X_tilde[i, j] = 0.0
                else
                    X_tilde[i, j] = - C_tilde[i, j] / (lam[i] + lam[j])
                end
            end
        end
        Sigma_II = U * X_tilde * U'
        
        # 计算短路残差 R_E 和 ρ_E
        F2 = eigen(Symmetric(Sigma_II))
        lam2 = F2.values
        V2 = F2.vectors
        
        inv_lam2 = zeros(BigFloat, N)
        for i in 1:N
            if abs(lam2[i]) > 1e-20
                inv_lam2[i] = 1 / lam2[i]
            end
        end
        
        R_E = dot(Sigma_IE, V2 * (inv_lam2 .* (V2' * Sigma_IE)))
        rho_E = R_E / Sigma_EE
        
        val = lambda^4 * rho_E
        rel_err = abs(val - C_theory) / C_theory
        
        @printf(" %.1e      | %.10e    | %.2e\n", 
                Float64(lambda), Float64(val), Float64(rel_err))
        flush(stdout)
    end
end

run_experiment()

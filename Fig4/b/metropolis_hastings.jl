using Printf
include("likelihood.jl")
function MH(e,N_range,p_range,q_range,sigma_range,tauD_range,delta_t,nb_it)

    # Sample vectors (final results)
    N       = zeros(nb_it)
    p       = zeros(nb_it)
    q       = zeros(nb_it)
    sigma   = zeros(nb_it)
    tau_D   = zeros(nb_it)
    L = 0
    L_max = 0
    idx = zeros(Int32,5)
    N_MLE       = 0
    p_MLE       = 0
    q_MLE       = 0
    sigma_MLE   = 0
    tau_D_MLE   = 0
    u = 0

    # When the initial parameters are very different from the ground truth parameters
    # (and especially when sigma is low), their likelihood will be zero.
    # To avoid running unnecessarily the algo, this loop ensures that the
    # initial parameters are valid.
    allowed_init = false
    while allowed_init == false
        # Initializations
        # Start MH at a random value
        idx_init = [rand(1:length(N_range)),rand(1:length(p_range)),rand(1:length(q_range)),rand(1:length(sigma_range)),rand(1:length(tauD_range))]
        N[1]        = N_range[idx_init[1]]
        p[1]        = p_range[idx_init[2]]
        q[1]        = q_range[idx_init[3]]
        sigma[1]    = sigma_range[idx_init[4]]
        tau_D[1]    = tauD_range[idx_init[5]]
        # Initialization of the likelihood at the MLE
        L_max = -Inf
        idx = idx_init

        # Current value of the likelihood
        L = ll_binomial_tau_D(e,Int32(N[1]),p[1],q[1],sigma[1],tau_D[1],delta_t)
        if !isinf(L)
            allowed_init = true
        end
    end

    for i in 2:nb_it
        @printf "Progress = %0.3f\n" (i/nb_it)*100

        # Updates values for the MLE
        if L > L_max
            N_MLE       = N[i-1]
            p_MLE       = p[i-1]
            q_MLE       = q[i-1]
            sigma_MLE   = sigma[i-1]
            tau_D_MLE   = tau_D[i-1]
            L_max = L
        end

        # Draw the next index

        idx_next    = zeros(Int32, length(idx))
        allowed = false
        while allowed == false

            if idx[1] == length(N_range)
                idx_next[1] = idx[1]+[-1,0][rand(1:2)]
            elseif idx[1] == 1
                idx_next[1] = idx[1]+[1,0][rand(1:2)]
            else
                idx_next[1] = idx[1]+[-1,0,1][rand(1:3)]
            end

            if idx[2] == length(p_range)
                idx_next[2] = idx[2]+[-1,0][rand(1:2)]
            elseif idx[2] == 1
                idx_next[2] = idx[2]+[1,0][rand(1:2)]
            else
                idx_next[2] = idx[2]+[-1,0,1][rand(1:3)]
            end

            if idx[3] == length(q_range)
                idx_next[3] = idx[3]+[-1,0][rand(1:2)]
            elseif idx[3] == 1
                idx_next[3] = idx[3]+[1,0][rand(1:2)]
            else
                idx_next[3] = idx[3]+[-1,0,1][rand(1:3)]
            end

            if idx[4] == length(sigma_range)
                idx_next[4] = idx[4]+[-1,0][rand(1:2)]
            elseif idx[4] == 1
                idx_next[4] = idx[4]+[1,0][rand(1:2)]
            else
                idx_next[4] = idx[4]+[-1,0,1][rand(1:3)]
            end

            if idx[5] == length(tauD_range)
                idx_next[5] = idx[5]+[-1,0][rand(1:2)]
            elseif idx[5] == 1
                idx_next[5] = idx[5]+[1,0][rand(1:2)]
            else
                idx_next[5] = idx[5]+[-1,0,1][rand(1:3)]
            end


            if (idx_next[1] in 1:length(N_range)) && (idx_next[2] in 1:length(p_range)) && (idx_next[3] in 1:length(q_range)) && (idx_next[4] in 1:length(sigma_range)) && (idx_next[5] in 1:length(tauD_range))
                allowed = true
            end
        end

        N_next      = N_range[idx_next[1]]
        p_next      = p_range[idx_next[2]]
        q_next      = q_range[idx_next[3]]
        sigma_next  = sigma_range[idx_next[4]]
        tau_D_next  = tauD_range[idx_next[5]]

        L_next = ll_binomial_tau_D(e,Int32(N_next),p_next,q_next,sigma_next,tau_D_next,delta_t)

        if isinf(L_next) || isnan(L_next)
            N[i]        = N[i-1]
            p[i]        = p[i-1]
            q[i]        = q[i-1]
            sigma[i]    = sigma[i-1]
            tau_D[i]    = tau_D[i-1]
        else
            alpha = L_next-L
            u = log(rand())

        # Updates the values of the likelihood and the index
            if u <= alpha
                N[i]        = N_next
                p[i]        = p_next
                q[i]        = q_next
                sigma[i]    = sigma_next
                tau_D[i]    = tau_D_next
                idx         = idx_next
                L           = L_next
            else
                N[i]        = N[i-1]
                p[i]        = p[i-1]
                q[i]        = q[i-1]
                sigma[i]    = sigma[i-1]
                tau_D[i]    = tau_D[i-1]
            end
        end
    end

    return N,p,q,sigma,tau_D,N_MLE,p_MLE,q_MLE,sigma_MLE,tau_D_MLE,L_max
end

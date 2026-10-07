# Calibration of the climate-only prototype in the emulator (Nelder-Mead)

function nelder_mead(f, x0, step; maxevals = 800, tol = 1e-5)
    n = length(x0)
    X = [copy(x0)];
    for j in 1:n
        ;
        x = copy(x0);
        x[j] += step[j];
        push!(X, x);
    end
    F = [f(x) for x in X];
    evals = n + 1
    while evals < maxevals
        o = sortperm(F);
        X, F = X[o], F[o]
        (F[end] - F[1]) < tol * (abs(F[1]) + tol) && break
        c = sum(X[1:n]) / n
        xr = c + (c - X[end]);
        fr = f(xr);
        evals += 1
        if fr < F[1]
            xe = c + 2 * (c - X[end]);
            fe = f(xe);
            evals += 1
            X[end], F[end] = fe < fr ? (xe, fe) : (xr, fr)
        elseif fr < F[n]
            X[end], F[end] = xr, fr
        else
            xc = fr < F[end] ? c + 0.5 * (xr - c) : c + 0.5 * (X[end] - c);
            fc = f(xc);
            evals += 1
            if fc < min(fr, F[end])
                X[end], F[end] = xc, fc
            else
                for j in 2:(n + 1)
                    ;
                    X[j] = X[1] + 0.5 * (X[j] - X[1]);
                    F[j] = f(X[j]);
                end;
                evals += n
            end
        end
    end
    o = sortperm(F);
    return X[o[1]], F[o[1]], evals
end

const NTREE = 4
cfg_of(θ) = Config(
    z_tree = exp(θ[1]),
    z_grass = exp(θ[2]),
    sigma = max(θ[3], 0.3),
    f0_max = clamp(θ[4], 0.3, 1.0),
    w_humid = max(θ[5], 0.0),
    ever_floor = clamp(θ[6], 0.0, 1.0),
)
function tree_of(θ, cfg)
    β = θ[7:(6 + NTREE)];
    L = LtreeW(cfg)
    return [
        logistic(
            β[1] + β[2] * L[i] + β[3] * FEAT[i].dry_all + β[4] * FEAT[i].Tgs,
        ) for i in eachindex(PTS)
    ]
end
function simulate(θ, idx)
    cfg = cfg_of(θ);
    tree = tree_of(θ, cfg)
    return [emulate(CW[i], FW[i], tree[i], everE[i], 1.0, cfg) for i in idx],
    tree
end
function losses(θ, idx)
    sims, tree = simulate(θ, idx)
    S = reduce(vcat, (s' for s in sims));
    w = WGT[idx]
    mse_lai = wmean(vec(mean((S .- OBS[idx, :]) .^ 2; dims = 2)), w)
    v = filter(i -> isfinite(CLMTREE[i]), idx)
    mse_tree = wmean((tree[v] .- CLMTREE[v]) .^ 2, WGT[v])
    return mse_lai, mse_tree
end

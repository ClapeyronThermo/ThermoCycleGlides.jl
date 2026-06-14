module ThermoCycleGlidesMetaheuristicsExt

using ThermoCycleGlides
using Metaheuristics
using ThermoCycleGlides.CommonSolve

function CommonSolve.solve(prob::ThermoCycleGlides.ThermoCycleProblem,
    alg::Metaheuristics.AbstractAlgorithm;
    N = 20)

    ℓ = ThermoCycleGlides._build_least_squares_objective(prob,N)
    lb,ub = ThermoCycleGlides.generate_box_solve_bounds(prob)
    x01 = ThermoCycleGlides.generate_initial_point(prob,lb,ub,:average)
    x02 = ThermoCycleGlides.generate_initial_point(prob,lb,ub,:default)
    Metaheuristics.set_user_solutions!(alg, x01, ℓ)
    Metaheuristics.set_user_solutions!(alg, x02, ℓ)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    opt_result = Metaheuristics.optimize(ℓ,bounds,alg)
    
    x_best = Metaheuristics.minimizer(opt_result)

    loss_opt_M = Metaheuristics.minimum(opt_result)
    status_code = Symbol(Metaheuristics.TerminationStatusCode(opt_result))
    f_calls = Metaheuristics.nfes(opt_result)
    residuals = NaN .* x_best #don't calculate residuals
    result = SolutionState(
        x_best,                 #x::Vector{T}
        f_calls,                #f_calls::I
        opt_result.iteration,   #iterations::I
        residuals,              #residuals::Vector{T}
        lb,                     #lb::Vector{T}
        ub,                     #ub::Vector{T}
        false,                  #autodiff::Bool
        0,                      #fd_order::I
        0.0,                    #lenx::T
        loss_opt_M,             #lenf::T
        status_code             #soltype::Symbol
    )
end

end #module
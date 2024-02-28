function _net_present_value(r::Float64, model::MyUSTreasuryZeroCouponBondModel, compounding::DiscreteCompoundingModel)

    # get data from the model -
    T = model.T;
    price = model.price; # the price is set, we are looking for the interest rate
    Vₚ = model.par

    # we are passing in the rate -
    rate = r;

    # compute the discount with this rate -
    𝒟 = (1+rate)^(T)
    future_payout = (1/𝒟)*Vₚ
   
    # compute the npv value -
    npv_value = (future_payout - price)

    # return -
    return npv_value
end

function _net_present_value(r::Float64, model::MyUSTreasuryZeroCouponBondModel, compounding::ContinuousCompoundingModel)

    # get data from the model -
    T = model.T;
    price = model.price; # the price is set, we are looking for the interest rate
    Vₚ = model.par

    # we are passing in the rate -
    rate = r;

    # compute the discount with this rate -
    𝒟 = exp(rate*T);
    future_payout = (1/𝒟)*Vₚ
   
    # compute the npv value -
    npv_value = (future_payout - price)

    # return -
    return npv_value
end


function _fitness(κ, model::MyUSTreasuryZeroCouponBondModel, compounding::T) where T <: AbstractCompoundingModel

    # grab the discount rate from the κ array -
	discount_rate = κ[1]
	
	# we need to min the NPV - 
	npv_value = _net_present_value(discount_rate, model, compounding)

	# return the fitness -
	return (npv_value)^2
end

"""
    YTM(model::MyUSTreasuryZeroCouponBondModel, compounding::T; rₒ::Float64 = 0.01) where T <: AbstractCompoundingModel
"""
function YTM(model::MyUSTreasuryZeroCouponBondModel, 
    compounding::T; rₒ::Float64 = 0.01) where T <: AbstractCompoundingModel

    # initialize -    
    xinitial = [rₒ]
	
	# setup bounds -
	L = 0.00001
	U = 0.99999
	
	# setup the objective function -
	OF(p) = _fitness(p, model, compounding)
    
    # call the optimizer -
    opt_result = optimize(OF, L, U, xinitial, Fminbox(BFGS()))

    # grab the solution -
    bgfs_soln = Optim.minimizer(opt_result)[1]

    # return -
    return bgfs_soln
end
# using EmpiricalCDFs


# struct Copula <: MutualInformationEstimator

# end

# function mutualinfo(e::Shannon, est::Copula, x, y)
#     joint_vars = columns(Dataset(x, y))
#     cdfs = [EmpiricalCDF() for var in joint_vars]
#     for (i, var) in enumerate(joint_vars)
#         append!(cdfs[i], var)
#         sort!(cdfs[i])
#     end

#     copula_space = Dataset(cdfs)

# end

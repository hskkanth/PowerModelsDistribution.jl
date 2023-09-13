
"Solve multinetwork load shedding problem with storage and flexible loads"
function solve_mn_mc_mld_simple_water(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return solve_mc_model(data, model_type, solver, build_mn_mc_mld_simple_water; multinetwork=true, kwargs...)
end

"Multinetwork load shedding problem including storage"
function build_mn_mc_mld_simple_water(pm::AbstractUnbalancedPowerModel)
    for (n, network) in nws(pm)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_generator_power_ne(pm; nw=n)
        variable_mc_bus_voltage(pm; nw=n)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)
        variable_mc_storage_power_mi(pm; nw=n, relax=true)
        variable_mc_storage_power_mi_on_off_ne(pm; nw=n)
        variable_mc_generator_indicator_ne(pm; nw=n, relax=false)
    end

    for (n, network) in nws(pm)
        constraint_mc_model_voltage(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for i in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, i; nw=n)
        end

        for i in ids(pm, n, :gen_ne)
            constraint_mc_generator_power_ne(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance_shed(pm, i; nw=n)
        end

        # for i in ids(pm, n, :load)
        #     # Constrain load variables if they are not connected to a pump.
        #     constraint_fixed_load(pm, i; nw = n)
        # end

        for i in ids(pm, n, :storage)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
            constraint_storage_complementarity_mi(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage_ne)
            constraint_mc_storage_losses_ne(pm, i; nw=n)
            constraint_mc_storage_thermal_limit_ne(pm, i; nw=n)
            constraint_storage_complementarity_mi_ne(pm, i; nw=n)
            constraint_storage_indicator_expand_ne(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_ohms_yt_from(pm, i; nw=n)
            constraint_mc_ohms_yt_to(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
            constraint_mc_ampacity_from(pm, i; nw=n)
            constraint_mc_ampacity_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
            constraint_mc_switch_ampacity(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))
    n_0 = network_ids[1]
    n_1 = network_ids[1]

    for i in ids(pm, :storage; nw=n_1)
        constraint_storage_state(pm, i; nw=n_1)
    end

    for i in ids(pm, :storage_ne; nw=n_1)
        constraint_storage_state_ne(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        for i in ids(pm, :storage_ne; nw=n_2)
            constraint_storage_state_ne(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    for i in ids(pm, :storage, nw=n_0)
        constraint_round_trip_storage(pm, i, n_0, n_1)
    end

    for i in ids(pm, :storage_ne, nw=n_0)
        constraint_round_trip_storage_ne(pm, i, n_0, n_1)
    end

    # objective_mc_min_load_setpoint_delta_simple(pm)
    objective_ne(pm)
end


"Multinetwork load shedding problem for Branch Flow model"
function build_mn_mc_mld_simple_water(pm::AbstractUBFModels)
    for (n, network) in nws(pm)
        variable_mc_branch_current(pm; nw=n)
        variable_mc_branch_power(pm; nw=n)
        variable_mc_switch_power(pm; nw=n)
        variable_mc_transformer_power(pm; nw=n)
        variable_mc_generator_power(pm; nw=n)
        variable_mc_bus_voltage(pm; nw=n)

        variable_mc_load_indicator(pm; nw=n, relax=true)
        variable_mc_shunt_indicator(pm; nw=n, relax=true)
        variable_mc_storage_power_mi(pm; nw=n, relax=true)
        variable_mc_storage_power_mi_on_off_ne(pm; nw=n)
    end

    for (n, network) in nws(pm)
        constraint_mc_model_current(pm; nw=n)

        for i in ids(pm, n, :ref_buses)
            constraint_mc_theta_ref(pm, i; nw=n)
        end

        for i in ids(pm, n, :gen)
            constraint_mc_generator_power(pm, i; nw=n)
        end

        for i in ids(pm, n, :bus)
            constraint_mc_power_balance_shed(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage)
            constraint_mc_storage_losses(pm, i; nw=n)
            constraint_mc_storage_thermal_limit(pm, i; nw=n)
            constraint_storage_complementarity_mi(pm, i; nw=n)
        end

        for i in ids(pm, n, :storage_ne)
            constraint_mc_storage_losses_ne(pm, i; nw=n)
            constraint_mc_storage_thermal_limit_ne(pm, i; nw=n)
            constraint_storage_complementarity_mi_ne(pm, i; nw=n)
            constraint_storage_indicator_expand_ne(pm, i; nw=n)
        end

        for i in ids(pm, n, :branch)
            constraint_mc_power_losses(pm, i; nw=n)
            constraint_mc_model_voltage_magnitude_difference(pm, i; nw=n)
            constraint_mc_voltage_angle_difference(pm, i; nw=n)
            constraint_mc_thermal_limit_from(pm, i; nw=n)
            constraint_mc_thermal_limit_to(pm, i; nw=n)
            constraint_mc_ampacity_from(pm, i; nw=n)
            constraint_mc_ampacity_to(pm, i; nw=n)
        end

        for i in ids(pm, n, :switch)
            constraint_mc_switch_state(pm, i; nw=n)
            constraint_mc_switch_thermal_limit(pm, i; nw=n)
            constraint_mc_switch_ampacity(pm, i; nw=n)
        end

        for i in ids(pm, n, :transformer)
            constraint_mc_transformer_power(pm, i; nw=n)
        end
    end

    network_ids = sort(collect(nw_ids(pm)))

    n_1 = network_ids[1]

    for i in ids(pm, :storage; nw=n_1)
        constraint_storage_state(pm, i; nw=n_1)
    end

    for i in ids(pm, :storage_ne; nw=n_1)
        constraint_storage_state_ne(pm, i; nw=n_1)
    end

    for n_2 in network_ids[2:end]
        for i in ids(pm, :storage; nw=n_2)
            constraint_storage_state(pm, i, n_1, n_2)
        end

        for i in ids(pm, :storage_ne; nw=n_2)
            constraint_storage_state_ne(pm, i, n_1, n_2)
        end

        n_1 = n_2
    end

    objective_mc_min_load_setpoint_delta_simple(pm)
end

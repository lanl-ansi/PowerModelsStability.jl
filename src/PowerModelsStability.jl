module PowerModelsStability

    import PowerModelsDistribution
    import PowerModelsDistribution: ids, ref, var, sol, con, nw_id_default

    const _PMD = PowerModelsDistribution
    const  PMS = PowerModelsStability

    import JuMP
    import JSON

    import LinearAlgebra: eigvals, eigvecs, dot, inv, norm
    const LA = LinearAlgebra

    # Create our module level logger (this will get precompiled)
    const _LOGGER = Memento.getlogger(@__MODULE__)

    # Register the module level logger at runtime so that folks can access the logger via `getlogger(PowerModelsStability)`
    # NOTE: If this line is not included then the precompiled `PowerModelsStability._LOGGER` won't be registered at runtime.
    __init__() = Memento.register(_LOGGER)

    "Suppresses information and warning messages output by PowerModelsStability, for fine grained control use of the Memento package"
    function silence()
        Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
        Memento.setlevel!(Memento.getlogger(PowerModelsStability), "error")
    end

    "allows the user to set the logging level without the need to add Memento"
    function logger_config!(level)
        Memento.config!(Memento.getlogger("PowerModelsStability"), level)
    end

    include("core/inverter.jl")
    include("core/constraint.jl")
    include("core/data.jl")

    include("data_model/eng2math.jl")

    include("io/preprocessing.jl")
    include("io/json.jl")
    include("io/common.jl")

    include("prob/common.jl")

    include("core/export.jl")
end

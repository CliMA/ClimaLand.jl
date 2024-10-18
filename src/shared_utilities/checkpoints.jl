function save_checkpoint(u, t, output_path, model)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    output_file = joinpath(output_dir, "day$day.$sec.hdf5")
    comms_ctx = ClimaComms.context(model)
    hdfwriter = InputOutput.HDF5Writer(output_file, comms_ctx)
    InputOutput.write_attributes!(hdfwriter,
                                  "/"
                                  Dict(
                                      "time" => t,
                                      "land_model_hash" => hash(model)
                                  )
                                  )
    InputOutput.write!(hdfwriter, u, "u")
    Base.close(hdfwriter)
    return nothing
end

function read_checkpoint(file_path, model)
    hdfreader = InputOutput.HDF5Reader(output_file, comms_ctx)
    Y = InputOutput.read_field(restart_reader, "u")
    attributes = InputOutput.read_attributes(hdfreader, "/")
    if hash(model) != attributes["land_model_hash"]
        @warn "Restart file $(file_path) was constructed with a different land model"
    end
    t = attributes["time"]
    Base.close(hdfreader)
    return Y, t
end

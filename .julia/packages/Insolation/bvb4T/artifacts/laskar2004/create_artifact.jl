using Downloads, Pkg.Artifacts
using Printf

insoln_file = Downloads.download(
    "http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLN.LA2004.BTL.ASC",
)
insolp_file = Downloads.download(
    "http://vo.imcce.fr/insola/earth/online/earth/La2004/INSOLP.LA2004.BTL.ASC",
)

artifact_tree_sha1 = create_artifact() do dir

    insoln_lines = reverse(readlines(insoln_file)) # stored in reverse order
    insolp_lines = readlines(insolp_file)[2:end]   # 0 is duplicated in both files 

    open(joinpath(dir, "INSOL.LA2004.BTL.csv"); write = true) do f
        println(f, "# t (kyr from J2000), ecc, obliq (rad), varpi (rad)")
        varpi_offset = Float64(pi)
        varpi_rad_prev = 1.0 + varpi_offset
        for line in vcat(insoln_lines, insolp_lines)
            kyr_s, ecc_s, obliq_rad_s, varpi_rad_s = split(line)
            kyr = parse(Float64, kyr_s)
            ecc = parse(Float64, replace(ecc_s, 'D' => 'e'))
            obliq_rad = parse(Float64, replace(obliq_rad_s, 'D' => 'e'))
            varpi_rad = parse(Float64, replace(varpi_rad_s, 'D' => 'e'))

            # keep continuous so it can be interpolated
            if varpi_rad + varpi_offset < varpi_rad_prev - pi
                varpi_offset += 2pi
            elseif varpi_rad + varpi_offset > varpi_rad_prev + pi
                varpi_offset -= 2pi
            end
            varpi_rad += varpi_offset
            varpi_rad_prev = varpi_rad
            @printf(
                f,
                "%.3f,%.16e,%.16e,%.16e\n",
                kyr,
                ecc,
                obliq_rad,
                varpi_rad
            )
        end
    end

    # add README
    cp(joinpath(@__DIR__, "README.md"), joinpath(dir, "README.md"))
end

artifact_path = joinpath(@__DIR__, "laskar2004.tar.gz")

artifact_archive_sha256 = archive_artifact(artifact_tree_sha1, artifact_path)

@info "Created artifact" artifact_path artifact_tree_sha1 artifact_archive_sha256

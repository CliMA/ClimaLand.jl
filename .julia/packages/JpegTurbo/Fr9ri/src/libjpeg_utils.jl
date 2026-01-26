"""
    versioninfo(io=stdout)

Print information about the libjpeg-turbo in use.
"""
function versioninfo(io=stdout)
    println(io, "JpegTurbo.jl version: ", project_info["version"])
    println(io, "libjpeg version: ", LibJpeg.JPEG_LIB_VERSION)
    println(io, "libjpeg-turbo version: ", LibJpeg.LIBJPEG_TURBO_VERSION)
    println(io, "bit mode: ", LibJpeg.BITS_IN_JSAMPLE)
    println(io, "SIMD: ", LibJpeg.WITH_SIMD == 1 ? "enabled" : "disabled")
end

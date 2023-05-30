# Load packages
using ParamViz
using JSServe
using ClimaLSM
using ClimaLSM.Canopy
FT = Float64

# Create server
IP = "127.0.0.1"
port = 9385
server = Server(IP, port; proxy_url="https://clima.westus3.cloudapp.azure.com/jsserve/")

# Load and create apps
include("leaf_An.jl"); An_app = An_app_f()
An_app = An_app_f() # need to run twice for unicode character... (bug)

# Route apps
route!(server,"/leaf_An"=>An_app)


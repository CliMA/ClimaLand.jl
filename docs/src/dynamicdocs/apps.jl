# Load packages
using ParamViz
using JSServe
using ClimaLSM
using ClimaLSM.Canopy
using ClimaLSM.Soil.Biogeochemistry
FT = Float64

# Create server
IP = "127.0.0.1"
port = 9385
server = Server(IP, port; proxy_url="https://clima.westus3.cloudapp.azure.com/jsserve/")

# Load and create apps
include("leaf_An.jl"); An_app = An_app_f()
An_app = An_app_f() # need to run twice for unicode character... (bug)
include("Beer.jl"); beer_app = Beer_app_f()
beer_app = Beer_app_f()
include("hetero_resp.jl") # DAMM Rh
Rh_app = Rh_app_f(); Rh_app = Rh_app_f(); 
include("leaf_An_ci.jl") # An_ci
An_ci_app = An_ci_app_f(); An_ci_app = An_ci_app_f(); 

# Route apps
route!(server,"/leaf_An"=>An_app)
route!(server,"/beer_APAR"=>beer_app)
route!(server,"/leaf_An_ci"=>An_ci_app)
route!(server,"/Rh"=>Rh_app)


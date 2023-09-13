# Load packages
using ParamViz
using JSServe
using ClimaLSM
using ClimaLSM.Canopy
using ClimaLSM.Soil.Biogeochemistry
FT = Float64

# Create server
IPa = "127.0.0.1"
port = 9385
server = Server(IPa, port; proxy_url="https://clima.westus3.cloudapp.azure.com/jsserve/")

# Load and create apps
include("apps/leaf_An.jl"); An_app = An_app_f(); An_app = An_app_f(); # need to run twice for unicode character... (bug)
include("apps/Beer.jl"); beer_app = Beer_app_f(); beer_app = Beer_app_f();
include("apps/hetero_resp.jl"); Rh_app = Rh_app_f(); Rh_app = Rh_app_f(); 
include("apps/leaf_An_ci.jl"); An_ci_app = An_ci_app_f(); An_ci_app = An_ci_app_f(); 

# Route apps
route!(server,"/leaf_An"=>An_app)
route!(server,"/beer_APAR"=>beer_app)
route!(server,"/leaf_An_ci"=>An_ci_app)
route!(server,"/Rh"=>Rh_app)

wait()

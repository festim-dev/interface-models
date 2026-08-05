import festim as F
import numpy as np

H = F.Species("H")
D = F.Species("D")

H2 = F.GasSpecies("H2")
HD = F.GasSpecies("HD")
D2 = F.GasSpecies("D2")


model = F.HydrogenTransportProblemDiscontinuous()

model.mesh = F.Mesh1D(vertices=np.linspace(0, 1, 100))

mat = F.Material(D_0=1, E_D=0)
vol = F.VolumeSubdomain1D(id=1, material=mat, borders=[0, 1])
right = F.SurfaceSubdomain1D(id=2, x=1)

model.subdomains = [vol, right]


enclosure = F.Enclosure(
    volume=1e-3, species=[H2, HD, D2], temperature=1000, surfaces={right: 1}
)


model.species = [H, D]

for species in [H, D]:
    species.subdomains = model.volume_subdomains

model.enclosures = [enclosure]


k_plus, k_minus = 10, 20

hh_h2 = F.SurfaceReactionBC(
    reactant=[H, H],
    gas_pressure=H2,
    k_r0=k_plus,
    E_kr=0,
    k_d0=k_minus,
    E_kd=0,
    subdomain=right,
)
hd_hd = F.SurfaceReactionBC(
    reactant=[H, D],
    gas_pressure=HD,
    k_r0=2 * k_plus,
    E_kr=0,
    k_d0=2 * k_minus,
    E_kd=0,
    subdomain=right,
)
dd_d2 = F.SurfaceReactionBC(
    reactant=[D, D],
    gas_pressure=D2,
    k_r0=k_plus,
    E_kr=0,
    k_d0=k_minus,
    E_kd=0,
    subdomain=right,
)

model.boundary_conditions = [hh_h2, hd_hd, dd_d2]

model.temperature = 300

model.initial_conditions = [
    F.InitialConcentration(species=H, value=1e10, volume=vol),
    F.InitialConcentration(species=D, value=1e10, volume=vol),
]

model.settings = F.Settings(atol=1e-10, rtol=1e-10, final_time=10, stepsize=0.1)


model.exports = [F.GasPressure(field=gp) for gp in [H2, HD, D2]]


from dolfinx.log import set_log_level, LogLevel

set_log_level(LogLevel.INFO)

model.initialise()


model.run()
print("H", vol.u.x.array[:])

print("H2 pressure:", model.exports[0].data[-1])
print("HD pressure:", model.exports[1].data[-1])
print("D2 pressure:", model.exports[2].data[-1])
import matplotlib.pyplot as plt

# plt.plot(model.exports[0].t, model.exports[0].data, label="H2")
# plt.plot(model.exports[1].t, model.exports[1].data, label="HD")
# plt.plot(model.exports[2].t, model.exports[2].data, label="D2")

# stack area plot instead
plt.stackplot(
    model.exports[0].t,
    model.exports[0].data,
    model.exports[1].data,
    model.exports[2].data,
    labels=["H2", "HD", "D2"],
)

plt.xlabel("Time (s)")
plt.ylabel("Gas pressure (Pa)")
plt.legend()
plt.show()

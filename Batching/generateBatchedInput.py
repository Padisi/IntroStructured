import os
import json
import numpy as np
import pyUAMMD

with open("./parameters.json") as f:
    inputData = json.load(f)

def oneSimulation(inputData):
    particles = np.loadtxt("ftsz.org")
    N = particles.shape[0]
    bonds = np.loadtxt("ftsz_bonds.dat")

    sigma  = inputData["sigma"]
    mass   = inputData["mass"]

    k = inputData["k"]
    L = inputData["L"]
    T  = inputData["T"]
    dt = inputData["dt"]
    nSteps = inputData["nSteps"]
    nWrite = inputData["nWrite"]


    sim=pyUAMMD.simulation()

    sim["system"]={}
    sim["system"]["info"] = {}
    sim["system"]["info"]["type"] = ["Simulation","Information"]
    sim["system"]["info"]["parameters"] = {"name":"ftsz_molecule"}


    sim["global"]={}

    sim["global"]["fundamental"] = {"type":["Fundamental","Time"]}

    sim["global"]["ensemble"]={}
    sim["global"]["ensemble"]["labels"] = ["box","temperature"]
    sim["global"]["ensemble"]["data"] = [[[L,L,L],T]]
    sim["global"]["ensemble"]["type"] = ["Ensemble","NVT"]

    sim["global"]["types"] = {}
    sim["global"]["types"]["labels"] =  ["name", "mass", "radius", "charge"]
    sim["global"]["types"]["data"] = [["A",mass,sigma/2,0.0]]
    sim["global"]["types"]["type"] = ["Types","Basic"]

    sim["global"]["units"] = {}
    sim["global"]["units"]["type"] = ["Units","None"]

    sim["integrator"]={}

    sim["integrator"]["schedule"]={}
    sim["integrator"]["schedule"]["type"]=["Schedule","Integrator"]
    sim["integrator"]["schedule"]["data"]=[[1,"brownian",nSteps]]
    sim["integrator"]["schedule"]["labels"]=["order","integrator","steps"]

    sim["integrator"]["brownian"]={}
    sim["integrator"]["brownian"]["parameters"]={"timeStep":dt,"viscosity":inputData["viscosity"]}
    sim["integrator"]["brownian"]["type"]=["Brownian","EulerMaruyama"]

    sim["state"]={}
    sim["state"]["labels"]=["id","position"]
    sim["state"]["data"]=[]

    x0 = particles[0,0]
    y0 = particles[0,1]
    z0 = particles[0,2]

    for i in range(particles.shape[0]):
        sim["state"]["data"].append([i,[particles[i,0]-x0,particles[i,1]-y0,particles[i,2]-z0+sigma/2]])

    sim["topology"]={}

    sim["topology"]["structure"]={}
    sim["topology"]["structure"]["labels"]=["id","type"]
    sim["topology"]["structure"]["data"]=[]
    for i in range(particles.shape[0]):
        sim["topology"]["structure"]["data"].append([i,"A"])

    sim["topology"]["forceField"]={}

    eps = 1

    sim["topology"]["forceField"]["fixed_particles"]={}
    sim["topology"]["forceField"]["fixed_particles"]["type"]=["Bond1","FixedHarmonic"]
    sim["topology"]["forceField"]["fixed_particles"]["labels"]=["id_i","K","r0","position"]
    sim["topology"]["forceField"]["fixed_particles"]["data"]=[[0.0,k,sigma/2,[0.0,0.0,0.0]]]

    sim["topology"]["forceField"]["harmonicBonds"]={}
    sim["topology"]["forceField"]["harmonicBonds"]["type"]=["Bond2","Harmonic"]
    sim["topology"]["forceField"]["harmonicBonds"]["labels"]=["id_i","id_j","r0","K"]
    sim["topology"]["forceField"]["harmonicBonds"]["data"]=[]
    for i in range(bonds.shape[0]):
        sim["topology"]["forceField"]["harmonicBonds"]["data"].append([int(bonds[i,0])-1,int(bonds[i,1])-1,bonds[i,2],bonds[i,3]])

    sim["topology"]["forceField"]["floor"]={}
    sim["topology"]["forceField"]["floor"]["type"]=["Surface","SurfaceWCAType1"]
    sim["topology"]["forceField"]["floor"]["parameters"]={"surfacePosition":0.0}
    sim["topology"]["forceField"]["floor"]["labels"]=["name","epsilon","sigma"]
    sim["topology"]["forceField"]["floor"]["data"]=[["A",eps,sigma]]

    sim["topology"]["forceField"]["wca_particles"]={}
    sim["topology"]["forceField"]["wca_particles"]["type"]=["NonBonded","WCAType2"]
    sim["topology"]["forceField"]["wca_particles"]["parameters"]={"cutOffFactor":2.5,"condition":"all"}
    sim["topology"]["forceField"]["wca_particles"]["labels"]=["name_i","name_j","epsilon","sigma"]
    sim["topology"]["forceField"]["wca_particles"]["data"]=[["A","A",eps,sigma]]

    sim["topology"]["forceField"]["nl"]={}
    sim["topology"]["forceField"]["nl"]["type"]=["VerletConditionalListSet","all"]
    sim["topology"]["forceField"]["nl"]["parameters"]={"cutOffVerletFactor":1.2}

    sim["simulationStep"]={}

    sim["simulationStep"]["saveState"]={}
    sim["simulationStep"]["saveState"]["type"]=["WriteStep","WriteStep"]
    sim["simulationStep"]["saveState"]["parameters"]={
        "intervalStep":nWrite,
        "outputFilePath":"output",
        "outputFormat":"sp",
    }

    return sim

#try:
#    os.makedirs("results")
#except:
#    pass

#sim.write("results/simulation.json")

total_sim = pyUAMMD.simulation()
total_sim = oneSimulation(inputData)
for n in range(1,inputData["nBatch"]):
    print(n)
    sim_i = oneSimulation(inputData)
    total_sim.append(sim_i)

os.makedirs("results", exist_ok=True)
total_sim.write("results/simulation.json")

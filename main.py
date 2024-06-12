import math
import matplotlib.pyplot as plt
import cmath
import numpy as np
import time

def readcdf2(cdf_file_path):
    bus_data = []
    branch_data = []
    bus_numbersx = set()
    V = np.ones((300, 1))  # flat start


    with open(cdf_file_path, 'r') as file:
        lines = file.readlines()

        reading_bus_data = False
        reading_branch_data = False

        for line in lines:
            if line.startswith('BUS DATA FOLLOWS'):
                reading_bus_data = True
                continue
            elif line.startswith('BRANCH DATA FOLLOWS'):
                reading_bus_data = False
                reading_branch_data = True
                continue

            if reading_bus_data:

                fields = line.split()
                if len(fields) >= 3:
                    try:
                        bus_number = int(fields[0])
                        bus_base = float(fields[11])
                        bus_type = int(fields[4])
                        P_load = float(fields[7])
                        Q_load = float(fields[8])
                        P_gen = float(fields[9])
                        Q_gen = float(fields[10])
                        Shunt_conductance_G = float(fields[15])
                        Shunt_susceptance_B = float(fields[16])
                        bus_data.append(
                            [bus_number, bus_base, bus_type, P_load, Q_load, P_gen, Q_gen, Shunt_conductance_G,
                             Shunt_susceptance_B])
                        bus_numbersx.add(bus_number)
                    except (ValueError, IndexError):
                        continue

            elif reading_branch_data:

                fields = line.split()
                if len(fields) >= 18:
                    try:
                        from_bus = int(fields[0])
                        to_bus = int(fields[1])
                        branch_type = int(fields[5])
                        branch_resistance = float(fields[6])
                        branch_reactance = float(fields[7])
                        line_charging = float(fields[8])
                        tr_turnsratio = float(fields[14])
                        tr_angle = float(fields[15])

                        if len(fields[16]) > 10:

                            Min_tapphase = 0.9043
                            Max_tapphase = 1.10435


                        else:
                            Min_tapphase = float(fields[16])
                            Max_tapphase = float(fields[17])
                        branch_data.append(
                            [from_bus, to_bus, branch_type, branch_resistance, branch_reactance, line_charging,
                             tr_turnsratio, tr_angle, Min_tapphase, Max_tapphase])
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing branch data: {e}. Skipping line: {line}")
                        continue

    return bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load, bus_type


cdf_file_path = 'ieee300cdf.txt'

bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load, bus_type = readcdf2(cdf_file_path)


def build_ybus(branch_data, bus_numbersx):
    bus_numbers = sorted(bus_numbersx)

    num_buses = len(bus_numbers)

    bus_index_map = {bus: idx for idx, bus in enumerate(bus_numbers)}
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    print("this is bus index map: ", bus_index_map)

    for branch in branch_data:
        from_bus = branch[0]
        to_bus = branch[1]
        R = branch[3]
        X = branch[4]
        charge = branch[5]
        tapratio = branch[6]

        if tapratio == 0:
            tapratio = 1

        phaseangle = branch[7]

        admittance = 1 / (R + (1j * X))
        Ycharg = 1j * charge / 2
        if phaseangle == 0:

            Ybus[bus_index_map[from_bus], bus_index_map[to_bus]] -= admittance / tapratio
            Ybus[bus_index_map[to_bus], bus_index_map[from_bus]] -= admittance / tapratio
            Ybus[bus_index_map[from_bus], bus_index_map[from_bus]] += (admittance / (tapratio) ** 2) + Ycharg
            Ybus[bus_index_map[to_bus], bus_index_map[to_bus]] += admittance + Ycharg

        else:
            newphase = cmath.rect(1, (phaseangle / math.pi))
            Ybus[bus_index_map[from_bus], bus_index_map[to_bus]] -= admittance / np.conj(newphase)
            Ybus[bus_index_map[to_bus], bus_index_map[from_bus]] -= admittance / newphase
            Ybus[bus_index_map[from_bus], bus_index_map[from_bus]] += admittance / (newphase) ** 2
            Ybus[bus_index_map[to_bus], bus_index_map[to_bus]] += admittance

    for bus in bus_data:
        bus_number = bus[0]

        diag_admittance = bus[7] + 1j * bus[8]  # y = g +jb

        Ybus[bus_index_map[bus_number], bus_index_map[bus_number]] += diag_admittance

    return Ybus

start_time = time.time()

cdf_file_path = 'ieee300cdf.txt'

bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load, bus_type = readcdf2(cdf_file_path)

Ybus = build_ybus(branch_data, bus_numbersx)

end_time = time.time()
# Calculate execution time
execution_time = end_time - start_time
print(f"Execution time of Ybus formation: {execution_time} seconds")


# Plotting sparsity pattern
plt.figure(figsize=(10, 8))
plt.spy(Ybus, markersize=1)
plt.title('Sparsity Pattern of Ybus Matrix')
plt.xlabel('Bus Index')
plt.ylabel('Bus Index')
plt.show()

import networkx as nx

# Create a graph from Ybus matrix (assuming it represents a network)
G = nx.Graph()
edges = np.argwhere(np.abs(Ybus) > 0)  # Get nonzero elements as edges

for edge in edges:
    if edge[0] != edge[1]:  # Avoid self-loops
        G.add_edge(edge[0], edge[1])

# Visualize the graph
plt.figure(figsize=(8, 5))
pos = nx.spring_layout(G)  # Layout nodes using spring algorithm
nx.draw(G, pos, with_labels=True, node_size=200, node_color='skyblue', font_size=8)
plt.title('Power System Network Visualization')
plt.show()

plt.figure(figsize=(8, 5))  # BEST LAYOUT
pos = nx.kamada_kawai_layout(G)
nx.draw(G, pos, with_labels=True, node_size=200, node_color='skyblue', font_size=8)
plt.title('Kamada-Kawai Layout')
plt.show()

P = P_gen-P_load
Q = Q_gen-Q_load
theta = np.zeros((300, 1))


def power_flow(Ybus, P, Q, V, theta, bus_type, tol=1e-3, max_iter=20):
    num_buses = len(Ybus)
    PQ_buses = np.where(bus_type == 0 or bus_type == 1)[0]
    PV_buses = np.where(bus_type == 2)[0]

    for iteration in range(max_iter):
        P_calc = np.zeros(num_buses)
        Q_calc = np.zeros(num_buses)

        for i in range(num_buses):
            for j in range(num_buses):
                P_calc[i] += V[i] * V[j] * (Ybus[i, j].real * np.cos(theta[i] - theta[j]) + Ybus[i, j].imag * np.sin(
                    theta[i] - theta[j]))
                Q_calc[i] += V[i] * V[j] * (Ybus[i, j].real * np.sin(theta[i] - theta[j]) - Ybus[i, j].imag * np.cos(
                    theta[i] - theta[j]))

        dP = P - P_calc
        dQ = Q - Q_calc

        mismatch = np.concatenate((dP[PQ_buses], dP[PV_buses], dQ[PQ_buses]))

        if np.max(np.abs(mismatch)) < tol:
            print(f'Converged in {iteration + 1} iterations.')
            break

        J = np.zeros((len(PQ_buses) + len(PV_buses) + len(PQ_buses), len(PQ_buses) + len(PV_buses) + len(PQ_buses)))

        for i, k in enumerate(np.concatenate((PQ_buses, PV_buses))):
            for j, m in enumerate(np.concatenate((PQ_buses, PV_buses))):
                if i == j:
                    J[i, j] = -Q_calc[k] - (V[k] ** 2) * Ybus[k, k].imag
                    J[i, j + len(PQ_buses) + len(PV_buses)] = P_calc[k] + (V[k] ** 2) * Ybus[k, k].real
                else:
                    J[i, j] = V[k] * V[m] * (Ybus[k, m].real * np.sin(theta[k] - theta[m]) - Ybus[k, m].imag * np.cos(
                        theta[k] - theta[m]))
                    J[i, j + len(PQ_buses) + len(PV_buses)] = V[k] * (
                                Ybus[k, m].real * np.cos(theta[k] - theta[m]) + Ybus[k, m].imag * np.sin(
                            theta[k] - theta[m]))

        for i, k in enumerate(PQ_buses):
            for j, m in enumerate(np.concatenate((PQ_buses, PV_buses))):
                if i == j:
                    J[i + len(PQ_buses) + len(PV_buses), j] = P_calc[k] + (V[k] ** 2) * Ybus[k, k].real
                    J[i + len(PQ_buses) + len(PV_buses), j + len(PQ_buses) + len(PV_buses)] = Q_calc[k] - (V[k] ** 2) * \
                                                                                              Ybus[k, k].imag
                else:
                    J[i + len(PQ_buses) + len(PV_buses), j] = V[k] * V[m] * (
                                Ybus[k, m].real * np.cos(theta[k] - theta[m]) + Ybus[k, m].imag * np.sin(
                            theta[k] - theta[m]))
                    J[i + len(PQ_buses) + len(PV_buses), j + len(PQ_buses) + len(PV_buses)] = -V[k] * (
                                Ybus[k, m].real * np.sin(theta[k] - theta[m]) - Ybus[k, m].imag * np.cos(
                            theta[k] - theta[m]))

        dx = np.linalg.solve(J, mismatch)

        theta[np.concatenate((PQ_buses, PV_buses))] += dx[:len(PQ_buses) + len(PV_buses)]
        V[PQ_buses] += dx[len(PQ_buses) + len(PV_buses):]

    return V, theta

start_time2 = time.time()

V_final, theta_final = power_flow(Ybus, P, Q, V, theta, bus_type)

end_time2 = time.time()


execution_time2 = end_time2 - start_time2
print(f"Execution time of power flow: {execution_time2} seconds")



for i in range(len(V_final)):
    print(f"Bus {i + 1}: V = {float(V_final[i]):.4f} p.u., Angle = {float(theta_final[i]):.4f} degrees")

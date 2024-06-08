import math
import matplotlib.pyplot as plt
import cmath
import numpy as np
def readcdf2(cdf_file_path):
    bus_data = []
    branch_data = []
    bus_numbersx = set()
    V = np.ones((300, 1)) # flat start

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
                        bus_data.append( [bus_number,  bus_base, bus_type, P_load, Q_load,  P_gen,  Q_gen, Shunt_conductance_G, Shunt_susceptance_B])
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
                        branch_data.append([from_bus, to_bus, branch_type, branch_resistance, branch_reactance, line_charging, tr_turnsratio ,tr_angle, Min_tapphase, Max_tapphase])
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing branch data: {e}. Skipping line: {line}")
                        continue


    return bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load



cdf_file_path = 'ieee300cdf.txt'


bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load = readcdf2(cdf_file_path)


print("Bus Data:")
for bus in bus_data:
    print(bus)

print("\nBranch Data:")
for branch in branch_data:
    print(branch)


print(bus_data[0][:])
x= (bus_data[4][1])*(bus_data[4][3])
print(x)
print(bus_numbersx)
print(V)


def build_ybus(branch_data, bus_numbersx):
    bus_numbers = sorted(bus_numbersx)
    print(bus_numbers)

    num_buses = len(bus_numbers)
    print(num_buses)
    bus_index_map = {bus: idx for idx, bus in enumerate(bus_numbers)}
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    print("this is bus index map: ",bus_index_map)

    for branch in branch_data:
        from_bus = branch[0]
        to_bus = branch[1]
        R = branch[3]
        X = branch[4]
        charge = branch[5]
        tapratio = branch[6]

        if tapratio == 0:
            tapratio = 1

        print(tapratio)
        phaseangle = branch[7]

        admittance = 1 / (R + (1j * X))
        Ycharg = 1j * charge/2
        if phaseangle == 0:

            Ybus[bus_index_map[from_bus], bus_index_map[to_bus]] -= admittance /tapratio
            Ybus[bus_index_map[to_bus], bus_index_map[from_bus]] -= admittance /tapratio
            Ybus[bus_index_map[from_bus], bus_index_map[from_bus]] += (admittance / (tapratio)**2) + Ycharg
            Ybus[bus_index_map[to_bus], bus_index_map[to_bus]] += admittance +Ycharg

        else:
            newphase = cmath.rect(1, (phaseangle / math.pi))
            Ybus[bus_index_map[from_bus], bus_index_map[to_bus]] -= admittance / np.conj(newphase)
            Ybus[bus_index_map[to_bus], bus_index_map[from_bus]] -= admittance / newphase
            Ybus[bus_index_map[from_bus], bus_index_map[from_bus]] += admittance / (newphase)**2
            Ybus[bus_index_map[to_bus], bus_index_map[to_bus]] += admittance


    for bus in bus_data:
        bus_number = bus[0]

        diag_admittance = bus[7] + 1j*bus[8] # y = g +jb

        Ybus[bus_index_map[bus_number], bus_index_map[bus_number]] += diag_admittance




    return Ybus


cdf_file_path = 'ieee300cdf.txt'


bus_data, branch_data, bus_numbersx, V, P_gen, Q_gen, P_load, Q_load = readcdf2(cdf_file_path)


Ybus = build_ybus(branch_data, bus_numbersx)


print("Ybus matrix:")
#for yb in Ybus:
#   print(yb)
print(len(Ybus))
print(Ybus[0][0],Ybus[0][1])
print(Ybus[1][0],Ybus[1][1])
print(Ybus[111][136],Ybus[111][111])
print(Ybus[2][1]) # normalde bağlı değilken neden 20j gözüküyor???


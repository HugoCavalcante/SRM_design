# -*- coding: utf-8 -*-
"""
Conversor de dados para o formato RASP.ENG (Força de empuxo de motor para foguetes).
Este é um dos formatos usados pelo OpenRocket para modelar motores de foguete e fazer simulações de vôo.

Formato: as primeiras linhas são comentários começados com ponto-e-vírgula.
A primeira linha não comentada contém o cabeçalho de identificação, com os dados:
1 - nome do motor, 2 - diâmetro (mm), 3 - comprimento (mm), 4 - atrasos das cargas de ejeção, separados por "-" 
(a documentação diz que pode ser usado um caracter "P", para motor sem estágios, mas "0" funciona com certeza),
5 - peso do propelente (kg), 6 - peso total (kg), 7 - Cód. do fabricante (NAR) 
Finalmente, duas colunas (separadas por espaço) com os dados: tempo (s) e força (N), terminando com força = 0.0.
Exemplo:
; Esta é uma linha de comentário
; Esta também
Rj-0 25 100 0 0.023 0.045 UFPB
0.020	2.485
0.057	12.695
0.089	25.660
0.100   0.0

Essa foi a descrição de um foguete Rj-0, com diâmetro externo 25 mm, comprimento 100 mm, com uma única carga (atraso 0, "plugged"),
massa do propelente 23 g, massa total (carga + propelentes) 45 g, fabricado por "UFPB". Versões antigas de softwarede simulação 
aceitam no máximo 32 pontos na série temporal, incluindo o ponto final. Falta testar se o OpenRocket também tem esta limitação para este formato.

Created on Wed Mar 20 10:37:18 2024

@author: hugo
"""


from numpy import loadtxt, savetxt, linspace, stack
import numpy as np

Plot_thrust = True

Data_Folder = "Motor_Thrust_Data"
#input_filename = "Motor_22_03_24_Ensaio_2.csv"
#input_filename = "Motor_Comerc_4_18-10-2024.csv"
input_filename = "Carcara_teoria_03-02-2025.csv"
output_filename = "Carcara_teoria_030225.eng"
input_filename = Data_Folder + "\\" + input_filename
output_filename = Data_Folder + "\\" + output_filename

Code = "Car-0"
diameter = 48
lenght = 333 
stage_delays = 0
prop_mass = 0.643
total_mass = 1.383
#Footer = f"{t_new[-1]+0.001} 0.0"
manufacturer_name = "CnE"
Header = "; Carcara - Modelo teorico com Planilha Nakka - 03-02-2025 \n"


### Comma to point converter (so you don't need to change the 'locale' to portuguese)
def comma2point(s):
    s = s.replace(b',', b'.')
    return s
conv = {0:comma2point, 1:comma2point}

#%% 

print(f"Reading input file: {input_filename}")

#Raw_data = loadtxt(input_filename, delimiter = ',') ## Ver a documentação de loadtxt para ajustar os parâmetros de leitura. Cogite-se também usar pandas framework ou o módulo csv.
#t, thrust = Raw_data[:,0], Raw_data[:,1]
#Raw_data = loadtxt(input_filename, delimiter = ',', usecols=0)
Raw_data = loadtxt(input_filename, delimiter = '\t', usecols=(0,1), skiprows=1, converters = conv)
t, thrust = Raw_data.copy().T

N_raw = len(t)
if N_raw < 2 or N_raw != len(thrust):
    print(f"Number of points too small ({N_raw}) or inconsistent. Aborting")
    exit(1)
print(f"Read {N_raw} points from file {input_filename}")
dt = (t[-1]-t[0])/(N_raw-1)


g = 9.78

# ### Podemos fazer um "downsampling" dos pontos medidos 
# from scipy.signal import resample_poly
# ### o número de pontos vai ser multiplicado pelo fator resample_factor = upsampling / downsampling
# upsampling = 1 #32
# downsampling = int(N_raw / 1000) #int(N_t / 10) # Note que N_t/N_t significa "inalterado", N_t/10 significa "1 ponto a cada 1000", N_t/1000 significa "1 ponto a cada 10", etc.
# resample_factor = upsampling / downsampling
# print(f"resampling factor = {upsampling}/{downsampling} = {resample_factor}") 

# #Data = Raw_data.copy()
# #Data[:,1] *= VtoN_factor
# #Data_new = resample_poly(Data, upsampling, downsampling, axis = 0)
# #t_new, thrust_new = Data_new[:,0], Data_new[:,1]
# #for i in range(len(t_new)-1):
# #    Data_new[i, 1] = max(0.001, Data_new[i, 1])
# thrust_new = resample_poly(thrust, upsampling, downsampling)


# ### terminating point must have thrust = 0, required by syntax
# thrust_new[-1] = 0.0
# N_new = len(thrust_new)
# dt_new = (t[-1]-t[0])/(N_new-1)
# t_new = linspace(t[0], t[-1], num = N_new)
# ### Devemos evitar o registro de forças negativas
# for i in range(len(t_new)-1):
#     thrust_new[i] = max(0.001, thrust_new[i])

# print(f"New number of points = {N_new}")


### Outras operações podem ser feitas nos dados, como produzir um deslocamento temporal
### Se o evento de inicial acontece no tempo t_offset (deveria ser 0)
#t_offset = 0.0 #0.01 #0.02 #-0.01
#t_new -=  t_offset

### Se quiser descartar tempos negativos
# discard_negative_times = False
# if discard_negative_times:
#     start_index = 0
#     for i in range(len(t_new)):
#         if t_new[i] < 0:
#             start_index = i+1
#     t_new = t_new[start_index:]
#     thrust_new = thrust_new[start_index:]
#     N_new = len(t_new)

# Data_new = stack((t_new, thrust_new)).T

### Desenhar os graficos para teste/comparação
#from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, grid, show
#from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, grid
if Plot_thrust:
    # from matplotlib.pyplot import subplots
    # fig, axs = subplots(2, 1, sharex = True)
    # axs[0].plot(t, Raw_data[:,1])  # measured signal, in volts
    # axs[0].plot(t, Raw_data)  # measured signal, in volts
    # axs[0].set_ylabel("signal (V)")
    # axs[0].grid()
    # axs[1].plot(t, thrust)  # partially processed signal (N)
    # axs[1].set_ylabel("thrust (N)")
    # axs[1].set_xlabel("time (s)")
    # #axs[1].plot(t_new, thrust_new, "o-", ms = 3, label = "processed")  # processed signal
    # #axs[1].plot(t_new+t_offset, thrust_new, "o-", ms = 3, label = "resampled")  # resampled signal without time shift
    # axs[1].grid()
    # axs[1].legend()
    from matplotlib.pyplot import plot, xlabel, ylabel, grid, show
    plot(t, thrust)
    xlabel('tempo (s)')
    ylabel('empuxo (N)')
    grid()
    show()


### Calculando alguns dados úteis:
#zero_thrust = 0.8 ### Valor de empuxo que pode ser considerado desprezível (ruído), em N
zero_thrust = 0.10 ### Valor de empuxo que pode ser considerado desprezível (ruído), em N
pulse_started = False
i = 1 
while (i<N_raw) and not(pulse_started):
    if (thrust[i] > zero_thrust) and (thrust[i-1] > zero_thrust):
        start_index_new = i
        start_time_new = t[i]
        pulse_started = True
    i += 1
pulse_ended = False
i = N_raw-3 
while (i>0) and not(pulse_ended):
    if (thrust[i] > zero_thrust) and (thrust[i+1] > zero_thrust):
        end_index_new = i
        end_time_new = t[i]
        pulse_ended = True
    i -= 1
if not(pulse_ended):
    end_index_new = N_raw-1
    end_time_new = t[-1]
pulse_duration = end_time_new - start_time_new
F_max = np.amax(thrust)
F_ave = np.mean(thrust[start_index_new:end_index_new])
total_impulse = pulse_duration*F_ave
Isp = total_impulse/(g*prop_mass)
print(f"Pulse duration = {pulse_duration:.3} s")
print(f"F_max = {F_max:.3} N")
print(f"F_ave = {F_ave:.3} N")
print(f"I_ave = {total_impulse:.3} N.s")
print(f"Isp = {Isp:.3} (1/s)")


Data_new = stack((t, thrust)).T

###

### Salvar arquivo de saída
#Header += f"{Code} {diameter} {lenght} {stage_delays} {prop_mass} {total_mass} {manufacturer_name}\n"
Header += f"{Code} {diameter} {lenght} {stage_delays} {prop_mass} {total_mass} {manufacturer_name}"
#savetxt(output_filename, Data_new, fmt='%.4g', delimiter = ' ', newline = "\n ", header = Header, footer = Footer, comments=";")
savetxt(output_filename, Data_new, fmt='%.4g', delimiter = ' ', newline = "\n ", header = Header, comments="")

print(f"ENG file saved as: {output_filename}")




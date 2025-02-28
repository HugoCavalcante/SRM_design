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

Plot_cal = True
Plot_thrust = True

Data_Folder = "Motor_Thrust_Data"
#input_filename = "Motor_22_03_24_Ensaio_2.csv"
#input_filename = "Motor_Comerc_4_18-10-2024.csv"
input_filename = Data_Folder + "\\" + input_filename
output_filename = Data_Folder + "\\" + output_filename

Code = "Rj-0"
diameter = 25
lenght = 100 
stage_delays = 0
prop_mass = 0.023
total_mass = 0.045
#Footer = f"{t_new[-1]+0.001} 0.0"
manufacturer_name = "CnE"
Header = "; Motor Comercial extraido de rojao treme-terra \n"


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
time, thrust = Raw_data.copy().T



# N_raw = len(t)
# if N_raw < 2 or N_raw != len(thrust):
#     print(f"Number of points too small ({N_raw}) or inconsistent. Aborting")
#     exit(1)
# print(f"Read {N_raw} points from file {input_filename}")
# dt = (t[-1]-t[0])/(N_raw-1)

N_raw = len(thrust)  # Número de pontos lidos do arquivo
N_t = 10240          # Número de pontos fixo na série temporal (característica do osciloscópio. Espera-se que N_t = N_raw)
time_base = 0.2      # "Base de tempo" usada no osciloscópio (tamanho da divisão da escala). T_total é aproximadamente 10*time_base (a tela tem 10 divisões)
T_total = 10.239*time_base # A taxa de amostragem não foi informada, não pode ser calculada dos dados. Estimando que dt seja um número redondo.
#T_total = 10.0*time_base # A taxa de amostragem não foi informada, não pode ser calculada dos dados. Estimando que dt seja um número "quebrado".
dt = T_total/(N_t-1)
t = linspace(0, T_total, num = N_t)

if N_raw != N_t:
    print(f"Number of points ({N_raw}) is inconsistent with expected length ({N_t}). Aborting")
    exit(1)

### Fazer aqui processamento adequado dos dados. Por exemplo, conversão/calibração do sinal medido em V para força em N.
### curva de calibração
g = 9.80665 # Fator de conversão de kgf para N
### Curva de calibração em abril de 2024
#F_cal = g*np.array([0.0, 0.542, 1.795, 4.047])  # massas de calibração, em kg
#V_cal = np.array([-6.91, -6.48, -5.38, -4.10])  # tensões correspondentes 
### Curva de calibração em 18/10/2024 
F_cal = g*np.array([0.0, 0.35, 0.548, 1.548, 2.006, 2.62, 4.62])  # massas de calibração, em kg
V_cal = np.array([0, 0.345, 0.525, 1.523, 1.996, 2.492, 4.323])  # tensões correspondentes 
V_off = -0.856004
V_cal = V_cal + V_off
F_max_calibraton = F_cal[-1] 
V_max_calibration = V_cal[-1]-V_off
VtoN_factor = F_max_calibraton / V_max_calibration
if Plot_cal:
    from matplotlib.pyplot import figure, plot, xlabel, ylabel, grid, legend 
    y_cal = linspace(0.0, g*10, num = 2)
    x_cal = V_off +(1/VtoN_factor)*y_cal
    fig_cal = figure("calibration curve")
    plot(V_cal, F_cal, 'o', label = "Calibration points")
    if V_off <0.0: 
        plot(x_cal, y_cal, label = f"x/{VtoN_factor:.3g} {V_off:.3g}")
    else:
        plot(x_cal, y_cal, label = f"x/{VtoN_factor:.3g} -{V_off:.3g}")
    xlabel("measured tension (V)")
    ylabel("calibration force (N)")
    grid()
    legend()


print(f"Voltage to force conversion factor = {VtoN_factor:.3} N/V")
print(f"Calibration offset = {V_off} V")

### correcting V_off 
#V_off = -6.38
thrust = VtoN_factor * (thrust-V_off)

print(f"Adjusted calibration offset = {V_off} V")

### Podemos fazer um "downsampling" dos pontos medidos 
from scipy.signal import resample_poly
### o número de pontos vai ser multiplicado pelo fator resample_factor = upsampling / downsampling
upsampling = 1 #32
downsampling = int(N_t / 1000) #int(N_t / 10) # Note que N_t/N_t significa "inalterado", N_t/10 significa "1 ponto a cada 1000", N_t/1000 significa "1 ponto a cada 10", etc.
resample_factor = upsampling / downsampling
print(f"resampling factor = {upsampling}/{downsampling} = {resample_factor}") 

#Data = Raw_data.copy()
#Data[:,1] *= VtoN_factor
#Data_new = resample_poly(Data, upsampling, downsampling, axis = 0)
#t_new, thrust_new = Data_new[:,0], Data_new[:,1]
#for i in range(len(t_new)-1):
#    Data_new[i, 1] = max(0.001, Data_new[i, 1])
thrust_new = resample_poly(thrust, upsampling, downsampling)


### terminating point must have thrust = 0, required by syntax
thrust_new[-1] = 0.0
N_new = len(thrust_new)
dt_new = (t[-1]-t[0])/(N_new-1)
t_new = linspace(t[0], t[-1], num = N_new)
### Devemos evitar o registro de forças negativas
for i in range(len(t_new)-1):
    thrust_new[i] = max(0.001, thrust_new[i])

print(f"New number of points = {N_new}")


### Outras operações podem ser feitas nos dados, como produzir um deslocamento temporal
### Se o evento de inicial acontece no tempo t_offset (deveria ser 0)
#t_offset = 0.0 #0.01 #0.02 #-0.01
#t_new -=  t_offset

### Se quiser descartar tempos negativos
discard_negative_times = False
if discard_negative_times:
    start_index = 0
    for i in range(len(t_new)):
        if t_new[i] < 0:
            start_index = i+1
    t_new = t_new[start_index:]
    thrust_new = thrust_new[start_index:]
    N_new = len(t_new)

Data_new = stack((t_new, thrust_new)).T

### Desenhar os graficos para teste/comparação
#from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, grid, show
#from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, grid
if Plot_thrust:
    from matplotlib.pyplot import subplots
    fig, axs = subplots(2, 1, sharex = True)
    #axs[0].plot(t, Raw_data[:,1])  # measured signal, in volts
    axs[0].plot(t, Raw_data)  # measured signal, in volts
    axs[0].set_ylabel("signal (V)")
    axs[0].grid()
    axs[1].plot(t, thrust)  # partially processed signal (N)
    axs[1].set_ylabel("thrust (N)")
    axs[1].set_xlabel("time (s)")
    axs[1].plot(t_new, thrust_new, "o-", ms = 3, label = "processed")  # processed signal
    #axs[1].plot(t_new+t_offset, thrust_new, "o-", ms = 3, label = "resampled")  # resampled signal without time shift
    axs[1].grid()
    axs[1].legend()


### Calculando alguns dados úteis:
#zero_thrust = 0.8 ### Valor de empuxo que pode ser considerado desprezível (ruído), em N
zero_thrust = 3.0 ### Valor de empuxo que pode ser considerado desprezível (ruído), em N
pulse_started = False
i = 1 
while (i<N_new) and not(pulse_started):
    if (thrust_new[i] > zero_thrust) and (thrust_new[i-1] > zero_thrust):
        start_index_new = i
        start_time_new = t_new[i]
        pulse_started = True
    i += 1
pulse_ended = False
i = N_new-3 
while (i>0) and not(pulse_ended):
    if (thrust_new[i] > zero_thrust) and (thrust_new[i+1] > zero_thrust):
        end_index_new = i
        end_time_new = t_new[i]
        pulse_ended = True
    i -= 1
if not(pulse_ended):
    end_index_new = N_new-1
    end_time_new = t_new[-1]
pulse_duration = end_time_new - start_time_new
F_max = np.amax(thrust_new)
F_ave = np.mean(thrust_new[start_index_new:end_index_new])
total_impulse = pulse_duration*F_ave
Isp = total_impulse/(g*prop_mass)
print(f"Pulse duration = {pulse_duration:.3} s")
print(f"F_max = {F_max:.3} N")
print(f"F_ave = {F_ave:.3} N")
print(f"I_ave = {total_impulse:.3} N.s")
print(f"Isp = {Isp:.3} (1/s)")

###

### Salvar arquivo de saída
#Header += f"{Code} {diameter} {lenght} {stage_delays} {prop_mass} {total_mass} {manufacturer_name}\n"
Header += f"{Code} {diameter} {lenght} {stage_delays} {prop_mass} {total_mass} {manufacturer_name}"
#savetxt(output_filename, Data_new, fmt='%.4g', delimiter = ' ', newline = "\n ", header = Header, footer = Footer, comments=";")
savetxt(output_filename, Data_new, fmt='%.4g', delimiter = ' ', newline = "\n ", header = Header, comments="")

print(f"ENG file saved as: {output_filename}")



